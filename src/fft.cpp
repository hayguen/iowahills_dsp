/*
 This software is part of iowahills_dsp, a set of DSP routines under MIT License.
 2016 By Daniel Klostermann, Iowa Hills Software, LLC  IowaHills.com
 Copyright (c) 2021  Hayati Ayguen <h_ayguen@web.de>
 All rights reserved.
*/

#include <iowahills/fft.h>

#include "non_throwing_vector.h"
#include <cmath>

//---------------------------------------------------------------------------
#define M_2PI     6.28318530717958647692   // 2*Pi
#define M_SQRT_2  0.707106781186547524401  // sqrt(2)/2
#define MAXIMUM_FFT_SIZE  1048576
#define MINIMUM_FFT_SIZE  8

//---------------------------------------------------------------------------
static void ReArrangeInput(double *InputR, double *InputI, double *BufferR, double *BufferI, int *RevBits, int N);
static void FillTwiddleArray(double *TwiddleR, double *TwiddleI, int N, TTransFormType Type);
static void Transform(double *InputR, double *InputI, double *BufferR, double *BufferI, const double *TwiddleR, const double *TwiddleI, int N);

//---------------------------------------------------------------------------
// This calculates the required FFT size for a given number of points.
int RequiredFFTSize(int NumPts)
{
 int N = MINIMUM_FFT_SIZE;
 while(N < NumPts && N < MAXIMUM_FFT_SIZE)
  {
   N *= 2;
  }
 return N;
}

//---------------------------------------------------------------------------

// This verifies that the FFT Size N = 2^M.   exponent M is returned - or 0 on error
// N must be >= 8 for the Twiddle calculations
int IsValidFFTSize(int N)
{
 if(N < MINIMUM_FFT_SIZE || N > MAXIMUM_FFT_SIZE || (N & (N - 1)) != 0)
     return(0);   // N & (N - 1) ensures a power of 2
 return ( (int)( log((double)N) / M_LN2 + 0.5 ) );         // return M where N = 2^M
}

//---------------------------------------------------------------------------

struct FFTBuffers
{
    FFTBuffers(int _N, TTransFormType _Type);

    bool good() const {
        return (N && BufferR && BufferI && TwiddleR && TwiddleI && RevBits);
    }

    const int LogTwoOfN;
    const int N;
    const TTransFormType Type;
    NonThrowingVector<double> BufferR;
    NonThrowingVector<double> BufferI;
    NonThrowingVector<double> TwiddleR;
    NonThrowingVector<double> TwiddleI;
    NonThrowingVector<int> RevBits;
};

FFTBuffers::FFTBuffers(int _N, TTransFormType _Type)
    : LogTwoOfN( IsValidFFTSize(_N) )
    , N( ( LogTwoOfN && (_Type == FORWARD || _Type == INVERSE) ) ? _N : 0 )    // Verify the FFT size and type
    , Type(_Type)
    , BufferR(N)
    , BufferI(N)
    , TwiddleR(N/2)
    , TwiddleI(N/2)
    , RevBits(N)
{
    if (good())
        FillTwiddleArray(TwiddleR, TwiddleI, _N, _Type);
}

FFTBuffers * prepareFFT(int N, TTransFormType Type)
{
    FFTBuffers *p = new(std::nothrow)  FFTBuffers(N, Type);
    if (p && p->good())
        return p;
    if (p)
        delete p;
    return 0;
}

void disposeFFT(FFTBuffers *buffers)
{
    if (buffers)
        delete buffers;
}

//---------------------------------------------------------------------------

int FFTex(FFTBuffers *buffers, double *InputR, double *InputI, int skip_fwd_div_N)
{
    if (!buffers)
        return 1;

    FFTBuffers &b = *buffers;
    const int N = b.N;
    const int LogTwoOfN = b.LogTwoOfN;
    const TTransFormType Type = b.Type;

    ReArrangeInput(InputR, InputI, b.BufferR, b.BufferI, b.RevBits, N);
    Transform(InputR, InputI, b.BufferR, b.BufferI, b.TwiddleR, b.TwiddleI, N);

    // The ReArrangeInput function swapped Input[] and Buffer[]. Then Transform()
    // swapped them again, LogTwoOfN times. Ultimately, these swaps must be done
    // an even number of times, or the pointer to Buffer gets returned.
    // So we must do one more swap here, for N = 16, 64, 256, 1024, ...
    const double OneOverN = (Type == FORWARD) ? (1.0 / (double)N) : 1.0;
    int j;

    if (skip_fwd_div_N || Type != FORWARD)
    {
      if(LogTwoOfN % 2 == 0)
      {
        for(j=0; j<N; j++) InputR[j] = b.BufferR[j];
        for(j=0; j<N; j++) InputI[j] = b.BufferI[j];
      }
    }
    else
    {
      if(LogTwoOfN % 2 == 1)
       {
        for(j=0; j<N; j++) InputR[j] = InputR[j] * OneOverN;
        for(j=0; j<N; j++) InputI[j] = InputI[j] * OneOverN;
       }
      else // if(LogTwoOfN % 2 == 0) then the results are still in Buffer.
       {
        for(j=0; j<N; j++) InputR[j] = b.BufferR[j] * OneOverN;
        for(j=0; j<N; j++) InputI[j] = b.BufferI[j] * OneOverN;
       }
    }

    return 0;
}



//---------------------------------------------------------------------------

// Fast Fourier Transform
// This code puts DC in bin 0 and scales the output of a forward transform by 1/N.
// InputR and InputI are the real and imaginary input arrays of length N.
// The output values are returned in the Input arrays.
// TTransFormType is either FORWARD or INVERSE (defined in the header file)
// 256 pts in 50 us
int FFT(double *InputR, double *InputI, int N, TTransFormType Type)
{
 // Memory allocation for all the arrays.
 FFTBuffers b(N, Type);

 if (!b.good())
  {
   // ShowMessage("FFT Memory Allocation Error");
   return 1;
  }

 const int skip_fwd_div_N = 0;
 return FFTex(&b, InputR, InputI, skip_fwd_div_N);
}
//---------------------------------------------------------------------------

// This puts the input arrays in bit reversed order.
// The while loop generates an array of bit reversed numbers. e.g.
// e.g. N=8: RevBits = 0,4,2,6,1,5,3,7   N=16: RevBits = 0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15
static
void ReArrangeInput(double *InputR, double *InputI, double *BufferR, double *BufferI, int *RevBits, int N)
{
 int j, k, J, K;

 J = N/2;
 K = 1;
 RevBits[0] = 0;
 while(J >= 1)
  {
   for(k=0; k<K; k++)
	{
	 RevBits[k+K] = RevBits[k] + J;
	}
   K *= 2;
   J /= 2;
  }


 // Move the rearranged input values to Buffer.
 // Take note of the pointer swaps at the top of the transform algorithm.
 for(j=0; j<N; j++)
  {
   BufferR[j] = InputR[ RevBits[j] ];
   BufferI[j] = InputI[ RevBits[j] ];
  }

}

//---------------------------------------------------------------------------

/*
 The Pentium takes a surprising amount of time to calculate the sine and cosine.
 You may want to make the twiddle arrays static if doing repeated FFTs of the same size.
 This uses 4 fold symmetry to calculate the twiddle factors. As a result, this function
 requires a minimum FFT size of 8.
*/
static
void FillTwiddleArray(double *TwiddleR, double *TwiddleI, int N, TTransFormType Type)
{
 int j;
 double Theta, TwoPiOverN;

 TwoPiOverN = M_2PI / (double) N;

 if(Type == FORWARD)
  {
   TwiddleR[0] = 1.0;
   TwiddleI[0] = 0.0;
   TwiddleR[N/4] = 0.0;
   TwiddleI[N/4] = -1.0;
   TwiddleR[N/8] = M_SQRT_2;
   TwiddleI[N/8] = -M_SQRT_2;
   TwiddleR[3*N/8] = -M_SQRT_2;
   TwiddleI[3*N/8] = -M_SQRT_2;
   for(j=1; j<N/8; j++)
	{
	 Theta = (double)j * -TwoPiOverN;
	 TwiddleR[j] = cos(Theta);
	 TwiddleI[j] = sin(Theta);
	 TwiddleR[N/4-j] = -TwiddleI[j];
	 TwiddleI[N/4-j] = -TwiddleR[j];
	 TwiddleR[N/4+j] = TwiddleI[j];
	 TwiddleI[N/4+j] = -TwiddleR[j];
	 TwiddleR[N/2-j] = -TwiddleR[j];
	 TwiddleI[N/2-j] = TwiddleI[j];
	}
  }

 else
  {
   TwiddleR[0] = 1.0;
   TwiddleI[0] = 0.0;
   TwiddleR[N/4] = 0.0;
   TwiddleI[N/4] = 1.0;
   TwiddleR[N/8] = M_SQRT_2;
   TwiddleI[N/8] = M_SQRT_2;
   TwiddleR[3*N/8] = -M_SQRT_2;
   TwiddleI[3*N/8] = M_SQRT_2;
   for(j=1; j<N/8; j++)
	{
	 Theta = (double)j * TwoPiOverN;
	 TwiddleR[j] = cos(Theta);
	 TwiddleI[j] = sin(Theta);
	 TwiddleR[N/4-j] = TwiddleI[j];
	 TwiddleI[N/4-j] = TwiddleR[j];
	 TwiddleR[N/4+j] = -TwiddleI[j];
	 TwiddleI[N/4+j] = TwiddleR[j];
	 TwiddleR[N/2-j] = -TwiddleR[j];
	 TwiddleI[N/2-j] = TwiddleI[j];
	}
  }

}

//---------------------------------------------------------------------------

// The Fast Fourier Transform.
static
void Transform(double *InputR, double *InputI, double *BufferR, double *BufferI, const double *TwiddleR, const double *TwiddleI, int N)
{
 int j, k, J, K, I, T;
 double *TempPointer;
 double TempR, TempI;

 J = N/2;     // J increments down to 1
 K = 1;       // K increments up to N/2
 while(J > 0) // Loops Log2(N) times.
  {
   // Swap pointers, instead doing this: for(j=0; j<N; j++) Input[j] = Buffer[j];
   // We start with a swap because of the swap in ReArrangeInput.
   TempPointer = InputR;
   InputR = BufferR;
   BufferR = TempPointer;
   TempPointer = InputI;
   InputI = BufferI;
   BufferI = TempPointer;

   I = 0;
   for(j=0; j<J; j++)
	{
	 T = 0;
	 for(k=0; k<K; k++) // Loops N/2 times for every value of J and K
	  {
	   TempR = InputR[K+I] * TwiddleR[T] - InputI[K+I] * TwiddleI[T];
	   TempI = InputR[K+I] * TwiddleI[T] + InputI[K+I] * TwiddleR[T];
	   BufferR[I]   = InputR[I] + TempR;
	   BufferI[I]   = InputI[I] + TempI;
	   BufferR[I+K] = InputR[I] - TempR;
	   BufferI[I+K] = InputI[I] - TempI;
	   I++;
	   T += J;
	  }
	 I += K;
	}
   K *= 2;
   J /= 2;
  }

}

//-----------------------------------------------------------------------------------------------

/*
 The only difficulty in writing an FFT is figuring out how to calculate the various array indexes.
 This shows how the index values change when doing a 16 pt FFT.
 This print only has value if you compare it to a butterfly chart. Then you can
 see how an FFT works. Use a 16 point decimation in time butterfly chart. We have one here:
 http://www.iowahills.com/FFTCode.html

 Note: The code above uses real variables. This print out came from code using complex variables as shown here.
 Buffer[I]   = Input[I] + Input[K+I] * Twiddle[T];
 Buffer[I+K] = Input[I] - Input[K+I] * Twiddle[T];

 N = 16
 J = 8    K = 1
 Buffer[0]  = Input[0]  + Input[1]  * Twiddle[0]   I = 0
 Buffer[1]  = Input[0]  - Input[1]  * Twiddle[0]   I = 0
 Buffer[2]  = Input[2]  + Input[3]  * Twiddle[0]   I = 2
 Buffer[3]  = Input[2]  - Input[3]  * Twiddle[0]   I = 2
 Buffer[4]  = Input[4]  + Input[5]  * Twiddle[0]   etc.
 Buffer[5]  = Input[4]  - Input[5]  * Twiddle[0]
 Buffer[6]  = Input[6]  + Input[7]  * Twiddle[0]
 Buffer[7]  = Input[6]  - Input[7]  * Twiddle[0]
 Buffer[8]  = Input[8]  + Input[9]  * Twiddle[0]
 Buffer[9]  = Input[8]  - Input[9]  * Twiddle[0]
 Buffer[10] = Input[10] + Input[11] * Twiddle[0]
 Buffer[11] = Input[10] - Input[11] * Twiddle[0]
 Buffer[12] = Input[12] + Input[13] * Twiddle[0]
 Buffer[13] = Input[12] - Input[13] * Twiddle[0]
 Buffer[14] = Input[14] + Input[15] * Twiddle[0]
 Buffer[15] = Input[14] - Input[15] * Twiddle[0]

 J = 4    K = 2
 Buffer[0]  = Input[0]  + Input[2]  * Twiddle[0]
 Buffer[2]  = Input[0]  - Input[2]  * Twiddle[0]
 Buffer[1]  = Input[1]  + Input[3]  * Twiddle[4]
 Buffer[3]  = Input[1]  - Input[3]  * Twiddle[4]
 Buffer[4]  = Input[4]  + Input[6]  * Twiddle[0]
 Buffer[6]  = Input[4]  - Input[6]  * Twiddle[0]
 Buffer[5]  = Input[5]  + Input[7]  * Twiddle[4]
 Buffer[7]  = Input[5]  - Input[7]  * Twiddle[4]
 Buffer[8]  = Input[8]  + Input[10] * Twiddle[0]
 Buffer[10] = Input[8]  - Input[10] * Twiddle[0]
 Buffer[9]  = Input[9]  + Input[11] * Twiddle[4]
 Buffer[11] = Input[9]  - Input[11] * Twiddle[4]
 Buffer[12] = Input[12] + Input[14] * Twiddle[0]
 Buffer[14] = Input[12] - Input[14] * Twiddle[0]
 Buffer[13] = Input[13] + Input[15] * Twiddle[4]
 Buffer[15] = Input[13] - Input[15] * Twiddle[4]

 J = 2    K = 4
 Buffer[0]  = Input[0]  + Input[4]  * Twiddle[0]
 Buffer[4]  = Input[0]  - Input[4]  * Twiddle[0]
 Buffer[1]  = Input[1]  + Input[5]  * Twiddle[2]
 Buffer[5]  = Input[1]  - Input[5]  * Twiddle[2]
 Buffer[2]  = Input[2]  + Input[6]  * Twiddle[4]
 Buffer[6]  = Input[2]  - Input[6]  * Twiddle[4]
 Buffer[3]  = Input[3]  + Input[7]  * Twiddle[6]
 Buffer[7]  = Input[3]  - Input[7]  * Twiddle[6]
 Buffer[8]  = Input[8]  + Input[12] * Twiddle[0]
 Buffer[12] = Input[8]  - Input[12] * Twiddle[0]
 Buffer[9]  = Input[9]  + Input[13] * Twiddle[2]
 Buffer[13] = Input[9]  - Input[13] * Twiddle[2]
 Buffer[10] = Input[10] + Input[14] * Twiddle[4]
 Buffer[14] = Input[10] - Input[14] * Twiddle[4]
 Buffer[11] = Input[11] + Input[15] * Twiddle[6]
 Buffer[15] = Input[11] - Input[15] * Twiddle[6]

 J = 1    K = 8
 Buffer[0]  = Input[0]  + Input[8]  * Twiddle[0]
 Buffer[8]  = Input[0]  - Input[8]  * Twiddle[0]
 Buffer[1]  = Input[1]  + Input[9]  * Twiddle[1]
 Buffer[9]  = Input[1]  - Input[9]  * Twiddle[1]
 Buffer[2]  = Input[2]  + Input[10] * Twiddle[2]
 Buffer[10] = Input[2]  - Input[10] * Twiddle[2]
 Buffer[3]  = Input[3]  + Input[11] * Twiddle[3]
 Buffer[11] = Input[3]  - Input[11] * Twiddle[3]
 Buffer[4]  = Input[4]  + Input[12] * Twiddle[4]
 Buffer[12] = Input[4]  - Input[12] * Twiddle[4]
 Buffer[5]  = Input[5]  + Input[13] * Twiddle[5]
 Buffer[13] = Input[5]  - Input[13] * Twiddle[5]
 Buffer[6]  = Input[6]  + Input[14] * Twiddle[6]
 Buffer[14] = Input[6]  - Input[14] * Twiddle[6]
 Buffer[7]  = Input[7]  + Input[15] * Twiddle[7]
 Buffer[15] = Input[7]  - Input[15] * Twiddle[7]

*/

//-----------------------------------------------------------------------------------------------


// Discrete Fourier Transform ( textbook code )
// This takes the same arguments as the FFT function.
// 256 pts in 1.720 ms
int DFT(double *InputR, double *InputI, int N, TTransFormType Type)
{
 int j, k, n;
 double Sign, Arg;

 NonThrowingVector<double> SumR(N);
 NonThrowingVector<double> SumI(N);
 NonThrowingVector<double> TwiddleReal(N);
 NonThrowingVector<double> TwiddleImag(N);

 if(!SumR || !SumI ||
    !TwiddleReal || !TwiddleImag ||
    (Type != FORWARD && Type != INVERSE) )
  {
   // ShowMessage("Incorrect DFT Type or unable to allocate memory");
   return 1;
  }

 // Calculate the twiddle factors and initialize the Sum arrays.
 if(Type == FORWARD)Sign = -1.0;
 else               Sign =  1.0;

 for(j=0; j<N; j++)
  {
   Arg = M_2PI * (double)j / (double)N;
   TwiddleReal[j] = cos(Arg);
   TwiddleImag[j] = Sign * sin(Arg);
   SumR[j] = SumI[j] = 0.0;
  }


 //Calculate the DFT
 for(j=0; j<N; j++) // Sum index
 for(k=0; k<N; k++) // Input index
  {
   n = (j*k) % N;
   SumR[j] += TwiddleReal[n] * InputR[k] - TwiddleImag[n] * InputI[k];
   SumI[j] += TwiddleReal[n] * InputI[k] + TwiddleImag[n] * InputR[k];
  }

 // Scale the result if doing a forward DFT, and move the result to the input arrays.
 if(Type == FORWARD)
  {
   for(j=0; j<N; j++)InputR[j] = SumR[j] / (double)N;
   for(j=0; j<N; j++)InputI[j] = SumI[j] / (double)N;
  }
 else  // Inverse DFT
  {
   for(j=0; j<N; j++)InputR[j] = SumR[j];
   for(j=0; j<N; j++)InputI[j] = SumI[j];
  }

 return 0;
}

//-----------------------------------------------------------------------------------------------

// This is a DFT for real valued samples. Since Samples is real, it can only do a forward DFT.
// The results are returned in OutputR  OutputI
// 256 pts in 700 us    30 us to calc the twiddles
int RealSigDFT(const double *Samples, double *OutputR, double *OutputI, int N)
{
 int j, k;
 double Arg;

 NonThrowingVector<double> TwiddleReal(N);
 NonThrowingVector<double> TwiddleImag(N);
 if(!TwiddleReal || !TwiddleImag)
  {
   // ShowMessage("Failed to allocate memory in RealSigDFT");
   return 1;
  }

 for(j=0; j<N; j++)
  {
   Arg = M_2PI * (double)j / (double)N;
   TwiddleReal[j] = cos(Arg);
   TwiddleImag[j] = -sin(Arg);
  }


 // Compute the DFT.
 // We have a real input, so only do the pos frequencies. i.e. j<N/2
 for(j=0; j<=N/2; j++)
  {
   OutputR[j] = 0.0;
   OutputI[j] = 0.0;
   for(k=0; k<N; k++)
	{
     OutputR[j] += Samples[k] * TwiddleReal[(j*k) % N];
     OutputI[j] += Samples[k] * TwiddleImag[(j*k) % N];
	}

   // Scale the result
   OutputR[j] /= (double)N;
   OutputI[j] /= (double)N;
  }

 // The neg freq components are the conj of the pos components because the input signal is real.
 for(j=1; j<N/2; j++)
  {
   OutputR[N-j] = OutputR[j];
   OutputI[N-j] = -OutputI[j];
  }

 return 0;
}

//---------------------------------------------------------------------------

// This is a single frequency DFT.
// This code uses iteration to calculate the Twiddle factors.
// To evaluate the frequency response of an FIR filter at Omega, set
// Samples[] = FirCoeff[]   N = NumTaps  0.0 <= Omega <= 1.0
// 256 pts in 15.6 us
double SingleFreqDFT(const double *Samples, int N, double Omega)
{
 int k;
 double SumR, SumI, zR, zI, TwiddleR, TwiddleI, Temp;

 TwiddleR =  cos(Omega * M_PI);
 TwiddleI = -sin(Omega * M_PI);
 zR = 1.0;    // z, as in e^(j*omega)
 zI = 0.0;
 SumR = 0.0;
 SumI = 0.0;

 for(k=0; k<N; k++)
  {
   SumR += Samples[k] * zR;
   SumI += Samples[k] * zI;

   // Calculate the complex exponential z by taking it to the kth power.
   Temp = zR * TwiddleR - zI * TwiddleI;
   zI =   zR * TwiddleI + zI * TwiddleR;
   zR = Temp;
  }

 /*
 // This is the more conventional implementation of the loop above.
 // It is a bit more accurate, but slower.
 for(k=0; k<N; k++)
  {
   SumR += Samples[k] *  cos((double)k * Omega * M_PI);
   SumI += Samples[k] * -sin((double)k * Omega * M_PI);
  }
 */

 return( sqrt(SumR*SumR + SumI*SumI) );
 // return( ComplexD(SumR, SumI) );// if phase is needed.
}

//---------------------------------------------------------------------------

// This shows how to adjust the delay of an FIR by a fractional amount.
// We take the FFT of the FIR coefficients to get to the frequency domain,
// then apply the Laplace delay operator, and then do an inverse FFT.

// Use this function last. i.e. After the window was applied to the coefficients.
// The Delay value is in terms of a fraction of a sample (not in terms of sampling freq).
// Delay may be pos or neg. Typically a filter's delay can be adjusted by +/- NumTaps/20
// without affecting its performance significantly. A typical Delay value would be 0.75
int AdjustFIRDelay(double *FirCoeff, int NumTaps, double Delay)
{
 int j;
 double Arg, Temp;
 int FFTSize = RequiredFFTSize(NumTaps+(int)fabs(Delay)+1); // Zero pad by at least Delay + 1 to prevent the impulse response from wrapping around.

 NonThrowingVector<double> FFTInputR(FFTSize);  // Real part
 NonThrowingVector<double> FFTInputI(FFTSize);  // Imag part
 if(!FFTInputR || !FFTInputI)
  {
   //ShowMessage("Unable to allocate memory in AdjustDelay");
   return 1;
  }
 for(j=0; j<FFTSize; j++)FFTInputR[j] = FFTInputI[j] = 0.0; // A mandatory init.
 for(j=0; j<NumTaps; j++)FFTInputR[j] = FirCoeff[j];        // Fill the real part with the FIR coeff.

 FFT(FFTInputR, FFTInputI, FFTSize, FORWARD);              // Do an FFT
 for(j=0; j<=FFTSize/2; j++)                               // Apply the Laplace Delay operator e^(-j*omega*Delay).
  {
   Arg = -Delay * (double)j / (double)FFTSize * M_2PI;     // This is -Delay * (the FFT bin frequency).
   Temp =         cos(Arg)*FFTInputR[j] - sin(Arg)*FFTInputI[j];
   FFTInputI[j] = cos(Arg)*FFTInputI[j] + sin(Arg)*FFTInputR[j];
   FFTInputR[j] = Temp;
  }
 for(j=1; j<FFTSize/2; j++) // Fill the neg freq bins with the conjugate values.
  {
   FFTInputR[FFTSize-j] = FFTInputR[j];
   FFTInputI[FFTSize-j] = -FFTInputI[j];
  }

 FFT(FFTInputR, FFTInputI, FFTSize, INVERSE); // Inverse FFT
 for(j=0; j<NumTaps; j++)
  {
   FirCoeff[j] = FFTInputR[j];
  }

 return 0;
}
//-----------------------------------------------------------------------------

