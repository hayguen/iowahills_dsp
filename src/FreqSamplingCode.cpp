
#include <iowahills/FreqSamplingCode.h>
#include <iowahills/FIRFilterCode.h> // for the definition of TFIRPassTypes

#include <math.h>
#include <iowahills/CplxDMath.hpp>
#include <iowahills/FFTCode.h>
#include <iowahills/LowPassPrototypes.h>

/*
 By Daniel Klostermann
 Iowa Hills Software, LLC  IowaHills.com
 If you find a problem, please leave a note at:
 http://www.iowahills.com/feedbackcomments.html
 May 1, 2016

 This code generates an FIR filter with the frequency over-sampling method. By this we
 mean that we always sample the frequency domain at least 1024 times (a large power of
 2 suitable for an FFT) even though we typically want a relatively small number of FIR
 coefficients (< 50).

 Most authors use N frequency samples to generate N taps. While valid, this tends to loose a
 significant amount of frequency domain information if the tap count is small.

 Using a large number of samples will generate a filter equivalent to those generated with
 the equations for a classic Rectangular Windowed FIR filter. Those equations were derived
 by using an infinite number of frequency samples, which is to say, performing an integral.

 To see how well this works, use this code to sample a simple rectanglular low pass response
 and compare the results to the coefficients generated by Windowed FIR code used the
 BasicFIR() function in FIRFilterCode.cpp

 See the example code in FilterKitMain for an example of using SampledFreqFIR().
 This function is called after the HofSReal array is filled with the desired magnitude response.
 The only trick to this method is getting the phase set correctly. There are two aspects to
 setting the phase. Using the correct slope, which is the same for all filters, but does
 depend on even or odd tap counts. And we must ensure that the phase value is correct at
 Omega = 0.0 and Pi.

 If the filter has a low pass response (magnitude != 0.0 at DC), then the phase at Omega=0 must
 be zero. If the filter has a high pass response (magnitude != 0.0 at Pi), then the phase
 must be zero of Omega = Pi.

 A band pass filter, which has neither a low or high pass response, can have any phase value at
 Omega = 0 and Pi. But a Notch filter, which has both a low and high pass response, must have
 zero phase at both zero and Pi. The code below should make this more clear.

 NumTaps     Number of FIR taps
 FirCoeff    The output array for the FIR coefficients.
 HofSReal    The array containing the desired magnitude response.
 HofSImag    The imag part of the response (zero when this function is called).
 OmegaC      The center frequency. (Only needed for notch filters.)
 PassType    firLPF, firHPF, firBPF, firNOTCH  Needed to set the phase properly. (defined in FIRFilterCode.h)
*/

int SampledFreqFIR(int NumTaps, double *FirCoeff, double *HofSReal, double *HofSImag, double OmegaC, TFIRPassTypes PassType)
{
 int j, CenterJ, NumSamples, StartJ;
 double dNumSamples, RadPerSample, Arg;
 NumSamples = NUM_POS_FREQ_SAMPLES;
 dNumSamples = (double)NumSamples;

 // Limit test NumTaps
 if(NumTaps > MAX_NUMTAPS)NumTaps = MAX_NUMTAPS;
 if(NumTaps > 2*NUM_POS_FREQ_SAMPLES)NumTaps = 2*NUM_POS_FREQ_SAMPLES;


 // Set the slope of the phase.
 RadPerSample = -M_PI_2 * (2.0*dNumSamples - 1.0)/dNumSamples;    // Even tap count.
 if(NumTaps % 2 == 1)RadPerSample = -M_PI;                        // Odd tap count.

 // Set the phase according to the type of response.
 switch(PassType)
  {
   case firLPF:   // Low pass and band pass   phase = 0 at DC
   case firBPF:
	for(j=0; j<NumSamples; j++)
	 {
	  Arg = RadPerSample * (double)j;  // For band pass filters ONLY, an arbitrary amount of phase can be added to Arg. e.g. Add Pi/2 to generate a Hilbert filter, or +/- Pi/4 to 2 different filters to generate a pair of 45 degree Hilberts.
      HofSImag[j] = HofSReal[j] * sin(Arg);
      HofSReal[j] = HofSReal[j] * cos(Arg);
	 }
	break;

   case firHPF:   // High pass   phase = 0 at Pi
	for(j=NumSamples; j>=0; j--)
	 {
	  Arg = RadPerSample * (double)(j-NumSamples);
      HofSImag[j] = HofSReal[j] * sin(Arg);
      HofSReal[j] = HofSReal[j] * cos(Arg);
	 }
	break;

   case firNOTCH:  // Notch   phase = 0 at DC and Pi
    CenterJ = (int)(OmegaC * dNumSamples);
    for(j=0; j<=CenterJ; j++)
     {
      Arg = RadPerSample * (double)j;
      HofSImag[j] = HofSReal[j] * sin(Arg);
      HofSReal[j] = HofSReal[j] * cos(Arg);
     }
    for(j=NumSamples; j>=CenterJ; j--)
     {
      Arg = RadPerSample * (double)(j-NumSamples);
      HofSImag[j] = HofSReal[j] * sin(Arg);
      HofSReal[j] = HofSReal[j] * cos(Arg);
     }
   break;

   default:
   case firALLPASS:
   case firNOT_FIR:
    return 1;  // signal error
  }

 // Fill the negative frequency bins of HofS with the conjugate of HofS for the FFT.
 for(j=1; j<NumSamples; j++)HofSReal[2*NumSamples-j] = HofSReal[j];
 for(j=1; j<NumSamples; j++)HofSImag[2*NumSamples-j] = -HofSImag[j];


 // The Fourier Transform requires the center freq bins to be 0 for LPF and BPF, 1 for HPF and Notch
 if(PassType == firLPF || PassType == firBPF)
  {
   HofSReal[NumSamples] = 0.0;
   HofSImag[NumSamples] = 0.0;
  }
 else
  {
   HofSReal[NumSamples] = 1.0;
   HofSImag[NumSamples] = 0.0;
  }

 // Do an inverse FFT on HofS to generate the impulse response. On return, HofSImag will be zero.
 FFT(HofSReal, HofSImag, 2*NumSamples, INVERSE);

 // We just generated an impulse response that is 2*NumSamples long. Since we used linear phase
 // the response will be symmetric about the center. In general, we only need a small number
 // of these taps, so we use the taps from the center of HofSReal, starting at StartJ.
 // We also need to scale the FFT's output by the size of the FFT.
 StartJ = NumSamples - NumTaps/2;
 for(j=0; j<NumTaps; j++)FirCoeff[j] = HofSReal[StartJ+j] / (2.0*dNumSamples);

  return 0;
}
//---------------------------------------------------------------------------

// This function shows how to sample an analog transfer function H(s) to generate an FIR filter.
// We didn't put much effort into this example because, in general, an analog prototype generates
// a rather poor FIR filter in the sense that it requires such a large number of taps to realize
// the response. For a discussion on this, see this article:
// http://iowahills.com/B2PolynomialFIRFilters.html
// In this example, we generate an FIR low pass from an Inverse Chebyshev low pass prototype.
// We sample the low pass prototype H(s) at frequencies determined by the bilinear transform.
void SampledFreqAnalog(int NumTaps, double *FirCoeff, double *HofSReal, double *HofSImag, double OmegaC)
{
 int j, k, NumSamples;
 double dNumSamples, Omega, Omega0;
 CplxD s, H;
 TFIRPassTypes PassType = firLPF;
 NumSamples = NUM_POS_FREQ_SAMPLES;
 dNumSamples = (double)NumSamples;

  // Limit test NumTaps
 if(NumTaps > MAX_NUMTAPS)NumTaps = MAX_NUMTAPS;
 if(NumTaps > 2*NUM_POS_FREQ_SAMPLES)NumTaps = 2*NUM_POS_FREQ_SAMPLES;

 // Define the low pass filter prototype
 TLowPassParams LPFProto;            // defined in LowPassPrototypes.h
 TSPlaneCoeff SCoeff;                // defined in LowPassPrototypes.h
 LPFProto.ProtoType = INVERSE_CHEBY; // BUTTERWORTH, CHEBYSHEV, GAUSSIAN, BESSEL, ADJUSTABLE, INVERSE_CHEBY, PAPOULIS, ELLIPTIC  (defined in LowPassPrototypes.h)
 LPFProto.NumPoles = 8;              // 1 <= NumPoles <= 12, 15, 20 Depending on the filter.
 LPFProto.Ripple = 0.25;             // 0.0 <= Ripple <= 1.0 dB     Chebyshev and Elliptic (less for high order Chebyshev).
 LPFProto.StopBanddB = 60.0;         // 20 <= StopBand <= 120 dB    Inv Cheby and Elliptic
 LPFProto.Gamma = 0.0;               // -1.0 <= Gamma <= 1.0        Adjustable Gauss  Controls the transition BW.

 // Get the prototype filter's 2nd order s plane coefficients.
 SCoeff = CalcLowPassProtoCoeff(LPFProto);

 // Evaluate the prototype's H(s)
 Omega0 = 1.0 / tan(OmegaC * M_PI_2);  // This sets the corner frequency.
 for(j=0; j<NumSamples; j++)
  {
   Omega = Omega0 * tan(M_PI_2 * (double)j / dNumSamples); // Frequencies per the bilinear transform.
   s = CplxD(0.0, Omega);
   H = CplxD(1.0,0.0);
   for(k=0; k<SCoeff.NumSections; k++)
	{
	 H *=  SCoeff.N2[k] * s * s + SCoeff.N1[k] * s + SCoeff.N0[k]; // The numerator
     H /=  SCoeff.D2[k] * s * s + SCoeff.D1[k] * s + SCoeff.D0[k]; // The denominator
     H *= SCoeff.D0[k] / SCoeff.N0[k];                             // The gain constants.
	}
   HofSReal[j] = H.re;  // We need to do this for the FFT, which uses real arrays.
   HofSImag[j] = H.im;
  }

 // Fill the negative frequency bins of HofS with the conjugate of HofS for the FFT.
 for(j=1; j<NumSamples; j++)HofSReal[2*NumSamples-j] = HofSReal[j];
 for(j=1; j<NumSamples; j++)HofSImag[2*NumSamples-j] = -HofSImag[j];

 // The Fourier Transform requires the center freq bins to be 0 for LPF and BPF, 1 for HPF and Notch
 if(PassType == firLPF || PassType == firBPF)
  {
   HofSReal[NumSamples] = 0.0;
   HofSImag[NumSamples] = 0.0;
  }
 else
  {
   HofSReal[NumSamples] = 1.0;
   HofSImag[NumSamples] = 0.0;
  }

 // Do an inverse FFT on HofS to generate the impulse response. On return, HofSImag will be zero.
 FFT(HofSReal, HofSImag, 2*NumSamples, INVERSE);

 // We just generated an impulse response that is 2*NumSamples long. We can use as many of these
 // as we like, but if you look at HofSReal you will see that the tail of the response goes to
 // zero rather quickly. The desired part of impulse response is at the beginning of HofSReal
 // instead of the center because H(s) didn't have linear phase (more like minimum phase).

 // Take the taps from the start of HofSReal and scale by the size of the FFT.
 for(j=0; j<NumTaps; j++)FirCoeff[j] = HofSReal[j] / (2.0*dNumSamples);

 // Most FIR responses benefit from the application of a window to reduce the effects of
 // truncating the impulse response. But a typical window, such as the Kaiser, can't be used
 // because this impulse response isn't symmetric about the center tap.
 // A Half Cosine window works quite nicely with these types of responses.
 for(j=0; j<NumTaps; j++)FirCoeff[j] = FirCoeff[j] * cos(M_PI_2 * j / (double)NumTaps);
}

//---------------------------------------------------------------------------