/*
 This software is part of iowahills_dsp, a set of DSP routines under MIT License.
 2016 By Daniel Klostermann, Iowa Hills Software, LLC  IowaHills.com
 Copyright (c) 2021  Hayati Ayguen <h_ayguen@web.de>
 All rights reserved.


 These are the window definitions. These windows can be used for either
 FIR filter design or with an FFT for spectral analysis.
 For definitions, see this article:  http://en.wikipedia.org/wiki/Window_function

 This function has 6 inputs
 Data is the array, of length N, containing the data to to be windowed.
 This data is either an FIR filter sinc pulse, or the data to be analyzed by an fft.

 WindowType is an enum defined in the header file.
 e.g. wtKAISER, wtSINC, wtHANNING, wtHAMMING, wtBLACKMAN, ...

 Alpha sets the width of the flat top.
 Windows such as the Tukey and Trapezoid are defined to have a variably wide flat top.
 As can be seen by its definition, the Tukey is just a Hanning window with a flat top.
 Alpha can be used to give any of these windows a partial flat top, except the Flattop and Kaiser.
 Alpha = 0 gives the original window. (i.e. no flat top)
 To generate a Tukey window, use a Hanning with 0 < Alpha < 1
 To generate a Bartlett window (triangular), use a Trapezoid window with Alpha = 0.
 Alpha = 1 generates a rectangular window in all cases. (except the Flattop and Kaiser)


 Beta is used with the Kaiser, Sinc, and Sine windows only.
 These three windows are used primarily for FIR filter design. Then
 Beta controls the filter's transition bandwidth and the sidelobe levels.
 All other windows ignore Beta.

 UnityGain controls whether the gain of these windows is set to unity.
 Only the Flattop window has unity gain by design. The Hanning window, for example, has a gain
 of 1/2.  UnityGain = true  sets the gain to 1, which preserves the signal's energy level
 when these windows are used for spectral analysis.

 Don't use this with FIR filter design however. Since most of the enegy in an FIR sinc pulse
 is in the middle of the window, the window needs a peak amplitude of one, not unity gain.
 Setting UnityGain = true will simply cause the resulting FIR filter to have excess gain.

 If using these windows for FIR filters, start with the Kaiser, Sinc, or Sine windows and
 adjust Beta for the desired transition BW and sidelobe levels (set Alpha = 0).
 While the FlatTop is an excellent window for spectral analysis, don't use it for FIR filter design.
 It has a peak amplitude of ~ 4.7 which causes the resulting FIR filter to have about this much gain.
 It works poorly for FIR filters even if you adjust its peak amplitude.
 The Trapezoid also works poorly for FIR filter design.

 If using these windows with an fft for spectral analysis, start with the Hanning, Gauss, or Flattop.
 When choosing a window for spectral analysis, you must trade off between resolution and amplitude
 accuracy. The Hanning has the best resolution while the Flatop has the best amplitude accuracy.
 The Gauss is midway between these two for both accuracy and resolution. These three were
 the only windows available in the HP 89410A Vector Signal Analyzer. Which is to say, these three
 are the probably the best windows for general purpose signal analysis.
*/

#include <iowahills/windowing.h>
#include <iowahills/fft.h>

#include "non_throwing_vector.h"

#include <cmath>

//---------------------------------------------------------------------------
#define M_2PI     6.28318530717958647692   // 2*Pi

static double iowa_Sinc(double x);
double iowa_Bessel(double x);

//---------------------------------------------------------------------------

int WindowData(double *Data, int N, TWindowType WindowType, double Alpha, double Beta, bool UnityGain)
{
 if(WindowType == wtNONE) return 1;

 int j, M, TopWidth;
 double dM;

 if(WindowType == wtKAISER ||  WindowType == wtFLATTOP )Alpha = 0.0;

 if(Alpha < 0.0)Alpha = 0.0;
 if(Alpha > 1.0)Alpha = 1.0;

 if(Beta < 0.0)Beta = 0.0;
 if(Beta > 10.0)Beta = 10.0;

 NonThrowingVector<double> WinCoeff(N+2);
 if(!WinCoeff)
  {
   // ShowMessage("Failed to allocate memory in WindowData() ");
   return 1;
  }

 TopWidth = (int)( Alpha * (double)N );
 if(TopWidth%2 != 0)TopWidth++;
 if(TopWidth > N)TopWidth = N;
 M = N - TopWidth;
 dM = M + 1;


 // Calculate the window for N/2 points, then fold the window over (at the bottom).
 // TopWidth points will be set to 1.
 if(WindowType == wtKAISER)
  {
   double Arg;
   for(j=0; j<M; j++)
	{
	 Arg = Beta * sqrt(1.0 - pow( ((double)(2*j+2) - dM) / dM, 2.0) );
     WinCoeff[j] = iowa_Bessel(Arg) / iowa_Bessel(Beta);
	}
  }

 else if(WindowType == wtSINC)  // Lanczos
  {
   for(j=0; j<M; j++)   WinCoeff[j] = iowa_Sinc((double)(2*j+1-M)/dM * M_PI );
   for(j=0; j<M; j++)   WinCoeff[j] = pow(WinCoeff[j], Beta);
  }

 else if(WindowType == wtSINE)  // Hanning if Beta = 2
  {
   for(j=0; j<M/2; j++) WinCoeff[j] = sin((double)(j+1) * M_PI / dM);
   for(j=0; j<M/2; j++) WinCoeff[j] = pow(WinCoeff[j], Beta);
  }

 else if(WindowType == wtHANNING)
  {
   for(j=0; j<M/2; j++) WinCoeff[j] = 0.5 - 0.5 * cos((double)(j+1) * M_2PI / dM);
  }

 else if(WindowType == wtHAMMING)
  {
   for(j=0; j<M/2; j++)
   WinCoeff[j] = 0.54 - 0.46 * cos((double)(j+1) * M_2PI / dM);
  }

 else if(WindowType == wtBLACKMAN)
  {
   for(j=0; j<M/2; j++)
	{
	 WinCoeff[j] = 0.42
	 - 0.50 * cos((double)(j+1) * M_2PI / dM)
	 + 0.08 * cos((double)(j+1) * M_2PI * 2.0 / dM);
	}
  }


 // Defined at: http://www.bth.se/fou/forskinfo.nsf/0/130c0940c5e7ffcdc1256f7f0065ac60/$file/ICOTA_2004_ttr_icl_mdh.pdf
 else if(WindowType == wtFLATTOP)
  {
   for(j=0; j<=M/2; j++)
	{
	 WinCoeff[j] = 1.0
	 - 1.93293488969227 * cos((double)(j+1) * M_2PI / dM)
	 + 1.28349769674027 * cos((double)(j+1) * M_2PI * 2.0 / dM)
	 - 0.38130801681619 * cos((double)(j+1) * M_2PI * 3.0 / dM)
	 + 0.02929730258511 * cos((double)(j+1) * M_2PI * 4.0 / dM);
	}
  }


 else if(WindowType == wtBLACKMAN_HARRIS)
  {
   for(j=0; j<M/2; j++)
	{
	 WinCoeff[j] = 0.35875
	 - 0.48829 * cos((double)(j+1) * M_2PI / dM)
	 + 0.14128 * cos((double)(j+1) * M_2PI * 2.0 / dM)
	 - 0.01168 * cos((double)(j+1) * M_2PI * 3.0 / dM);
	}
  }

 else if(WindowType == wtBLACKMAN_NUTTALL)
  {
   for(j=0; j<M/2; j++)
	{
	 WinCoeff[j] = 0.3535819
	 - 0.4891775 * cos((double)(j+1) * M_2PI / dM)
	 + 0.1365995 * cos((double)(j+1) * M_2PI * 2.0 / dM)
	 - 0.0106411 * cos((double)(j+1) * M_2PI * 3.0 / dM);
	}
  }

 else if(WindowType == wtNUTTALL)
  {
   for(j=0; j<M/2; j++)
	{
	 WinCoeff[j] = 0.355768
	 - 0.487396 * cos((double)(j+1) * M_2PI / dM)
	 + 0.144232 * cos((double)(j+1) * M_2PI * 2.0 / dM)
	 - 0.012604 * cos((double)(j+1) * M_2PI * 3.0 / dM);
	}
  }

 else if(WindowType == wtKAISER_BESSEL)
  {
   for(j=0; j<=M/2; j++)
	{
	 WinCoeff[j] = 0.402
	 - 0.498 * cos(M_2PI * (double)(j+1) / dM)
	 + 0.098 * cos(2.0 * M_2PI * (double)(j+1) / dM)
	 + 0.001 * cos(3.0 * M_2PI * (double)(j+1) / dM);
	}
  }

 else if(WindowType == wtTRAPEZOID) // Rectangle for Alpha = 1  Triangle for Alpha = 0
  {
   int K = M/2;
   if(M%2)K++;
   for(j=0; j<K; j++)   WinCoeff[j] = (double)(j+1) / (double)K;
  }


 // This definition is from http://en.wikipedia.org/wiki/Window_function (Gauss Generalized normal window)
 // We set their p = 2, and use Alpha in the numerator, instead of Sigma in the denominator, as most others do.
 // Alpha = 2.718 puts the Gauss window response midway between the Hanning and the Flattop (basically what we want).
 // It also gives the same BW as the Gauss window used in the HP 89410A Vector Signal Analyzer.
 else if(WindowType == wtGAUSS)
  {
   for(j=0; j<M/2; j++)
    {
     WinCoeff[j] = ((double)(j+1) - dM/2.0) / (dM/2.0) * 2.7183;
     WinCoeff[j] *= WinCoeff[j];
     WinCoeff[j] = exp(-WinCoeff[j]);
    }
  }

 else // Error.
  {
   // ShowMessage("Incorrect window type in WindowFFTData");
   return 1;
  }

 // Fold the coefficients over.
 for(j=0; j<M/2; j++)   WinCoeff[N-j-1] = WinCoeff[j];

 // This is the flat top if Alpha > 0. Cannot be applied to a Kaiser or Flat Top.
 if(WindowType != wtKAISER &&  WindowType != wtFLATTOP)
  {
   for(j=M/2; j<N-M/2; j++) WinCoeff[j] = 1.0;
  }


 // UnityGain = true will set the gain of these windows to 1. Don't use this with FIR filter design.
 if(UnityGain)
  {
   double Sum = 0.0;
   for(j=0; j<N; j++)   Sum += WinCoeff[j];
   Sum /= (double)N;
   if(Sum != 0.0)for(j=0; j<N; j++) WinCoeff[j] /= Sum;
  }

 // Apply the window to the data.
 for(j=0; j<N; j++) Data[j] *= WinCoeff[j];

 return 0;
}

//---------------------------------------------------------------------------

// This gets used with the Kaiser window.
double iowa_Bessel(double x)
{
 double Sum=0.0, XtoIpower;
 int i, j, Factorial;
 for(i=1; i<10; i++)
  {
   XtoIpower = pow(x/2.0, (double)i);
   Factorial = 1;
   for(j=1; j<=i; j++)Factorial *= j;
   Sum += pow(XtoIpower / (double)Factorial, 2.0);
  }
 return(1.0 + Sum);
}

//-----------------------------------------------------------------------------

// This gets used with the Sinc window.
static double iowa_Sinc(double x)
{
 if(x > -1.0E-5 && x < 1.0E-5)return(1.0);
 return(sin(x)/x);
}

//---------------------------------------------------------------------------

