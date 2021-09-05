/*
 This software is part of iowahills_dsp, a set of DSP routines under MIT License.
 2016 By Daniel Klostermann, Iowa Hills Software, LLC  IowaHills.com
 Copyright (c) 2021  Hayati Ayguen <h_ayguen@web.de>
 All rights reserved.


 RectWinFIR() generates the impulse response for a rectangular windowed low pass, high pass,
 band pass, or notch filter. Then a window, such as the Kaiser, is applied to the FIR coefficients.
 See the FilterKitMain.cpp file for an example on how to use this code.

 double FirCoeff[MAXNUMTAPS];
 int NumTaps;                        NumTaps can be even or odd and < MAXNUMTAPS
 TPassTypeName PassType;             PassType is defined in the header file. firLPF, firHPF, firBPF, firNOTCH, firALLPASS
 double OmegaC  0.0 < OmegaC < 1.0   The corner freq, or center freq if BPF or NOTCH
 double BW      0.0 < BW < 1.0       The band width if BPF or NOTCH
*/

#include <iowahills/fir.h>
#include <iowahills/goertzel.h>

#include "non_throwing_vector.h"

#include <cmath>

//---------------------------------------------------------------------------

#define NUM_FREQ_ERR_PTS  1000    // these are only used in the FIRFreqError function.
#define dNUM_FREQ_ERR_PTS 1000.0

//---------------------------------------------------------------------------

// Rectangular Windowed FIR. The equations used here are developed in numerous textbooks.
int RectWinFIR(double *FirCoeff, int NumTaps, TFIRPassTypes PassType, double OmegaC, double BW)
{
 int j;
 double Arg, OmegaLow, OmegaHigh;

 switch(PassType)
  {
   case firLPF:    // Low Pass
	for(j=0; j<NumTaps; j++)
     {
      Arg = (double)j - (double)(NumTaps-1) / 2.0;
      FirCoeff[j] = OmegaC * iowa_Sinc(OmegaC * Arg * M_PI);
     }
    break;

   case firHPF:     // High Pass
    if(NumTaps % 2 == 1) // Odd tap counts
     {
      for(j=0; j<NumTaps; j++)
       {
        Arg = (double)j - (double)(NumTaps-1) / 2.0;
        FirCoeff[j] = iowa_Sinc(Arg * M_PI) - OmegaC * iowa_Sinc(OmegaC * Arg * M_PI);
       }
     }

    else  // Even tap counts
      {
       for(j=0; j<NumTaps; j++)
        {
         Arg = (double)j - (double)(NumTaps-1) / 2.0;
         if(Arg == 0.0)FirCoeff[j] = 0.0;
         else FirCoeff[j] = cos(OmegaC * Arg * M_PI) / M_PI / Arg  + cos(Arg * M_PI);
        }
      }
   break;

   case firBPF:   // Band Pass
    OmegaLow  = OmegaC - BW/2.0;
    OmegaHigh = OmegaC + BW/2.0;
	for(j=0; j<NumTaps; j++)
     {
      Arg = (double)j - (double)(NumTaps-1) / 2.0;
      if(Arg == 0.0)FirCoeff[j] = 0.0;
      else FirCoeff[j] =  ( cos(OmegaLow * Arg * M_PI) - cos(OmegaHigh * Arg * M_PI) ) / M_PI / Arg ;
     }
   break;

   case firNOTCH:  // Notch,  if NumTaps is even, the response at Pi is attenuated.
    OmegaLow  = OmegaC - BW/2.0;
    OmegaHigh = OmegaC + BW/2.0;
	for(j=0; j<NumTaps; j++)
     {
      Arg = (double)j - (double)(NumTaps-1) / 2.0;
      FirCoeff[j] =  iowa_Sinc(Arg * M_PI) - OmegaHigh * iowa_Sinc(OmegaHigh * Arg * M_PI) - OmegaLow * iowa_Sinc(OmegaLow * Arg * M_PI);
     }
   break;

   case firALLPASS: // All Pass, this is trivial, but it shows how an fir all pass (delay) can be done.
	for(j=0; j<NumTaps; j++)FirCoeff[j] = 0.0;
    FirCoeff[(NumTaps-1) / 2] = 1.0;
    break;

   default:
   case firNOT_FIR:
    return 1;  // signal error
  }
 // Now use the FIRFilterWindow() function to reduce the sinc(x) effects.

 return 0;
}

//---------------------------------------------------------------------------

// This gets used in numerous places above.
// If using Visual C, this (and the header file declaration) may need to be commented out.
double iowa_Sinc(double x)
{
 if(x > -1.0E-5 && x < 1.0E-5)return(1.0);
 return(sin(x)/x);
}

//---------------------------------------------------------------------------

// Used to reduce the sinc(x) effects on a set of FIR coefficients. This will, unfortunately,
// widen the filter's transition band, but the stop band attenuation will improve dramatically.
int FIRFilterWindow(double *FIRCoeff, int N, TWindowType WindowType, double Beta)
{
 if(WindowType == wtNONE) return 1;

 int j;
 double dN;

 if(Beta < 0.0)Beta = 0.0;
 if(Beta > 10.0)Beta = 10.0;

 NonThrowingVector<double> WinCoeff(N+2);
 if(!WinCoeff)
  {
   // ShowMessage("Failed to allocate memory in WindowData() ");
   return 1;
  }

 // Calculate the window for N/2 points, then fold the window over (at the bottom).
 dN = N + 1; // a double
 if(WindowType == wtKAISER)
  {
   double Arg;
   for(j=0; j<N; j++)
	{
	 Arg = Beta * sqrt(1.0 - pow( ((double)(2*j+2) - dN) / dN, 2.0) );
     WinCoeff[j] = iowa_Bessel(Arg) / iowa_Bessel(Beta);
	}
  }

 else if(WindowType == wtSINC)  // Lanczos
  {
   for(j=0; j<N; j++)   WinCoeff[j] = iowa_Sinc((double)(2*j+1-N)/dN * M_PI );
   for(j=0; j<N; j++)   WinCoeff[j] = pow(WinCoeff[j], Beta);
  }

 else if(WindowType == wtSINE)  // Hanning if Beta = 2
  {
   for(j=0; j<N/2; j++) WinCoeff[j] = sin((double)(j+1) * M_PI / dN);
   for(j=0; j<N/2; j++) WinCoeff[j] = pow(WinCoeff[j], Beta);
  }

 else // Error.
  {
   // ShowMessage("Incorrect window type in WindowFFTData");
   return 1;
  }

 // Fold the coefficients over.
 for(j=0; j<N/2; j++)   WinCoeff[N-j-1] = WinCoeff[j];

 // Apply the window to the FIR coefficients.
 for(j=0; j<N; j++) FIRCoeff[j] *= WinCoeff[j];

 return 0;
}

//---------------------------------------------------------------------------

// This implements an FIR filter. The register shifts are done by rotating the indexes.
void FilterWithFIR(double *FirCoeff, int NumTaps, double *Signal, double *FilteredSignal, int NumSigPts)
{
 int j, k, n, Top = 0;
 double y, Reg[MAX_NUMTAPS];

 for(j=0; j<NumTaps; j++)
   Reg[j] = 0.0;

 for(j=0; j<NumSigPts; j++)
  {
   Reg[Top] = Signal[j];
   y = 0.0;
   n = 0;

   // The FirCoeff index increases while the Reg index decreases.
   for(k=Top; k>=0; k--)
    {
     y += FirCoeff[n++] * Reg[k];
    }
   for(k=NumTaps-1; k>Top; k--)
    {
     y += FirCoeff[n++] * Reg[k];
    }
   FilteredSignal[j] = y;

   Top++;
   if(Top >= NumTaps)Top = 0;
  }

}

//---------------------------------------------------------------------------

// This code is equivalent to the code above. It uses register shifts, which makes it
// less efficient, but it is easier to follow (i.e. compare to a FIR flow chart).
void FilterWithFIR2(double *FirCoeff, int NumTaps, double *Signal, double *FilteredSignal, int NumSigPts)
{
 int j, k;
 double y, Reg[MAX_NUMTAPS];

 for(j=0; j<NumTaps; j++)
   Reg[j] = 0.0; // Init the delay registers.

 for(j=0; j<NumSigPts; j++)
 {
  // Shift the register values down and set Reg[0].
  for(k=NumTaps; k>1; k--)
    Reg[k-1] = Reg[k-2];
  Reg[0] = Signal[j];

  y = 0.0;
  for(k=0; k<NumTaps; k++)
    y += FirCoeff[k] * Reg[k];
  FilteredSignal[j] = y;
 }

}

//---------------------------------------------------------------------------

// This function is used to correct the corner frequency values on FIR filters.
// We normally specify the 3 dB frequencies when specifing a filter. The Parks McClellan routine
// uses OmegaC and BW to set the 0 dB band edges, so its final OmegaC and BW values are not close
// to -3 dB. The Rectangular Windowed filters are better for high tap counts, but for low tap counts,
// their 3 dB frequencies are also well off the mark.

// To use this function, first calculate a set of FIR coefficients, then pass them here, along with
// OmegaC and BW. This calculates a corrected OmegaC for low and high pass filters. It calcultes a
// corrected BW for band pass and notch filters. Use these corrected values to recalculate the FIR filter.

// The Goertzel algorithm is used to calculate the filter's magnitude response at the single
// frequency defined in the loop. We start in the pass band and work out to the -20dB freq.

void FIRFreqError(double *Coeff, int NumTaps, int PassType, double *OmegaC, double *BW)
{
 int j, J3dB, CenterJ;
 double Omega, CorrectedOmega, CorrectedBW, Omega1, Omega2, Mag;

 // In these loops, we break at -20 dB to ensure that large ripple is ignored.
 if(PassType == firLPF)
  {
   J3dB = 10;
   for(j=0; j<NUM_FREQ_ERR_PTS; j++)
	{
     Omega = (double)j / dNUM_FREQ_ERR_PTS;
     Mag = Goertzel(Coeff, NumTaps, Omega);
	 if(Mag > 0.707)J3dB = j;   // J3dB will be the last j where the response was > -3 dB
	 if(Mag < 0.1)break;        // Stop when the response is down to -20 dB.
	}
   Omega = (double)J3dB / dNUM_FREQ_ERR_PTS;
  }

 else if(PassType == firHPF)
  {
   J3dB = NUM_FREQ_ERR_PTS - 10;
   for(j=NUM_FREQ_ERR_PTS-1; j>=0; j--)
	{
     Omega = (double)j / dNUM_FREQ_ERR_PTS;
     Mag = Goertzel(Coeff, NumTaps, Omega);
	 if(Mag > 0.707)J3dB = j;  // J3dB will be the last j where the response was > -3 dB
	 if(Mag < 0.1)break;       // Stop when the response is down to -20 dB.
	}
   Omega = (double)J3dB / dNUM_FREQ_ERR_PTS;
  }

 else if(PassType == firBPF)
  {
   CenterJ = (int)(dNUM_FREQ_ERR_PTS * *OmegaC);
   J3dB = CenterJ;
   for(j=CenterJ; j>=0; j--)
	{
     Omega = (double)j / dNUM_FREQ_ERR_PTS;
     Mag = Goertzel(Coeff, NumTaps, Omega);
	 if(Mag > 0.707)J3dB = j;
	 if(Mag < 0.1)break;
	}
   Omega1 = (double)J3dB / dNUM_FREQ_ERR_PTS;

   J3dB = CenterJ;
   for(j=CenterJ; j<NUM_FREQ_ERR_PTS; j++)
	{
     Omega = (double)j / dNUM_FREQ_ERR_PTS;
     Mag = Goertzel(Coeff, NumTaps, Omega);
	 if(Mag > 0.707)J3dB = j;
	 if(Mag < 0.1)break;
	}
   Omega2 = (double)J3dB / dNUM_FREQ_ERR_PTS;
  }

 // The code above starts in the pass band. This starts in the stop band.
 else // PassType == firNOTCH
  {
   CenterJ = (int)(dNUM_FREQ_ERR_PTS * *OmegaC);
   J3dB = CenterJ;
   for(j=CenterJ; j>=0; j--)
	{
     Omega = (double)j / dNUM_FREQ_ERR_PTS;
     Mag = Goertzel(Coeff, NumTaps, Omega);
	 if(Mag <= 0.707)J3dB = j;
	 if(Mag > 0.99)break;
	}
   Omega1 = (double)J3dB/dNUM_FREQ_ERR_PTS;

   J3dB = CenterJ;
   for(j=CenterJ; j<NUM_FREQ_ERR_PTS; j++)
	{
     Omega = (double)j / dNUM_FREQ_ERR_PTS;
     Mag = Goertzel(Coeff, NumTaps, Omega);
	 if(Mag <= 0.707)J3dB = j;
	 if(Mag > 0.99)break;
	}
   Omega2 = (double)J3dB / dNUM_FREQ_ERR_PTS;
  }


 // This calculates the corrected OmegaC and BW and error checks the values.
 if(PassType == firLPF || PassType == firHPF )
  {
   CorrectedOmega = *OmegaC * 2.0 - Omega;  // This is usually OK.
   if(CorrectedOmega < 0.001)CorrectedOmega = 0.001;
   if(CorrectedOmega > 0.99)CorrectedOmega = 0.99;
   *OmegaC = CorrectedOmega;
  }

 else // PassType == firBPF || PassType == firNOTCH
  {
   CorrectedBW = *BW * 2.0 - (Omega2 - Omega1);  // This routinely goes neg with Notch.
   if(CorrectedBW < 0.01)CorrectedBW = 0.01;
   if(CorrectedBW > *BW * 2.0)CorrectedBW = *BW * 2.0;
   if(CorrectedBW > 0.98)CorrectedBW = 0.98;
   *BW = CorrectedBW;
  }

}

//-----------------------------------------------------------------------------

