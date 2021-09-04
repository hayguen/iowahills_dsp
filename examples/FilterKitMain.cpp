
/*
 By Daniel Klostermann
 Iowa Hills Software, LLC  IowaHills.com
 If you find a problem, please leave a note at:
 http://www.iowahills.com/feedbackcomments.html
 May 1, 2016

 This file is an example of how to use the code kit.

 The code in this kit is essentially the same code used in the Iowa Hills IIR and FIR programs,
 but without the GUI. If you compare the output from this code to the output from the Iowa Hills
 programs, you may find small differences in the coefficients. The reasons for this are numerous,
 and of little concern so long as the filters generated here behave in the desired manner.
*/

#include <stdio.h>
#include <string.h>

#include "FilterKitMain.h"

#include <iowahills/FFTCode.h>             // The FFT code and the windowing code (Kaiser and Sinc windows, etc.)
#include <iowahills/IIRFilterCode.h>       // The IIR filter code.
#include <iowahills/FIRFilterCode.h>       // The FIR filter code.
#include <iowahills/NewParksMcClellan.h>   // The Parks McClellan FIR algorithm code.
#include <iowahills/FreqSamplingCode.h>  //

//---------------------------------------------------------------------------

// ExampleFIRCall shows how to use the functions in the FIRFilterCode and NewParksMcClellan files.
// Select between RectWinFIR and NewParksMcClellan by commenting out the undesired type in 2 places.

// Both of these algorithms are defined by band edges, but these edges don't correspond
// to the filter's -3 dB point. This uses the FIRFreqError function to set the -3 dB at OmegaC
// This also shows how to use FIRFilterWindow to apply a window, such as the Kaiser to the coefficients.

// Frequency values (OmegaC, BW, ParksWidth) are in terms of Nyquist.
// For example, if the sampling freq = 20 kHz then Nyquist = 10 kHz.
// For a low pass filter with a 3 dB corner at 1 kHz set OmegaC = 0.1 (1 kHz / 10 kHz)
void ExampleFIRCall(FILE *OutputFileParam)
{
 int k;
 bool FreqCorrection = true;          // Frequency correction is usually desired, but don't use it with an AllPass filter.
 double RealHofZ[NUM_SAMPLES];       // Real and imag parts of H(z). Used with the FFT.
 double ImagHofZ[NUM_SAMPLES];
 int NumTaps = 40;                    // 4 <= NumTaps < 256 for windowed FIR
                                      // Should be odd for windowed high pass.
                                      // 9 <= NumTaps < 127 for Parks McClellan
                                      // Must be odd for Parks high Pass and notch (checked in the PM code)
 double FIRCoeff[MAX_NUMTAPS];        // FIR filter coefficients.  MAX_NUMTAPS = 256
 double OmegaC = 0.2;                 // 0.0 < OmegaC < 1.0
 double BW = 0.1;                     // 0.0 < BandWidth < 1.0
 TFIRPassTypes PassType = firLPF;     // firLPF, firHPF, firBPF, firNOTCH, firALLPASS  See FIRFilterCode.h
 TWindowType WindowType = wtKAISER;   // wtNONE, wtKAISER, wtSINC, and others.   See the FFT header file.
 double WinBeta = 4.0;                // 0 <= WinBeta <= 10.0  This controls the Kaiser and Sinc windows.
 //double ParksWidth = 0.15;            // 0.01 <= ParksWidth <= 0.3 The transition bandwidth.
                                      // 0.01 <= ParksWidth <= 0.15 if BPF or NOTCH or if NumTaps > 70


 // Use either RectWinFIR or NewParksMcClellan to calculate the FIR Coefficients.
 RectWinFIR(FIRCoeff, NumTaps, PassType, OmegaC, BW);
 FIRFilterWindow(FIRCoeff, NumTaps, WindowType, WinBeta); // Use a window with RectWinFIR.
 //NewParksMcClellan(FIRCoeff, NumTaps, PassType, OmegaC, BW, ParksWidth); // This doesn't require a window, but one may be used.

 if(FreqCorrection && PassType != firALLPASS)
  {
   double OrigOmega = OmegaC;
   double OrigBW = BW;

   // This function corrects OmegaC for LPF and HPF. It corrects BW for BPF and Notch.
   FIRFreqError(FIRCoeff, NumTaps, PassType, &OmegaC, &BW);

   // Recalculate the filter with the corrected OmegaC and BW values.
   RectWinFIR(FIRCoeff, NumTaps, PassType, OmegaC, BW);
   FIRFilterWindow(FIRCoeff, NumTaps, WindowType, WinBeta);  // Use a window with RectWinFIR.
   //NewParksMcClellan(FIRCoeff, NumTaps, PassType, OmegaC, BW, ParksWidth); // This doesn't require a window, but one may be used.

   OmegaC = OrigOmega; // Restore these in case they are needed.
   BW = OrigBW;
  }

 // We can use AdjustDelay to make fractional delay adjustments to the filter. Typically, the
 // delay can be adjusted by +/- NumTaps/20 without affecting the filter's performance significantly.
 // This is done to align signals and is done as the last step (after the application of the window).
 // AdjustDelay(FIRCoeff, NumTaps, 0.5);

 // Calculate the frequency response of the filter with an FFT.
 for(k=0; k<NUM_SAMPLES; k++)RealHofZ[k] = ImagHofZ[k] = 0.0;             // Init the arrays
 for(k=0; k<NumTaps; k++)RealHofZ[k] = FIRCoeff[k] * (double)NUM_SAMPLES; // Need to do this scaling to account for the 1/N scaling done by the forward FFT.
 FFT(RealHofZ, ImagHofZ, NUM_SAMPLES, FORWARD);                           // The FFT's results are returned in input arrays, RealHofZ and ImagHofZ.


 //Print the FIR coefficients to a text file.
 FILE * OutputFile = OutputFileParam;
 if (!OutputFileParam)
   OutputFile = fopen("FIR Filter Coeff.txt","w");
 fprintf(OutputFile,"\noutput of ExampleFIRCall(): FIR Filter %d Coeff:\n", NumTaps);
 for(int j=0; j<NumTaps; j++)
  {
   fprintf(OutputFile,"\n %9.9f ", FIRCoeff[j]);
  }
 fprintf(OutputFile,"\n");
 if (!OutputFileParam)
   fclose(OutputFile);
}

//---------------------------------------------------------------------------

// This shows how to use the functions in the IIRFilterCode file.
// Frequency values (OmegaC, BW, ParksWidth) are in terms of Nyquist.
void ExampleIIRCall(FILE *OutputFileParam)
{
 int j, N;
 TIIRFilterParams IIRFilt;  // Defined in IIRFilterCode.h
 TIIRCoeff IIRCoeff;

 // This structure must be filled before calling CalcIIRFilterCoeff().
 IIRFilt.IIRPassType = iirLPF;        // iirLPF, iirHPF, iirBPF, iirNOTCH, iirALLPASS  (defined in IIRFilterCode.h)
 IIRFilt.OmegaC = 0.2;                // 0.0 < OmegaC < 1.0        3 dB freq for low and high pass filters, center freq for band pass and notch filters.
 IIRFilt.BW = 0.1;                    // 0.0 < BandWidth < 1.0     3 dB bandwidth for bandpass and notch filters
 IIRFilt.dBGain = 0.0;                // -60.0 < dBGain < 60.0     All filters

 // Define the low pass prototype. These are range checked in LowPassPrototypes.cpp
 IIRFilt.ProtoType = INVERSE_CHEBY;   // BUTTERWORTH, CHEBYSHEV, GAUSSIAN, BESSEL, ADJUSTABLE, INVERSE_CHEBY, PAPOULIS, ELLIPTIC  (defined in LowPassPrototypes.h)
 IIRFilt.NumPoles = 6;                // 1 <= NumPoles <= 12, 15, 20 Depending on the filter.
 IIRFilt.Ripple = 0.1;                // 0.0 <= Ripple <= 1.0 dB     Chebyshev and Elliptic (less for high order Chebyshev).
 IIRFilt.StopBanddB = 60.0;           // 20 <= StopBand <= 120 dB    Inv Cheby and Elliptic
 IIRFilt.Gamma = 0.0;                 // -1.0 <= Gamma <= 1.0        Adjustable Gauss  Controls the transition BW.


 // This will fill the IIRCoeff struct with the 2nd order IIR coefficients.
 IIRCoeff = CalcIIRFilterCoeff(IIRFilt);

 // If desired, this will create an Nth order poly from the 2nd order polys in IIRCoeff.
 double DenomCoeff[25], NumerCoeff[25];
 N = RebuildPoly(IIRCoeff.NumSections, DenomCoeff, IIRCoeff.a2, IIRCoeff.a1, IIRCoeff.a0 );
 N = RebuildPoly(IIRCoeff.NumSections, NumerCoeff, IIRCoeff.b2, IIRCoeff.b1, IIRCoeff.b0 );

  // This calculates the frequency response of the filter by doing a DFT of the IIR coefficients.
 double RealHofZ[NUM_SAMPLES];   // Real and imag parts of H(z). Used with the function IIRFreqResponse.
 double ImagHofZ[NUM_SAMPLES];
 IIRFreqResponse(IIRCoeff, IIRCoeff.NumSections, RealHofZ, ImagHofZ, NUM_SAMPLES);

 // This is an alternative way to calculate the filter's frequency response using the FFT.
 // We send an impulse through the filter, and calc the FFT of the filters output.
 // Since the FFT scales the output of a forward transform by 1/N, we use N = NUM_SAMPLES instead of 1 for the impulse.
 double Samples[NUM_SAMPLES];
 for(j=0; j<NUM_SAMPLES; j++)Samples[j] = RealHofZ[j] = ImagHofZ[j] = 0.0;
 Samples[0] = NUM_SAMPLES;                                 // The impulse.
 FilterWithIIR(IIRCoeff, Samples, RealHofZ, NUM_SAMPLES);  // Filter the impulse. RealHofZ is used to store the filtered output.
 FFT(RealHofZ, ImagHofZ, NUM_SAMPLES, FORWARD);            // The FFT's results are returned in input arrays, RealHofZ and ImagHofZ.


 // Print the IIR coefficients to a text file in 3 formats.
 FILE * OutputFile = OutputFileParam;
 if (!OutputFileParam)
   OutputFile = fopen("IIR Filter Coeff.txt","w");
 fprintf(OutputFile,"\noutput of ExampleIIRCall(): IIR Filter Coeff %d sections:\n", IIRCoeff.NumSections);
 for(j=0; j<IIRCoeff.NumSections; j++)
  {
   fprintf(OutputFile,"\n Section %d", j);
   fprintf(OutputFile,"\n a0= %9.9f  a1= %9.9f  a2= %9.9f", IIRCoeff.a0[j], IIRCoeff.a1[j], IIRCoeff.a2[j]);
   fprintf(OutputFile,"\n b0= %9.9f  b1= %9.9f  b2= %9.9f ",IIRCoeff.b0[j], IIRCoeff.b1[j], IIRCoeff.b2[j]);
   fprintf(OutputFile,"\n");
  }
 for(j=0; j<IIRCoeff.NumSections; j++)
  {
   fprintf(OutputFile,"\n  %9.9f \n  %9.9f \n  %9.9f", IIRCoeff.a0[j], IIRCoeff.a1[j], IIRCoeff.a2[j]);
   fprintf(OutputFile,"\n  %9.9f \n  %9.9f \n  %9.9f ",IIRCoeff.b0[j], IIRCoeff.b1[j], IIRCoeff.b2[j]);
   fprintf(OutputFile,"\n");
  }

 fprintf(OutputFile,"\n Nth Order Coeff (N=%d).\n b's are numerator, a's are denominator, each %d coeffs:", N, N+1);
 for(j=N; j>=0; j--)
  {
   fprintf(OutputFile,"\n b%d %9.9f ", j, NumerCoeff[j]);
  }
 fprintf(OutputFile,"\n");
 for(j=N; j>=0; j--)
  {
   fprintf(OutputFile,"\n a%d %9.9f ", j, DenomCoeff[j]);
  }

 fprintf(OutputFile,"\n");
 if (!OutputFileParam)
   fclose(OutputFile);
}

//---------------------------------------------------------------------------

// This example shows how to generate a table of 2nd order filter coefficients
// for low pass prototype filters.
void ExampleLowPassProtoType(FILE *OutputFileParam)
{
 int j, NumPoles;
 double Ripple, StopBanddB;  // Gamma looks unused
 TLowPassParams LPFProto;   // defined in LowPassPrototypes.h
 TSPlaneCoeff SPlaneCoeff;  // defined in LowPassPrototypes.h

 // These are the 5 variables that define all the available low pass prototypes.
 // These are range checked in the function CalcLowPassProtoCoeff().
 LPFProto.ProtoType = BUTTERWORTH;   // BUTTERWORTH, CHEBYSHEV, GAUSSIAN, BESSEL, ADJUSTABLE, INVERSE_CHEBY, PAPOULIS, ELLIPTIC  (defined in LowPassPrototypes.h)
 LPFProto.NumPoles = 6;              // 1 <= NumPoles <= 12, 15, 20 Depending on the filter.
 LPFProto.Ripple = 0.1;              // 0.0 <= Ripple <= 1.0 dB     Chebyshev and Elliptic (less for high order Chebyshev).
 LPFProto.StopBanddB = 60.0;         // 20 <= StopBand <= 120 dB    Inv Cheby and Elliptic
 LPFProto.Gamma = 0.0;               // -1.0 <= Gamma <= 1.0        Adjustable Gauss  Controls the transition BW.


 // This example generates a table of 2nd order S plane coefficients for an Elliptc low pass prototype.
 FILE * OutputFile = OutputFileParam;
 if (!OutputFileParam)
   OutputFile = fopen("Elliptic S Plane Coeff.txt","w");
 fprintf(OutputFile,"\noutput of ExampleLowPassProtoType(): Elliptic S Plane Coeff:\n");
 fprintf(OutputFile,"The format is: H(s) = (N2*s^2 + N1*s + N0) / (D2*s^2 + D1*s + D0) \n");
 fprintf(OutputFile,"The 1st N2=0 for odd pole counts (no zero) else N2=1\n");
 fprintf(OutputFile,"The 1st D2=0 for odd pole counts (a real pole) else D2=1\n");
 fprintf(OutputFile,"Numer N2 N1 N0 \n");
 fprintf(OutputFile,"Denom D2 D1 D0 \n");
 LPFProto.ProtoType = ELLIPTIC;
 for(NumPoles=2; NumPoles<=10; NumPoles++)
  {
   for(Ripple=0.0; Ripple<0.501; Ripple+=0.25)
    {
     for(StopBanddB=40.0; StopBanddB<100.01; StopBanddB+=20.0)
      {
       LPFProto.NumPoles = NumPoles;
       LPFProto.Ripple = Ripple;
       LPFProto.StopBanddB = StopBanddB;
       SPlaneCoeff = CalcLowPassProtoCoeff(LPFProto);
       fprintf(OutputFile,"\n\n Poles=%d  Ripple=%1.2f  Atten=%1.2fdB  %d Sections", NumPoles, Ripple, StopBanddB, SPlaneCoeff.NumSections);
       for(j=0; j<SPlaneCoeff.NumSections; j++)
        {
         fprintf(OutputFile,"\n Numer %d  %1.2f  %1.15f  %1.15f ",j+1, SPlaneCoeff.N2[j], SPlaneCoeff.N1[j], SPlaneCoeff.N0[j]);
         fprintf(OutputFile,"\n Denom %d  %1.2f  %1.15f  %1.15f ",j+1, SPlaneCoeff.D2[j], SPlaneCoeff.D1[j], SPlaneCoeff.D0[j]);
        }
      }
    }
  }
 fprintf(OutputFile,"\n");
 if (!OutputFileParam)
   fclose(OutputFile);

 // This generates a table of Butterworth coefficients, which will match what you find in any text.
 // The other filter coefficients will not match most filter tables because we frequency scale the
 // polynomials so that their 3 dB corner is at 1 rad/sec, which some tabulated filters aren't close to.
 FILE *OutPutFile2 = OutputFileParam;
 if (!OutputFileParam)
   OutPutFile2 = fopen("Butterworth S Plane Coeff.txt","w");
 fprintf(OutputFile,"\noutput of ExampleLowPassProtoType(): Butterworth S Plane Coeff:\n");
 fprintf(OutPutFile2,"The format is: H(s) = (N2*s^2 + N1*s + N0) / (D2*s^2 + D1*s + D0) \n");
 fprintf(OutPutFile2,"The 1st N2=0 for odd pole counts (no zero) else N2=1\n");
 fprintf(OutPutFile2,"The 1st D2=0 for odd pole counts (a real pole) else D2=1\n");
 fprintf(OutPutFile2,"Numer N2 N1 N0\n");
 fprintf(OutPutFile2,"Denom D2 D1 D0\n");
 LPFProto.ProtoType = BUTTERWORTH;
 for(NumPoles=1; NumPoles<=20; NumPoles++)
  {
   LPFProto.NumPoles = NumPoles;
   SPlaneCoeff = CalcLowPassProtoCoeff(LPFProto);
   fprintf(OutPutFile2,"\n\n Poles=%d    %d Sections", NumPoles, SPlaneCoeff.NumSections);
   for(j=0; j<SPlaneCoeff.NumSections; j++)
    {
     fprintf(OutPutFile2,"\n Denom %d  %1.2f  %1.15f  %1.15f ",j+1, SPlaneCoeff.D2[j], SPlaneCoeff.D1[j], SPlaneCoeff.D0[j]);
    }
  }
 fprintf(OutPutFile2,"\n");
 if (!OutputFileParam)
   fclose(OutPutFile2);
}

//---------------------------------------------------------------------------
/*
 This shows how to use frequency sampling to generate an FIR filter. We sample the
 positive frequency domain of the desired magnitude response NUM_POS_FREQ_SAMPLES
 times, regardless the number of taps needed. NUM_POS_FREQ_SAMPLES is defined in
 FIRFreqSampligCode.h and needs to be a power of 2 for the FFT.  1024 works fine for most filters.

 This example generates a conventional low pass, high pass, band pass, or notch filter.
*/
void ExampleFreqSampling(FILE *OutputFileParam)
{
 int j, NumTaps, RightJ, LeftJ;
 double HofSReal[2*NUM_POS_FREQ_SAMPLES], HofSImag[2*NUM_POS_FREQ_SAMPLES];
 double dNumSamples, OmegaC, BW, WinBeta, FIRCoeff[MAX_NUMTAPS];

 dNumSamples = (double)NUM_POS_FREQ_SAMPLES;
 for(j=0; j<MAX_NUMTAPS; j++)FIRCoeff[j] = 0.0;   // Init these arrays
 for(j=0; j<2*NUM_POS_FREQ_SAMPLES; j++)HofSReal[j] = HofSImag[j] = 0.0;

 TFIRPassTypes PassType = firLPF;  // firLPF firBPF firHPF firNOTCH  (Defined in FIRFilterCode.h)
 TWindowType WindowType= wtKAISER; // wtKAISER wtSINC wtSINE (Defined in FFTCode.h)

 NumTaps = 101;  // 3 < NumTaps < MAX_NUMTAPS which is defined in FIRFilterCode.h
 WinBeta = 6.0;  // 0.0 < WinBeta < 10.0  Used with the Kaiser and Sinc Windows.
 OmegaC = 0.2;   // The band edge for low and high pass filters. The center freq for notch and bandpass.
 BW = 0.1;       // The bandwidth for notch an bandpass filters.


 // Use OmegaC and BW to define the filter.
 switch(PassType)
  {
   case firLPF:
	RightJ = (int)(OmegaC * dNumSamples);
	for(j=0; j<RightJ; j++)HofSReal[j] = 1.0;
   break;

   case firHPF:
	LeftJ = (int)(OmegaC * dNumSamples);
	for(j=LeftJ; j<NUM_POS_FREQ_SAMPLES; j++)HofSReal[j] = 1.0;
   break;

   case firBPF:
	RightJ = (int)((OmegaC + BW/2.0)  * dNumSamples);
	LeftJ  = (int)((OmegaC - BW/2.0)  * dNumSamples);
    for(j=LeftJ; j<RightJ; j++)HofSReal[j] = 1.0;
    break;

   case firNOTCH:
    for(j=0; j<NUM_POS_FREQ_SAMPLES; j++)HofSReal[j] = 1.0;
	RightJ = (int)((OmegaC + BW/2.0)  * dNumSamples);
	LeftJ  = (int)((OmegaC - BW/2.0)  * dNumSamples);
    for(j=LeftJ; j<RightJ; j++)HofSReal[j] = 0.0;
    break;

   default:
   case firALLPASS:
   case firNOT_FIR:
    return;
  }

 // Generate the FIR taps.
 SampledFreqFIR(NumTaps, FIRCoeff, HofSReal, HofSImag, OmegaC, PassType);

 // Apply a window to the coefficients to reduce the sidelobe levels.
 FIRFilterWindow(FIRCoeff, NumTaps, WindowType, WinBeta);

  //Print the FIR coefficients to a text file.
 FILE * OutputFile = OutputFileParam;
 if (!OutputFileParam)
   OutputFile = fopen("Sampled Freq FIR.txt","w");
 fprintf(OutputFile,"\noutput of ExampleFreqSampling(): Sampled Freq FIR: %d coeffs:\n", NumTaps);
 for(int j=0; j<NumTaps; j++)
  {
   fprintf(OutputFile,"\n %9.9f ", FIRCoeff[j]);
  }
 fprintf(OutputFile,"\n");
 if (!OutputFileParam)
   fclose(OutputFile);
}
//---------------------------------------------------------------------------

// This example generates a custom FIR filter with 2 band pass responses and a low pass response.
void ExampleFreqSampling2(FILE *OutputFileParam)
{
 int j, NumTaps, RightJ, LeftJ;
 double HofSReal[2*NUM_POS_FREQ_SAMPLES], HofSImag[2*NUM_POS_FREQ_SAMPLES];
 double dNumSamples, OmegaC, BW, WinBeta, FIRCoeff[MAX_NUMTAPS];

 dNumSamples = (double)NUM_POS_FREQ_SAMPLES;
 for(j=0; j<MAX_NUMTAPS; j++)FIRCoeff[j] = 0.0;   // Init these arrays
 for(j=0; j<2*NUM_POS_FREQ_SAMPLES; j++)HofSReal[j] = HofSImag[j] = 0.0;

 TFIRPassTypes PassType = firLPF;  // firLPF firBPF firHPF firNOTCH
 TWindowType WindowType= wtKAISER; // wtKAISER wtSINC wtSINE

 NumTaps = 150;
 WinBeta = 4.0;

 OmegaC = 0.1;
 RightJ = (int)(OmegaC * dNumSamples);
 for(j=0; j<RightJ; j++)HofSReal[j] = 1.0;   // Low pass Gain = 0 dB

 OmegaC = 0.3;
 BW = 0.1;
 RightJ = (int)((OmegaC + BW/2.0)  * dNumSamples);
 LeftJ  = (int)((OmegaC - BW/2.0)  * dNumSamples);
 for(j=LeftJ; j<RightJ; j++)HofSReal[j] = 2.0;  // Bandpass at Omega = .3  Gain = +6 dB

 OmegaC = 0.5;
 BW = 0.1;
 RightJ = (int)((OmegaC + BW/2.0)  * dNumSamples);
 LeftJ  = (int)((OmegaC - BW/2.0)  * dNumSamples);
 for(j=LeftJ; j<RightJ; j++)HofSReal[j] = 0.5;  // Bandpass at Omega = .5  Gain = -6 dB


 // Generate the FIR taps.
 SampledFreqFIR(NumTaps, FIRCoeff, HofSReal, HofSImag, OmegaC, PassType);

 // Apply a window to the coefficients to reduce the sidelobe levels.
 FIRFilterWindow(FIRCoeff, NumTaps, WindowType, WinBeta);

  //Print the FIR coefficients to a text file.
 FILE * OutputFile = OutputFileParam;
 if (!OutputFileParam)
   OutputFile = fopen("Sampled Freq FIR 2.txt","w");
 fprintf(OutputFile,"\noutput of ExampleFreqSampling2(): Sampled Freq FIR 2: %d coeffs\n", NumTaps);
 for(int j=0; j<NumTaps; j++)
  {
   fprintf(OutputFile,"\n %9.9f ", FIRCoeff[j]);
  }
 fprintf(OutputFile,"\n");
 if (!OutputFileParam)
   fclose(OutputFile);
}

//---------------------------------------------------------------------------

// This example shows how to generate an FIR filter from an analog response.
// This function doesn't do much, all the work is done in SampledFreqAnalog.
void AnalogFreqSampling(FILE *OutputFileParam)
{
 int j, NumTaps;
 double HofSReal[2*NUM_POS_FREQ_SAMPLES], HofSImag[2*NUM_POS_FREQ_SAMPLES];
 double OmegaC, FIRCoeff[MAX_NUMTAPS];

 for(j=0; j<MAX_NUMTAPS; j++)FIRCoeff[j] = 0.0;   // Init these arrays
 for(j=0; j<2*NUM_POS_FREQ_SAMPLES; j++)HofSReal[j] = HofSImag[j] = 0.0;

 NumTaps = 150;
 OmegaC = 0.1;

 // Generate the FIR taps.
 SampledFreqAnalog(NumTaps, FIRCoeff, HofSReal, HofSImag, OmegaC);

 //Print the FIR coefficients to a text file.
 FILE * OutputFile = OutputFileParam;
 if (!OutputFileParam)
   OutputFile = fopen("Sampled Freq FIR Analog.txt","w");
 fprintf(OutputFile,"\noutput of AnalogFreqSampling(): Sampled Freq FIR Analog: %d coeffs\n", NumTaps);
 for(int j=0; j<NumTaps; j++)
  {
   fprintf(OutputFile,"\n %9.9f ", FIRCoeff[j]);
  }
 fprintf(OutputFile,"\n");
 if (!OutputFileParam)
   fclose(OutputFile);
}

//---------------------------------------------------------------------------

int main(int argc, char* argv[])
{
 FILE *OutputFile = stdout;
 if (1 < argc && !strcmp(argv[1], "-"))
   OutputFile = stdout;
 if (1 < argc && !strcmp(argv[1], "f"))
   OutputFile = NULL;
 if (1 >= argc)
   fprintf(stderr,
   "usage: %s [-|f]\n"
   "  option - writes to stdout\n"
   "  option f writes to files\n"
   "  default is stdout.\n\n", argv[0]);

 ExampleFIRCall(OutputFile);          // Geneartes either a Rectangular Windowed FIR or a Parks McClellan FIR.
 ExampleIIRCall(OutputFile);          // IIR Filter coefficients
 ExampleLowPassProtoType(OutputFile); // 2nd order s plane coefficients for a Butterworth, Cheby, Bessel, etc.
 ExampleFreqSampling(OutputFile);     // An FIR filter using frequency sampling.
 ExampleFreqSampling2(OutputFile);    // A 2nd example of FIR using freq sampling.
 AnalogFreqSampling(OutputFile);      // An analog type FIR filter using freq sampling.
 return 0;
}

