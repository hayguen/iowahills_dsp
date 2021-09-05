/*
 This software is part of iowahills_dsp, a set of DSP routines under MIT License.
 2016 By Daniel Klostermann, Iowa Hills Software, LLC  IowaHills.com
 Copyright (c) 2021  Hayati Ayguen <h_ayguen@web.de>
 All rights reserved.
*/

#ifndef IIRFilterCodeH
#define IIRFilterCodeH

#include "lowpass_prototypes.h"    // defines TFilterPoly and IOWA_HILLS_ARRAY_DIM


enum TIIRPassTypes {iirLPF, iirHPF, iirBPF, iirNOTCH, iirALLPASS};

struct TIIRCoeff {
    double a0[IOWA_HILLS_ARRAY_DIM]; double a1[IOWA_HILLS_ARRAY_DIM]; double a2[IOWA_HILLS_ARRAY_DIM]; double a3[IOWA_HILLS_ARRAY_DIM]; double a4[IOWA_HILLS_ARRAY_DIM];
    double b0[IOWA_HILLS_ARRAY_DIM]; double b1[IOWA_HILLS_ARRAY_DIM]; double b2[IOWA_HILLS_ARRAY_DIM]; double b3[IOWA_HILLS_ARRAY_DIM]; double b4[IOWA_HILLS_ARRAY_DIM];
    int NumSections;
};

struct TIIRFilterParams {
    TIIRPassTypes IIRPassType;     // Defined above: Low pass, High Pass, etc.
    double OmegaC;                 // The IIR filter's 3 dB corner freq for low pass and high pass, the center freq for band pass and notch.
    double BW;                     // The IIR filter's 3 dB bandwidth for band pass and notch filters.
    double dBGain;                 // Sets the Gain of the filter

    // These define the low pass prototype to be used
    TFilterPoly ProtoType;  // Butterworth, Cheby, etc.
    int NumPoles;           // Pole count
    double Ripple;          // Passband Ripple for the Elliptic and Chebyshev
    double StopBanddB;      // Stop Band Attenuation in dB for the Elliptic and Inverse Chebyshev
    double Gamma;           // Controls the transition bandwidth on the Adjustable Gauss. -1 <= Gamma <= 1
};

TIIRCoeff CalcIIRFilterCoeff(TIIRFilterParams IIRFilt);
void FilterWithIIR(TIIRCoeff IIRCoeff, double *Signal, double *FilteredSignal, int NumSigPts);
double SectCalc(int j, int k, double x, TIIRCoeff IIRCoeff);
void IIRFreqResponse(TIIRCoeff IIRCoeff, int NumSections, double *RealHofZ, double *ImagHofZ, int NumPts);

#endif

