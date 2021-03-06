/*
 This software is part of iowahills_dsp, a set of DSP routines under MIT License.
 2016 By Daniel Klostermann, Iowa Hills Software, LLC  IowaHills.com
 Copyright (c) 2021  Hayati Ayguen <h_ayguen@web.de>
 All rights reserved.
*/

#ifndef LowPassPrototypesH
#define LowPassPrototypesH

#include "CplxDMath.hpp"

// This MUST be at least 2*MAX_POLE_COUNT because some filter polys are defined in terms of 2 * NumPoles
#define IOWA_HILLS_ARRAY_DIM 50

// These coeff form H(s) = (N2*s^2 + N1*s + N0) / (D2*s^2 + D1*s + D0)
// NumSections is the number of 1st and 2nd order polynomial factors .
struct TSPlaneCoeff {
    double N2[IOWA_HILLS_ARRAY_DIM]; double N1[IOWA_HILLS_ARRAY_DIM]; double N0[IOWA_HILLS_ARRAY_DIM];
    double D2[IOWA_HILLS_ARRAY_DIM]; double D1[IOWA_HILLS_ARRAY_DIM]; double D0[IOWA_HILLS_ARRAY_DIM];
    int NumSections;
};

// These are the available filter polynomials. NOT_IIR is for code testing.
enum TFilterPoly {
    BUTTERWORTH, GAUSSIAN, BESSEL, ADJUSTABLE, CHEBYSHEV,
    INVERSE_CHEBY, PAPOULIS, ELLIPTIC, NOT_IIR
};

// This structure defines the low pass filter prototype.
// The 3 dB corner frequency is 1 rad/sec for all filters.
struct TLowPassParams {
    TFilterPoly ProtoType;  // Butterworth, Cheby, etc.
    int NumPoles;           // Pole count
    double Ripple;          // Passband Ripple for the Elliptic and Chebyshev
    double StopBanddB;      // Stop Band Attenuation in dB for the Elliptic and Inverse Cheby
    double Gamma;           // Controls the transition bandwidth on the Adjustable Gauss. -1 <= Gamma <= 1
};


TSPlaneCoeff CalcLowPassProtoCoeff(TLowPassParams Filt);
void SetCornerFreq(int Count, double *D2, double *D1, double *D0, double *N2, double *N1, double *N0);
int GetFilterCoeff(int RootCount, CplxD *Roots, double *A2, double *A1, double *A0);
int RebuildPoly(int PolyCount, double *PolyCoeff, double *A2, double *A1, double *A0 );

#endif

