/*
 This software is part of iowahills_dsp, a set of DSP routines under MIT License.
 2016 By Daniel Klostermann, Iowa Hills Software, LLC  IowaHills.com
 Copyright (c) 2021  Hayati Ayguen <h_ayguen@web.de>
 All rights reserved.
*/

#ifndef ParksMcClellanH
#define ParksMcClellanH

// If the fir.cpp file is in the project along with this ParksMcClellan file,
// we need to include fir.h for the TFIRPassTypes enum.
#include "fir.h"
#include "limits.h"

//---------------------------------------------------------------------------

// previously Global variables.
struct ParksGlobals {
    int HalfTapCount;
    int ExchangeIndex[PARKS_SMALL];
    double LeGrangeD[PARKS_SMALL], Alpha[PARKS_SMALL], CosOfGrid[PARKS_SMALL], DesPlus[PARKS_SMALL];
    double Coeff[PARKS_SMALL], Edge[PARKS_SMALL], BandMag[PARKS_SMALL], InitWeight[PARKS_SMALL];

    double DesiredMag[PARKS_BIG], Grid[PARKS_BIG], Weight[PARKS_BIG];
};

//---------------------------------------------------------------------------

// BW is only for PassType in [firBPF, firNOTCH]
// for PassType == firLPF: OmegaC is Pass band edge, (OmegaC + ParksWidth) is Stop band edge
int ParksMcClellan(double *FirCoeff, int NumTaps, TFIRPassTypes PassType, double OmegaC, double BW, double ParksWidth);
int ParksMcClellanEx(ParksGlobals *PtrGlobals, double *FirCoeff, int NumTaps, TFIRPassTypes PassType, double OmegaC, double BW, double ParksWidth);

#endif

