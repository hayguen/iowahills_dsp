/*
 This software is part of iowahills_dsp, a set of DSP routines under MIT License.
 2016 By Daniel Klostermann, Iowa Hills Software, LLC  IowaHills.com
 Copyright (c) 2021  Hayati Ayguen <h_ayguen@web.de>
 All rights reserved.
*/

#ifndef FIRFilterCodeH
#define FIRFilterCodeH

#include "windowing.h"   // For the definition of TWindowType
#include "limits.h"
//---------------------------------------------------------------------------

enum TFIRPassTypes {firLPF, firHPF, firBPF, firNOTCH, firALLPASS, firNOT_FIR};

void FilterWithFIR(double *FirCoeff, int NumTaps, double *Signal, double *FilteredSignal, int NumSigPts);
void FilterWithFIR2(double *FirCoeff, int NumTaps, double *Signal, double *FilteredSignal, int NumSigPts);
int  RectWinFIR(double *FirCoeff, int NumTaps, TFIRPassTypes PassType, double OmegaC, double BW);
double iowa_Sinc(double x);
void FIRFreqError(double *Coeff, int NumTaps, int PassType, double *OmegaC, double *BW);
int FIRFilterWindow(double *FIRCoeff, int N, TWindowType WindowType, double Beta);

#endif

