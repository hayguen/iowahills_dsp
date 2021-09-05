/*
 This software is part of iowahills_dsp, a set of DSP routines under MIT License.
 2016 By Daniel Klostermann, Iowa Hills Software, LLC  IowaHills.com
 Copyright (c) 2021  Hayati Ayguen <h_ayguen@web.de>
 All rights reserved.
*/

#ifndef FreqSamplingCodeH
#define FreqSamplingCodeH

#include "fir.h" // for the definition of TFIRPassTypes
#include "limits.h"


int SampledFreqFIR(int NumTaps, double *FirCoeff, double *HofSReal, double *HofSImag, double OmegaC, TFIRPassTypes PassType);
void SampledFreqAnalog(int NumTaps, double *FirCoeff, double *HofSReal, double *HofSImag, double OmegaC);

#endif

