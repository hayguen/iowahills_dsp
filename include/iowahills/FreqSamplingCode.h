//---------------------------------------------------------------------------

#ifndef FreqSamplingCodeH
#define FreqSamplingCodeH
#include "CplxDMath.hpp"
#include "FIRFilterCode.h" // for the definition of TFIRPassTypes

 #define NUM_POS_FREQ_SAMPLES 1024  // needs to be 1024 for the kit test app

int SampledFreqFIR(int NumTaps, double *FirCoeff, double *HofSReal, double *HofSImag, double OmegaC, TFIRPassTypes PassType);
void SampledFreqAnalog(int NumTaps, double *FirCoeff, double *HofSReal, double *HofSImag, double OmegaC);
#endif
