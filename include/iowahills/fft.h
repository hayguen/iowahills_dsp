/*
 This software is part of iowahills_dsp, a set of DSP routines under MIT License.
 2016 By Daniel Klostermann, Iowa Hills Software, LLC  IowaHills.com
 Copyright (c) 2021  Hayati Ayguen <h_ayguen@web.de>
 All rights reserved.
*/

#ifndef FFTCodeH
#define FFTCodeH

//---------------------------------------------------------------------------
enum TTransFormType {FORWARD, INVERSE};

int RequiredFFTSize(int NumPts);
int IsValidFFTSize(int x);  // valid is only 8 .. 1048576

// FFT() calculation is always in-place: output values overwrite the input
// FORWARD transform does always normalize the output by 1/N.
// InputR array is for the real components of the complex data
// InputI array is for the imag components of the complex data
int FFT(double *InputR, double *InputI, int N, TTransFormType Type);

struct FFTBuffers;
FFTBuffers * prepareFFT(int N, TTransFormType Type);
void disposeFFT(FFTBuffers *buffers);
// the FFTBuffers can be pre-allocated and reused for multiple calls of FFTex()
//   - except when transformation length N or type Type changes
// skip_fwd_div_N is a bool flag (zero - or not) to do what it's name says
int FFTex(FFTBuffers *buffers, double *InputR, double *InputI, int skip_fwd_div_N);


// DFT() and RealSigDFT() are NOT Fast: they consume O(N^2) !
int DFT(double *InputR, double *InputI, int N, TTransFormType Type);
int RealSigDFT(const double *Samples, double *OutputR, double *OutputI, int N);

// Omega is frequency, normalized to half samplerate. Thus in [0 .. 1]
double SingleFreqDFT(const double *Samples, int N, double Omega);

int AdjustFIRDelay(double *FirCoeff, int NumTaps, double Delay);

#endif

