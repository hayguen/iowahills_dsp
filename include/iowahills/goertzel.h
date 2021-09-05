/*
 This software is part of iowahills_dsp, a set of DSP routines under MIT License.
 2016 By Daniel Klostermann, Iowa Hills Software, LLC  IowaHills.com
 Copyright (c) 2021  Hayati Ayguen <h_ayguen@web.de>
 All rights reserved.
*/

#ifndef GoertzelH
#define GoertzelH

// Omega is frequency, normalized to half samplerate. Thus in [0 .. 1]
double Goertzel(const double *Samples, int N, double Omega);

#endif

