/*
 This software is part of iowahills_dsp, a set of DSP routines under MIT License.
 2016 By Daniel Klostermann, Iowa Hills Software, LLC  IowaHills.com
 Copyright (c) 2021  Hayati Ayguen <h_ayguen@web.de>
 All rights reserved.
*/

#ifndef IOWA_HILLS_DSP_LIMITS_H
#define IOWA_HILLS_DSP_LIMITS_H

//---------------------------------------------------------------------------
// fir
#define MAX_NUMTAPS 1024            // limit was 256

//---------------------------------------------------------------------------
// freq_sampling
#define NUM_POS_FREQ_SAMPLES (4 * MAX_NUMTAPS)      // limit was 1024

//---------------------------------------------------------------------------
// parks_mcclellan
#define PARKS_MAX_NUM_TAPS  511     // This was the limit set in the original code: 127
#define PARKS_BIG   ( (8 * PARKS_MAX_NUM_TAPS * 11) / 10 )  // Per the original code, this must be > 8 * MaxNumTaps = 8 * 128 = 1024 -> 1100
#define PARKS_SMALL (  2 * PARKS_MAX_NUM_TAPS )             // This needs to be greater than or equal to PARKS_MAX_NUM_TAPS

#endif

