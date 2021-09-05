/*
 This software is part of iowahills_dsp, a set of DSP routines under MIT License.
 2016 By Daniel Klostermann, Iowa Hills Software, LLC  IowaHills.com
 Copyright (c) 2021  Hayati Ayguen <h_ayguen@web.de>
 All rights reserved.
*/

#ifndef QuadRootsRevHH
#define QuadRootsRevHH
//---------------------------------------------------------------------------

void QuadRoots(long double *P, long double *RealPart, long double *ImagPart);
void CubicRoots(long double *P, long double *RealPart, long double *ImagPart);
void BiQuadRoots(long double *P, long double *RealPart, long double *ImagPart);
void ReversePoly(long double *P, int N);
void InvertRoots(int N, long double *RealRoot, long double *ImagRoot);

#endif

