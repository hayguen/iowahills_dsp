/*
 This software is part of iowahills_dsp, a set of DSP routines under MIT License.
 2016 By Daniel Klostermann, Iowa Hills Software, LLC  IowaHills.com
 Copyright (c) 2021  Hayati Ayguen <h_ayguen@web.de>
 All rights reserved.
*/

#ifndef PFiftyOneRevEH
#define PFiftyOneRevEH
//---------------------------------------------------------------------------

#include "CplxDMath.hpp"

// TUpdateStatus represent the status of the UpdateTUV function.
// The first 3 are returned from UpdateTUV. The last two are sent to UpdateTUV.
enum TUpdateStatus {UPDATED, BAD_ANGLE, ZERO_DEL, DAMPER_ON, DAMPER_OFF};

#define P51_MAXDEGREE   100       // The max poly order allowed. Used at the top of P51. This was set arbitrarily.
#define P51_ARRAY_SIZE  102       // P51 uses the new operator. P51 arrays must be MaxDegree + 2

int FindRoots(int N, double *Coeff, CplxD *Roots);
int PFiftyOne(long double *Coeff, int Degree, long double *RealRoot, long double *ImagRoot);
int QuadIterate(int &FirstDamperlIter, int Iter, long double *P, long double *QP, long double *K, long double *QK, int N, long double *TUV, TUpdateStatus *UpdateStatus);
void UpdateTUV(int &FirstDamperlIter, int Iter, long double *P, int N, long double *QP, long double *K, long double *QK, long double *TUV, TUpdateStatus *UpdateStatus);
int RealIterate(long double &PrevQPN, int P51_Iter, long double *P, long double *QP, long double *K, long double *QK, int N, long double *RealZero);
void QuadraticFormula(long double *TUV, long double *RealRoot, long double *ImagRoot);
void QuadSynDiv(long double *P, int N, long double *TUV, long double *Q);
void DerivOfP(long double *P, int N, long double *dP);
int SetTUVandK(long double *P, int N, long double *TUV, long double *RealK, long double *QuadK, long double X, int AngleNumber, int TypeOfQuadK);

#endif

