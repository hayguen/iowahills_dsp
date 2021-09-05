/*
 This software is part of iowahills_dsp, a set of DSP routines under MIT License.
 2016 By Daniel Klostermann, Iowa Hills Software, LLC  IowaHills.com
 Copyright (c) 2021  Hayati Ayguen <h_ayguen@web.de>
 All rights reserved.
*/

#ifndef LowPassRootsH
#define LowPassRootsH

#include "CplxDMath.hpp"

void ReverseCoeff(double *P, int N);
int ButterworthPoly(int NumPoles, CplxD *Roots);
int GaussianPoly(int NumPoles, CplxD *Roots);
int AdjustablePoly(int NumPoles, CplxD *Roots, double Gamma);
int ChebyshevPoly(int NumPoles, double Ripple, CplxD *Roots);
int BesselPoly(int NumPoles, CplxD *Roots);
int InvChebyPoly(int NumPoles, double StopBanddB, CplxD *ChebyPoles, CplxD *ChebyZeros, int *ZeroCount);
int PapoulisPoly(int NumPoles, CplxD *Roots);
int EllipticPoly(int FiltOrder, double Ripple, double DesiredSBdB, CplxD *EllipPoles, CplxD *EllipZeros, int *ZeroCount);

#endif

