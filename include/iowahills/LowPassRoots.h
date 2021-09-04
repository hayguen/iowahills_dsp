//---------------------------------------------------------------------------

#ifndef LowPassRootsH
#define LowPassRootsH

#include "CplxDMath.hpp"
#define MAX_ELLIP_ITER 15
#define ELLIPARRAYSIZE 20  // needs to be > 10 and >= Max Num Poles + 1


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


