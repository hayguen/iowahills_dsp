/*
 This software is part of iowahills_dsp, a set of DSP routines under MIT License.
 2016 By Daniel Klostermann, Iowa Hills Software, LLC  IowaHills.com
 Copyright (c) 2021  Hayati Ayguen <h_ayguen@web.de>
 All rights reserved.
*/

#include <iowahills/goertzel.h>

#include <cmath>

//---------------------------------------------------------------------------

// Goertzel is essentially a single frequency DFT, but without phase information.
// Its simplicity allows it to run about 3 times faster than a single frequency DFT.
// It is typically used to find a tone embedded in a signal. A DTMF tone for example.
// 256 pts in 6 us
double Goertzel(const double *Samples, int N, double Omega)
{
 if (!N)
     return 0.0;
 int j;
 double Reg0, Reg1, Reg2;        // 3 shift registers
 double CosVal, Mag;
 Reg1 = Reg2 = 0.0;

 CosVal = 2.0 * cos(M_PI * Omega );
 for (j=0; j<N; j++)
  {
   Reg0 = Samples[j] + CosVal * Reg1 - Reg2;
   Reg2 = Reg1;  // Shift the values.
   Reg1 = Reg0;
  }
 Mag = Reg2 * Reg2 + Reg1 * Reg1 - CosVal * Reg1 * Reg2;

 if(Mag > 0.0)Mag = sqrt(Mag);
 else Mag = 1.0E-12;

 return(Mag);
}

//---------------------------------------------------------------------------

