/*
 This software is part of iowahills_dsp, a set of DSP routines under MIT License.
 2016 By Daniel Klostermann, Iowa Hills Software, LLC  IowaHills.com
 Copyright (c) 2021  Hayati Ayguen <h_ayguen@web.de>
 All rights reserved.
*/

#ifndef WindowingH
#define WindowingH

//---------------------------------------------------------------------------
// Must retain the order on the 1st line for legacy FIR code.
enum TWindowType {
    wtFIRSTWINDOW, wtNONE, wtKAISER, wtSINC, wtHANNING,
    wtHAMMING, wtBLACKMAN, wtFLATTOP, wtBLACKMAN_HARRIS,
    wtBLACKMAN_NUTTALL, wtNUTTALL, wtKAISER_BESSEL, wtTRAPEZOID,
    wtGAUSS, wtSINE, wtTEST
};

int WindowData(double *Data, int N, TWindowType WindowType, double Alpha, double Beta, bool UnityGain);

double iowa_Bessel(double x);

#endif

