
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>

// https://www.sciencedirect.com/topics/engineering/mcclellan-algorithm
// https://dsp.stackexchange.com/questions/37646/filter-order-rule-of-thumb
// https://dsp.stackexchange.com/questions/20514/what-is-the-relation-of-the-transition-bands-width-and-the-filter-order-for-the/20515


#define M_2PI     6.28318530717958647692   // 2*Pi


int main(int argc, char* argv[])
{
    double OmegaPass = -1.0;
    double OmegaStop = -1.0;
    double decim = -1.0;
    double ripple = 0.1;    // in dB
    double atten = 80.0;    // attenuation
    int NumTaps = -1;
    bool printUsage = false;

    for (int idx = 1; idx < argc; ++idx)
    {
        if (!strcmp(argv[idx], "-h") || !strcmp(argv[idx], "--help"))
        {
            printUsage = true;
            break;
        }

        else if ( (!strcmp(argv[idx], "-n") || !strcmp(argv[idx], "--numtaps") ) && idx +1 < argc )
        {
            NumTaps = atoi(argv[idx+1]);
            ++idx;
        }

        else if ( (!strcmp(argv[idx], "-s") || !strcmp(argv[idx], "--stop") ) && idx +1 < argc )
        {
            OmegaStop = atof(argv[idx+1]);
            ++idx;
            if ( OmegaStop < 0.0 )
                OmegaStop = 0.0;
            else if (OmegaStop > 1.0)
                OmegaStop = 1.0;
        }

        else if ( (!strcmp(argv[idx], "-p") || !strcmp(argv[idx], "--pass") ) && idx +1 < argc )
        {
            OmegaPass = atof(argv[idx+1]);
            ++idx;
            if ( OmegaPass < 0.0 )
                OmegaPass = 0.0;
            else if (OmegaPass > 1.0)
                OmegaPass = 1.0;
        }

        else if ( (!strcmp(argv[idx], "-a") || !strcmp(argv[idx], "--atten") ) && idx +1 < argc )
        {
            atten = atof(argv[idx+1]);
            ++idx;
            if ( atten < 0.0 )
                atten = 0.0;
            else if (atten > 200.0)
                atten = 200.0;
        }

        else if ( (!strcmp(argv[idx], "-r") || !strcmp(argv[idx], "--ripple") ) && idx +1 < argc )
        {
            ripple = atof(argv[idx+1]);
            ++idx;
            if ( ripple < 0.0 )
                ripple = 0.0;
            else if (ripple > 20.0)
                ripple = 20.0;
        }

        else
            fprintf(stderr, "ignoring unknown option '%s' at idx %d\n", argv[idx], idx);
    }

    if (1 >= argc || printUsage)
    {
        fprintf(stderr,
                "usage:\n"
                "  %s [-p|--pass <omega_pass>] [-s|--stop <omega_stop>] [-a|--atten <a>] [-r|--ripple <r>]\n"
                "estimates orders for Parks McClellan filter design - and Kaiser Window  and ..\n"
                , argv[0]);
        return 0;
    }

    if (OmegaPass < 0.0)
        OmegaPass = 0.2;
    if (OmegaStop < 0.0)
        OmegaStop = 1.3 * OmegaPass;
    if (OmegaStop > 1.0)
        OmegaStop = 1.0;

    double df = (OmegaStop - OmegaPass) / 2.0;
    double d1 = std::pow(10.0, ripple / 20.0) - 1.0;
    double d2 = std::pow(10.0, -atten / 20.0);

    if (decim < 0.0 && NumTaps < 0)
    {
        // formula identical to Kaiser's Formula from
        // https://dsp.stackexchange.com/questions/37646/filter-order-rule-of-thumb
        double nom = -10.0 * std::log10( d1 * d2) - 13.0;
        double den = 14.6 * df;
        double flt_M = std::ceil( 1.0 + nom / den );
        int M = int(flt_M + 0.5);
        printf("filter spec: %f .. %f omega => df %f, ripple %f, atten %f\n"
               , OmegaPass, OmegaStop, df, ripple, atten );
        printf("firpm filter len M = %d\n", M);
    }

    if (decim < 0.0 && NumTaps < 0)
    {
        // see https://dsp.stackexchange.com/questions/20514/what-is-the-relation-of-the-transition-bands-width-and-the-filter-order-for-the/20515
        double normFstop = OmegaStop / 2.0;
        double normFpass = OmegaPass / 2.0;
        //double normFcut = (normFstop + normFpass) / 2.0f;   //low pass filter 6dB cutoff

        //calculate Kaiser-Bessel window shape factor, Beta, from stopband attenuation
        double beta;
        if (atten < 20.96)
            beta = 0.0;
        else if (atten >= 50.0)
            beta = 0.1102 * (atten - 8.71);
        else
            beta = 0.5842 * std::pow((atten - 20.96), 0.4) + 0.07886 * (atten - 20.96);

        int NumTaps = (atten - 8.0) / (2.285 * M_2PI * (normFstop - normFpass) ) + 1;
        printf("kaiser design filter len M = %d  with beta = %f\n", NumTaps, beta);
    }

    if (decim < 0.0 && NumTaps < 0)
    {
        // https://dsp.stackexchange.com/questions/37646/filter-order-rule-of-thumb
        double flt_M = ( 2.0 / df ) * (atten / 22.0);
        int NumTaps = int( ceil(flt_M + 0.5) );
        printf("Fred Harris 'rule of thumb' filter len = %d\n", NumTaps);
    }


    return 0;
}

