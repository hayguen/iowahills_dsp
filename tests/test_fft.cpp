
#include <iowahills/fft.h>

#include <stdio.h>
#include <cmath>
#include <vector>

#define M_2PI     6.28318530717958647692   // 2*Pi


typedef std::vector<double> vecT;

void copy_input(int N, const vecT& inp_i, const vecT& inp_q, vecT& out_i, vecT& out_q)
{
    for (int n = 0; n < N; ++n) {
        out_i[n] = inp_i[n];
        out_q[n] = inp_q[n];
    }
}

int compare_result(double eps, int N, const vecT& ref_i, const vecT& ref_q, const vecT& out_i, const vecT& out_q)
{
    for (int k = 0; k < N; ++k) {
        double d_i = out_i[k] - ref_i[k];
        double d_q = out_q[k] - ref_q[k];
        if ( std::abs(d_i) >= eps || std::abs(d_q) >= eps )
        {
            fprintf(stderr, "difference in output to reference at index k = %d: %f + i* %f\n", k, d_i, d_q);
            fprintf(stderr, "  reference: %f + i* %f\n", ref_i[k], ref_q[k]);
            fprintf(stderr, "  output:    %f + i* %f\n", out_i[k], out_q[k]);
            return 1;
        }
    }
    return 0;
}



int main(int, char* [])
{
    const int N = RequiredFFTSize(16);
    std::vector<double> inp_i(N), inp_q(N);
    std::vector<double> wrk_i(N), wrk_q(N);
    std::vector<double> out_i(N), out_q(N);

    FFTBuffers * fft_handle = prepareFFT(N, FORWARD);
    if (!fft_handle)
        return 1;

    for (int k = 0; k < N; ++k )
    {
        // generate carrier with frequency for k : omega = (2*k) / N : 0 .. 1 = positive; 1 .. 2 => negative!
        const double k_dphi = k * M_2PI / N;
        for (int n = 0; n < N; ++n ) {
            double n_phi = n * k_dphi;
            inp_i[n] = cos(n_phi);
            inp_q[n] = sin(n_phi);
        }

        // run stateless FFT
        copy_input(N, inp_i, inp_q,  out_i, out_q);
        FFT(out_i.data(), out_q.data(), N, FORWARD);

        // run/compare stateless FFT again
        copy_input(N, inp_i, inp_q,  wrk_i, wrk_q);
        FFT(wrk_i.data(), wrk_q.data(), N, FORWARD);
        int r = compare_result(1E-8, N, out_i, out_q, wrk_i, wrk_q);
        if (r) {
            fprintf(stderr, "error in result comparison 1 for generated carrier frequency k = %d\n", k);
            return 1;
        }

        // run/compare 1st FFT with state
        copy_input(N, inp_i, inp_q,  wrk_i, wrk_q);
        const int skip_fwd_div_N = 0;
        FFTex(fft_handle, wrk_i.data(), wrk_q.data(), skip_fwd_div_N);
        r = compare_result(1E-8, N, out_i, out_q, wrk_i, wrk_q);
        if (r) {
            fprintf(stderr, "error in result comparison 1 for generated carrier frequency k = %d\n", k);
            return 1;
        }

        // run/compare 2nd FFT with state
        copy_input(N, inp_i, inp_q,  wrk_i, wrk_q);
        FFTex(fft_handle, wrk_i.data(), wrk_q.data(), skip_fwd_div_N);
        r = compare_result(1E-8, N, out_i, out_q, wrk_i, wrk_q);
        if (r) {
            fprintf(stderr, "error in result comparison 2 for generated carrier frequency k = %d\n", k);
            return 1;
        }


        for (int kk = 0; kk < N; ++kk )
        {
            double pwr = out_i[kk] * out_i[kk] + out_q[kk] * out_q[kk];
            double mag = std::sqrt(pwr);
            if ( k == kk )
            {
                if ( std::abs(mag - 1.0) > 1E-8 )
                {
                    fprintf(stderr, "error in result: generated carrier frequency k = %d  should produce magnitude ~= 1.0 at kk = %d\n", k, kk);
                    fprintf(stderr, "  produced magnitude is %f\n", mag);
                    return 1;
                }
            }
            else
            {
                if ( std::abs(mag) > 1E-8 )
                {
                    fprintf(stderr, "error in result: generated carrier frequency k = %d  should produce magnitude ~= 0.0 at kk = %d\n", k, kk);
                    fprintf(stderr, "  produced magnitude is %f\n", mag);
                    return 1;
                }
            }
        }

    }

    fprintf(stderr, "fft tests passed\n");
    disposeFFT(fft_handle);
    return 0;
}

