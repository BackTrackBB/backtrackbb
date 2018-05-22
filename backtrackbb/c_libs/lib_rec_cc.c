#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void initlib_rec_cc() {}
void PyInit_lib_rec_cc() {}

/* C implementation of the algorithm by Young&Vliet, 1995 */
void _gausscoeff(double sigma, double *A, int *nA, double *B, int *nB)
{
    int i;
    double q;
    double *b;

    if (sigma > 0.5) {
        q = 0.98711*sigma - 0.96330;
    } else if (sigma == 0.5) {
        q = 3.97156 - 4.14554 * sqrt(1.0 - 0.26891*sigma);
    } else {
        perror("Sigma for Gaussian filter must be >=0.5 samples.\n");
        exit(1);
    }

    b = (double *) malloc(4 * sizeof(double));
    b[0] = 1.57825 + 2.44413*q + 1.4281*pow(q, 2) + 0.422205*pow(q, 3);
    b[1] = 2.44413*q + 2.85619*pow(q, 2) + 1.26661*pow(q, 3);
    b[2] = -(1.4281*pow(q, 2) + 1.26661*pow(q, 3));
    b[3] = 0.422205*pow(q, 3);

    B[0] = 1.0 - ((b[1] + b[2] + b[3])/b[0]);

    A[0] = 1;
    for (i=1; i<4; i++) {
        A[i] = -b[i]/b[0];
    }

    *nA = 4;
    *nB = 1;
    free(b);
}


void _lfilter(const double *signal, double *filt_signal, int npts,
             const double *A, int nA, const double *B, int nB)
{
    //a[0]*y[n] = b[0]*x[n] + b[1]*x[n-1] + ... + b[nb]*x[n-nb]
    //                    - a[1]*y[n-1] - ... - a[na]*y[n-na]
    int n, na, nb;

    for (n=0; n<npts; n++) {
        filt_signal[n] = 0;
        for (nb=0; nb<nB; nb++) {
            if (nb > n) break;
            filt_signal[n] += B[nb] * signal[n-nb];
        }
        for (na=1; na<nA; na++) {
            if (na > n) break;
            filt_signal[n] -= A[na] * filt_signal[n-na];
        }
        filt_signal[n] /= A[0];
    }
}


void _reverse(const double *signal, double *rev_signal, int npts)
{
    int n, end;
    double tmp;

    for (n=0; n<npts; n++)
        rev_signal[n] = signal[n];

    end = npts - 1;
    for (n = 0; n < npts/2; n++) {
        tmp = rev_signal[n];
        rev_signal[n] = rev_signal[end];
        rev_signal[end] = tmp;
        end--;
    }
}

#ifdef _MSC_VER
__declspec(dllexport)
#endif
void _Gaussian1D(double *signal, int npts, double sigma)
{
    // signal is overwritten
    double *A=NULL, *B=NULL;
    int n, nA, nB;
    double *rev_filt_signal, *filt_signal;

    if (npts < 4) {
        perror("Signal too short\n");
        exit(1);
    }

    A = (double *) malloc(4 * sizeof(double));
    B = (double *) malloc(4 * sizeof(double));
    _gausscoeff(sigma, A, &nA, B, &nB);

    rev_filt_signal = (double *) malloc(npts * sizeof(double));
    filt_signal = (double *) malloc(npts * sizeof(double));

    _lfilter(signal, filt_signal, npts, A, nA, B, nB);
    _reverse(filt_signal, rev_filt_signal, npts);
    _lfilter(rev_filt_signal, filt_signal, npts, A, nA, B, nB);
    _reverse(filt_signal, rev_filt_signal, npts);
    for (n=0; n<npts; n++)
        signal[n] = rev_filt_signal[n];

    free(A);
    free(B);
    free(rev_filt_signal);
    free(filt_signal);
}


#ifdef _MSC_VER
__declspec(dllexport)
#endif
void _local_CCr(const double *signal1, const double *signal2, int npts,
                double *cc, int lmax, double sigma)
{
    int n;
    int l, l_f, l_g;
    double *_cc;
    double shift1_1, shift2_2, shift1_2, shift2_1;

    for (l=-lmax; l<lmax; l++) {
        l_f = (int) floor(l/2.);
        l_g = (int) ceil(l/2.);

        _cc = &cc[npts*(l+lmax)];

        for(n=0; n<npts; n++) {
            if ((n - l_f) >= 0 && (n - l_f) < npts)
               shift1_1 = signal1[n - l_f];
            else
               shift1_1 = 0;
            if ((n + l_g) >= 0 && (n + l_g) < npts)
               shift2_2 = signal2[n + l_g];
            else
               shift2_2 = 0;
            if ((n - l_g) >= 0 && (n - l_g) < npts)
               shift1_2 = signal1[n - l_g];
            else
               shift1_2 = 0;
            if ((n + l_f) >= 0 && (n + l_f) < npts)
               shift2_1 = signal2[n + l_f];
            else
               shift2_1 = 0;
            _cc[n] = 0.5 * (shift1_1 * shift2_2 + shift1_2 * shift2_1);
        }
        if (sigma > 0)
            _Gaussian1D(_cc, npts, sigma);
    }
}
