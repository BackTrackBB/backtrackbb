/*
 * lib_rec_filter.c
 *
 * Recursive bandpass (or highpass) filtering
 *
 * (c) 2013-2018 - Natalia Poiata <poiata@ipgp.fr>,
 *                 Claudio Satriano <satriano@ipgp.fr>
 */
#include <stdlib.h>

void initlib_rec_filter() {}
void PyInit_lib_rec_filter() {}

/*
 * Bandpass (or highpass) filtering by cascade of simple, first-order
 * recursive highpass and lowpass filters. The number of poles gives the
 * number of filter stages.
 */
#ifdef _MSC_VER
__declspec(dllexport)
#endif
void _recursive_filter(const double *signal, double *filt_signal, int npts,
                       float C_HP, float C_LP, int npoles,
                       double *filterH, double *filterL,
                       double *prev_sample_value, int memory_sample)
{
    int i, n;
    double *_filterH;
    double *_filterH0;
    double *_filterL;
    double _prev_sample_value = *prev_sample_value;
    double s0, s1;

    _filterH = (double *) malloc(npoles * sizeof(double));
    _filterH0 = (double *) malloc(npoles * sizeof(double));
    _filterL = (double *) malloc(npoles * sizeof(double));
    for (n=0; n<npoles; n++) {
        _filterH[n] = filterH[n];
        _filterH0[n] = 0;
        _filterL[n] = filterL[n];
    }

    if (memory_sample < 0 || memory_sample >= npts) {
        memory_sample = npts-1;
    }

    for (i=0; i<npts; i++) {
        for (n=0; n<npoles; n++) {
            if (n == 0) {
                s0 = _prev_sample_value;
                s1 = signal[i];
            } else {
                s0 = _filterH0[n-1];
                s1 = _filterH[n-1];
            }
            _filterH0[n] = _filterH[n];
            _filterH[n] = C_HP * (_filterH[n] + s1 - s0);
        }
        if (C_LP < 0) {
            /* high-pass filter */
            filt_signal[i] = _filterH[npoles-1];
        } else {
            for (n=0; n<npoles; n++) {
                if (n == 0) {
                    s0 = _filterH[npoles-1];
                } else {
                    s0 = _filterL[n-1];
                }
                _filterL[n] = _filterL[n] + C_LP * (s0 - _filterL[n]);
            }
            filt_signal[i] = _filterL[npoles-1];
        }
        _prev_sample_value = signal[i];
        /* save memory values */
        if (i == memory_sample) {
            for (n=0; n<npoles; n++) {
                filterH[n] = _filterH[n];
                filterL[n] = _filterL[n];
            }
            *prev_sample_value = _prev_sample_value;
        }
    }

    free(_filterH);
    free(_filterL);
}
