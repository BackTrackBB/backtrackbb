#include <stdlib.h>
#include <math.h>

void initlib_rec_rms() {}
void PyInit_lib_rec_rms() {}

#ifdef _MSC_VER
__declspec(dllexport)
#endif
void _recursive_rms(const double *signal, double *rms_signal, int npts, float C_WIN,
        double *mean_sq, int memory_sample, int initialize)
{
    int i,j;
    int n_win = (int) 1/C_WIN;
    double _mean_sq = *mean_sq;

    if (memory_sample < 0 || memory_sample >= npts) {
        memory_sample = npts-1;
    }

    if (initialize) {
        for (j=0; j<n_win; j++) {
            _mean_sq = _mean_sq + pow(signal[j], 2);
        }
        _mean_sq = sqrt(_mean_sq/n_win);
    }

    for (i=0; i<npts; i++) {
        _mean_sq = sqrt(C_WIN * pow(signal[i], 2.) + (1 - C_WIN) * pow(_mean_sq, 2.));
        rms_signal[i] = _mean_sq;
        /* save memory values */
        if (i == memory_sample) {
            *mean_sq = _mean_sq;
        }
    }
}
