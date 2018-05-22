/*
 * lib_rec_hos.c
 *
 * Recursive higher-order statistics.
 *
 * (c) 2015-2018 - Natalia Poiata <poiata@ipgp.fr>,
 *                 Claudio Satriano <satriano@ipgp.fr>
 * (c) 2013-2014 - Natalia Poiata <poiata@ipgp.fr>,
 *                 Claudio Satriano <satriano@ipgp.fr>,
 *                 Pierre Romanet <romanet@ipgp.fr>
 */
#include <stdlib.h>
#include <math.h>

void initlib_rec_hos() {}
void PyInit_lib_rec_hos() {}

#ifdef _MSC_VER
__declspec(dllexport)
#endif
void _recursive_hos(const double *signal, double *hos_signal, int npts,
        float sigma_min, float C_WIN, int order,
        double *mean, double *var, double *hos, int memory_sample, int initialize)
{
    int i;
    int n_win = (int) 1/C_WIN;
    double power = order/2;
    double var_temp;
    double _mean = *mean;
    double _var = *var;
    double _hos = *hos;

    if (memory_sample < 0 || memory_sample >= npts) {
        memory_sample = npts-1;
    }

    if (initialize) {
        for (i=0; i<n_win; i++) {
            _mean = C_WIN * signal[i] + (1 - C_WIN) * _mean;
            _var = C_WIN * pow((signal[i] - _mean), 2.0) + (1 - C_WIN) * _var;
        }
    }

    for (i=0; i<npts; i++) {
        _mean = C_WIN * signal[i] + (1 - C_WIN) * _mean;
        var_temp = C_WIN * pow((signal[i] - _mean), 2.0) + (1 - C_WIN) * _var;
        if (var_temp > sigma_min) {
            /* if sigma_min < 0, this is always true */
            _var = var_temp;
        } else {
            _var = sigma_min;
        }
        _hos = C_WIN * (pow((signal[i] - _mean), order) / pow(_var, power)) + (1 - C_WIN) * _hos;
        hos_signal[i] = _hos;
        /* save memory values */
        if (i == memory_sample) {
            *mean = _mean;
            *var = _var;
            *hos = _hos;
        }
    }
}
