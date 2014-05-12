/*
 * lib_rec_hos.c
 *
 * Recursive higher-order statistics.
 *
 * (c) 2013-2014 - Natalia Poiata <poiata@ipgp.fr>,
 *                 Claudio Satriano <satriano@ipgp.fr>,
 *                 Pierre Romanet <romanet@ipgp.fr>
 */
#include <stdlib.h>
#include <math.h>

void _recursive_hos(const double *signal, double *hos_signal, int npts,
        float sigma_min, float C_WIN, int order1, int order2, float power2)
{
    int i;
    double mean = 0;
    double var_temp = 1.0;
    double var = 1;
    double hos = 0;

    for (i=0; i<npts; i++) {
        mean = C_WIN * signal[i] + (1 - C_WIN) * mean;

        var_temp = C_WIN * pow((signal[i] - mean),order2) + (1 - C_WIN) * var;
        if (var_temp > sigma_min) {
            var = var_temp;
        } else {
            var = sigma_min;
        }

        hos = C_WIN * (pow((signal[i]-mean), order1) / pow((var), power2)) + (1 - C_WIN) * hos;
        hos_signal[i] = hos;
    }
}
