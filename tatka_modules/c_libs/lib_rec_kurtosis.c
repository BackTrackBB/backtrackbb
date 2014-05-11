#include <stdlib.h>
#include <math.h>

void _recursive_kurtosis(const double *signal, double *kurt_signal, int npts,
        float sigma_min, float C_WIN, int order1, int order2, float power2)
{
    int i;
    double mean = 0;
    double var_temp = 1.0;
    double var =1;
    double kurt = 0;

    for (i=0; i<npts; i++) {
        mean = C_WIN * signal[i] + (1 - C_WIN) * mean;

        var_temp = C_WIN * pow((signal[i] - mean),order2) + (1 - C_WIN) * var;
        if (var_temp > sigma_min) {
            var = var_temp;
        } else {
            var = sigma_min;
        }

        kurt = C_WIN * (pow((signal[i]-mean), order1) / pow((var), power2)) + (1 - C_WIN) * kurt;
        kurt_signal[i] = kurt;
    }
}
