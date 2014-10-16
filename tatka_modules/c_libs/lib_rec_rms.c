#include <stdlib.h>
#include <math.h>

void _recursive_rms(const double *signal, double *rms_signal, int npts, float C_WIN)
{
    int i,j;
    int n_win;
    double mean_sq = 0;

    n_win = (int) 1/C_WIN;
    for(j=0;j<n_win;j++){
	mean_sq = mean_sq + pow(signal[j],2);
    }
    mean_sq = sqrt(mean_sq/n_win);

    for (i=0; i<npts; i++) {
        mean_sq = sqrt(C_WIN * pow(signal[i], 2.) + (1 - C_WIN) * pow(mean_sq, 2.));
        rms_signal[i] = mean_sq;
    }
}
