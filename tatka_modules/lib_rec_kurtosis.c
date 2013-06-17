#include <stdlib.h>
#include <math.h> 

void _recursive_kurtosis(const double *signal, double *kurt_signal, int npts, float C_WIN)
{
    int i;
    
    double mean = 0;
    double var = 1.0;
    double kurt = 0;
    
    for (i=0; i<npts; i++) {
		mean = C_WIN * signal[i] + (1 - C_WIN) * mean;
		var = C_WIN * pow((signal[i] - mean),2.0) + (1 - C_WIN) * var;
		kurt = C_WIN * (pow((signal[i]-mean),4.0)/pow((var),2.0)) + (1 - C_WIN) * kurt;
		kurt_signal[i] = kurt;
    }
}
