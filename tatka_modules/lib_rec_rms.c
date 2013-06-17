#include <stdlib.h>
#include <math.h> 

void _recursive_rms(const double *signal, double *rms_signal, int npts, float C_WIN)
{
    int i;

    double mean_sq = 0;
    
    for (i=0; i<npts; i++) {
	mean_sq = sqrt(C_WIN*pow(signal[i],2.)+(1-C_WIN)*pow(mean_sq,2.));	
	rms_signal[i] = mean_sq;
    }
}
