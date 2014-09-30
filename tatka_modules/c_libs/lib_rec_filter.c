#include <stdlib.h>

void _recursive_filter(const double *signal, double *filt_signal, int npts, float C_HP, float C_LP)
{
    int i;

    double filterH1 = 0;
    double filterH1_0 = 0;
    double filterH2 = 0;
    double filterL1 = 0;
    double filterL2 = 0;

    filt_signal[0] = signal[0];

    for (i=0; i<npts; i++) {
        filterH1_0 = filterH1;
        filterH1 = C_HP * (filterH1 + signal[i] - signal[i-1]);
        filterH2 = C_HP * (filterH2 + filterH1 - filterH1_0);
        filterL1 = filterL1 + C_LP * (filterH2 - filterL1);
        filterL2 = filterL2 + C_LP * (filterL1 - filterL2);
        filt_signal[i] = filterL2;
    }
}
