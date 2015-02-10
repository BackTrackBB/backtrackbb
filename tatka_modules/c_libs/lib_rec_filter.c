#include <stdlib.h>

void _recursive_filter_BP(const double *signal, double *filt_signal, int npts, float C_HP, float C_LP,
        double *filterH1, double *filterH2, double *filterL1, double *filterL2, double *previous_sample)
{
    int i;
    double filterH1_0 = 0;

    for (i=0; i<npts; i++) {
        filterH1_0 = *filterH1;
        *filterH1 = C_HP * (*filterH1 + signal[i] - *previous_sample);
        *filterH2 = C_HP * (*filterH2 + *filterH1 - filterH1_0);
        *filterL1 = *filterL1 + C_LP * (*filterH2 - *filterL1);
        *filterL2 = *filterL2 + C_LP * (*filterL1 - *filterL2);
        filt_signal[i] = *filterL2;
        *previous_sample = signal[i];
    }
}

void _recursive_filter_HP(const double *signal, double *filt_signal, int npts, float C_HP,
        double *filterH1, double *filterH2, double *previous_sample)
{
    int i;
    double filterH1_0 = 0;

    for (i=0; i<npts; i++) {
        filterH1_0 = *filterH1;
        *filterH1 = C_HP * (*filterH1 + signal[i] - *previous_sample);
        *filterH2 = C_HP * (*filterH2 + *filterH1 - filterH1_0);
        filt_signal[i] = *filterH2;
        *previous_sample = signal[i];
    }
}
