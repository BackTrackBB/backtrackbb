#include <stdlib.h>

void _recursive_filter_BP(const double *signal, double *filt_signal, int npts, float C_HP, float C_LP,
        double *filterH1, double *filterH2, double *filterL1, double *filterL2, double *prev_sample_value,
        int memory_sample)
{
    int i;
    double filterH1_0 = 0;
    double _filterH1 = *filterH1;
    double _filterH2 = *filterH2;
    double _filterL1 = *filterL1;
    double _filterL2 = *filterL2;
    double _prev_sample_value = *prev_sample_value;

    if (memory_sample < 0 || memory_sample >= npts) {
        memory_sample = npts-1;
    }

    for (i=0; i<npts; i++) {
        filterH1_0 = _filterH1;
        _filterH1 = C_HP * (_filterH1 + signal[i] - _prev_sample_value);
        _filterH2 = C_HP * (_filterH2 + _filterH1 - filterH1_0);
        _filterL1 = _filterL1 + C_LP * (_filterH2 - _filterL1);
        _filterL2 = _filterL2 + C_LP * (_filterL1 - _filterL2);
        filt_signal[i] = _filterL2;
        _prev_sample_value = signal[i];
        /* save memory values */
        if (i == memory_sample) {
            *filterH1 = _filterH1;
            *filterH2 = _filterH2;
            *filterL1 = _filterL1;
            *filterL2 = _filterL2;
            *prev_sample_value = _prev_sample_value;
        }
    }
}

void _recursive_filter_HP(const double *signal, double *filt_signal, int npts, float C_HP,
        double *filterH1, double *filterH2, double *prev_sample_value, int memory_sample)
{
    int i;
    double filterH1_0 = 0;
    double _filterH1 = *filterH1;
    double _filterH2 = *filterH2;
    double _prev_sample_value = *prev_sample_value;

    if (memory_sample < 0 || memory_sample >= npts) {
        memory_sample = npts-1;
    }

    for (i=0; i<npts; i++) {
        filterH1_0 = _filterH1;
        _filterH1 = C_HP * (_filterH1 + signal[i] - _prev_sample_value);
        _filterH2 = C_HP * (_filterH2 + _filterH1 - filterH1_0);
        filt_signal[i] = _filterH2;
        _prev_sample_value = signal[i];
        /* save memory values */
        if (i == memory_sample) {
            *filterH1 = _filterH1;
            *filterH2 = _filterH2;
            *prev_sample_value = _prev_sample_value;
        }
    }
}
