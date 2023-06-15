# -*- coding: utf8 -*-
import numpy as np
import sys
try:
    from .rec_filter import recursive_filter
    from .rec_rms import recursive_rms
    from .rec_hos import recursive_hos
    from .rec_gauss_filter import recursive_gauss_filter
    from .rosenberger import rosenberger
except ImportError:
    from rec_filter import recursive_filter
    from rec_rms import recursive_rms
    from rec_hos import recursive_hos
    from rec_gauss_filter import recursive_gauss_filter
    from rosenberger import rosenberger


def make_LinFq(f_min, f_max, delta, nfreq):
    """Calculate linearly spaced frequency array for MB filtering."""
    f_ny = float(1./(2*delta))
    if f_max > f_ny:
        f_max = f_ny
    freq = np.linspace(f_min, f_max, nfreq)
    return freq


def make_LogFq(f_min, f_max, delta, nfreq):
    """Calculate log spaced frequency array for MB filtering."""
    f_ny = float(1./(2*delta))
    if f_max > f_ny:
        f_max = f_ny
    freq = np.logspace(np.log2(f_min), np.log2(f_max), nfreq, base=2)
    return freq


def MBfilter_CF(st, frequencies,
                CN_HP, CN_LP,
                filter_norm, filter_npoles=2,
                var_w=True,
                CF_type='envelope', CF_decay_win=1.0,
                hos_order=4,
                rosenberger_decay_win=1.0,
                rosenberger_filter_power=1.0,
                rosenberger_filter_threshold=None,
                rosenberger_normalize_each=False,
                wave_type='P',
                hos_sigma=None,
                rec_memory=None,
                full_output=False):
    """
    Multiband filtering.

    Performs MBfiltering using 2HP+2LP recursive filter
    and calculates the characteristic function (CF)
    for each band.
    """
    delta = st[0].stats.delta
    Tn = 1. / frequencies
    Nb = len(frequencies)
    CF_decay_nsmps = CF_decay_win / delta
    rosenberger_decay_nsmps = rosenberger_decay_win / delta

    if hos_sigma is None:
        hos_sigma = -1.

    # Single component analysis
    if len(st) < 2:
        # Use just the first trace in stream
        tr = st[0]
        y = tr.data

        YN1 = np.zeros((Nb, len(y)), float)
        CF1 = np.zeros((Nb, len(y)), float)

        for n in range(Nb):
            if rec_memory is not None:
                rmem = rec_memory[(tr.id, wave_type)][n]
            else:
                rmem = None

            YN1[n] = recursive_filter(y, CN_HP[n], CN_LP[n],
                                      filter_npoles, rmem)
            YN1[n] /= filter_norm[n]

            if var_w and CF_type == 'envelope':
                CF_decay_nsmps_mb = (Tn[n]/delta) * CF_decay_nsmps
            else:
                CF_decay_nsmps_mb = CF_decay_nsmps

            # Define the decay constant
            CF_decay_constant = 1 / CF_decay_nsmps_mb

            # Calculates CF for each MBF signal
            if CF_type == 'envelope':
                CF1[n] = recursive_rms(YN1[n], CF_decay_constant, rmem)

            if CF_type == 'kurtosis':
                CF1[n] = recursive_hos(YN1[n], CF_decay_constant,
                                       hos_order, hos_sigma, rmem)

    # 2 (horizontal) components analysis
    elif len(st) == 2:
        # Assumes that 2 horizontal components are used
        tr1 = st.select(channel='*[E,W,1]')[0]
        tr2 = st.select(channel='*[N,S,2]')[0]

        y1 = tr1.data
        y2 = tr2.data

        # Initializing arrays
        YN_E = np.zeros((Nb, len(y1)), float)
        YN_N = np.zeros((Nb, len(y1)), float)
        YN1 = np.zeros((Nb, len(y1)), float)
        CF1 = np.zeros((Nb, len(y1)), float)

        for n in range(Nb):
            if rec_memory is not None:
                rmem1 = rec_memory[(tr1.id, wave_type)][n]
                rmem2 = rec_memory[(tr2.id, wave_type)][n]
            else:
                rmem1 = None
                rmem2 = None

            YN_E[n] = recursive_filter(y1, CN_HP[n], CN_LP[n],
                                       filter_npoles, rmem1)
            YN_E[n] /= filter_norm[n]
            YN_N[n] = recursive_filter(y2, CN_HP[n], CN_LP[n],
                                       filter_npoles, rmem2)
            YN_N[n] /= filter_norm[n]
            # Combining horizontal components
            YN1[n] = np.sqrt(np.power(YN_E[n], 2) + np.power(YN_N[n], 2))

            if var_w and CF_type == 'envelope':
                CF_decay_nsmps_mb = (Tn[n] / delta) * CF_decay_nsmps
            else:
                CF_decay_nsmps_mb = CF_decay_nsmps

            # Define the decay constant
            CF_decay_constant = 1 / CF_decay_nsmps_mb

            # Calculates CF for each MBF signal
            if CF_type == 'envelope':
                CF1[n] = recursive_rms(YN1[n], CF_decay_constant, rmem1)

            if CF_type == 'kurtosis':
                CF1[n] = recursive_hos(YN1[n], CF_decay_constant,
                                       hos_order, hos_sigma, rmem1)

    # 3 components analysis, includes polarization P and S decomposition
    else:
        # Vertical
        tr1 = st.select(channel='*[Z,U,D]')[0]
        # Horizontals
        tr2 = st.select(channel='*[E,W,1]')[0]
        tr3 = st.select(channel='*[N,S,2]')[0]

        y1 = tr1.data
        y2 = tr2.data
        y3 = tr3.data

        # Initializing arrays
        YN1 = np.zeros((Nb, len(y1)), float)
        YN2 = np.zeros((Nb, len(y1)), float)
        YN3 = np.zeros((Nb, len(y1)), float)
        CF1 = np.zeros((Nb, len(y1)), float)
        filteredDataP = np.zeros((Nb, len(y1)), float)
        filteredDataS = np.zeros((Nb, len(y1)), float)
        if full_output:
            CF2 = np.zeros((Nb, len(y1)), float)

        for n in range(Nb):
            if rec_memory is not None:
                rmem1 = rec_memory[(tr1.id, wave_type)][n]
                rmem2 = rec_memory[(tr2.id, wave_type)][n]
                rmem3 = rec_memory[(tr3.id, wave_type)][n]
            else:
                rmem1 = None
                rmem2 = None
                rmem3 = None

            YN1[n] = recursive_filter(y1, CN_HP[n], CN_LP[n],
                                      filter_npoles, rmem1)
            YN1[n] /= filter_norm[n]
            YN2[n] = recursive_filter(y2, CN_HP[n], CN_LP[n],
                                      filter_npoles, rmem2)
            YN2[n] /= filter_norm[n]
            YN3[n] = recursive_filter(y3, CN_HP[n], CN_LP[n],
                                      filter_npoles, rmem3)
            YN3[n] /= filter_norm[n]

            # Define the decay constant
            rosenberger_decay_constant = 1 / rosenberger_decay_nsmps

            # print('Rosenberger in process {}/{}\r'.format(n+1, Nb),
            #       sys.stdout.flush())

            # third value returned by rosenberger() is the polarizaion filter,
            # which we do not use here
            filt_dataP, filt_dataS, _ =\
                rosenberger(YN2[n], YN3[n], YN1[n],
                            rosenberger_decay_constant,
                            pol_filter_power=rosenberger_filter_power,
                            pol_filter_threshold=rosenberger_filter_threshold,
                            normalize_each=rosenberger_normalize_each)

            # Use vertical component for P data
            filteredDataP[n] = filt_dataP[0, :]
            # Use vector composition of the two horizontal component for S data
            filteredDataS[n] = np.sqrt(np.power(filt_dataS[1, :], 2) +
                                       np.power(filt_dataS[2, :], 2))

            if var_w and CF_type == 'envelope':
                CF_decay_nsmps_mb = (Tn[n]/delta) * CF_decay_nsmps
            else:
                CF_decay_nsmps_mb = CF_decay_nsmps

            # Define the decay constant
            CF_decay_constant = 1 / CF_decay_nsmps_mb

            if CF_type == 'envelope':
                if wave_type == 'P':
                    CF1[n] = recursive_rms(filteredDataP[n],
                                           CF_decay_constant, rmem1)
                    if full_output:
                        CF2[n] = recursive_rms(filteredDataS[n],
                                               CF_decay_constant, rmem2)
                else:
                    CF1[n] = recursive_rms(filteredDataS[n],
                                           CF_decay_constant, rmem1)
                    if full_output:
                        CF2[n] = recursive_rms(filteredDataP[n],
                                               CF_decay_constant, rmem2)

            if CF_type == 'kurtosis':
                if wave_type == 'P':
                    CF1[n] = recursive_hos(filteredDataP[n],
                                           CF_decay_constant,
                                           hos_order, hos_sigma, rmem1)
                    if full_output:
                        CF2[n] = recursive_hos(filteredDataS[n],
                                               CF_decay_constant,
                                               hos_order, hos_sigma, rmem2)
                else:
                    CF1[n] = recursive_hos(filteredDataS[n],
                                           CF_decay_constant,
                                           hos_order, hos_sigma, rmem1)
                    if full_output:
                        CF2[n] = recursive_hos(filteredDataP[n],
                                               CF_decay_constant,
                                               hos_order, hos_sigma, rmem2)

    if full_output:
        return YN1, CF1, CF2, Tn, Nb, filteredDataP, filteredDataS
    else:
        return YN1, CF1, Tn, Nb


def GaussConv(data_in, sigma):
    derivative_data = np.gradient(data_in)
    derivative_data[derivative_data < 0] = 0
    CF_gaussian = recursive_gauss_filter(derivative_data, sigma)
    return CF_gaussian


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from generate_signal import generate_signal_noise2, generate_signal_expSin
    from rec_filter import rec_filter_coeff, rec_filter_norm
    from obspy.core import read, Trace, Stream

    # if arguments, read the file
    if len(sys.argv) == 2:
        filename = sys.argv[1]
        print(filename)
        data = read(filename)
        signal = np.array(data[0].data, dtype='float')
        sys.exitfunc()
    elif len(sys.argv) >= 2:
        # to be completed to account for extra arguments
        pass
    else:
        noise = generate_signal_noise2(1000, 0.05)
        signal = generate_signal_expSin(300, 0.005, 0.5, noise,
                                        0.5, 500, 0.05, 1)

    # sampling interval
    # TODO: fix code below
    # if 'data' in locals():
    #    delta = data[0].stats.sampling_rate
    # else:
    delta = 0.01

    tr = Trace(signal)
    tr.stats.delta = delta
    tr.stats.channel = 'HHZ'
    st = Stream(tr)

    frequencies = make_LinFq(1, 40, delta, 8)
    CN_HP, CN_LP = rec_filter_coeff(frequencies, delta)
    filter_norm = rec_filter_norm(frequencies, delta, CN_HP, CN_LP, 2)
    YN, CF, Tn, Nb = MBfilter_CF(st, frequencies, CN_HP, CN_LP, filter_norm,
                                 var_w=False, CF_type='kurtosis',
                                 CF_decay_win=0.1)

    fig1 = plt.figure()
    fig2 = plt.figure()
    max2 = CF.max()

    for n in range(Nb):
        ax1 = fig1.add_subplot(Nb+1, 1, n+1)
        ax2 = fig2.add_subplot(Nb+1, 1, n+1)
        # ax2.set_ylim((0, max2))
        ax1.plot(YN[n])
        ax2.plot(CF[n])

    ax1.plot(signal, 'g')
    ax2.plot(recursive_hos(signal, 0.1, order=4), 'g')

    plt.show()
