import numpy as np
import scipy as sp
import sys
from obspy.signal.util import smooth
from scipy.signal import gaussian
from rec_filter import recursive_filter
from rec_rms import recursive_rms
from rec_hos import recursive_hos
#from RosenbergerAlgorithm import rosenberger
from rosenberger import rosenberger


def make_LinFq(f_min, f_max, delta, nfreq):
    """
    Calculates linearly spaced frequency array for MB filtering
    """
    f_ny = float(1./(2*delta))
    if f_max > f_ny:
        f_max = f_ny
    freq = np.linspace(f_min, f_max, nfreq)
    return freq[::-1]


def make_LogFq(f_min, f_max, delta, nfreq):
    """
    Calculates log spaced frequency array for MB filtering
    """
    f_ny = float(1./(2*delta))
    if f_max > f_ny:
        f_max = f_ny
    freq = np.logspace(np.log2(f_min), np.log2(f_max), nfreq, base=2)
    return freq[::-1]


def MBfilter_CF(st, frequencies, var_w=True,
                CF_type='envelope', CF_decay_win=1.0,
                filter_type='bandpass',
                order=4, rosenberger_decay_win=1.0,
                wave_type='P', hos_sigma=None,
                rec_memory=None,
                full_output=False):
    """
    Performs MBfiltering using 2HP+2LP recursive filter
    and calculates the characteristic function (CF)
    for each band.
    """
    delta = st[0].stats.delta
    Nb = len(frequencies)
    Tn = 1./frequencies
    wn = Tn/(2*np.pi)
    CN_HP = wn/(wn+delta)        # high-pass filter constant
    if filter_type == 'bandpass':
        CN_LP = delta/(wn+delta)        # low-pass filter constant
    elif filter_type == 'highpass':
        CN_LP = [None,] * len(CN_HP)
    else:
        raise ValueError, 'Wrong filter type: %s' % filter_type
    CF_decay_nsmps = int(CF_decay_win / delta)
    rosenberger_decay_nsmps = int(rosenberger_decay_win / delta)

    if hos_sigma is None:
        hos_sigma = -1.

    # Less than 3 components
    if len(st) < 3:
        # Use just the first trace in stream
        tr = st[0]
        y = tr.data
        y = y - y.mean()
        y[0] = np.mean(y[0:CF_decay_nsmps])

        YN1 = np.zeros((Nb, len(y)), float)
        CF1 = np.zeros((Nb, len(y)), float)

        for n in xrange(Nb):
            if rec_memory is not None:
                rmem = rec_memory[n]
            else:
                rmem = None
            YN1[n] = recursive_filter(y, CN_HP[n], CN_LP[n], rmem)

            if var_w and CF_type == 'envelope':
                CF_decay_nsmps_mb = (Tn[n]/delta)*CF_decay_nsmps
            else:
                CF_decay_nsmps_mb = CF_decay_nsmps

            if CF_type == 'envelope':
                CF1[n] = recursive_rms(YN1[n], 1./CF_decay_nsmps_mb, rmem)

            if CF_type == 'hilbert':
                # This does not support recursive memory!
                CF1[n] = smooth(abs(sp.signal.hilbert(YN1[n])), CF_decay_nsmps_mb)

            if CF_type == 'kurtosis':
                CF1[n] = recursive_hos(YN1[n], 1./CF_decay_nsmps_mb,
                                      hos_sigma, order, rmem)

    # More than 3 components
    else:
        # We assume that component names are: Z, E, N
        #TODO: generalize this
        tr1 = st.select(component='Z')[0]
        tr2 = st.select(component='E')[0]
        tr3 = st.select(component='N')[0]

        y1 = tr1.data
        y2 = tr2.data
        y3 = tr3.data

        y1 = y1 - y1.mean()
        y1[0] = np.mean(y1[0:CF_decay_nsmps])
        y2 = y2 - y2.mean()
        y2[0] = np.mean(y2[0:CF_decay_nsmps])
        y3 = y3 - y3.mean()
        y3[0] = np.mean(y3[0:CF_decay_nsmps])

        # Initializing arrays
        YN1 = np.zeros((Nb, len(y1)), float)
        YN2 = np.zeros((Nb, len(y1)), float)
        YN3 = np.zeros((Nb, len(y1)), float)
        CF1 = np.zeros((Nb, len(y1)), float)
        filteredDataP = np.zeros((Nb, len(y1)), float)
        filteredDataS = np.zeros((Nb, len(y1)), float)
        if full_output:
            CF2 = np.zeros((Nb, len(y1)), float)

        for n in xrange(Nb):
            YN1[n] = recursive_filter(y1, CN_HP[n], CN_LP[n])
            YN2[n] = recursive_filter(y2, CN_HP[n], CN_LP[n])
            YN3[n] = recursive_filter(y3, CN_HP[n], CN_LP[n])

            print 'Rosenberger in process {}/{}\r'.format(n+1, Nb),
            sys.stdout.flush()

            filt_dataP, filt_dataS, U =\
                    rosenberger(YN2[n], YN3[n], YN1[n], 1./rosenberger_decay_nsmps)

            # Use vertical component for P data
            filteredDataP[n] = filt_dataP[0,:]
            # Use vector composition of the two horizontal component for S data
            filteredDataS[n] = np.sqrt(np.power(filt_dataS[1,:], 2) +
                                       np.power(filt_dataS[2,:], 2))

            if var_w:
                CF_decay_nsmps_mb = b * Tn[n]/delta
            else:
                CF_decay_nsmps_mb = CF_decay_nsmps

            if CF_type == 'envelope':
                if wave_type == 'P':
                    CF1[n] = recursive_rms(filteredDataP[n], 1./CF_decay_nsmps_mb)
                    if full_output:
                        CF2[n] = recursive_rms(filteredDataS[n], 1./CF_decay_nsmps_mb)
                else:
                    CF1[n] = recursive_rms(filteredDataS[n], 1./CF_decay_nsmps_mb)
                    if full_output:
                        CF2[n] = recursive_rms(filteredDataP[n], 1./CF_decay_nsmps_mb)

            if CF_type == 'kurtosis':
                if wave_type == 'P':
                    CF1[n] = recursive_hos(filteredDataP[n], 1./CF_decay_nsmps_mb,
                                          hos_sigma, order)
                    if full_output:
                        CF2[n] = recursive_hos(filteredDataS[n], 1./CF_decay_nsmps_mb,
                                          hos_sigma, order)
                else:
                    CF1[n] = recursive_hos(filteredDataS[n], 1./CF_decay_nsmps_mb,
                                          hos_sigma, order1, order2, power2)
                    if full_output:
                        CF2[n] = recursive_hos(filteredDataP[n], 1./CF_decay_nsmps_mb,
                                          hos_sigma, order)

    if full_output:
        return YN1, CF1, CF2, Tn, Nb, filteredDataP, filteredDataS
    else:
        return YN1, CF1, Tn, Nb


def GaussConv(data_in, sigma):
    derivative_data = np.gradient(data_in)
    derivative_data[derivative_data < 0] = 0
    gauss_window = gaussian(len(derivative_data), sigma)
    CF_gaussian = np.convolve(derivative_data, gauss_window, mode='same')
    return CF_gaussian


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from generate_signal import generate_signal_noise2, generate_signal_expSin
    from obspy.core import read, Trace, Stream

    #if arguments, read the file
    if len(sys.argv) == 2:
        filename = sys.argv[1]
        print filename
        data = read(filename)
        signal = np.array(data[0].data, dtype='float')
        sys.exitfunc()
    elif len(sys.argv) >= 2:
        #to be completed to account for extra arguments
        pass
    else:
        noise = generate_signal_noise2(1000, 0.05)
        signal = generate_signal_expSin(300, 0.005, 0.5, noise, 0.5, 500, 0.05, 1)

    # sampling interval
    #TODO: fix code below
    #if 'data' in locals():
    #    delta = data[0].stats.sampling_rate
    #else:
    delta = 0.01
    print 'Sampling interval = {}'.format(delta)

    tr = Trace(signal)
    tr.stats.delta = delta
    tr.stats.channel = 'HHZ'
    st = Stream(tr)

    frequencies = make_LinFq(1, 40, delta, 8)
    YN, CF, Tn, Nb = MBfilter_CF(st, frequencies, var_w=False,
                                 CF_type='kurtosis', CF_decay_win=0.1)

    fig1 = plt.figure()
    fig2 = plt.figure()
    max2 = CF.max()

    for n in range(Nb):
        ax1 = fig1.add_subplot(Nb+1, 1, n+1)
        ax2 = fig2.add_subplot(Nb+1, 1, n+1)
        #ax2.set_ylim((0, max2))
        ax1.plot(YN[n])
        ax2.plot(CF[n])

    ax1.plot(signal, 'g')
    ax2.plot(recursive_hos(signal, 0.1, 0.001, 4, 2, 2),'g')

    plt.show()
