import os
from ctypes import CDLL, c_int, c_float, c_double, c_void_p, POINTER, byref
from numpy.ctypeslib import ndpointer
import numpy as np


libpath = os.path.join(os.path.dirname(__file__), os.pardir, 'lib', 'lib_rec_filter.so')
lib_rec_filter = CDLL(libpath)

lib_rec_filter._recursive_filter.argtypes = [
        ndpointer(dtype=np.float64), #signal
        ndpointer(dtype=np.float64), #filt_signal
        c_int, #npts
        c_float, #C_HP
        c_float, #C_LP
        c_int, #npoles
        ndpointer(dtype=np.float64), #filterH
        ndpointer(dtype=np.float64), #filterL
        POINTER(c_double), #prev_sample_value
        c_int #memory_sample
        ]
lib_rec_filter._recursive_filter.restype = c_void_p


def recursive_filter(signal, C_HP, C_LP=None, npoles=2, rec_memory=None):
    signal = np.array(signal, dtype=np.float64)
    filt_signal = np.zeros(len(signal))

    if rec_memory is not None:
        filterH = rec_memory.filterH
        filterL = rec_memory.filterL
        prev_sample_value = rec_memory.prev_sample_value
        memory_sample = rec_memory.memory_sample
    else:
        filterH = np.zeros(npoles)
        filterL = np.zeros(npoles)
        prev_sample_value = c_double(0)
        memory_sample = -1

    if C_LP is None:
        C_LP = -1

    lib_rec_filter._recursive_filter(
        signal, filt_signal, signal.size, C_HP, C_LP,
        npoles, filterH, filterL,
        byref(prev_sample_value), memory_sample)

    if rec_memory is not None:
        rec_memory.filterH = filterH
        rec_memory.filterL = filterL
        rec_memory.prev_sample_value = prev_sample_value

    return filt_signal


def rec_filter_coeff(freq, delta):
    if isinstance(freq, tuple) or isinstance(freq, list):
        freq = np.array(freq)
    T = 1. / freq
    w = T / (2*np.pi)
    C_HP = w / (w + delta)      # high-pass filter constant
    C_LP = delta / (w + delta)  # low-pass filter constant
    return C_HP, C_LP


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    def do_fft(signal):
        npts = len(signal)
        if not npts % 2:
            npts -= 1
        fft = np.fft.rfft(signal, n=npts) * delta
        fftfreq = np.fft.fftfreq(len(signal), d=delta)
        fftfreq = fftfreq[0:fft.size]
        return fftfreq, fft

    # Generate a delta function
    npts = 5001
    signal = np.zeros(npts)
    signal[int(npts/2.)] += 1
    delta = 0.01

    fig, ax = plt.subplots(1)
    ax.set_xlim((1./(npts*delta), 1./(2*delta)))
    ax.set_ylim((1e-7, 1e-1))
    ax.grid(True, which='both')
    ax.set_xlabel('Frequency (Hz)')
    for freq in (0.1, 1, 10, 30, 40):
        C_HP, C_LP = rec_filter_coeff(freq, delta)
        signal_filt_BP_2 = recursive_filter(signal, C_HP, C_LP, npoles=2)
        signal_filt_BP_4 = recursive_filter(signal, C_HP, C_LP, npoles=4)
        signal_filt_BP_6 = recursive_filter(signal, C_HP, C_LP, npoles=6)

        ax.axvline(freq, color='gray', lw=2)
        fftfreq, fft = do_fft(signal_filt_BP_2)
        h1, = ax.loglog(fftfreq, np.abs(fft), color='red', lw=2)
        fftfreq, fft = do_fft(signal_filt_BP_4)
        h2, = ax.loglog(fftfreq, np.abs(fft), color= 'blue', lw=2)
        fftfreq, fft = do_fft(signal_filt_BP_6)
        h3, = ax.loglog(fftfreq, np.abs(fft), color= 'green', lw=2)
    ax.legend([h1, h2, h3], ['2-poles', '4-poles', '6-poles'])
    plt.show()
