# -*- coding: utf8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from ctypes import CDLL, c_int, c_float, c_double, c_void_p, POINTER, byref
from numpy.ctypeslib import ndpointer
try:
    from .lib_names import get_lib_path
except (ImportError, ValueError):
    from lib_names import get_lib_path


lib_rec_filter = CDLL(get_lib_path('lib_rec_filter'))
lib_rec_filter._recursive_filter.argtypes = [
        ndpointer(dtype=np.float64),  # signal
        ndpointer(dtype=np.float64),  # filt_signal
        c_int,  # npts
        c_float,  # C_HP
        c_float,  # C_LP
        c_int,  # npoles
        ndpointer(dtype=np.float64),  # filterH
        ndpointer(dtype=np.float64),  # filterL
        POINTER(c_double),  # prev_sample_value
        c_int  # memory_sample
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


def rec_filter_coeff(freqs, delta):
    freqs = np.array(freqs, ndmin=1)
    nyq = 1. / (2*delta)
    T = 1. / freqs
    w = T / (2*np.pi)
    rel_freqs = freqs / nyq
    # empirical approach to fix filter shape close to the Nyquist:
    w[rel_freqs >= 0.2] /= (rel_freqs[rel_freqs >= 0.2] * 7)
    C_HP = w / (w + delta)      # high-pass filter constant
    C_LP = delta / (w + delta)  # low-pass filter constant
    return C_HP, C_LP


def rec_filter_norm(freqs, delta, C_HP, C_LP, npoles):
    """Empirical approach for computing filter normalization coefficients."""
    freqs = np.array(freqs, ndmin=1)
    norm = np.zeros(len(freqs))
    for n, freq in enumerate(freqs):
        length = 4. / freq
        time = np.arange(0, length+delta, delta)
        signal = np.sin(freq * 2 * np.pi * time)
        signal_filt = recursive_filter(signal, C_HP[n], C_LP[n], npoles)
        norm[n] = signal_filt.max()
    return norm


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    def do_fft(signal):
        npts = len(signal)
        if not npts % 2:
            npts -= 1
        fft = np.fft.rfft(signal, n=npts)
        fftfreq = np.fft.fftfreq(len(signal), d=delta)
        fftfreq = fftfreq[0:fft.size]
        return fftfreq, fft

    # Generate a delta function
    npts = 5001
    signal = np.zeros(npts)
    signal[int(npts/2.)] += 1
    delta = 0.01

    freqs = (0.1, 1, 10, 20, 40)
    C_HP, C_LP = rec_filter_coeff(freqs, delta)
    norm_2 = rec_filter_norm(freqs, delta, C_HP, C_LP, npoles=2)
    norm_4 = rec_filter_norm(freqs, delta, C_HP, C_LP, npoles=4)
    norm_6 = rec_filter_norm(freqs, delta, C_HP, C_LP, npoles=6)

    fig, ax = plt.subplots(1)
    ax.set_xlim((1./(npts*delta), 1./(2*delta)))
    ax.set_ylim((1e-6, 100))
    ax.grid(True, which='both')
    ax.set_xlabel('Frequency (Hz)')
    for n, freq in enumerate(freqs):
        signal_filt_BP_2 = recursive_filter(signal, C_HP[n], C_LP[n], npoles=2)
        signal_filt_BP_2 /= norm_2[n]
        signal_filt_BP_4 = recursive_filter(signal, C_HP[n], C_LP[n], npoles=4)
        signal_filt_BP_4 /= norm_4[n]
        signal_filt_BP_6 = recursive_filter(signal, C_HP[n], C_LP[n], npoles=6)
        signal_filt_BP_6 /= norm_6[n]

        ax.axvline(freq, color='gray', lw=2)
        fftfreq, fft = do_fft(signal_filt_BP_2)
        h1, = ax.loglog(fftfreq, np.abs(fft), color='red', lw=2)
        fftfreq, fft = do_fft(signal_filt_BP_4)
        h2, = ax.loglog(fftfreq, np.abs(fft), color='blue', lw=2)
        fftfreq, fft = do_fft(signal_filt_BP_6)
        h3, = ax.loglog(fftfreq, np.abs(fft), color='green', lw=2)
    ax.legend([h1, h2, h3], ['2-poles', '4-poles', '6-poles'])
    plt.show()
