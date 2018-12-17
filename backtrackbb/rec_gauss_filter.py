# -*- coding: utf8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from ctypes import CDLL, c_int, c_double, c_void_p
from numpy.ctypeslib import ndpointer
try:
    from .lib_names import get_lib_path
except (ImportError, ValueError):
    from lib_names import get_lib_path


lib_rec_cc = CDLL(get_lib_path('lib_rec_cc'))
lib_rec_cc._Gaussian1D.argtypes = [
        ndpointer(dtype=np.float64),  # signal
        c_int,                        # npts
        c_double                      # sigma
        ]
lib_rec_cc._Gaussian1D.restype = c_void_p


def recursive_gauss_filter(signal, sigma):
    filt_signal = np.array(signal, dtype=np.float64)
    lib_rec_cc._Gaussian1D(filt_signal, filt_signal.size, sigma)
    return filt_signal


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    signal = np.zeros(1001)
    signal[500] = 1
    filt_signal = recursive_gauss_filter(signal, 10)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(signal)
    ax.plot(filt_signal/max(filt_signal))
    plt.show()
