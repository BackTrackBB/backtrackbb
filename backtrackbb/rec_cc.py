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
lib_rec_cc._local_CCr.argtypes = [
        ndpointer(dtype=np.float64),  # signal1
        ndpointer(dtype=np.float64),  # signal2
        c_int,                        # npts
        ndpointer(dtype=np.float64),  # cc
        c_int,                        # lmax
        c_double                      # sigma
        ]
lib_rec_cc._local_CCr.restype = c_void_p


def local_CCr(signal1, signal2, t_lag, fs, sigma=None):
    if not isinstance(signal1, np.ndarray):
        signal1 = np.array(signal1, dtype=np.float64)
    if not isinstance(signal2, np.ndarray):
        signal2 = np.array(signal2, dtype=np.float64)
    if signal1.size != signal2.size:
        raise RuntimeError('Signals must have the same size.')

    if sigma is not None:
        # Compute sigma in samples
        sigma *= fs
        if sigma < 0.5:
            raise RuntimeError(
                'Sigma for Gaussian filter must be >= 0.5 samples.')
    else:
        # Do not smooth
        sigma = -999.
    lmax = int(t_lag * fs)
    cc = np.zeros(2*lmax * len(signal1))
    lib_rec_cc._local_CCr(
            signal1, signal2, signal1.size,
            cc, lmax, sigma)
    return cc.reshape(2*lmax, len(signal1))


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    signal1 = np.zeros(6001)
    signal2 = np.zeros(6001)
    signal1[1000:1200] = 1
    signal2[1200:1300] = 1
    cc = local_CCr(signal1, signal2, 300, 1, 100)
    cc_no_filt = local_CCr(signal1, signal2, 300, 1, None)
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.imshow(cc_no_filt)
    ax2 = fig.add_subplot(212, sharex=ax1, sharey=ax1)
    ax2.imshow(cc)
    plt.show()
