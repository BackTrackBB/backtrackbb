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


lib_rec_rms = CDLL(get_lib_path('lib_rec_rms'))
lib_rec_rms._recursive_rms.argtypes = [
        ndpointer(dtype=np.float64),  # signal
        ndpointer(dtype=np.float64),  # rms_signal
        c_int,  # npts
        c_float,  # C_WIN
        POINTER(c_double),  # mean_sq
        c_int,  # memory_sample
        c_int  # initialize
        ]
lib_rec_rms._recursive_rms.restype = c_void_p


def recursive_rms(signal, C_WIN, rec_memory=None):
    signal = np.array(signal, dtype=np.float64)
    rms_signal = np.zeros(len(signal))

    if rec_memory is not None:
        mean_sq = rec_memory.mean_sq
        memory_sample = rec_memory.memory_sample
        initialize = int(rec_memory.initialize)
    else:
        mean_sq = c_double(0)
        memory_sample = -1
        initialize = 1

    lib_rec_rms._recursive_rms(
            signal, rms_signal, signal.size, C_WIN,
            byref(mean_sq), memory_sample, initialize)

    if rec_memory is not None:
        rec_memory.mean_sq = mean_sq
        rec_memory.initialize = False

    return rms_signal


if __name__ == '__main__':
    signal = np.ones(60000)
    rms_signal = recursive_rms(signal, 1)
    print(rms_signal)
