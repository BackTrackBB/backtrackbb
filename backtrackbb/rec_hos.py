# rec_hos.py
#
# Recursive higher-order statistics.
#
# (c) 2015-2018 - Natalia Poiata <poiata@ipgp.fr>,
#                 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013-2014 - Natalia Poiata <poiata@ipgp.fr>,
#                 Claudio Satriano <satriano@ipgp.fr>,
#                 Pierre Romanet <romanet@ipgp.fr>
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from ctypes import CDLL, c_int, c_float, c_double, c_void_p, POINTER, byref
from numpy.ctypeslib import ndpointer
try:
    from .lib_names import get_lib_path
except (ImportError, ValueError):
    from lib_names import get_lib_path


lib_rec_hos = CDLL(get_lib_path('lib_rec_hos'))
lib_rec_hos._recursive_hos.argtypes = [
        ndpointer(dtype=np.float64),  # signal
        ndpointer(dtype=np.float64),  # hos_signal
        c_int,  # npts
        c_float,  # sigma_min
        c_float,  # C_WIN
        c_int,  # order
        POINTER(c_double),  # mean
        POINTER(c_double),  # var
        POINTER(c_double),  # hos
        c_int,  # memory_sample
        c_int  # initialize
        ]
lib_rec_hos._recursive_hos.restype = c_void_p


def recursive_hos(signal, C_WIN, order=4, sigma_min=-1., rec_memory=None):
    """Recursive computation of higher-horder statistics."""
    try:
        C_WIN = float(C_WIN)
    except ValueError:
        print('C_WIN should be a double')

    signal = np.array(signal, dtype=np.float64)
    hos_signal = np.zeros(len(signal))

    if rec_memory is not None:
        mean = rec_memory.mean
        var = rec_memory.var
        hos = rec_memory.hos
        memory_sample = rec_memory.memory_sample
        initialize = int(rec_memory.initialize)
    else:
        mean = c_double(0)
        var = c_double(1)
        hos = c_double(0)
        memory_sample = -1
        initialize = 1

    lib_rec_hos._recursive_hos(signal, hos_signal, signal.size,
                               sigma_min, C_WIN, order,
                               byref(mean), byref(var), byref(hos),
                               memory_sample, initialize)

    if rec_memory is not None:
        rec_memory.mean = mean
        rec_memory.var = var
        rec_memory.hos = hos
        rec_memory.initialize = False

    return hos_signal


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    signal = np.ones(601)
    signal[300] += 1
    kurt_signal = recursive_hos(signal, 0.5, sigma_min=0.001, order=4)
    plt.plot(kurt_signal)
    plt.show()
