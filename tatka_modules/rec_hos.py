# rec_hos.py
#
# Recursive higher-order statistics.
#
# (c) 2015      - Natalia Poiata <poiata@ipgp.fr>,
#                 Claudio Satriano <satriano@ipgp.fr>
# (c) 2013-2014 - Natalia Poiata <poiata@ipgp.fr>,
#                 Claudio Satriano <satriano@ipgp.fr>,
#                 Pierre Romanet <romanet@ipgp.fr>
import os
from ctypes import CDLL, c_int, c_float, c_double, c_void_p, POINTER, byref
from numpy.ctypeslib import ndpointer
import numpy as np


libpath = os.path.join(os.path.dirname(__file__), os.pardir, 'lib', 'lib_rec_hos.so')
lib_rec_hos = CDLL(libpath)

lib_rec_hos._recursive_hos.argtypes = [
        ndpointer(dtype=np.float64), #signal
        ndpointer(dtype=np.float64), #hos_signal
        c_int, #npts
        c_float, #sigma_min
        c_float, #C_WIN
        c_int, #order
        POINTER(c_double), #mean
        POINTER(c_double), #var
        POINTER(c_double) #hos
        ]
lib_rec_hos._recursive_hos.restype = c_void_p


def recursive_hos(signal, C_WIN, sigma_min, order, rec_memory=None):
    """
        Recursive computation of higher-horder statistics.
    """

    try:
        C_WIN = float(C_WIN)
    except ValueError:
        print 'C_rosenberger should be a double'

    signal = np.array(signal, dtype=np.float64)
    hos_signal = np.zeros(len(signal))

    if rec_memory is not None:
        mean = c_double(rec_memory.mean)
        var = c_double(rec_memory.var)
        hos = c_double(rec_memory.hos)
    else:
        mean = c_double(0)
        var = c_double(1)
        hos = c_double(0)

    lib_rec_hos._recursive_hos(
            signal, hos_signal, signal.size, sigma_min, C_WIN, order,
            byref(mean), byref(var), byref(hos))

    if rec_memory is not None:
        rec_memory.mean = mean.value
        rec_memory.var = var.value
        rec_memory.hos = hos.value

    return hos_signal


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    signal = np.ones(60000)
    signal[30000] += 1
    kurt_signal = recursive_hos(signal, 0.5, 0.001, 4)
    plt.plot(kurt_signal)
    plt.show()
