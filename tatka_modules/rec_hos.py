# rec_hos.py
#
# Recursive higher-order statistics.
#
# (c) 2013-2014 - Natalia Poiata <poiata@ipgp.fr>,
#                 Claudio Satriano <satriano@ipgp.fr>,
#                 Pierre Romanet <romanet@ipgp.fr>
import os
import ctypes
from numpy.ctypeslib import ndpointer
import numpy as np


libpath = os.path.join(os.path.dirname(__file__), os.pardir, 'lib', 'lib_rec_hos.so')
lib_rec_hos = ctypes.CDLL(libpath)

lib_rec_hos._recursive_hos.argtypes = [
        ndpointer(ctypes.c_double),
        ndpointer(ctypes.c_double),
        ctypes.c_int,
        ctypes.c_float,
        ctypes.c_float,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_float
        ]
lib_rec_hos._recursive_hos.restype = ctypes.c_void_p


def recursive_hos(signal, C_WIN, sigma_min, order1, order2, power2):
    """
        Recursive computation of higher-horder statistics.
    """

    try:
        C_WIN = float(C_WIN)
    except ValueError:
        print 'C_rosenberger should be a double'

    hos_signal = np.zeros(len(signal))
    lib_rec_hos._recursive_hos(
            signal, hos_signal, signal.size, sigma_min, C_WIN, order1, order2, power2)
    return hos_signal


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    signal = np.ones(60000)
    signal[30000]+=1
    kurt_signal = recursive_hos(signal, 0.5, 0.001, 4, 2, 2)
    plt.plot(kurt_signal)
    plt.show()
