import os
import ctypes
from numpy.ctypeslib import ndpointer
import numpy as np


libpath = os.path.join(os.path.dirname(__file__), os.pardir, 'lib', 'lib_rec_filter.so')
lib_rec_filter = ctypes.CDLL(libpath)

lib_rec_filter._recursive_filter_BP.argtypes = [
        ndpointer(dtype=np.float64),
        ndpointer(dtype=np.float64),
        ctypes.c_int,
        ctypes.c_float,
        ctypes.c_float
        ]
lib_rec_filter._recursive_filter_BP.restype = ctypes.c_void_p

lib_rec_filter._recursive_filter_HP.argtypes = [
        ndpointer(dtype=np.float64),
        ndpointer(dtype=np.float64),
        ctypes.c_int,
        ctypes.c_float
        ]
lib_rec_filter._recursive_filter_HP.restype = ctypes.c_void_p


def recursive_filter(signal, C_HP, C_LP=None):
    signal = np.array(signal, dtype=np.float64)
    filt_signal = np.zeros(len(signal))
    if C_LP is not None:
        lib_rec_filter._recursive_filter_BP(
            signal, filt_signal, signal.size, C_HP, C_LP)
    else:
        lib_rec_filter._recursive_filter_HP(
            signal, filt_signal, signal.size, C_HP)
    return filt_signal


if __name__ == '__main__':
    signal = np.ones(60000)
    filt_signal = recursive_filter(signal, 0.1, 0.1)
    print filt_signal
