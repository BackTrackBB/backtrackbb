import os
import ctypes
from numpy.ctypeslib import ndpointer
import numpy as np


libpath = os.path.join(os.path.dirname(__file__), os.pardir, 'lib', 'lib_rec_kurtosis.so')
lib_rec_kurtosis = ctypes.CDLL(libpath)

lib_rec_kurtosis._recursive_kurtosis.argtypes = [
        ndpointer(ctypes.c_double),
        ndpointer(ctypes.c_double),
        ctypes.c_int,
        ctypes.c_float,
        ]
lib_rec_kurtosis._recursive_kurtosis.restype = ctypes.c_void_p

def recursive_kurtosis(signal, C_WIN):
    kurt_signal = np.zeros(len(signal))
    lib_rec_kurtosis._recursive_kurtosis(
            signal, kurt_signal, signal.size, C_WIN)
    return kurt_signal

if __name__ == '__main__':
    signal = np.ones(60000)
    kurt_signal = recursive_kurtosis(signal, 0.1)
    print kurt_signal
