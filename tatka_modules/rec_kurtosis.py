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
        ctypes.c_float,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_float
        ]
lib_rec_kurtosis._recursive_kurtosis.restype = ctypes.c_void_p

def recursive_kurtosis(signal, C_WIN, sigma_min, order1, order2, power2):

    kurt_signal = np.zeros(len(signal))
    lib_rec_kurtosis._recursive_kurtosis(
            signal, kurt_signal, signal.size, sigma_min, C_WIN, order1, order2, power2)
    return kurt_signal

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    signal = np.ones(60000)
    signal[30000]+=1
    kurt_signal = recursive_kurtosis(signal, 0.5, 0.001, 4, 2, 2)
    plt.plot(kurt_signal)
    plt.show()
