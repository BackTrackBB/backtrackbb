# -*-coding:UTF-8 -*
"""Generate different signals with gaussian noise."""

from math import ceil
from math import exp
from math import sin
import numpy as np


#Definition d'une gausienne
def gaussian(x, mu, sigma):
    return exp(-0.5*((x-mu)/sigma)**2)


#Definition d'une sin(x)exp(-x)
def function_sinExp(x, decreasingFactor, pulse):
    return exp(-decreasingFactor*x)*sin(pulse*x)


# Definition du signal gaussien
def generate_signalG(position, sigma, noise):
    signalG = noise.copy()
    for i in range(int(ceil(position-3*sigma)), int(ceil(position+3*sigma))):
        signalG[i-1] = signalG[i-1] + gaussian(i, position, sigma)
    return signalG


# signal dirac
def generate_signalD(position, noise):
    signalD = noise.copy()
    signalD[position-1] = signalD[position-1]+1
    return signalD


# signal cacher dans le coda
def generate_signal_expSin(position1, decreasingFactor1, pulse1, noise,
                           amplitude2=1, position2=0, decreasingFactor2=0,
                           pulse2=0):
    signal = noise.copy()
    for i in range(position1, len(noise)-1):
        signal[i-1] = signal[i-1] + \
                      function_sinExp(i-position1, decreasingFactor1, pulse1)
    for i in range(position2, len(noise)-1):
        signal[i-1] = signal[i-1] + \
                      amplitude2*function_sinExp(i-position2,
                                                 decreasingFactor2, pulse2)
    return signal


# Noise (unniform distribution)
def generate_signal_noise(lengthSignal):
    noise = np.random.random_sample((lengthSignal,))-0.5
    return noise


# Noise (Gaussian distribution)
def generate_signal_noise2(lengthSignal, sigmaNoise):
    noise = np.random.randn(lengthSignal,) * sigmaNoise
    return noise


# it does not work like other function,indeed it has no noise
def generate_signal_triangle(position, slope, maxValue, lengthSignal,
                             sigmaNoise=0):
    signal = np.random.randn(lengthSignal,) * sigmaNoise
    for i in range(lengthSignal):
        if i > (position-maxValue/slope) and i < position:
            signal[i] = signal[i]+(slope*(i-position)+maxValue)
        elif i >= position and i < (position+maxValue/slope):
            signal[i] = signal[i]+(-slope*(i-position)+maxValue)
        else:
            signal[i] = signal[i]
    return signal
