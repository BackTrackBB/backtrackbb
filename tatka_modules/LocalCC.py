import numpy as np
import scipy as sp
#from recursive_cc import local_CCr          #Python version
from rec_cc import local_CCr                 #C version

class LocalCC():
    '''
    '''
    len_in_signal = int(0)
    n_lags = int(0)
    time_in_signal = None
    cc_time_lags = None
    smoothed_cc = None
    cc = None
    def __init__(self, in_signal1=None, in_signal2=None, samp_rate=None,
                 max_time_lag=None, zero_time=None, sigma=None):
        self.in_signal1 = in_signal1
        self.in_signal2 = in_signal2
        self.samp_rate = samp_rate
        self.max_time_lag = max_time_lag
        self.zero_time = zero_time
        self.sigma = sigma
        self.define_param()
        self.local_cc()
        

    def define_param(self):
        self.len_in_signal = len(self.in_signal1)
        self.time_in_signal = np.arange(self.len_in_signal) / self.samp_rate
        self.n_lags = 2*int(self.max_time_lag * self.samp_rate)
        self.cc_time_lags = sp.linspace(-self.n_lags/2/self.samp_rate,
                                         self.n_lags/2/self.samp_rate,
                                         self.n_lags, endpoint=True)

    def local_cc(self):
        self.smoothed_cc, self.cc  = local_CCr(self.in_signal1, self.in_signal2,
                                               self.max_time_lag,
                                               self.samp_rate, self.sigma)
        lag, time = [int(n) for n in np.where(self.smoothed_cc == self.smoothed_cc.max())]

        time = self.zero_time + time/self.samp_rate
        lag = lag/self.samp_rate - self.max_time_lag

        self.arrival1 = time - lag/2.
        self.arrival2 = time + lag/2.
