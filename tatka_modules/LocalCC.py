import numpy as np
import scipy as sp
#from recursive_cc import local_CCr          #Python vertion
from rec_cc import local_CCr                 #C vertion

class LocalCC():
    '''
    '''
    len_in_signal = int(0)
    n_lags = int(0)
    time_in_signal = None
    cc_time_lags = None
    smoothed_cc = None
    cc = None
    def __init__(self,in_signal1=None,in_signal2=None,samp_rate=None
                 ,max_time_lag=None,zero_time=None,sms=None):
        self.in_signal1 = in_signal1
        self.in_signal2 = in_signal2
        self.samp_rate = samp_rate
        self.max_time_lag = max_time_lag
        self.zero_time = zero_time
        self.sms = sms
        self.define_param(in_signal1,in_signal2,samp_rate,max_time_lag)
        self.local_cc(in_signal1,in_signal2,samp_rate,max_time_lag,sms)
        

    def define_param(self,in_signal1,in_signal2,samp_rate,max_time_lag):
        self.len_in_signal = len(in_signal1)
        self.time_in_signal = np.arange(self.len_in_signal)/ samp_rate
        self.n_lags = 2*int(max_time_lag*samp_rate)
        self.cc_time_lags = sp.linspace(-self.n_lags/2/samp_rate,
                                        self.n_lags/2/samp_rate,
                                    self.n_lags, endpoint=True)

    def local_cc(self,in_signal1,in_signal2,samp_rate,max_time_lag,sms):
        self.smoothed_cc, self.cc  = local_CCr(in_signal1,in_signal2,
                                               max_time_lag,samp_rate,sms)
