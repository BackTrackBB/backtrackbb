import numpy as np
import scipy as sp
from tatka_modules.LocalCC import LocalCC


def sta_GRD_Proj(stream, ttime_GRIDS, sta1, sta2, beg, end, shift,
                 fs_sampling, max_lag, start_time, sigma,
                 nnx, nny, nnz):

    trace1 = stream.select(station=sta1)[0]
    trace2 = stream.select(station=sta2)[0]
    sig1 = trace1.data[beg-shift:end+shift]/max(abs(trace1.data[beg-shift:end+shift]))
    sig2 = trace2.data[beg-shift:end+shift]/max(abs(trace2.data[beg-shift:end+shift]))

    corr = LocalCC(sig1, sig2, fs_sampling, max_lag, start_time, sms=sigma)
    a = corr.smoothed_cc
    t_lag = corr.cc_time_lags

    ## Max value of local_cc in given window
    local_cc = np.amax(a[:,shift+1:shift+1+end],axis=1)
##    local_cc = np.amax(a,axis=1)

    ## Projecting LCC for the station pair on the grid of theoretical t_times
    ## Check which of the functions for the interpolations is faster?
    ##function = sp.interpolate.interp1d(t_lag,local_cc)
    function = sp.interpolate.UnivariateSpline(t_lag,local_cc,k=1,s=0)

    ## Output_Grid
    tt_array1 = ttime_GRIDS[sta1].array
    tt_array2 = ttime_GRIDS[sta2].array
    stationPair_projGrid = function((tt_array2 - tt_array1).flatten()).reshape(nnx, nny, nnz)

    return stationPair_projGrid
    
