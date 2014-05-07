import numpy as np
import scipy as sp
from tatka_modules.LocalCC import LocalCC


def sta_GRD_Proj(stream, ttime_GRIDS, sta1, sta2, t_b, t_e, shift,
                 fs_sampling, start_time, config,
                 nnx, nny, nnz, arrival_times,tau_max=0.):

    max_lag = config.time_lag
    if tau_max!=0.:
        max_lag = tau_max
##        print 'max_lag', max_lag
    sigma = config.smooth_lcc

    beg = int(t_b * fs_sampling)
    end = int(t_e * fs_sampling)

    trace1 = stream.select(station=sta1)[0]
    trace2 = stream.select(station=sta2)[0]
    sig1 = trace1.data[beg-shift:end+shift]/max(abs(trace1.data[beg-shift:end+shift]))
    sig2 = trace2.data[beg-shift:end+shift]/max(abs(trace2.data[beg-shift:end+shift]))

    corr = LocalCC(sig1, sig2, fs_sampling, max_lag, start_time+t_b, sigma)

    if config.do_smooth_lcc:
        a = corr.smoothed_cc
    else:
        a = corr.cc
    
    t_lag = corr.cc_time_lags

    arrival_times[sta1].append(corr.arrival1)
    arrival_times[sta2].append(corr.arrival2)

    ## Max value of local_cc in given window
    local_cc = np.amax(a[:,shift+1:shift+1+end], axis=1)
##    local_cc = np.amax(a,axis=1)

    ## Projecting LCC for the station pair on the grid of theoretical t_times
    ## Check which of the functions for the interpolations is faster?
    ##function = sp.interpolate.interp1d(t_lag,local_cc)
    function = sp.interpolate.UnivariateSpline(t_lag,local_cc, k=1, s=0)

    ## Output_Grid
    tt_array1 = ttime_GRIDS[sta1].array
    tt_array2 = ttime_GRIDS[sta2].array
    stationPair_projGrid = function((tt_array2 - tt_array1).flatten()).reshape(nnx, nny, nnz)

    return stationPair_projGrid
    
