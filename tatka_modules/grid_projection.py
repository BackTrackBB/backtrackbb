import numpy as np
import scipy as sp
from LocalCC import LocalCC


def sta_GRD_Proj(config, stream, ttime_GRIDS, sta1, sta2, wave1, wave2, t_b, t_e, tau_max=None):

    max_lag = config.time_lag
    if tau_max is not None:
        max_lag = tau_max

    fs_sampling = stream[0].stats.sampling_rate
    start_time = stream[0].stats.starttime

    beg = int(t_b * fs_sampling)
    end = int(t_e * fs_sampling)

    trace1 = stream.select(station=sta1, channel=wave1)[0]
    trace2 = stream.select(station=sta2, channel=wave2)[0]
    sig1 = trace1.data[beg:end]/max(abs(trace1.data[beg:end]))
    sig2 = trace2.data[beg:end]/max(abs(trace2.data[beg:end]))

    t_lag, local_cc, arrival1, arrival2 =\
        LocalCC(sig1, sig2, fs_sampling, max_lag, start_time+t_b, config)

    ## Max value of local_cc in given window
    local_cc_1d = np.amax(local_cc, axis=1)

    ## Projecting LCC for the station pair on the grid of theoretical t_times
    ## Check which of the functions for the interpolations is faster?
    ##function = sp.interpolate.interp1d(t_lag, local_cc_1d)
    function = sp.interpolate.UnivariateSpline(t_lag, local_cc_1d, k=1, s=0)

    ## Output_Grid
    tt_array1 = ttime_GRIDS[sta1][wave1].array
    tt_array2 = ttime_GRIDS[sta2][wave2].array
    stationPair_projGrid = function((tt_array2 - tt_array1).flatten()).reshape(tt_array1.shape)

    return stationPair_projGrid, arrival1, arrival2
