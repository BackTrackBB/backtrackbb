import numpy as np
import scipy as sp
from tatka_modules.LocalCC import LocalCC


def sta_GRD_Proj(traces, ttime_GRIDS,sta_list, sta_pairs, beg, end, shift,\
                 fs_sampling, max_lag,start_time,sigma,nnx,nny,nnz,sta_pairs_id):
##   for l in xrange(len(comb_sta)):
    ind1 = sta_list.index(sta_pairs[sta_pairs_id][0])
    ind2 = sta_list.index(sta_pairs[sta_pairs_id][1])
    sig1 = traces[ind1].data[beg-shift:end+shift]/max(abs(traces[ind1].data[beg-shift:end+shift]))
    sig2 = traces[ind2].data[beg-shift:end+shift]/max(abs(traces[ind2].data[beg-shift:end+shift]))


##    sig1 = traces[ind1].data[beg:end]/max(abs(traces[ind1].data[beg:end]))
##    sig2 = traces[ind2].data[beg:end]/max(abs(traces[ind2].data[beg:end]))


    corr = LocalCC(sig1,sig2,fs_sampling,max_lag,start_time,sms=sigma)
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
    ##stationPair_projGrid = function(ttime_GRIDS[ind2]-ttime_GRIDS[ind1])
    stationPair_projGrid = function((ttime_GRIDS[ind2] -
                                     ttime_GRIDS[ind1]).flatten()).reshape(nnx, nny, nnz)
    return stationPair_projGrid
    
