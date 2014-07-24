import numpy as np
from obspy.signal.util import utlGeoKm

def dist_to_evt(trg0,trg):
    dist_to_evt0 = np.sqrt(np.power((trg0.x-trg.x),2)+\
                   np.power((trg0.y-trg.y),2)+\
                   np.power((trg0.z-trg.z),2))
    time_shift_evt0 = abs(trg0.origin_time-trg.origin_time)
    return dist_to_evt0,time_shift_evt0

def groupe_triggers(triggers,config):
    sorted_trg = sorted(triggers, key=lambda x: x.max_grid,
                          reverse=True)
    n=0
    while True:
        for trg in sorted_trg[n+1:]:
            dist,time_shift = dist_to_evt(sorted_trg[n],trg)
            if dist <= config.max_dist and time_shift <= config.max_time_shift:
                sorted_trg.remove(trg)
        n+=1
        if n == len(sorted_trg):
            break
    return sorted(sorted_trg, key=lambda x: x.beg_win)


