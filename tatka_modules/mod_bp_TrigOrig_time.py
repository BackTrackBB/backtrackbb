import numpy as np
from tatka_modules.NLLGrid import NLLGrid

def TrOrig_time(config,stations,GRD_sta,xx_trig, yy_trig, zz_trig,
                rec_start_time,arrival_times,trig_time):
    ##-------------------------------------        
    dt_min = config.dt_min
    time_dev = []
    ##-------------------------------------
    orig_time = 0

    for sta in stations:
        tr_time = 0
        for pick_times in arrival_times[sta]:
            tr_time += pick_times-rec_start_time+config.cut_start

        trig_time[sta].append(tr_time/len(arrival_times[sta]))
        
        tt_time = GRD_sta[sta].get_value(xx_trig,yy_trig,zz_trig)
        trig_time[sta].append(tt_time)
        
        orig_time += trig_time[sta][0] - trig_time[sta][1]
    bp_origin_time = rec_start_time+orig_time/len(stations)-config.cut_start

    for sta in stations:
        trig_time[sta].append(orig_time/len(stations)+trig_time[sta][1])
        trig_time[sta].append(abs(trig_time[sta][0]-trig_time[sta][2]))
        time_dev.append(abs(trig_time[sta][0]-trig_time[sta][2]))

    if len(filter(lambda x: x>= dt_min,time_dev))> 1:
        bp_origin_time = 0
        tt_time = 0
        orig_time = 0
        mm=0        
        for sta in stations:
            tr_time = 0
            if trig_time[sta][3] <= dt_min:
                mm += 1
                orig_time += trig_time[sta][0] - trig_time[sta][1]
        bp_origin_time = rec_start_time + orig_time/mm - config.cut_start
        
        for sta in stations:
            trig_time[sta][2] = 0.
            trig_time[sta][2] = orig_time/mm+trig_time[sta][1]

    return bp_origin_time, trig_time
