from collections import defaultdict

def TrOrig_time(config, stations, GRD_sta, xx_trig, yy_trig, zz_trig,
                rec_start_time, arrival_times):
    ##-------------------------------------
    dt_min = config.dt_min
    time_dev = []
    ##-------------------------------------
    orig_time = 0

    # mean picked time relative to the 'rec_start_time' for each station - trig_time[sta][0]
    trig_time = defaultdict(list)

    wave = 'P'
    for sta in stations:
        tr_time = 0
        for pick_times in arrival_times[sta][wave]:
            tr_time += pick_times - rec_start_time + config.cut_start

        trig_time[sta].append(tr_time/len(arrival_times[sta][wave]))

        # theoretical travel times for each station - trig_time[sta][1]
        tt_time = GRD_sta[sta][wave].get_value(xx_trig, yy_trig, zz_trig)
        trig_time[sta].append(tt_time)

        # initial estimation of relative & absolute origin times
        orig_time += trig_time[sta][0] - trig_time[sta][1]

    bp_origin_time = rec_start_time + orig_time/len(stations) - config.cut_start

    for sta in stations:
        # theoretical arrival time relative to 'rec_start_time' - trig_time[sta][2]
        trig_time[sta].append(orig_time/len(stations) + trig_time[sta][1])

        # difference btw theoretical and picked arrival times - trig_time[sta][3]
        trig_time[sta].append(abs(trig_time[sta][0] - trig_time[sta][2]))
        time_dev.append(abs(trig_time[sta][0] - trig_time[sta][2]))

        # picked arrival time relative to origin time - trig_time[sta][4]
        trig_time[sta].append(rec_start_time - bp_origin_time + trig_time[sta][0]- config.cut_start)

    # loop for recalculating the origin time by disregardin picks with time_dev > dt_min
    if len(filter(lambda x: x >= dt_min, time_dev)) > 1:
        bp_origin_time = 0
        tt_time = 0
        orig_time = 0
        mm = 0
        for sta in stations:
            tr_time = 0
            if trig_time[sta][3] <= dt_min:
                mm += 1
                orig_time += trig_time[sta][0] - trig_time[sta][1]
        if mm!=0:
            bp_origin_time = rec_start_time + orig_time/mm - config.cut_start

            for sta in stations:
                trig_time[sta][2] = 0.
                trig_time[sta][2] = orig_time/mm + trig_time[sta][1]
        else:
            bp_origin_time = None
            trig_time = defaultdict(list)

### different formulation for providing picked arrival time and its standard deviation
### need checks to see if can be used
##    if bp_origin_time:
##        for sta in stations:
##            tr_time = []
##            for pick_times in arrival_times[sta]:
##                tr = st.select(station=sta)[0]
##                tr_time.append(pick_times - bp_origin_time)
####            print sta,np.mean(tr_time,axis=0), np.std(tr_time,axis=0)

    return bp_origin_time, trig_time
