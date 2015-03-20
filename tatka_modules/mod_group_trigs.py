import math

def trig_dist(trg1, trg2):
    dist = math.sqrt((trg1.x - trg2.x)**2 + (trg1.y - trg2.y)**2)
    time_shift = abs(trg1.origin_time - trg2.origin_time)
    return dist, time_shift

def group_triggers(config, triggers):
    sorted_trg = sorted(triggers, key=lambda x: x.max_grid,
                        reverse=True)
    n = 0
    while True:
        for trg in sorted_trg[n+1:]:
            dist, time_shift = trig_dist(sorted_trg[n], trg)
            if (dist <= config.max_dist and
                time_shift <= config.max_time_shift):
                sorted_trg.remove(trg)
        n += 1
        if n == len(sorted_trg):
            break

    return sorted(sorted_trg, key=lambda x: x.beg_win)
