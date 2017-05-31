# -*- coding: utf8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import math


def trig_dist(trg1, trg2):
    dist = math.sqrt((trg1.x - trg2.x)**2 + (trg1.y - trg2.y)**2)
    time_diff = abs(trg1.origin_time - trg2.origin_time)
    return dist, time_diff


def group_triggers(config, triggers):
    # Sorting triggers in the 
    sorted_trg = sorted(triggers, key=lambda x: x.max_grid,
                        reverse=True)

    for n, trg1 in enumerate(sorted_trg):
        # Currently, we compare all the triggers, 
        # because trigger list was sorted.
        #
        # For speed, number of triggers to be grouped can be reduced 
        # to a number N (ex N=50) 
        # for trg2 in sorted_trg[n+1:n+50]:
        
        for trg2 in sorted_trg[n+1:]:
            
            dist, time_diff = trig_dist(trg1, trg2)
            
            if (dist <= config.group_min_dist and
                    time_diff <= config.group_min_time_diff):
                
                if trg2.max_grid <= trg1.max_grid:
                    sorted_trg.remove(trg2)
                else:
                    sorted_trg.remove(trg1)
                    break

    return sorted(sorted_trg, key=lambda x: x.eventid)
