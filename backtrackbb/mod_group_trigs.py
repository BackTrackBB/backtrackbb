# -*- coding: utf8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import math


def trig_dist(trg1, trg2):
    dist = math.sqrt((trg1.x - trg2.x)**2 + (trg1.y - trg2.y)**2)
    time_diff = abs(trg1.origin_time - trg2.origin_time)
    return dist, time_diff


def group_triggers(config, triggers):

    for n, trg1 in enumerate(triggers):
        # we assume (for speed) than no more than
        # 10 triggers need to be grouped
        for trg2 in triggers[n+1:n+10]:
            dist, time_diff = trig_dist(trg1, trg2)
            if (dist <= config.group_min_dist and
                    time_diff <= config.group_min_time_diff):
                if trg2.max_grid <= trg1.max_grid:
                    triggers.remove(trg2)
                else:
                    triggers.remove(trg1)
                    break

    return triggers
