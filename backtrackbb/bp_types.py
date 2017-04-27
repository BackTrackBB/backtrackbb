# bp_types.py
# Data types for BackTrackBB
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from obspy import UTCDateTime
from ctypes import c_double
from itertools import product
from scipy.stats import trim_mean


def _time_average(times, trimfraction=None):
    times = np.array(list(times))
    if not len(times):
        raise ValueError('times array is empty')
    ref_time = times[0]
    times -= ref_time
    if trimfraction is not None:
        mean = trim_mean(times, trimfraction)
    else:
        mean = np.mean(times)
    return ref_time + mean


class Trigger():
    def __init__(self,
                 eventid=None, list_picks=None,
                 x=None, y=None, z=None,
                 i=None, j=None, k=None,
                 max_grid=None,
                 ntraces=None,
                 beg_win=None, end_win=None,
                 center_win=None):
        self.eventid = None
        self.picks = []
        self.x = x
        self.y = y
        self.z = z
        self.i = i
        self.j = j
        self.k = k
        self.max_grid = max_grid
        self.ntraces = ntraces
        self.beg_win = beg_win
        self.end_win = end_win
        self.center_win = center_win
        self.trigger_level = None
        self.lat = None
        self.lon = None
        self.origin_time = None
        self.valid = True

    def __str__(self):
        s = '%s ' % self.eventid
        s += 'X %s ' % self.x
        s += 'Y %s ' % self.y
        s += 'Z %s ' % self.z
        s += 'MaxStack '
        fmt = '%.1e ' if self.max_grid < 0.01 else '%.3f '
        s += fmt % (self.max_grid)
        s += 'Ntraces %s ' % self.ntraces
        s += 'BEG %s ' % self.beg_win
        s += 'END %s ' % self.end_win
        if (self.lat is not None and
                self.lon is not None):
            s += ' LAT %.5f LON %.5f' % (self.lat, self.lon)
        if (self.origin_time is not None):
            s += ' T_ORIG %s' % (self.origin_time)
        return s

    def from_str(self, string):
        word = string.split()
        # sanity check
        try:
            if word[1] != 'X' or word[19] != 'T_ORIG':
                raise ValueError('Not a valid trigger string')
        except IndexError:
            raise ValueError('Not a valid trigger string')
        self.eventid = word[0]
        self.x = float(word[2])
        self.y = float(word[4])
        self.z = float(word[6])
        self.max_grid = float(word[8])
        self.ntraces = int(word[10])
        self.beg_win = float(word[12])
        self.end_win = float(word[14])
        self.lat = float(word[16])
        self.lon = float(word[18])
        self.origin_time = UTCDateTime(word[20])

    def add_pick(self, pick):
        self.picks.append(pick)

    def make_picks(self, stations, arrival_types,
                   arrival_times=None, grids=None):
        for sta, arr in product(stations, arrival_types):
            pick = Pick()
            pick.station = sta
            pick.arrival_type = arr
            # pick is validated by Pick.from_arrival_times()
            pick.valid = False
            if arrival_times is not None:
                pick.from_arrival_times(arrival_times)
            if grids is not None:
                pick.theor_time =\
                    grids[sta][arr].get_value(self.x, self.y, self.z)
            self.add_pick(pick)

    def get_picks(self, station=None, arrival_type=None):
        return [pick for pick in self.picks
                if (station is not None and pick.station == station)
                or (arrival_type is not None
                    and pick.arrival_type == arrival_type)]

    def compute_origin_time(self, dt_min, update_evid=False):
        # initial estimation of origin time
        # avoid utliers by using a 20% trimmed mean
        otime = _time_average((_pick.pick_time_absolute - _pick.theor_time
                              for _pick in self.picks if _pick.valid),
                              trimfraction=0.2)
        # invalidate picks too distant from theoretical time
        for pick in self.picks:
            pick.theor_time_absolute = otime + pick.theor_time
            if pick.pick_time_absolute is not None:
                pick.time_dev =\
                    abs(pick.pick_time_absolute - pick.theor_time_absolute)
                if pick.time_dev >= dt_min:
                    pick.valid = False
        # recalculate origin time by disregarding new invalid picks
        try:
            otime = _time_average(_pick.pick_time_absolute - _pick.theor_time
                                  for _pick in self.picks if _pick.valid)
            # recompute theoretical time with updated origin_time
            for pick in self.picks:
                pick.theor_time_absolute = otime + pick.theor_time
        except ValueError:
            pass
        # compute pick_time and theor_time respect to origin_time
        for pick in self.picks:
            if pick.pick_time_absolute is not None:
                pick.pick_time = pick.pick_time_absolute - otime
            else:
                pick.pick_time = -10.0
            pick.theor_time = pick.theor_time_absolute - otime
        self.origin_time = otime

    def set_eventid(self, eventid=None):
        if eventid is None:
            try:
                eventid = self.origin_time.strftime('%Y%m%d_%H%M') + 'A'
            except Exception:
                raise ValueError('unable to set eventid form origin time')
        self.eventid = eventid
        for pick in self.picks:
            pick.eventid = self.eventid

    def check_validity(self):
        if self.origin_time is None:
            self.valid = False
            return
        # check if there is at least a valid pick
        for pick in self.picks:
            if pick.valid:
                self.valid = True
                return
        self.valid = False


class Pick():
    def __init__(self, eventid=None, station=None, arrival_type=None):
        self.eventid = eventid
        self.station = station
        self.arrival_type = arrival_type
        self.theor_time = None
        self.theor_time_absolute = None
        self.pick_time = None
        self.pick_time_absolute = None
        self.time_dev = None
        self.valid = True

    def __str__(self):
        s = ' sta %s ' % self.station
        s += ' Ph %s ' % self.arrival_type
        s += ' TT %.2f ' % self.theor_time
        s += ' PT %.2f ' % self.pick_time
        return s

    def from_str(self, string):
        word = string.split()
        # sanity check
        try:
            if word[0] != 'sta' or word[6] != 'PT':
                raise ValueError('Not a valid pick string')
        except IndexError:
            raise ValueError('Not a valid pick string')
        self.station = word[1]
        self.arrival_type = word[3]
        self.theor_time = float(word[5])
        self.pick_time = float(word[7])
        #TODO: read and write these fields?
        #self.travel_time = None
        #self.time_dev = None

    def from_arrival_times(self, arrival_times):
        station = self.station
        arrival_type = self.arrival_type
        try:
            _arrival_times = arrival_times[station][arrival_type]
        except KeyError:
            self.valid = False
            return
        try:
            self.pick_time_absolute = _time_average(_arrival_times)
            self.valid = True
        except ValueError:
            self.valid = False


class RecursiveMemory():
    def __init__(self, trid=None, wave=None, band=None,
                 nsamples=0, overlap=0, filter_npoles=2):
        self.trid = trid
        self.wave = wave
        self.band = band
        self.nsamples = int(nsamples)
        self.overlap = int(overlap)
        self.filter_npoles = filter_npoles
        self.filterH = np.zeros(self.filter_npoles)
        self.filterL = np.zeros(self.filter_npoles)
        self.prev_sample_value = c_double(0)
        self.memory_sample = self.nsamples - self.overlap - 1
        self.mean_sq = c_double(0)
        self.mean = c_double(0)
        self.var = c_double(1)
        self.hos = c_double(0)
        self.initialize = True

    def __str__(self):
        s = '%s %s %s\n' % (self.trid, self.wave, self.band)
        s += 'filterH1: %f\n' % self.filterH1.value
        s += 'filterH2: %f\n' % self.filterH2.value
        s += 'filterL1: %f\n' % self.filterL1.value
        s += 'filterL2: %f\n' % self.filterL2.value
        s += 'prev_sample_value: %f\n' % self.prev_sample_value.value
        s += 'memory_sample: %d\n' % self.memory_sample
        s += 'mean_sq: %f\n' % self.mean_sq.value
        s += 'mean: %f\n' % self.mean.value
        s += 'var: %f\n' % self.var.value
        s += 'hos: %f\n' % self.hos.value
        s += 'initialize: %s\n' % self.initialize
        return s
