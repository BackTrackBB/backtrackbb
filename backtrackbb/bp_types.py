# bp_types.py
# Data types for BackTrackBB
from obspy import UTCDateTime
from ctypes import c_double
import numpy as np


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

    def get_picks(self, station=None, arrival_type=None):
        return [pick for pick in self.picks
                if (station is not None and pick.station == station)
                or (arrival_type is not None
                    and pick.arrival_type == arrival_type)]


class Pick():
    def __init__(self,
                 eventid=None, station=None,
                 arrival_type=None):
        self.eventid = eventid
        self.station = station
        self.arrival_type = arrival_type
        self.theor_time = None
        self.pick_time = None
        self.travel_time = None
        self.time_dev = None

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
