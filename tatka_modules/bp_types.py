# bp_types.py
# Data types for BackProj
from obspy import UTCDateTime

class Trigger():
    def __init__(self,
                 eventid=None,list_picks=None,
                 x=None, y=None, z=None,
                 i=None, j=None, k=None,
                 max_grid=None,
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
        self.beg_win = beg_win
        self.end_win = end_win
        self.center_win = center_win
        self.lat = None
        self.lon = None
        self.origin_time = None

    def __str__(self):
        s = '%s ' % self.eventid
        s += 'X %s ' % self.x
        s += 'Y %s ' % self.y
        s += 'Z %s ' % self.z
        s += 'MaxStack %.3f ' % (self.max_grid)
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
            if word[1] != 'X' or word[17] != 'T_ORIG':
                raise ValueError, 'Not a valid trigger string'
        except IndexError:
            raise ValueError, 'Not a valid trigger string'
        self.eventid = word[0]
        self.x = float(word[2])
        self.y = float(word[4])
        self.z = float(word[6])
        self.max_grid = float(word[8])
        self.beg_win = float(word[10])
        self.end_win = float(word[12])
        self.lat = float(word[14])
        self.lon = float(word[16])
        self.origin_time = UTCDateTime(word[18])

    def add_pick(self, pick):
        self.picks.append(pick)

    def get_picks(self, station=None, arrival_type=None):
        return [pick for pick in self.picks
                if (station is not None and pick.station == station)
                or (arrival_type is not None and pick.arrival_type == arrival_type)]


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
                raise ValueError, 'Not a valid pick string'
        except IndexError:
            raise ValueError, 'Not a valid pick string'
        self.station = word[1]
        self.arrival_type = word[3]
        self.theor_time = float(word[5])
        self.pick_time = float(word[7])
        #TODO: read and write these fields?
        #self.travel_time = None
        #self.time_dev = None


class RecursiveMemory():
    def __init__(self):
        self.filterH1 = 0.
        self.filterH2 = 0.
        self.filterL1 = 0.
        self.filterL2 = 0.
        self.previous_sample = 0.
        self.mean_sq = 0.
        self.mean = 0.
        self.var = 1.
        self.hos = 0.

    def __str__(self):
        s = 'filterH1: %f\n' % self.filterH1
        s += 'filterH2: %f\n' % self.filterH2
        s += 'filterL1: %f\n' % self.filterL1
        s += 'filterL2: %f\n' % self.filterL2
        s += 'previous_sample: %f\n' % self.previous_sample
        s += 'mean_sq: %f\n' % self.mean_sq
        s += 'mean: %f\n' % self.mean
        s += 'var: %f\n' % self.var
        s += 'hos: %f\n' % self.hos
        return s
