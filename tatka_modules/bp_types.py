# bp_types.py
# Data types for BackProj

class Trigger():
    def __init__(self,
                 x=None, y=None, z=None,
                 beg_win=None, end_win=None,
                 center_win=None):
        self.x = x
        self.y = y
        self.z = z
        self.beg_win = beg_win
        self.end_win = end_win
        self.center_win = center_win

    def __str__(self):
        s = 'X %s ' % self.x
        s += 'Y %s ' % self.y
        s += 'Z %s ' % self.z
        s += 'BEG %s ' % self.beg_win
        s += 'END %s' % self.end_win
        return s
