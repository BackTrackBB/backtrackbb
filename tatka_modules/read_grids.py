# -*- coding: utf8 -*-
import os
from collections import defaultdict
from NLLGrid import NLLGrid

def read_grids(config):
    GRD_sta = defaultdict(dict)
    coord_sta = {}

    for grid_type in config.grid_type:
        for station in config.stations:
            grid_bname = '.'.join(('layer', grid_type, station, 'time'))
            grid_bname = os.path.join(config.grid_dir, grid_bname)
            grid = NLLGrid(grid_bname)
            coord_sta[station] = (grid.sta_x, grid.sta_y)
            GRD_sta[station][grid_type] = grid

    return GRD_sta, coord_sta
