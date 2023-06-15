# -*- coding: utf8 -*-
from obspy.signal.util import util_geo_km


def read_locationTremor(infile, hour, lat_or, lon_or, depth=35.):
    f = open(infile, 'r')
    vardict = {'yyyy': 0, 'mm': 1, 'dd': 2, 'hh': 3, 'lat_e': 4,
               'lon_e': 5, 'energy': 6, 'n_eve': 7}

    epicenter = dict.fromkeys(vardict.keys())

    for key in epicenter.keys():
        epicenter[key] = []

    for line in f:
        data = [float(x) for x in line.split()]
        for var, index in vardict.items():
            epicenter[var].append(data[index])

    XX = []
    YY = []
    ZZ = []
    hours_dic = epicenter.get('hh')
    for i, hh in enumerate(hours_dic):
        if hh == float(hour):
            xeq, yeq = util_geo_km(lon_or, lat_or, epicenter['lon_e'][i],
                                   epicenter['lat_e'][i])
            XX.append(xeq)
            YY.append(yeq)
            ZZ.append(depth)
    return XX, YY, ZZ


def read_locationEQ(infile, day, hours, lat_zero, lon_zero, depth=35.):
    f = open(infile, 'r')
    vardict = {'yyyy': 0, 'mm': 1, 'dd': 2, 'hh': 3, 'min': 4, 'sec': 5,
               'es': 6, 'lat_e': 7, 'e_lat': 8, 'lon_e': 9, 'e_lon': 10,
               'depth_e': 11, 'e_depth': 12, 'mag': 13}
    epicenter = dict.fromkeys(vardict.keys())

    for key in epicenter.keys():
        epicenter[key] = []

    for line in f:
        data = [float(x) for x in line.split()]
        for var, index in vardict.items():
            epicenter[var].append(data[index])

    xx = []
    yy = []
    zz = []

    hours_dic = epicenter.get('hh')
    for i, hh in enumerate(hours_dic):
        year = str(epicenter['yyyy'][i])[2:4] + '0' + \
               str(int(epicenter['mm'][i])) + str(int(epicenter['dd'][i]))
        if year == day and str(int(hh)) == hours:
            xeq, yeq = util_geo_km(lon_zero, lat_zero, epicenter['lon_e'][i],
                                   epicenter['lat_e'][i])
            xx.append(xeq)
            yy.append(yeq)
            zz.append(epicenter['depth_e'][i])
    return xx, yy, zz
