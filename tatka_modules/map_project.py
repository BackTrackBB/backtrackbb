# -*- coding: utf8 -*-
import os
import ctypes

libpath = os.path.join(os.path.dirname(__file__), os.pardir,
                       'lib', 'lib_map_project.so')
lib_map_project = ctypes.CDLL(libpath)

lib_map_project.get_transform.argtypes = [
        ctypes.c_int,
        ctypes.c_char_p
        ]

lib_map_project.latlon2rect.argtypes = [
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ]

lib_map_project.rect2latlon.argtypes = [
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ]


def get_transform(trans_type, orig_lat, orig_lon,
                  first_std_paral=None,
                  second_std_paral=None,
                  map_rot=0.,
                  ellipsoid=None):
    #TODO: add other trans_types
    if trans_type == 'LAMBERT':
        s = 'LAMBERT %s %s %s %s %s %s' %\
                (ellipsoid, orig_lat, orig_lon,
                 first_std_paral, second_std_paral, map_rot)
    if trans_type == 'SIMPLE':
        s = 'SIMPLE %s %s %s' %\
                (orig_lat, orig_lon, map_rot)

    lib_map_project.get_transform(0, s)


def latlon2rect(lat, lon):
    x = ctypes.c_double()
    y = ctypes.c_double()
    lib_map_project.latlon2rect(0, lat, lon, ctypes.byref(x), ctypes.byref(y))
    return x.value, y.value


def rect2latlon(x, y):
    lat = ctypes.c_double()
    lon = ctypes.c_double()
    lib_map_project.rect2latlon(0, x, y, ctypes.byref(lat), ctypes.byref(lon))
    return lat.value, lon.value


if __name__ == '__main__':
    lat0 = 51.656300
    lon0 = 7.742580
    std_par1 = 50.
    std_par2 = 52.
    rot = 0.
    get_transform('LAMBERT', lat0, lon0, std_par1, std_par2,
                  rot, 'Clarke-1880')
    print latlon2rect(lat0, lon0)
    x, y = latlon2rect(50.1, 7.2)
    print rect2latlon(x, y)
    get_transform('SIMPLE', lat0, lon0, None, None, rot, None)
    print latlon2rect(lat0, lon0)
    x, y = latlon2rect(50.1, 7.2)
    print rect2latlon(x, y)
