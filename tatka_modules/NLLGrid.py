# -*- coding: utf8 -*-
# NLLGrid.py
#
# Reading and writing of NLL grid files
#
# (c) 2013-2014 - Natalia Poiata <poiata@ipgp.fr>,
#                 Claudio Satriano <satriano@ipgp.fr>
import numpy as np
from array import array


class NLLGrid():
    '''
    Class for reading and writing NLL grid files
    '''

    def __init__(self,
                 basename=None,
                 nx=0, ny=0, nz=0,
                 x_orig=0., y_orig=0., z_orig=0.,
                 dx=0., dy=0., dz=0.):
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.x_orig = x_orig
        self.y_orig = y_orig
        self.z_orig = z_orig
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.type = None
        self.proj_name = None
        self.orig_lat = float(0)
        self.orig_lon = float(0)
        self.first_std_paral = None
        self.second_std_paral = None
        self.map_rot = float(0)
        self.station = None
        self.sta_x = float(0)
        self.sta_y = float(0)
        self.sta_z = float(0)
        self.array = None
        self.basename = basename
        if basename:
            self.read_hdr_file()
            self.read_buf_file()

    def __getitem__(self, key):
        if self.array is not None:
            return self.array[key]

    def init_array(self):
        self.array = np.zeros((self.nx, self.ny, self.nz), float)

    def read_hdr_file(self, basename=None):
        '''Reads header file of NLL grid format'''
        if basename is not None:
            self.basename = basename
        filename = self.basename + '.hdr'

        # read header file
        with open(filename, 'r') as fp:
            lines = fp.readlines()

        # extract information
        vals = lines[0].split()
        self.nx = int(vals[0])
        self.ny = int(vals[1])
        self.nz = int(vals[2])
        self.x_orig = float(vals[3])
        self.y_orig = float(vals[4])
        self.z_orig = float(vals[5])
        self.dx = float(vals[6])
        self.dy = float(vals[7])
        self.dz = float(vals[8])
        self.type = vals[9]

        lines.pop(0)

        for line in lines:
            vals = line.split()
            if vals[0] == 'TRANSFORM':
                if vals[1] == 'NONE':
                    self.proj_name = 'NONE'
                if vals[1] == 'SIMPLE':
                    self.proj_name = 'SIMPLE'
                    self.orig_lat = float(vals[3])
                    self.orig_lon = float(vals[5])
                    self.map_rot = float(vals[7])
                if vals[1] == 'LAMBERT':
                    self.proj_name = 'LAMBERT'
                    self.ellipsoid = vals[3]
                    self.orig_lat = float(vals[5])
                    self.orig_lon = float(vals[7])
                    self.first_std_paral = float(vals[9])
                    self.second_std_paral = float(vals[11])
                    self.map_rot = float(vals[13])
            else:
                self.station = vals[0]
                self.sta_x = float(vals[1])
                self.sta_y = float(vals[2])
                self.sta_z = float(vals[3])

    def read_buf_file(self, basename=None):
        '''Reads buf file as a 3d array'''
        if basename is not None:
            self.basename = basename
        filename = self.basename + '.buf'

        with open(filename, 'rb') as fp:
            buf = array('f')
            buf.fromfile(fp, self.nx * self.ny * self.nz)
        self.array = np.array(buf).reshape(self.nx, self.ny, self.nz)

    def write_hdr_file(self, basename=None):
        '''Writes header file of NLL grid format'''
        if basename is not None:
            self.basename = basename
        filename = self.basename + '.hdr'

        lines = []
        lines.append('%d %d %d  %.6f %.6f %.6f  %.6f %.6f %.6f %s\n' %
                (self.nx, self.ny, self.nz,
                 self.x_orig, self.y_orig, self.z_orig,
                 self.dx, self.dy, self.dz,
                 self.type))
        if self.station is not None:
            lines.append('%s %.6f %.6f %.6f\n' %
                    (self.station, self.sta_x, self.sta_y, self.sta_z))
        if self.proj_name == 'NONE':
            lines.append('TRANSFORM  NONE\n')
        if self.proj_name == 'SIMPLE':
            lines.append('TRANSFORM  SIMPLE  LatOrig %.6f  LongOrig %.6f  RotCW %.6f\n' %
                    (self.orig_lat, self.orig_lon, self.map_rot))
        if self.proj_name == 'LAMBERT':
            line = 'TRANSFORM  LAMBERT RefEllipsoid %s  ' % self.ellipsoid
            line += 'LatOrig %.6f  LongOrig %.6f  ' % (self.orig_lat, self.orig_lon)
            line += 'FirstStdParal %.6f  SecondStdParal %.6f  RotCW %.6f\n' %\
                    (self.first_std_paral, self.second_std_paral, self.map_rot)
            lines.append(line)

        with open(filename, 'w') as fp:
            for line in lines:
                fp.write(line)

    def write_buf_file(self, basename=None):
        '''Writes buf file as a 3d array'''
        if self.array is None:
            return

        if basename is not None:
            self.basename = basename
        filename = self.basename + '.buf'

        with open(filename, 'wb') as fp:
            self.array.astype(np.float32).tofile(fp)

    def get_xyz(self, i, j, k):
        x = i * self.dx + self.x_orig
        y = j * self.dy + self.y_orig
        z = k * self.dz + self.z_orig
        return x, y, z

    def get_ijk(self, x, y, z):
        i = int((x - self.x_orig) / self.dx)
        j = int((y - self.y_orig) / self.dy)
        k = int((z - self.z_orig) / self.dz)
        return i, j, k

    def get_value(self, x, y, z):
        if self.array is None:
            return
        i, j, k = self.get_ijk(x, y, z)
        return self.array[i, j, k]

    def get_extent(self):
        extent = (self.x_orig - self.dx / 2,
                  self.x_orig + self.nx * self.dx + self.dx / 2,
                  self.y_orig - self.dy / 2,
                  self.y_orig + self.ny * self.dy + self.dy / 2,
                  self.z_orig - self.dz / 2,
                  self.z_orig + self.nz * self.dz + self.dz / 2
                  )
        return extent

    def get_xy_extent(self):
        return self.get_extent()[0:4]

    def get_xz_extent(self):
        return self.get_extent()[0:2] + self.get_extent()[4:]

    def get_zx_extent(self):
        return self.get_extent()[4:] + self.get_extent()[0:2]

    def get_yz_extent(self):
        return self.get_extent()[2:]

    def get_zy_extent(self):
        return self.get_extent()[4:] + self.get_extent()[2:4]

    def max(self):
        if self.array is not None:
            return self.array.max()
