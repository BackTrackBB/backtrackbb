# -*- coding: utf8 -*-

import numpy as np
from array import array

class NLLGrid():
    '''
    Class for reading NLL grid files
    '''

    def __init__(self, basename=''):
        self.nx = int(0)
        self.ny = int(0)
        self.nz = int(0)
        self.x_orig = float(0)
        self.y_orig = float(0)
        self.z_orig = float(0)
        self.dx = float(0)
        self.dy = float(0)
        self.dz = float(0)
        self.type = ''
        self.proj_name = ''
        self.orig_lat = float(0)
        self.orig_lon = float(0)
        self.map_rot = float(0)
        self.station = ''
        self.sta_x = float(0)
        self.sta_y = float(0)
        self.sta_z = float(0)
        self.array = None
        self.basename = basename
        if basename:
            self.read_hdr_file(basename)
            self.read_buf_file(basename)

    def __getitem__(self, key):
        if self.array is not None:
            return self.array[key]

    def read_hdr_file(self, basename):
        '''Reads header file of NLLoc output'''
        self.basename = basename
        filename = basename + '.hdr'

        # read header file
        f = open(filename)
        lines = f.readlines()
        f.close()

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

    def read_buf_file(self, basename):
        '''reads buf file as a 3d array'''
        self.basename = basename
        buffilename = basename + '.buf'

        f1 = open(buffilename, 'rb')
        buf = array('f')
        buf.fromfile(f1,
                     self.nx *
                     self.ny *
                     self.nz)
        f1.close()
        self.array = np.array(buf).reshape(self.nx, self.ny, self.nz)

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
