# -*- coding: utf8 -*-
"""
Reading and writing of NonLinLoc grid files.

:copyright:
    2013-2018 Claudio Satriano <satriano@ipgp.fr>,
              Natalia Poiata <poiata@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import math
import numpy as np
from scipy.ndimage import zoom
from ctypes import Union, c_float, c_ushort
from copy import deepcopy
from pyproj import Proj
try:
    from collections.abc import Iterable
except ImportError:
    from collections import Iterable


valid_grid_types = (
    'VELOCITY',
    'VELOCITY_METERS',
    'SLOWNESS',
    'VEL2',
    'SLOW2',
    'SLOW2_METERS',
    'SLOW_LEN',
    'STACK',
    'TIME',
    'TIME2D',
    'PROB_DENSITY',
    'MISFIT',
    'ANGLE',
    'ANGLE2D'
)

valid_float_types = {
    # NLL_type: numpy_type
    'FLOAT': 'float32',
    'DOUBLE': 'float64'
}

valid_projections = (
    'NONE',
    'SIMPLE',
    'LAMBERT'
)

valid_ellipsoids = (
    'WGS-84',
    'GRS-80',
    'WGS-72',
    'Australian',
    'Krasovsky',
    'International',
    'Hayford-1909',
    'Clarke-1880',
    'Clarke-1866',
    'Airy',
    'Bessel',
    'Hayford-1830',
    'Sphere'
)


class TakeOffAngles(Union):
    """Union-style class for decoding take off angles."""

    _fields_ = [('fval', c_float),
                ('ival', c_ushort*2)]


class NLLGrid(object):
    """Class for manipulating NLL grid files.

    It has methods to read and write grid files,
    compute statistics and plot.
    """

    def __init__(self,
                 basename=None,
                 nx=1, ny=1, nz=1,
                 x_orig=0., y_orig=0., z_orig=0.,
                 dx=1., dy=1., dz=1.):
        """Init a NLLGrid object."""
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.x_orig = x_orig
        self.y_orig = y_orig
        self.z_orig = z_orig
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.__type = None
        self.__proj_name = None
        self.__proj_ellipsoid = None
        self.orig_lat = float(0)
        self.orig_lon = float(0)
        self.first_std_paral = 0
        self.second_std_paral = 0
        self.map_rot = float(0)
        self.station = None
        self.sta_x = float(0)
        self.sta_y = float(0)
        self.sta_z = float(0)
        self.float_type = 'FLOAT'
        self.__array = None
        self.xyz_mean = None
        self.xyz_cov = None
        self.ellipsoid = None
        if basename is not None:
            self.basename = self.remove_extension(basename)
            self.read_hdr_file()
            self.read_buf_file()
        else:
            self.basename = None

    def __str__(self):
        """Info string."""
        s = 'basename: {}\n'.format(self.basename)
        s += 'nx: {} ny: {} nz: {}\n'.format(self.nx, self.ny, self.nz)
        s += 'x_orig: {} y_orig: {} z_orig: {}\n'.format(
            self.x_orig, self.y_orig, self.z_orig)
        s += 'dx: {} dy: {} dz: {}\n'.format(self.dx, self.dy, self.dz)
        s += 'grid_type: {}\n'.format(self.type)
        s += 'float_type: {}\n'.format(self.float_type)
        if self.station is not None:
            s += 'station: {} sta_x: {} sta_y: {} sta_z: {}\n'.format(
                self.station, self.sta_x, self.sta_y, self.sta_z)
        s += 'transform: {}'.format(self.get_transform_line())
        return s

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__())

    def __getitem__(self, key):
        """Make the grid object array-like."""
        if self.type in ['ANGLE', 'ANGLE2D']:
            return self.dip[key]
        elif self.array is not None:
            return self.array[key]

    @property
    def array(self):
        """Property getter."""
        return self.__array

    @array.setter
    def array(self, array_data):
        array_data = np.asarray(array_data)
        if array_data.ndim != 3:
            raise ValueError('Only 3D arrays are supported')
        self.nx, self.ny, self.nz = array_data.shape
        self.__array = array_data

    @property
    def type(self):
        """Property getter."""
        return self.__type

    @type.setter
    def type(self, grid_type):
        try:
            grid_type = grid_type.upper()
        except AttributeError:
            raise ValueError('Grid type must be a string')
        if grid_type not in valid_grid_types:
            msg = 'Invalid grid type: {}\n'.format(grid_type)
            msg += 'Valid grid types are: {}'.format(valid_grid_types)
            raise ValueError(msg)
        self.__type = grid_type

    @property
    def float_type(self):
        """Property getter."""
        return self.__float_type

    @float_type.setter
    def float_type(self, float_type):
        try:
            float_type = float_type.upper()
        except AttributeError:
            raise ValueError('Float type must be a string')
        if float_type not in valid_float_types:
            msg = 'Invalid float type: {}\n'.format(float_type)
            msg += 'Valid grid types are: {}'.format(
                tuple(valid_float_types.keys()))
            raise ValueError(msg)
        self.__np_float_type = valid_float_types[float_type]
        self.__float_type = float_type

    @property
    def proj_name(self):
        """Property getter."""
        return self.__proj_name

    @proj_name.setter
    def proj_name(self, pname):
        try:
            pname = pname.upper()
        except AttributeError:
            raise ValueError('Projection name must be a string')
        if pname not in valid_projections:
            msg = 'Invalid projection name: {}\n'.format(pname)
            msg += 'Valid projection names: {}'.format(valid_projections)
            raise ValueError(msg)
        self.__proj_name = pname

    @property
    def proj_ellipsoid(self):
        """Property getter."""
        return self.__proj_ellipsoid

    @proj_ellipsoid.setter
    def proj_ellipsoid(self, ellipsoid):
        try:
            # here we just use upper() to check if ellipsoid is a string
            ellipsoid.upper()
        except AttributeError:
            raise ValueError('Ellipsoid must be a string')
        if ellipsoid not in valid_ellipsoids:
            msg = 'Invalid ellipsoid: {}\n'.format(ellipsoid)
            msg += 'Valid ellipsoids: {}'.format(valid_ellipsoids)
            raise ValueError(msg)
        self.__proj_ellipsoid = ellipsoid

    def remove_extension(self, basename):
        """Remove '.hdr' or '.buf' suffixes, if there."""
        bntmp = basename.rsplit('.hdr', 1)[0]
        return bntmp.rsplit('.buf', 1)[0]

    def init_array(self):
        """Init the array to zeros.

        Example of usage:
        >>> grd = NLLGrid(nx=20, ny=20, nz=30, dx=2., dy=2., dz=2.)
        >>> grd.init_array()
        >>> grd.array[2, 4, 10] = 3.
        """
        self.array = np.zeros((self.nx, self.ny, self.nz), float)

    def read_hdr_file(self, basename=None):
        """Read header file of NLL grid format."""
        if basename is not None:
            self.basename = self.remove_extension(basename)
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
        try:
            self.float_type = vals[10]
        except IndexError:
            self.float_type = 'FLOAT'

        lines.pop(0)

        for line in lines:
            vals = line.split()
            if not vals:
                # skip empty lines
                continue
            if vals[0] in ['TRANS', 'TRANSFORM']:
                if vals[1] == 'NONE':
                    self.proj_name = 'NONE'
                if vals[1] == 'SIMPLE':
                    self.proj_name = 'SIMPLE'
                    self.orig_lat = float(vals[3])
                    self.orig_lon = float(vals[5])
                    self.map_rot = float(vals[7])
                if vals[1] == 'LAMBERT':
                    self.proj_name = 'LAMBERT'
                    self.proj_ellipsoid = vals[3]
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
        """Read buf file as a 3d array."""
        if basename is not None:
            self.basename = self.remove_extension(basename)
        filename = self.basename + '.buf'

        with open(filename, 'rb') as fp:
            nitems = self.nx * self.ny * self.nz
            buf = np.fromfile(fp, dtype=self.__np_float_type, count=nitems)
            if len(buf) < nitems:
                raise ValueError(
                    'Not enough data values in buf file! '
                    '({} < {})'.format(len(buf), nitems))
        if self.type in ['ANGLE', 'ANGLE2D']:
            take_off_angles = (TakeOffAngles * len(buf))()
            for _i, _val in enumerate(buf):
                take_off_angles[_i].fval = _val
            self.azimuth = np.array(
                [t.ival[1]/10. for t in take_off_angles]
                ).reshape(self.nx, self.ny, self.nz)
            self.dip = np.array(
                [(t.ival[0]//16)/10. for t in take_off_angles]
                ).reshape(self.nx, self.ny, self.nz)
            self.quality = np.array(
                [t.ival[0] % 16 for t in take_off_angles]
                ).reshape(self.nx, self.ny, self.nz)
            self.azimuth[self.quality == 0] = np.nan
            self.dip[self.quality == 0] = np.nan
        else:
            self.array = np.array(buf).reshape(self.nx, self.ny, self.nz)

    def write_hdr_file(self, basename=None):
        """Write header file of NLL grid format."""
        if basename is not None:
            self.basename = basename
        filename = self.basename + '.hdr'

        lines = []
        lines.append('{} {} {}  {:.6f} {:.6f} {:.6f}  '
                     '{:.6f} {:.6f} {:.6f} {} {}\n'.format(
                        self.nx, self.ny, self.nz,
                        self.x_orig, self.y_orig, self.z_orig,
                        self.dx, self.dy, self.dz,
                        self.type, self.float_type))
        if self.station is not None:
            lines.append('{} {:.6f} {:.6f} {:.6f}\n'.format(
                self.station, self.sta_x, self.sta_y, self.sta_z))
        line = self.get_transform_line()
        if line is not None:
            lines.append('{}\n'.format(line))

        with open(filename, 'w') as fp:
            for line in lines:
                fp.write(line)

    def write_buf_file(self, basename=None):
        """Write buf file as a 3d array."""
        if self.type in ['ANGLE', 'ANGLE2D']:
            raise NotImplementedError(
                'Writing buf file not implemented for ANGLE grid.')
        if self.array is None:
            return
        if basename is not None:
            self.basename = basename
        filename = self.basename + '.buf'
        with open(filename, 'wb') as fp:
            self.array.astype(self.__np_float_type).tofile(fp)

    def get_transform_line(self):
        """Get the transform line in NLL hdr format."""
        if self.proj_name == 'NONE':
            return 'TRANSFORM  NONE'
        if self.proj_name == 'SIMPLE':
            line = 'TRANSFORM  SIMPLE  '
            line += 'LatOrig {:.6f}  LongOrig {:.6f}  RotCW {:.6f}'.format(
                self.orig_lat, self.orig_lon, self.map_rot)
            return line
        if self.proj_name == 'LAMBERT':
            line = 'TRANSFORM  LAMBERT RefEllipsoid {}  '.format(
                self.proj_ellipsoid)
            line += 'LatOrig {:.6f}  LongOrig {:.6f}  '.format(
                self.orig_lat, self.orig_lon)
            line += 'FirstStdParal {:.6f}  SecondStdParal {:.6f}  '.format(
                self.first_std_paral, self.second_std_paral)
            line += 'RotCW {:.6f}'.format(self.map_rot)
            return line

    def get_xyz(self, i, j, k):
        """Get cartesian coordinates (x, y, z) for grid indexes (i, j, k)."""
        x = i * self.dx + self.x_orig
        y = j * self.dy + self.y_orig
        z = k * self.dz + self.z_orig
        return x, y, z

    def get_ijk(self, x, y, z):
        """Get grid indexes (i, j, k) for cartesian coordinates (x, y, z)."""
        i = int((x - self.x_orig) / self.dx)
        j = int((y - self.y_orig) / self.dy)
        k = int((z - self.z_orig) / self.dz)
        return i, j, k

    def get_ijk_max(self):
        """Return the indexes (i, j, k) of the grid max point."""
        if self.array is None:
            return None
        return np.unravel_index(self.array.argmax(), self.array.shape)

    def get_ijk_min(self):
        """Return the indexes (i,j,k) of the grid min point."""
        if self.array is None:
            return None
        return np.unravel_index(self.array.argmin(), self.array.shape)

    def get_xyz_max(self):
        """Return the coordinates (x,y,z) of the grid max point."""
        ijk_max = self.get_ijk_max()
        if ijk_max is None:
            return None
        return self.get_xyz(*ijk_max)

    def get_xyz_min(self):
        """Return the coordinates (x,y,z) of the grid min point."""
        ijk_min = self.get_ijk_min()
        if ijk_min is None:
            return None
        return self.get_xyz(*ijk_min)

    def get_ijk_mean(self):
        """Return the indexes (i,j,k) of the grid mean point."""
        xyz_mean = self.get_xyz_mean()
        if xyz_mean is None:
            return None
        return self.get_ijk(*xyz_mean)

    def get_xyz_mean(self):
        """Return the coordinates (x,y,z) of the grid mean point."""
        if self.array is None:
            return None
        xx = np.arange(0, self.nx) * self.dx + self.x_orig
        yy = np.arange(0, self.ny) * self.dy + self.y_orig
        zz = np.arange(0, self.nz) * self.dz + self.z_orig
        yarray, xarray, zarray = np.meshgrid(yy, xx, zz)
        array_sum = self.array.sum()
        xmean = (xarray * self.array).sum()/array_sum
        ymean = (yarray * self.array).sum()/array_sum
        zmean = (zarray * self.array).sum()/array_sum
        self.xyz_mean = (xmean, ymean, zmean)
        return (xmean, ymean, zmean)

    def get_xyz_cov(self):
        """Return the grid covariance respect to the (x,y,z) mean point."""
        if self.array is None:
            return None
        xyz_mean = self.get_xyz_mean()
        xx = np.arange(0, self.nx) * self.dx + self.x_orig
        yy = np.arange(0, self.ny) * self.dy + self.y_orig
        zz = np.arange(0, self.nz) * self.dz + self.z_orig
        yarray, xarray, zarray = np.meshgrid(yy, xx, zz)
        array_sum = self.array.sum()
        cov = np.zeros((3, 3))
        cov[0, 0] = (np.power(xarray, 2) * self.array).sum()/array_sum \
            - (xyz_mean[0] * xyz_mean[0])
        cov[0, 1] = cov[1, 0] =\
            (xarray * yarray * self.array).sum()/array_sum \
            - (xyz_mean[0] * xyz_mean[1])
        cov[0, 2] = cov[2, 0] = \
            (xarray * zarray * self.array).sum()/array_sum \
            - (xyz_mean[0] * xyz_mean[2])
        cov[1, 1] = (np.power(yarray, 2) * self.array).sum()/array_sum \
            - (xyz_mean[1] * xyz_mean[1])
        cov[1, 2] = cov[2, 1] = \
            (yarray * zarray * self.array).sum()/array_sum \
            - (xyz_mean[1] * xyz_mean[2])
        cov[2, 2] = (np.power(zarray, 2) * self.array).sum()/array_sum \
            - (xyz_mean[2] * xyz_mean[2])
        self.xyz_cov = cov
        return cov

    def get_xyz_ellipsoid(self):
        """Return the 68% confidence ellipsoid."""
        from ellipsoid import Ellipsoid3D
        # The following code is a python translation of the
        # CalcErrorEllipsoid() c-function from the NonLinLoc package,
        # written by Anthony Lomax
        cov = self.get_xyz_cov()
        if cov is None:
            return None

        u, s, v = np.linalg.svd(cov)

        del_chi_2 = 3.53  # 3.53: value for 68% conf
        ell = Ellipsoid3D()
        ell.az1 = math.degrees(math.atan2(u[0, 0], u[1, 0]))
        if ell.az1 < 0.0:
            ell.az1 += 360.0
        ell.dip1 = math.degrees(math.asin(u[2, 0]))
        ell.len1 = math.sqrt(del_chi_2) / math.sqrt(1.0 / s[0])
        ell.az2 = math.degrees(math.atan2(u[0, 1], u[1, 1]))
        if ell.az2 < 0.0:
            ell.az2 += 360.0
        ell.dip2 = math.degrees(math.asin(u[2, 1]))
        ell.len2 = math.sqrt(del_chi_2) / math.sqrt(1.0 / s[1])
        ell.len3 = math.sqrt(del_chi_2) / math.sqrt(1.0 / s[2])

        self.ellipsoid = ell
        return ell

    def get_value(self, x, y, z, array=None):
        """Get the array value at specified cartesian coordinates (x, y, z)."""
        if array is None:
            if self.array is None:
                return
            else:
                array = self.array
        min_x, max_x, min_y, max_y, min_z, max_z = self.get_extent()
        if not (min_x <= x <= max_x and min_y <= y <= max_y and
                min_z <= z <= max_z):
            raise ValueError('point {} outside the grid.'.format((x, y, z)))
        i, j, k = self.get_ijk(x, y, z)
        return array[i, j, k]

    def get_extent(self):
        """Get the grid extent in cartesian units (generally km)."""
        extent = (self.x_orig - self.dx / 2,
                  self.x_orig + self.nx * self.dx + self.dx / 2,
                  self.y_orig - self.dy / 2,
                  self.y_orig + self.ny * self.dy + self.dy / 2,
                  self.z_orig - self.dz / 2,
                  self.z_orig + self.nz * self.dz + self.dz / 2
                  )
        return extent

    def get_xy_extent(self):
        """Get the grid xy extent in cartesian units (generally km)."""
        return self.get_extent()[0:4]

    def get_xz_extent(self):
        """Get the grid xz extent in cartesian units (generally km)."""
        return self.get_extent()[0:2] + self.get_extent()[4:]

    def get_zx_extent(self):
        """Get the grid zx extent in cartesian units (generally km)."""
        return self.get_extent()[4:] + self.get_extent()[0:2]

    def get_yz_extent(self):
        """Get the grid yz extent in cartesian units (generally km)."""
        return self.get_extent()[2:]

    def get_zy_extent(self):
        """Get the grid zy extent in cartesian units (generally km)."""
        return self.get_extent()[4:] + self.get_extent()[2:4]

    def max(self):
        """Get the grid max value."""
        if self.type in ['ANGLE', 'ANGLE2D']:
            return np.nanmax(self.dip)
        if self.array is not None:
            return np.nanmax(self.array)

    def resample(self, dx, dy, dz):
        """Resample grid to (dx, dy, dz)."""
        if self.type in ['ANGLE', 'ANGLE2D']:
            raise NotImplementedError(
                'Resample not implemented for ANGLE grid.')
        zoom_x = self.dx / dx
        zoom_y = self.dy / dy
        zoom_z = self.dz / dz
        self.array = zoom(self.array, (zoom_x, zoom_y, zoom_z))
        self.nx, self.ny, self.nz = self.array.shape
        if self.type == 'SLOW_LEN':
            self.array *= dx / self.dx
        self.dx = dx
        self.dy = dy
        self.dz = dz

    def get_plot_axes(self, figure=None, ax_xy=None):
        """Get the axes for the three projections, plus the colorbar axis.

        Requires Matplotlib.
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        xmin, xmax, ymin, ymax, zmin, zmax = self.get_extent()
        if figure is None and ax_xy is None:
            figure = plt.figure()
        if figure is None and ax_xy is not None:
            figure = ax_xy.get_figure()
        if ax_xy is None:
            ax_xy = figure.add_subplot(111)

        ratio = float(xmax - xmin) / (ymax - ymin)
        plot_xz_size = ((zmax - zmin)/(xmax - xmin))*100
        plot_yz_size = plot_xz_size / ratio
        plot_cbar_size = 5  # percent
        xz_size = '%f %%' % plot_xz_size
        yz_size = '%f %%' % plot_yz_size
        cb_size = '%f %%' % plot_cbar_size

        # ax_xy
        divider = make_axes_locatable(ax_xy)
        plt.setp(ax_xy.get_xticklabels(), visible=False)
        ax_xy.set_xlim(xmin, xmax)
        ax_xy.set_ylim(ymin, ymax)
        ax_xy.set_aspect('equal', 'datalim')
        plt.setp(ax_xy.get_yticklabels(), rotation=90, fontsize=12)

        # ax_yz
        ax_yz = divider.append_axes(
            'right', size=yz_size, pad=0.05, sharey=ax_xy)
        plt.setp(ax_yz.get_yticklabels(), visible=False)
        ax_yz.set_xlim(zmin, zmax)
        ax_yz.set_ylim(ymin, ymax)
        plt.setp(ax_yz.get_xticklabels(), rotation=90, fontsize=12)
        plt.setp(ax_yz.get_yticklabels(), rotation=90, fontsize=12)

        # ax_xz
        ax_xz = divider.append_axes(
            'bottom', size=xz_size, pad=0.05, sharex=ax_xy)
        ax_xz.set_xlim(xmin, xmax)
        ax_xz.set_ylim(zmax, zmin)

        # color-bar
        ax_cb = divider.append_axes('bottom', size=cb_size, pad=0.5)

        return ax_xy, ax_xz, ax_yz, ax_cb

    def plot(self, slice_index=None, handle=False, figure=None, ax_xy=None,
             vmin=None, vmax=None, cmap=None, array=None):
        """Plot the grid using three orthogonal projections.

        Requires Matplotlib.
        """
        import matplotlib.pyplot as plt
        from matplotlib import ticker

        if array is None:
            if self.array is None:
                return
            else:
                array = self.array

        ax_xy, ax_xz, ax_yz, ax_cb = self.get_plot_axes(figure, ax_xy)
        if figure is None:
            figure = ax_xy.get_figure()

        if slice_index is None:
            slice_index = list(map(int, (self.nx/2, self.ny/2, self.nz/2)))
        if slice_index == 'max':
            slice_index = self.get_ijk_max()
        if slice_index == 'min':
            slice_index = self.get_ijk_min()

        if vmin is None:
            vmin = np.nanmin(array)
        if vmax is None:
            vmax = np.nanmax(array)

        hnd = ax_xy.imshow(np.transpose(array[:, :, slice_index[2]]),
                           vmin=vmin, vmax=vmax, cmap=cmap,
                           origin='lower', extent=self.get_xy_extent(),
                           zorder=-10)
        ax_xz.imshow(np.transpose(array[:, slice_index[1], :]),
                     vmin=vmin, vmax=vmax, cmap=cmap,
                     origin='lower', extent=self.get_xz_extent(),
                     aspect='auto', zorder=-10)
        ax_yz.imshow(array[slice_index[0], :, :],
                     vmin=vmin, vmax=vmax, cmap=cmap,
                     origin='lower', extent=self.get_zy_extent(),
                     aspect='auto', zorder=-10)

        x_slice, y_slice, z_slice = self.get_xyz(*slice_index)
        ax_xy.axhline(y_slice, color='w', linestyle='dashed', zorder=-1)
        ax_xy.axvline(x_slice, color='w', linestyle='dashed', zorder=-1)
        ax_xz.axhline(z_slice, color='w', linestyle='dashed', zorder=-1)
        ax_yz.axvline(z_slice, color='w', linestyle='dashed', zorder=-1)

        fmt = '%.1e' if np.nanmax(array) <= 0.01 else '%.2f'
        cb = figure.colorbar(
            hnd, cax=ax_cb, orientation='horizontal', format=fmt)
        cb.locator = ticker.LinearLocator(numticks=3)
        cb.update_ticks()

        if handle:
            return (ax_xy, ax_xz, ax_yz), cb
        else:
            plt.show()

    def plot_3D_point(self, axes, point, color='r'):
        """Plot a point (i, j, k) on the grid."""
        ax_xy, ax_xz, ax_yz = axes
        ax_xy.scatter(point[0], point[1], color=color)
        ax_xz.scatter(point[0], point[2], color=color)
        ax_yz.scatter(point[2], point[1], color=color)

    def plot_ellipsoid(self, axes, ellipsoid=None, mean_xyz=None):
        """Plot an ellipsoid on the grid."""
        from ellipsoid import Vect3D, ellipsiod2Axes, toEllipsoid3D
        ax_xy, ax_xz, ax_yz = axes

        if ellipsoid is None:
            ellipsoid = self.get_xyz_ellipsoid()
        expect = Vect3D()
        if mean_xyz is None:
            mean_xyz = self.get_xyz_mean()
        expect.x, expect.y, expect.z = mean_xyz

        pax1, pax2, pax3 = ellipsiod2Axes(ellipsoid)

        ellArray12 = toEllipsoid3D(pax1, pax2, expect, 100)
        ellArray13 = toEllipsoid3D(pax1, pax3, expect, 100)
        ellArray23 = toEllipsoid3D(pax2, pax3, expect, 100)
        ell12 = np.array([(vect.x, vect.y) for vect in ellArray12])
        ell13 = np.array([(vect.x, vect.y) for vect in ellArray13])
        ell23 = np.array([(vect.x, vect.y) for vect in ellArray23])

        ax_xy.plot(ell12[:, 0], ell12[:, 1])
        ax_xy.plot(ell13[:, 0], ell13[:, 1])
        ax_xy.plot(ell23[:, 0], ell23[:, 1])

        ell12 = np.array([(vect.x, vect.z) for vect in ellArray12])
        ell13 = np.array([(vect.x, vect.z) for vect in ellArray13])
        ell23 = np.array([(vect.x, vect.z) for vect in ellArray23])
        ax_xz.plot(ell12[:, 0], ell12[:, 1])
        ax_xz.plot(ell13[:, 0], ell13[:, 1])
        ax_xz.plot(ell23[:, 0], ell23[:, 1])

        ell12 = np.array([(vect.y, vect.z) for vect in ellArray12])
        ell13 = np.array([(vect.y, vect.z) for vect in ellArray13])
        ell23 = np.array([(vect.y, vect.z) for vect in ellArray23])
        ax_yz.plot(ell12[:, 1], ell12[:, 0])
        ax_yz.plot(ell13[:, 1], ell13[:, 0])
        ax_yz.plot(ell23[:, 1], ell23[:, 0])

    def project(self, lon, lat):
        """Project lon, lat into grid coordinates."""
        if self.proj_name == 'LAMBERT':
            ellipsoids = {
                'WGS-84': 'WGS84',
                'GRS-80': 'GRS80',
                'WGS-72': 'WGS72',
                'Australian': 'aust_SA',
                'Krasovsky': 'krass',
                'International': 'new_intl',
                'Hayford-1909': 'intl',
                'Clarke-1880': 'clrk80',
                'Clarke-1866': 'clrk66',
                'Airy': 'airy',
                'Bessel': 'bessel',
                # 'Hayford-1830':
                'Sphere': 'sphere'
                }
            try:
                ellps = ellipsoids[self.proj_ellipsoid]
            except KeyError:
                raise ValueError(
                    'Ellipsoid not supported: {}'.format(self.proj_ellipsoid))
            p = Proj(proj='lcc', lat_0=self.orig_lat, lon_0=self.orig_lon,
                     lat_1=self.first_std_paral, lat_2=self.second_std_paral,
                     ellps=ellps)
        elif self.proj_name == 'SIMPLE':
            p = Proj(proj='eqc', lat_0=self.orig_lat, lon_0=self.orig_lon)
        else:
            raise ValueError('Projection not supported: {}'.format(
                self.proj_name))
        x, y = p(lon, lat)
        x = np.array(x)
        y = np.array(y)
        x[np.isnan(lon)] = np.nan
        y[np.isnan(lat)] = np.nan
        x /= 1000.
        y /= 1000.
        # inverse rotation
        theta = np.radians(self.map_rot)
        x1 = x*np.cos(theta) + y*np.sin(theta)
        y1 = -x*np.sin(theta) + y*np.cos(theta)
        x = x1
        y = y1
        # Try to return the same type of lon, lat
        if not isinstance(lon, np.ndarray):
            if isinstance(lon, Iterable):
                x = type(lon)(x)
            else:
                x = float(x)
        if not isinstance(lat, np.ndarray):
            if isinstance(lat, Iterable):
                y = type(lat)(y)
            else:
                y = float(y)
        return x, y

    def copy(self):
        """Get a deep copy of the grid object."""
        return deepcopy(self)


def main():
    """Test code.

    It generates a gaussian grid and computes the 3D ellipsoid
    around the grid mean.
    """
    import matplotlib.pyplot as plt

    # http://stackoverflow.com/q/17190649
    def gauss3D(shape=(3, 3, 3), sigmax=0.5, sigmay=0.5, sigmaz=0.5, theta=0):
        m, n, k = [(ss-1.)/2. for ss in shape]
        y, x, z = np.ogrid[-m:m+1, -n:n+1, -k:k+1]
        # xy rotation
        theta = np.radians(theta)
        x2 = math.cos(theta)*x - math.sin(theta)*y
        y2 = math.sin(theta)*x + math.cos(theta)*y
        h = np.exp(-(x2*x2) / (2.*sigmax*sigmax)) *\
            np.exp(-(y2*y2) / (2.*sigmay*sigmay)) *\
            np.exp(-(z*z) / (2.*sigmaz*sigmaz))
        h[h < np.finfo(h.dtype).eps*h.max()] = 0
        sumh = h.sum()
        if sumh != 0:
            h /= sumh
        return h

    # Generate the grid
    nx = 101
    ny = 201
    nz = 11
    x_orig = -50
    y_orig = -100
    grd = NLLGrid(nx=nx, ny=ny, nz=nz,
                  dx=1, dy=1, dz=1,
                  x_orig=x_orig, y_orig=y_orig)
    print(grd)
    grd.array = gauss3D((nx, ny, nz), 20, 10, 2, 30)

    # Compute statistics
    mean_xyz = grd.get_xyz_mean()
    max_ijk = grd.get_ijk_max()
    max_xyz = grd.get_xyz_max()

    # Plotting
    axes, cb = grd.plot(max_ijk, handle=True)
    grd.plot_3D_point(axes, mean_xyz, color='g')
    grd.plot_3D_point(axes, max_xyz, color='r')
    grd.plot_ellipsoid(axes, mean_xyz=mean_xyz)
    plt.show()


if __name__ == '__main__':
    main()
