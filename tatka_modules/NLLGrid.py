# -*- coding: utf8 -*-
# NLLGrid.py
#
# Reading and writing of NLL grid files
#
# (c) 2013-2015 - Natalia Poiata <poiata@ipgp.fr>,
#                 Claudio Satriano <satriano@ipgp.fr>
import math
import numpy as np
from array import array
from copy import deepcopy


class NLLGrid():
    '''
    Class for manipulating NLL grid files.
    It has methods to read and write grid files,
    compute statistics and plot.
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
        self.xyz_mean = None
        self.xyz_cov = None
        self.ellipsoid = None
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
        '''
        Reads header file of NLL grid format
        '''
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
        '''
        Reads buf file as a 3d array
        '''
        if basename is not None:
            self.basename = basename
        filename = self.basename + '.buf'

        with open(filename, 'rb') as fp:
            buf = array('f')
            buf.fromfile(fp, self.nx * self.ny * self.nz)
        self.array = np.array(buf).reshape(self.nx, self.ny, self.nz)

    def write_hdr_file(self, basename=None):
        '''
        Writes header file of NLL grid format
        '''
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
        '''
        Writes buf file as a 3d array
        '''
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

    def get_ijk_max(self):
        '''
        Returns the indexes (i,j,k) of the grid max point
        '''
        if self.array is None:
            return None
        return np.unravel_index(self.array.argmax(), self.array.shape)

    def get_xyz_max(self):
        '''
        Returns the coordinates (x,y,z) of the grid max point
        '''
        ijk_max = self.get_ijk_max()
        if ijk_max is None:
            return None
        return self.get_xyz(*ijk_max)

    def get_ijk_mean(self):
        '''
        Returns the indexes (i,j,k) of the grid mean point
        '''
        xyz_mean = self.get_xyz_mean()
        if xyz_mean is None:
            return None
        return self.get_ijk(*xyz_mean)

    def get_xyz_mean(self):
        '''
        Returns the coordinates (x,y,z) of the grid mean point
        '''
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
        '''
        Returns the grid covariance with respect to the (x,y,z) mean point
        '''
        if self.array is None:
            return None
        xyz_mean = self.get_xyz_mean()
        xx = np.arange(0, self.nx) * self.dx + self.x_orig
        yy = np.arange(0, self.ny) * self.dy + self.y_orig
        zz = np.arange(0, self.nz) * self.dz + self.z_orig
        yarray, xarray, zarray = np.meshgrid(yy, xx, zz)
        array_sum = self.array.sum()
        cov = np.zeros((3,3))
        cov[0,0] = (np.power(xarray, 2) * self.array).sum()/array_sum - (xyz_mean[0] * xyz_mean[0])
        cov[0,1] = cov[1,0] = (xarray * yarray * self.array).sum()/array_sum - (xyz_mean[0] * xyz_mean[1])
        cov[0,2] = cov[2,0] = (xarray * zarray * self.array).sum()/array_sum - (xyz_mean[0] * xyz_mean[2])
        cov[1,1] = (np.power(yarray, 2) * self.array).sum()/array_sum - (xyz_mean[1] * xyz_mean[1])
        cov[1,2] = cov[2,1] = (yarray * zarray * self.array).sum()/array_sum - (xyz_mean[1] * xyz_mean[2])
        cov[2,2] = (np.power(zarray, 2) * self.array).sum()/array_sum - (xyz_mean[2] * xyz_mean[2])
        self.xyz_cov = cov
        return cov

    def get_xyz_ellipsoid(self):
        '''
        Returns the 68% confidence ellipsoid with respect to the (x,y,z) mean point
        '''
        from ellipsoid import Ellipsoid3D
        # The following code is a python translation of the CalcErrorEllipsoid()
        # c-function from the NonLinLoc package, written by Anthony Lomax
        cov = self.get_xyz_cov()
        if cov is None:
            return None

        u, s, v = np.linalg.svd(cov)

        del_chi_2 = 3.53 #3.53: value for 68% conf
        ell = Ellipsoid3D()
        ell.az1 = math.degrees(math.atan2(u[0,0], u[1,0]))
        if ell.az1 < 0.0:
            ell.az1 += 360.0
        ell.dip1 = math.degrees(math.asin(u[2,0]))
        ell.len1 = math.sqrt(del_chi_2) / math.sqrt(1.0 / s[0])
        ell.az2 = math.degrees(math.atan2(u[0,1], u[1,1]))
        if ell.az2 < 0.0:
            ell.az2 += 360.0
        ell.dip2 = math.degrees(math.asin(u[2,1]))
        ell.len2 = math.sqrt(del_chi_2) / math.sqrt(1.0 / s[1]);
        ell.len3 = math.sqrt(del_chi_2) / math.sqrt(1.0 / s[2]);

        self.ellipsoid = ell
        return ell

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

    def get_plot_axes(self, figure=None, ax_xy=None):
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
        plot_cbar_size = 5 #percent
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

        # ax_xz
        ax_xz = divider.append_axes('bottom', size=xz_size, pad=0.05, sharex=ax_xy)
        ax_xz.set_xlim(xmin, xmax)
        ax_xz.set_ylim(zmax, zmin)

        # ax_yz
        ax_yz = divider.append_axes('right', size=yz_size, pad=0.05, sharey=ax_xy)
        plt.setp(ax_yz.get_yticklabels(), visible=False)
        ax_yz.set_xlim(zmin, zmax)
        ax_yz.set_ylim(ymin, ymax)
        plt.setp(ax_yz.get_xticklabels(), rotation=90, fontsize=12)
        plt.setp(ax_yz.get_yticklabels(), rotation=90, fontsize=12)

        # color-bar
        ax_cb = divider.append_axes('bottom', size=cb_size, pad=0.5)

        return ax_xy, ax_xz, ax_yz, ax_cb

    def plot(self, slice_index, handle=False, figure=None, ax_xy=None,
             vmin=None, vmax=None, cmap=None):
        import matplotlib.pyplot as plt
        from matplotlib import ticker

        ax_xy, ax_xz, ax_yz, ax_cb = self.get_plot_axes(figure, ax_xy)
        if figure is None:
            figure = ax_xy.get_figure()

        hnd = ax_xy.imshow(np.transpose(self.array[:, :, slice_index[2]]),
                           vmin=vmin, vmax=vmax, cmap=cmap,
                           origin='lower', extent=self.get_xy_extent())
        ax_xy.set_adjustable('box-forced')
        ax_xz.imshow(np.transpose(self.array[:, slice_index[1], :]),
                     vmin=vmin, vmax=vmax, cmap=cmap,
                     origin='lower', extent=self.get_xz_extent(), aspect='auto')
        ax_xz.set_adjustable('box-forced')
        ax_yz.imshow(self.array[slice_index[0], :, :],
                     vmin=vmin, vmax=vmax, cmap=cmap,
                     origin='lower', extent=self.get_zy_extent(), aspect='auto')
        ax_yz.set_adjustable('box-forced')

        fmt = '%.1e' if self.max() <= 0.01 else '%.2f'
        cb = figure.colorbar(hnd, cax=ax_cb, orientation='horizontal', format=fmt)
        cb.locator = ticker.LinearLocator(numticks=3)
        cb.update_ticks()

        if handle:
            return (ax_xy, ax_xz, ax_yz), cb
        else:
            plt.show()

    def plot_3D_point(self, axes, point, color='r'):
        ax_xy, ax_xz, ax_yz = axes
        ax_xy.scatter(point[0], point[1], color=color)
        ax_xz.scatter(point[0], point[2], color=color)
        ax_yz.scatter(point[2], point[1], color=color)

    def plot_ellipsoid(self, axes, ellipsoid=None, mean_xyz=None):
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
        ax_xz.plot(ell12[:,0], ell12[:,1])
        ax_xz.plot(ell13[:,0], ell13[:,1])
        ax_xz.plot(ell23[:,0], ell23[:,1])

        ell12 = np.array([(vect.y, vect.z) for vect in ellArray12])
        ell13 = np.array([(vect.y, vect.z) for vect in ellArray13])
        ell23 = np.array([(vect.y, vect.z) for vect in ellArray23])
        ax_yz.plot(ell12[:,1], ell12[:,0])
        ax_yz.plot(ell13[:,1], ell13[:,0])
        ax_yz.plot(ell23[:,1], ell23[:,0])

    def copy(self):
        return deepcopy(self)


def main():
    '''
    Test code generating a gaussian grid and computing
    the 3D ellipsoid around the grid mean
    '''
    import matplotlib.pyplot as plt

    #http://stackoverflow.com/questions/17190649/how-to-obtain-a-gaussian-filter-in-python
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
    grd.init_array()
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
