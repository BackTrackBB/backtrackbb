#!/usr/bin/python
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
import matplotlib.pyplot as plt
#
from NLLGrid import NLLGrid


bname1 = 'time/layer.P.QF14.time'
#bname1 = 'time/layer.P.U66B.time'
bname2 = 'time/layer.P.QF16.time'

grid1 = NLLGrid(bname1)
grid2 = NLLGrid(bname2)

extent = grid1.get_xy_extent()
grid1.init_xyz_arrays()
grid2.init_xyz_arrays()
dif_grid1_grid2 = grid1.array-grid2.array

##
fig = plt.figure(figsize=(15, 10))

cmap = 'Spectral'

nz_grid = 10


ax1 = fig.add_subplot(131)

hnd = ax1.imshow(np.flipud(np.transpose(grid1.array[:, :, nz_grid])),
                 extent=extent, vmin=0., vmax=100., cmap=cmap)
cb1 = plt.colorbar(hnd, orientation='horizontal')

hnd = plt.contour(grid1.x_array, grid1.y_array,
                  np.transpose(grid1.array[:, :, nz_grid]), colors='k')
ax1.scatter(grid1.sta_x, grid1.sta_y, marker='v', s=100, linewidths=1, c='r')

ax1.set_xlim(min(grid1.x_array), max(grid1.x_array))
ax1.set_ylim(min(grid1.y_array), max(grid1.y_array))
cb1.set_label('Time[sec]')
##

ax2 = fig.add_subplot(132)
hnd = ax2.imshow(np.flipud(np.transpose(grid2.array[:, :, nz_grid])),
                 extent=extent, vmin=0., vmax=100., cmap=cmap)
cb2 = plt.colorbar(hnd, orientation='horizontal')
##
hnd = plt.contour(grid2.x_array, grid2.y_array,
                  np.transpose(grid2.array[:, :, nz_grid]), colors='k')
ax2.scatter(grid2.sta_x, grid2.sta_y, marker='v', s=100, linewidths=1, c='r')

ax2.set_xlim(min(grid2.x_array), max(grid2.x_array))
ax2.set_ylim(min(grid2.y_array), max(grid2.y_array))
cb2.set_label('Time[sec]')

##

ax3 = fig.add_subplot(133)
hnd = ax3.imshow(np.flipud(np.transpose(dif_grid1_grid2[:, :, nz_grid])),
                 extent=extent, cmap=cmap)
cb3 = plt.colorbar(hnd, orientation='horizontal')
hnd = plt.contour(grid2.x_array, grid2.y_array,
                  np.transpose(dif_grid1_grid2[:, :, nz_grid]), colors='k')
cb3.set_label('Time_difference[sec]')
##
ax3.scatter(grid2.sta_x, grid2.sta_y, marker='v', s=100, linewidths=1, c='r')
ax3.scatter(grid1.sta_x, grid1.sta_y, marker='v', s=100, linewidths=1, c='r')
ax3.set_xlim(min(grid2.x_array), max(grid2.x_array))
ax3.set_ylim(min(grid2.y_array), max(grid2.y_array))
##
##

ax1.set_xlabel('X[km]')
ax2.set_xlabel('X[km]')
ax3.set_xlabel('X[km]')
##
ax1.set_ylabel('Y[km]')
##
ax1.set_title('time grid: ' + grid1.station)
ax2.set_title('time grid: ' + grid2.station)
ax3.set_title('dif. time grid: ' + grid1.station+'-'+grid2.station)
plt.show()
