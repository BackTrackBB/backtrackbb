import numpy as np
from array import array
#

class NLLGrid():
  '''
    Class for reading NLL grid files
  '''
  nx=int(0)
  ny=int(0)
  nz=int(0)
  x_orig=float(0)
  y_orig=float(0)
  z_orig=float(0)
  dx=float(0)
  dy=float(0)
  dz=float(0)
  proj_name=''
  orig_lat=float(0)
  orig_lon=float(0)
  map_rot=float(0)
  station=''
  sta_x=float(0)
  sta_y=float(0)
  sta_z=float(0)
  array=None
  x_array=None
  y_array=None
  z_array=None

  def __init__(self, basename=''):
    self.basename = basename
    if basename:
      self.read_hdr_file(basename)
      self.read_buf_file(basename)

  def read_hdr_file(self, basename):
    '''Reads header file of NLLoc output'''
    self.basename = basename
    filename = basename + '.hdr'

    # read header file
    f=open(filename)
    lines=f.readlines()
    f.close()
  
    # extract information
    vals=lines[0].split()
    self.nx=int(vals[0])
    self.ny=int(vals[1])
    self.nz=int(vals[2])
    self.x_orig=float(vals[3])
    self.y_orig=float(vals[4])
    self.z_orig=float(vals[5])
    self.dx=float(vals[6])
    self.dy=float(vals[7])
    self.dz=float(vals[8])

    for line in lines:
      if line.split()[0]=='TRANSFORM':
        if line.split()[1]=='NONE': info['proj_name']='TRANS_NONE'
        if line.split()[1]=='SIMPLE': 
          self.proj_name='TRANS_SIMPLE'
          self.orig_lat=float(line.split()[3])
          self.orig_lon=float(line.split()[5])
          self.map_rot=float(line.split()[7])

      if len(line.split())==4:
          self.station=line.split()[0]
          self.sta_x=float(line.split()[1])
          self.sta_y=float(line.split()[2])
          self.sta_z=float(line.split()[3])

  def read_buf_file(self, basename):
    '''reads buf file as a 3d array'''
    self.basename = basename
    buffilename = basename + '.buf'

    f1=open(buffilename,'rb')
    buf=array('f')
    buf.fromfile(f1,
                 self.nx*\
                 self.ny*
                 self.nz)
    f1.close()
    self.array=np.array(buf).reshape(self.nx, self.ny, self.nz)

  def init_xyz_arrays(self):
    self.x_array = np.arange(self.nx)*self.dx+self.x_orig
    self.y_array = np.arange(self.ny)*self.dy+self.y_orig
    self.z_array = np.arange(self.nz)*self.dz+self.z_orig

  def get_extent(self):
      extent = (self.x_orig - self.dx/2,
                self.x_orig + self.nx*self.dx + self.dx/2,
                self.y_orig - self.dy/2,
                self.y_orig + self.ny*self.dy + self.dy/2,
                self.z_orig - self.dz/2,
                self.z_orig + self.nz*self.dz + self.dz/2
                )
      return extent

  def get_xy_extent(self):
      return self.get_extent()[0:4]

  def get_yz_extent(self):
      return self.get_extent()[2:]


      
