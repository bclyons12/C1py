# -*- coding: utf-8 -*-
"""
plot_field

Plot fields from a C1.h5 file

M3D-C1 suite

Brendan Carrick Lyons
12/23/2015
"""
import numpy as np
import matplotlib.pyplot as mpl
import fio_py
from read_field import read_field
#from C1field import C1field


def plot_field(field,range=None,xrange=None,yrange=None,cmap='jet',lcfs=None,
               slice=0, filename='C1.h5', points=200, phi=0.,
               rrange=None, zrange=None, linfac=1., iabs=None, iphase=None,
               isum=None, iavg=None, idiff=None,ilinear=None, iequil=None, 
               icomplex=None, ntor=None):

    if isinstance(field,basestring):
        # Read this field
        field = read_field(field, slice=slice, filename=filename, phi=phi,
                           points=points, rrange=rrange, zrange=zrange,
                           linfac=linfac, iabs=iabs, iphase=iphase, 
                           isum=isum, iavg=iavg, idiff=idiff, ilinear=ilinear, 
                           iequil=iequil, icomplex=icomplex, ntor=ntor)
                        
    data = np.real(field.data)
        
    if range is None:
        vmin = data.min()
        vmax = data.max()
    else:
        vmin = range[0]
        vmax = range[1]
    
    extent = [field.R.data[0],field.R.data[-1],
              field.Z.data[0],field.Z.data[-1]]

    f, ax = mpl.subplots()
    
    im=ax.imshow(data.T, origin='lower', vmin=vmin, vmax=vmax,
                 extent=extent, cmap=cmap)
                     
    cb = f.colorbar(im,format='%1.3g')
    
    if xrange is not None:
        ax.set_xlim(xrange)
    if yrange is not None:
        ax.set_ylim(yrange)
        
    if lcfs is not None:
        if not isinstance(filename,basestring):
            filename = filename[0]
            
        isrc = fio_py.open_source(fio_py.FIO_M3DC1_SOURCE,filename)
        fio_py.get_options(isrc)
        
        psi = read_field('psi', slice=-1, filename=filename, points=points,
                         rrange=rrange, zrange=zrange, iequil=1)
                         
        ipsi_lcfs = fio_py.get_series(isrc, fio_py.FIO_LCFS_PSI)

        psi_lcfs = fio_py.eval_series(ipsi_lcfs, 0.)

        fio_py.close_series(ipsi_lcfs)
        
        ax.contour(psi.data.T,[psi_lcfs],hold='on',origin='lower',
                   extent=extent,colors='w',linewidths=1)

        
    
    return (f, ax, im, cb)

