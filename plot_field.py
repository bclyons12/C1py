# -*- coding: utf-8 -*-
"""
plot_field

Plot fields from a C1.h5 file

M3D-C1 suite

Brendan Carrick Lyons
12/23/2015
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import fio_py
from read_field import read_field
from get_wall import get_wall
import seaborn as sns
sns.set_style('white')

def plot_field(field,filename='C1.h5', points=200,  slice=0, 
               range=None,xrange=None,yrange=None, rrange=None, zrange=None, 
               palette=None, lcfs=None, bound=None, linfac=1., phi=0.,
               iabs=None, iphase=None, isum=None, iavg=None, idiff=None,
               ilinear=None, iequil=None, icomplex=None, ntor=None):

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
        

    if palette is None:
        if vmin >= 0.:
            palette = 'inferno'
            col_lcfs = 'w'
        else:
            palette = 'RdBu_r'
            col_lcfs = 'k'
    else:
        col_lcfs = 'w'
    
    cmap = mpl.colors.ListedColormap(sns.color_palette(palette, 256))
    
    extent = [field.R.data[0],field.R.data[-1],
              field.Z.data[0],field.Z.data[-1]]

    f, ax = plt.subplots(figsize=[8,12])
    
    im=ax.imshow(data.T, origin='lower', vmin=vmin, vmax=vmax,
                 extent=extent, cmap=cmap)
                     
    cb = f.colorbar(im,format='%1.3g')
    
    if xrange is None:
        xrange = [field.R.data[0],field.R.data[-1]]
    if yrange is None:
        yrange = [field.Z.data[0],field.Z.data[-1]]
        
    ax.set_xlim(xrange)
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
                   extent=extent,colors=col_lcfs,linewidths=1)

    if bound is not None:
        (Wi,Wo) = get_wall()
        ax.plot(Wi[:,0],Wi[:,1],'--',color=col_lcfs,linewidth=1)
        ax.plot(Wo[:,0],Wo[:,1],'--',color=col_lcfs,linewidth=1)
    
    return (f, ax, im, cb)

