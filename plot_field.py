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
from mpl_toolkits.axes_grid1 import make_axes_locatable
import fio_py
from C1py.read_field import read_field
from C1py.get_wall import get_wall
import seaborn as sns
sns.set_style('white')

def plot_field(field,filename='C1.h5', points=200,  slice=0, 
               range=None,xrange=None,yrange=None, rrange=None, zrange=None, 
               palette=None, lcfs=None, bound=None, fac=1., linfac=1., phi=0.,
               iabs=None, iphase=None, isum=None, iavg=None, idiff=None,
               ilinear=None, iequil=None, icomplex=None, ntor=None,
               title=None,fs=1.0,ax=None,symrange=False,cb_label=None,
               minval=None, maxval=None, nimrod=False,make_cb=True):

    if isinstance(field,str):
        # Read this field
        if title is None:
            title = field
        field = read_field(field, slice=slice, filename=filename, phi=phi,
                           points=points, rrange=rrange, zrange=zrange,
                           linfac=linfac, iabs=iabs, iphase=iphase, 
                           isum=isum, iavg=iavg, idiff=idiff, ilinear=ilinear, 
                           iequil=iequil, icomplex=icomplex, ntor=ntor,nimrod=nimrod)
                        
    data = np.nan_to_num(fac*np.real(field.data))
        
    if range is None:
        if symrange:
            vmin = -abs(data).max()
            vmax =  abs(data).max()
        else:
            vmin = data.min()
            vmax = data.max()
    else:
        vmin = range[0]
        vmax = range[1]
    if minval is not None:
        vmin = minval
    if maxval is not None:
        vmax = maxval

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

    if ax is None:
        f, ax = plt.subplots(figsize=[fs*9,fs*12])
    else:
        f = None
    
    im=ax.imshow(data.T, origin='lower', vmin=vmin, vmax=vmax,
                 extent=extent, cmap=cmap)
                     
    if xrange is None:
        xrange = [field.R.data[0],field.R.data[-1]]
    if yrange is None:
        yrange = [field.Z.data[0],field.Z.data[-1]]
        
    ax.set_xlim(xrange)
    ax.set_xlabel(r'$R$ (m)',fontsize=fs*28)
    ax.set_ylim(yrange)
    ax.set_ylabel(r'$Z$ (m)',fontsize=fs*28)
    if title is not None:
        ax.set_title(title,fontsize=fs*32)
    ax.tick_params(labelsize=fs*24)
    
    if make_cb:
        if f is not None:
            cb = f.colorbar(im,format='%1.3g')
            cb.ax.tick_params(labelsize=fs*24)
        else:
            div = make_axes_locatable(ax)
            cax = div.append_axes("right",size="10%",pad=0.02)
            cb = plt.colorbar(im,cax=cax,format='%1.3g')
            cb.ax.tick_params(labelsize=fs*24)
        if cb_label is not None:
            cb.ax.get_yaxis().labelpad=fs*36
            cb.ax.set_ylabel(cb_label,rotation=270,fontsize=fs*24)
    else:
        cb = None

    if lcfs is not None:
        if not isinstance(filename,str):
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

    if isinstance(bound,tuple):
        R0, a, delta, Z0, b = bound
        theta = np.linspace(-np.pi/2.,3*np.pi/2,1000)
        R = R0 + a*np.cos(theta + delta*np.sin(theta))
        Z = Z0 + b*np.sin(theta)
        ax.plot(R,Z,'-',color=col_lcfs,linewidth=fs*3)
    elif bound is not None:
        (Wi,Wo) = get_wall()
        ax.plot(Wi[:,0],Wi[:,1],'--',color=col_lcfs,linewidth=1)
        ax.plot(Wo[:,0],Wo[:,1],'--',color=col_lcfs,linewidth=1)
    
    if f is not None:
        f.tight_layout()
    else:
        plt.tight_layout()

    return (f, ax, im, cb)

