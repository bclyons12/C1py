# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 16:02:41 2016

@author: lyonsbc
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import fio_py
from read_field import read_field
from get_wall import get_wall
import seaborn as sns
sns.set_style('white')

def plot_shape(folder='./', rrange=None, zrange=None,bound=None):

    f, ax = plt.subplots(figsize=[8,12])
    
    if isinstance(folder,basestring):
           folder = [folder]
        
    Nf = len(folder)
    if Nf <=6:
        cols = sns.color_palette(n_colors=Nf,palette='colorblind')
    else:
        cols = sns.color_palette(n_colors=Nf,palette='inferno')
    
    cols[-1] = (0.,0.,0.)    
    
    for i in range(Nf):
        
        fn = folder[i]+'/C1.h5'

        isrc = fio_py.open_source(fio_py.FIO_M3DC1_SOURCE,fn)
        fio_py.get_options(isrc)
        
        ipsi_axis = fio_py.get_series(isrc, fio_py.FIO_MAGAXIS_PSI)
        ipsi_lcfs = fio_py.get_series(isrc, fio_py.FIO_LCFS_PSI)
    
        psi_axis = fio_py.eval_series(ipsi_axis, 0.)
        psi_lcfs = fio_py.eval_series(ipsi_lcfs, 0.)
    
        fio_py.close_series(ipsi_axis)
        fio_py.close_series(ipsi_lcfs)

#        print([psi_axis,psi_lcfs])
        
        levels = (np.arange(12.)/10.)*(psi_lcfs - psi_axis) + psi_axis
        if psi_lcfs < psi_axis:
            levels = np.flipud(levels)
        
        psi = read_field('psi',slice=-1,filename=fn,rrange=rrange,
                         zrange=zrange)

        col = mpl.colors.rgb2hex(cols[i-1])
        if i==0:
            ax.contour(psi.R,psi.Z,psi.data.T,levels,colors=col,
                            linewidths=3,linestyles='solid')
            h, = ax.plot([np.inf,np.inf],[np.inf,np.inf],'-',color=col,
                         linewidth=3)
        else:
            ax.contour(psi.R,psi.Z,psi.data.T,levels,hold='on',colors=col,
                            linewidths=3,linestyles='solid')
            h, = ax.plot([np.inf,np.inf],[np.inf,np.inf],'-',color=col,
                         linewidth=3)
        
        h.set_label(folder[i])
    
    if bound is not None:
        (Wi,Wo) = get_wall()
        ax.plot(Wi[:,0],Wi[:,1],'k--',linewidth=1)
        ax.plot(Wo[:,0],Wo[:,1],'k--',linewidth=1)
        
    
    if rrange is None:
        rrange = [min(psi.R),max(psi.R)]
    if zrange is  None:
        zrange = [min(psi.Z),max(psi.Z)]
        
    ax.set_xlim(rrange)
    ax.set_ylim(zrange)
    ax.set_xlabel('R',fontsize=24)
    ax.set_ylabel('Z',fontsize=24)
    ax.tick_params(labelsize=20)
    ax.legend(fontsize=16)
    plt.tight_layout()
    
    return
