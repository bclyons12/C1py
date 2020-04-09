# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 16:02:41 2016

@author: lyonsbc
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import fio_py
from C1py.read_field import read_field
from C1py.get_wall import get_wall
import seaborn as sns

def plot_shape(folder='./', rrange=None, zrange=None, bound=False, ax=None,
               legend=True, fs=1.0, Nlvl_in=10, Nlvl_out=1, linewidth=3, 
               title=None):
    
    if ax is None:
        sns.set_style('white')
        f, ax = plt.subplots(figsize=[fs*8,fs*12])
    else:
        f = None
    
    if isinstance(folder,str):
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
        
        levels = (np.arange(Nlvl_in+Nlvl_out+1.)/Nlvl_in)
        levels = levels*(psi_lcfs-psi_axis)+psi_axis
        if psi_lcfs < psi_axis:
            levels = np.flipud(levels)
        
        psi = read_field('psi',slice=-1,filename=fn,rrange=rrange,
                         zrange=zrange)

        col = mpl.colors.rgb2hex(cols[i-1])
        if i==0:
            hold='off'
        else:
            hold='on'
            
        ax.contour(psi.R,psi.Z,psi.data.T,levels,hold=hold,colors=col,
                        linewidths=linewidth,linestyles='solid')
                        
        if legend:
            h, = ax.plot([np.inf,np.inf],[np.inf,np.inf],'-',color=col,
                         linewidth=linewidth)
            
            h.set_label(folder[i])
    
    if bound:
        (Wi,Wo) = get_wall()
        ax.plot(Wi[:,0],Wi[:,1],'r-',linewidth=1)
        ax.plot(Wo[:,0],Wo[:,1],'r-',linewidth=1)
        
    
    if rrange is None:
        rrange = [min(psi.R),max(psi.R)]
    if zrange is  None:
        zrange = [min(psi.Z),max(psi.Z)]
        
    ax.set_xlim(rrange)
    ax.set_ylim(zrange)
    ax.set_xlabel(r'$R$ (m)',fontsize=28*fs)
    ax.set_ylabel(r'$Z$ (m)',fontsize=28*fs)
    ax.tick_params(labelsize=24*fs)
    if title is not None:
        ax.set_title(title,fontsize=32*fs)
    if legend:
        ax.legend(fontsize=24*fs)
    if f is not None:
        plt.tight_layout()
    
    return
