#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 18:32:50 2019

@author: blyons
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import griddata

def plot_contour(x, y, V, vmin=None, vmax=None, fs=1.0,
                 xrange=None, yrange=None, xlabel=None, ylabel=None):
    
    if vmin is None:
        vmin = 0.
    if vmax is None:
        vmax = V.max()
    cmap = mpl.colors.ListedColormap(sns.color_palette("inferno", 256))
        
    xl = x[0]
    xu = x[-1]
    yl = y[0]
    yu = y[-1]    
    
    X0,Y0 = np.meshgrid(x,y)
    X0 = X0.ravel()
    Y0 = Y0.ravel()
    V2 = V.ravel()
    
    xi = np.linspace(xl,xu,1025)
    yi = np.linspace(yl,yu,1026)
    
    Xi,Yi = np.meshgrid(xi,yi)
    
    V2 = griddata(X0,Y0,V2,Xi,Yi,method='linear')
    
    f, ax = plt.subplots(figsize=[fs*11.,fs*9.])    
    
    extent = [xl,xu,yl,yu]
    if xrange is None:
        xrange = [xl,xu]
    if yrange is None:
        yrange = [yl,yu]
    aspect = 1.0*(xrange[1]-xrange[0])/(yrange[1]-yrange[0])
    
    
    im=ax.imshow(V2, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax,
                 extent=extent, aspect=aspect, interpolation='nearest')
       
    ax.set_xlim(xrange)
    ax.set_ylim(yrange)
    
    ax.tick_params(labelsize=fs*28)
    if xlabel is not None:
        ax.set_xlabel(xlabel,fontsize=fs*32)
    if ylabel is not None:
        ax.set_ylabel(ylabel,fontsize=fs*32)
    cb = plt.colorbar(im,ax=ax,format='%1.3g')#,ticks=[0.,0.25,0.5,0.75,1.,1.25])
    cb.ax.tick_params(labelsize=fs*28)
    plt.tight_layout()

    return (f,ax)
