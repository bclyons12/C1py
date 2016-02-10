# -*- coding: utf-8 -*-
"""
plot_Bwall

Author:       Brendan Carrick Lyons
Date created: Thu Jan 28 16:25:19 2016
Date edited:  
"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

def plot_Bwall(B1,B2,range=None,title=None,legs=None):
    
    f, ax = plt.subplots(figsize=[12,8])
    
    Z1 = B1.Z.data
    Bw1 = B1.data
    
    Z2 = B2.Z.data
    Bw2 = B2.data
    
    b1r, = plt.plot(Z1,np.real(Bw1),'b',linewidth=3)
    b1i, = plt.plot(Z1,np.imag(Bw1),'b-.',linewidth=3)
    b2r, = plt.plot(Z2,np.real(Bw2),'r',linewidth=3)
    b2i, = plt.plot(Z2,np.imag(Bw2),'r-.',linewidth=3)
    
    if legs is not None:
        b1r.set_label(legs[0]+' - Re{}')
        b1i.set_label(legs[0]+' - Im{}')
        b2r.set_label(legs[1]+' - Re{}')
        b2i.set_label(legs[1]+' - Im{}')
        ax.legend(fontsize=18)
    
    ax.set_xlim([-1.25,1.25])
    if range is not None:
        ax.set_ylim(range)
    
    if title is not None:
        f.suptitle(title,fontsize=22)
    
    ax.set_xlabel('Z',fontsize=20)
    ax.set_ylabel('$B_Z$',fontsize=20)
    ax.tick_params(labelsize=18)
    
    return (f,ax)