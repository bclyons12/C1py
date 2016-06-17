# -*- coding: utf-8 -*-
"""
plot_Bwall

Author:       Brendan Carrick Lyons
Date created: Thu Jan 28 16:25:19 2016
Date edited:  
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_Bwall(Bs,yrange=None,title=None,legs=None,palette=None):
    
    f, ax = plt.subplots(figsize=[12,8])
    sns.set_style('white')
    N = len(Bs)

    if palette is not None:
        cols = sns.color_palette(n_colors=N,palette=palette)
    else:
        cols = sns.husl_palette(n_colors=N,s=1.,l=0.55)
    
    hs = []
    
    for i in range(N):
    
        Z  = Bs[i].Z.data
        Bw = Bs[i].data
    
        hr, = plt.plot(Z,np.real(Bw),'-', color=cols[i],linewidth=4)
        #hi, = plt.plot(Z,np.imag(Bw),'-.',color=cols[i],linewidth=3)
              
        if legs is not None:
            hr.set_label(legs[i])#+' - Re{}')
            #hi.set_label(legs[i]+' - Im{}')
            ax.legend(fontsize=22,loc='best')
            
        hs.append(hr)
        #hs.append(hi)
    
    ax.set_xlim([-1.25,1.25])
    if yrange is not None:
        ax.set_ylim(yrange)
    
    if title is not None:
        ax.set_title(title,fontsize=32)
    
    ax.set_xlabel('Z',fontsize=28)
    ax.set_ylabel(r'HFS Re{$\delta B_{Z,pl}$} (G/kA)',fontsize=28)
    ax.tick_params(labelsize=24)
    #plt.tight_layout()
    
    return (f,ax,hs)
