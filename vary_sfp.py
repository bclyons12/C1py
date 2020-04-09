#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 17:50:50 2017

@author: blyons
"""

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def vary_sfp(a1=1.0, a2=2, a3=1, a4p=0.02, a4v=0.15, a5p=0.02, a5v=0.15,
             a6=0.03, a7=0.003, lc1=100., lc2=100., Wc=0.07, Psic=0.16):

    Psi = np.linspace(0,2,2e3)
    
    h1_i = lambda x: a4p*(1.-np.exp(-np.abs(x/a1-1.)**a2)) + a7
    h2_i = lambda x: a5p*(1.-np.exp(-np.abs(x/a1-1.)**a2)) + a6
    h1_e = lambda x: a4v*(1.-np.exp(-np.abs(x/a1-1.)**a3)) + a7
    h2_e = lambda x: a5v*(1.-np.exp(-np.abs(x/a1-1.)**a3)) + a6
    
    h1 = np.piecewise(Psi,[Psi<=a1,Psi>a1],[h1_i,h1_e])
    h2 = np.piecewise(Psi,[Psi<=a1,Psi>a1],[h2_i,h2_e])
    
    o1 = lc1*(1.+((Psi-Psic)/Wc)**2)
    o2 = lc2*(1.+((Psi-Psic)/Wc)**2)
    
    h1 = (h1**-1. + o1**-1.)**-1.
    h2 = (h2**-1. + o2**-1.)**-1.
    
    # calc defaults
    
    a1, a2, a3, a4p, a4v, a5p, a5v = (1.0, 2, 1, 0.02, 0.15, 0.02, 0.15)
    a6, a7, lc1, lc2, Wc, Psic     = (0.03, 0.003, 100., 100., 0.07, 0.16)
    
    d1_i = lambda x: a4p*(1.-np.exp(-np.abs(x/a1-1.)**a2)) + a7
    d2_i = lambda x: a5p*(1.-np.exp(-np.abs(x/a1-1.)**a2)) + a6
    d1_e = lambda x: a4v*(1.-np.exp(-np.abs(x/a1-1.)**a3)) + a7
    d2_e = lambda x: a5v*(1.-np.exp(-np.abs(x/a1-1.)**a3)) + a6
    
    d1 = np.piecewise(Psi,[Psi<=a1,Psi>a1],[d1_i,d1_e])
    d2 = np.piecewise(Psi,[Psi<=a1,Psi>a1],[d2_i,d2_e])
    
    o1 = lc1*(1.+((Psi-Psic)/Wc)**2)
    o2 = lc2*(1.+((Psi-Psic)/Wc)**2)
    
    d1 = (d1**-1. + o1**-1.)**-1.
    d2 = (d2**-1. + o2**-1.)**-1.
    
    cols = sns.color_palette(n_colors=6,palette='colorblind')
    
    f,ax = plt.subplots()
    d, = ax.plot(Psi,1./d1,'k--',linewidth=3)
    h, = ax.plot(Psi,1./h1,'-',linewidth=3,color=cols[2])
    d.set_label('Default')
    h.set_label('Modified')
    ax.legend(loc='upper center')
    ax.set_xlabel(r'$\Psi$',fontsize=16)
    ax.set_ylabel(r'$h_n$',fontsize=16)
#    ax.set_ylim([0,a4v+a7])
    ax.set_title('Normal mesh size')
    f.tight_layout()
    
    f,ax = plt.subplots()
    d, = ax.plot(Psi,1./d2,'k--',linewidth=3)
    h, = ax.plot(Psi,1./h2,'-',linewidth=3,color=cols[2])
    d.set_label('Default')
    h.set_label('Modified')
    ax.legend(loc='upper center')
    ax.set_xlabel(r'$\Psi$',fontsize=16)
    ax.set_ylabel(r'$h_ts$',fontsize=16)
#    ax.set_ylim([0,a5v+a6])
    ax.set_title('Tangential mesh size')
    f.tight_layout()
    
    return