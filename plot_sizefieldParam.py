# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 16:02:41 2016

@author: lyonsbc
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import copy

def plot_sizefieldParam(a1=1.0, a2=2., a3=1., a4p=0.02, a4v=0.15, a5p=0.02,
                        a5v=0.15, a6=0.03, a7=0.003, 
                        lc1=100., lc2=100.0, Wc=0.07, Pc=0.16, 
                        qfile=None, over=False, ax = None):
    
    psii = np.linspace(0,a1,1024)
    psie = np.linspace(a1,2.0,1024)
    
    psi = copy.deepcopy(psii)
    
    if qfile is not None:
        q = np.loadtxt(qfile)
        f = interp1d(q[:,0],q[:,1])
        q = f(psii)
        ms = np.array([1.,2.,3.])
        for m in ms:
            tmp =  (q - m)/(q*n*0.05)
            psi = psi + 1.0/(tmp**2 + 1)
            psi[where(psi>1.)] = 1.
    
    h1 = a4p*(1-np.exp(-np.abs(psi/a1 - 1)**a2)) + a7
    h1e = a4v*(1-np.exp(-np.abs(psie/a1 - 1)**a3)) + a7
    h2 = a5p*(1-np.exp(-np.abs(psi/a1 - 1)**a2)) + a6
    h2e = a5v*(1-np.exp(-np.abs(psie/a1 - 1)**a3)) + a6
    
    h1 = 1./(1./h1 + (1./lc1)*(1./(1.+((psi-Pc)/Wc)**2)))
    h1e = 1./(1./h1e + (1./lc1)*(1./(1.+((psie-Pc)/Wc)**2)))
    h2 = 1./(1./h2 + (1./lc2)*(1./(1.+((psi-Pc)/Wc)**2)))
    h2e = 1./(1./h2e + (1./lc2)*(1./(1.+((psie-Pc)/Wc)**2)))
    
    if ax is None:
        if not over:
            f,ax = plt.subplots()
        else:
            ax = plt.gca()
    cp = {False:'r', True:'m'}
    ct = {False:'b', True:'c'}
    ax.plot(psi, h1, cp[over], lw=3)
    ax.plot(psie, h1e, cp[over], lw=3)
    ax.plot(psi, h2, ct[over], lw=3)
    ax.plot(psie, h2e, ct[over], lw=3)
    if not over:
        ax.set_ylim([0,np.max(h2e)])
    print(a1, a2, a3, a4p, a4v, a5p, a5v, a6, a7, lc1, lc2, Wc, Pc)
    return ax
    
