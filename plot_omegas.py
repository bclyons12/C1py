# -*- coding: utf-8 -*-
"""
plot_omegas

Created on Wed Aug 31 05:48:57 2016

@author: blyons
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

def plot_omegas(file='omegas.txt',Zimp=None,xrange=None,yrange=None,
                title=r'Rotation profiles',qs=None,fluids=2):
    
    oms = np.loadtxt(file,skiprows=1)
    psi = oms[:,0]
    
    if Zimp is None:
        
        # Rotations all correspond to main ion
        
        om_ExB = oms[:,1]
        om_i   = oms[:,2]
        om_e   = oms[:,3]
        om_di  = oms[:,4]
        
        om_de  = om_e - om_ExB
        om_z   = 0.*psi
        om_dz  = 0.*psi
        
    else:
        
        # Rotations all correspond to impurity ion
        
        om_z   = oms[:,2]
        om_dz  = oms[:,4]/Zimp
        om_di  = oms[:,4] 
        om_de  = oms[:,3] - oms[:,1]
        
        om_ExB = om_z - om_dz
        om_i   = om_ExB + om_di
        om_e   = om_ExB + om_de
        
    cols = sns.color_palette(n_colors=6,palette='colorblind')
    mpl.rcParams.update({'font.family': 'serif'})
    f,ax = plt.subplots()
    ax.plot([0,1],[0,0],'k:',linewidth=4)
    if qs is not None:
        for P in qs:
            ax.plot([P,P],yrange,'k--',linewidth=2)
    if fluids == 2:
        h_i,   = ax.plot(psi,om_i  ,'-',linewidth=4,color=cols[2])
        h_di,  = ax.plot(psi,om_di ,'-',linewidth=4,color=cols[3])
        h_ExB, = ax.plot(psi,om_ExB,'-',linewidth=4,color=cols[0])
        h_e,   = ax.plot(psi,om_e  ,'-',linewidth=4,color=cols[1])

        h_i.set_label(r'$\omega_i$')
        h_di.set_label(r'$\omega_{*i}$')
        h_ExB.set_label(r'$\omega_{E\times B}$')
        h_e.set_label(r'$\omega_e$')    
    
    elif fluids == 1:
        
        h_ExB, = ax.plot(psi,om_ExB,'-',linewidth=4,color=cols[0])
        h_ExB.set_label(r'$\omega_{E\times B}$')
        
    if Zimp is None:
        ax.legend(fontsize=18,ncol=2,loc='lower left')
    else:
        if fluids == 2:
            h_z,   = ax.plot(psi,om_z  ,'--',linewidth=4,color=cols[2])
            h_dz,  = ax.plot(psi,om_dz ,'--',linewidth=4,color=cols[3])
            h_z.set_label(r'$\omega_z$')
            h_dz.set_label(r'$\omega_{*z}$')
        
        ax.legend(fontsize=20,ncol=3,loc='lower left')
        
    if xrange is not None:
        ax.set_xlim(xrange)
    if yrange is not None:
        ax.set_ylim(yrange)
    ax.set_xlabel(r'$\Psi$',fontsize=24)
    ax.set_ylabel(r'krad/s',fontsize=24)
    ax.tick_params(labelsize=20)
    ax.set_title(title,fontsize=28)
    plt.tight_layout()
    
    return f,ax