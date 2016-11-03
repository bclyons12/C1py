# -*- coding: utf-8 -*-
"""
calc_ExB

Created on Wed Aug 31 05:48:57 2016

@author: blyons
"""

import numpy as np
from scipy.interpolate import interp1d
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

def calc_ExB(imp='profile_omega',omegas='omegas.txt',Zimp=1,
             output='profile_omega.ExB'):
    
    om_z = np.loadtxt(imp)
    psi_z = om_z[:,0]
    om_z = om_z[:,1]
    
    oms = np.loadtxt(omegas,skiprows=1)
    psi = oms[:,0]
    om_ExB = oms[:,1]
    om_dz  = oms[:,4]/Zimp
    
    f = interp1d(psi,om_dz,bounds_error=False)
    om_dz = f(psi_z)
    
    new_ExB = om_z - om_dz
    
    np.savetxt(output,np.column_stack((psi_z,new_ExB)),delimiter="    ",
               fmt='%1.6f')
               
    
    cols = sns.color_palette(n_colors=6,palette='colorblind')
    mpl.rcParams.update({'font.family': 'serif'})
    f,ax = plt.subplots()
    ax.plot([0,1],[0,0],'k:',linewidth=4)
    h_ExB, = ax.plot(psi,om_ExB,'-',linewidth=4,color=cols[0])
    h_new,   = ax.plot(psi_z,new_ExB  ,'--',linewidth=4,color=cols[2])

    h_ExB.set_label(r'Old')
    h_new.set_label(r'New')    
    
    ax.legend(fontsize=18,ncol=2,loc='lower left')
       
    ax.set_ylim([-200,200])
    ax.set_xlabel(r'$\Psi$',fontsize=24)
    ax.set_ylabel(r'krad/s',fontsize=24)
    ax.tick_params(labelsize=20)
    plt.tight_layout()