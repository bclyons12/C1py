# -*- coding: utf-8 -*-
"""
plot_Bmns

Plot Bmn vs Psi and m, overlayed with the m=nq resonance line from a Bmns Dataset

Author:       Brendan Carrick Lyons
Date created: Thu Jan 21 12:53:28 2016
Date edited:  
"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from load_Bmns import load_Bmns

def plot_Bmns(folder='./', phasing=0., slice=0, cur_up=1., cur_low=1.,
              machine=None, iplasma=None,  ntor=None, cmax=None, 
              interp='nearest'):
    
    (Bc1, Bmars, Bipec) = load_Bmns(folder=folder,phasing=phasing,slice=slice, 
                                    cur_up=cur_up,cur_low=cur_low,code='all',
                                    machine=machine,iplasma=iplasma,ntor=ntor)
                                    
    f, axes = plt.subplots(2,2,figsize=(13,12))
    
    
    #f.suptitle(folder)
    
    ax = axes[0,0]
    im = plot_Bmn(Bc1, phasing=phasing, cmax=cmax, interp=interp, ax=ax,
                  uniform=True, solo=False, title=r'M3D-$C^1$')
    
    ax = axes[0,1]
    im = plot_Bmn(Bmars, phasing=phasing, cmax=cmax, interp=interp, ax=ax,
                  solo=False, title='MARS')

    ax = axes[1,0]
    im = plot_Bmn(Bipec, phasing=phasing, cmax=cmax, interp=interp, ax=ax,
                  solo=False, title='IPEC')
    
    ax = axes[1,1]
    f.colorbar(im,cax=ax,orientation='horizontal')
    ax.set_aspect(0.1)
    
    for (i,j), ax in np.ndenumerate(axes):
        ax.tick_params(labelsize=18)
        if ax is not axes[1,1]:
            ax.set_xlabel('m',fontsize=20)
            ax.set_ylabel(r'$\Psi$',fontsize=20)
        else:
            ax.set_title(r'|$B_{mn}$| (G)',fontsize=22)
            
    plt.tight_layout()
    
    
    
    return f, axes
                  
def plot_Bmn(Bmns,phasing=0.,cmax=None,interp='nearest',ax=None,uniform=False,
             solo=True,title=r'|$B_{mn}$| (G)',xrange=None,yrange=None):
    
    p0 = Bmns.Psi.data
    m0 = Bmns.m.data
    
    if Bmns.phasing.size > 1:
        B0 = Bmns.sel(phasing=phasing).Bmn.data
    elif Bmns.phasing.data == phasing:
        B0 = Bmns.Bmn.data
    else:
        print("Error: phasing does not correspond to Bmns.phasing")
        return None
    
    q    = Bmns.q.data
    ntor = Bmns.attrs['ntor']
    
    
    vmin = 0.
    if cmax is None:
        vmax = B0.max()
    else:
        vmax = cmax
        
    if xrange is None:
        xrange = [-10*abs(ntor),10*abs(ntor)]
    if yrange is None:
        yrange = [0,1]
    
    ml = m0[0]
    mu = m0[-1]
    pl = p0[0]
    pu = p0[-1]    
    
    if not uniform:
        # interpolate onto uniform grid
            
        M0,P0 = np.meshgrid(m0,p0)
        M0 = M0.ravel()
        P0 = P0.ravel()
        B0 = B0.ravel()
        
        mi = np.linspace(ml,mu,len(m0))
        pi = np.linspace(pl,pu,len(p0))
        
        Mi,Pi = np.meshgrid(mi,pi)
    
        Bmn = griddata(M0,P0,B0,Mi,Pi,interp='linear')
    
    else:
        Bmn = B0
    
    if ax is None:
        f, ax = plt.subplots(figsize=[12,9])
    
    extent = [ml,mu,pl,pu]
    
    aspect = 0.9*(xrange[1]-xrange[0])/(yrange[1]-yrange[0])
    
    im=ax.imshow(Bmn, origin='lower', vmin=vmin, vmax=vmax,
                 extent=extent, cmap='jet',aspect=aspect,
                 interpolation=interp)
   
    ax.set_xlim(xrange)
    ax.set_ylim(yrange)
    
    ax.set_title(title,fontsize=24)
    if solo:
        ax.tick_params(labelsize=18)
        ax.set_xlabel('m',fontsize=22)
        ax.set_ylabel(r'$\Psi$',fontsize=22)
        cb = plt.colorbar(im,ax=ax,format='%1.3g')
        cb.ax.tick_params(labelsize=18)
        plt.tight_layout()
    
    ax.plot(ntor*q,p0,'w--',linewidth=3)
    
    return im
