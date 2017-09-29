# -*- coding: utf-8 -*-
"""


Author:       Brendan Carrick Lyons
Date created: Thu Feb 25 18:17:00 2016
Date edited:  
"""

import numpy as np
import xarray as xr
from numpy import pi
from geofac import geofac

def load_Bres(folder='./',slice=0,phasing=0.,machine='diiid',cur_up=1.,
               cur_low=1.,cur_mid=None,ntor=None,phase=False):
                   
    fup  = folder+'/Bres_upper-'+str(slice)+'.txt'
    flow = folder+'/Bres_lower-'+str(slice)+'.txt'    
    
    if machine is 'diiid':
        conv = 1.0
    elif machine is 'aug':
        conv = -1.0
    elif machine is 'kstar':
        conv = 1.0
        
    fac = geofac(machine=machine,ntor=ntor)
    
#    print fac
        
    cur_up  = fac*cur_up
    cur_low = fac*cur_low
        
    A_up   = np.loadtxt(fup)
    m      = A_up[:,0]
    mag_up = A_up[:,1]
    ph_up  = A_up[:,2]
    Psi    = A_up[:,3]
    q      = A_up[:,4]
    Bru    = mag_up*(np.cos(ph_up) + 1j*np.sin(ph_up))

    A_low   = np.loadtxt(flow)
    mag_low = A_low[:,1]
    ph_low  = A_low[:,2]
    Brl     = mag_low*(np.cos(ph_low) + 1j*np.sin(ph_low))
    
    cur  = xr.DataArray(np.cos(pi*phasing/180.)+conv*1j*np.sin(pi*phasing/180.),
                             {'phasing':phasing})
                             
    Bres_up  = cur_up*xr.DataArray(Bru,[('Psi',Psi)])
    Bres_low = cur_low*xr.DataArray(Brl,[('Psi',Psi)])
    Bres = cur*Bres_up + Bres_low
    
    if cur_mid is not None:
        
        cur_mid = fac*cur_mid
        fmid = folder+'/Bres_middle-'+str(slice)+'.txt'
        A_mid   = np.loadtxt(fmid)
        mag_mid = A_mid[:,1]
        ph_mid  = A_mid[:,2]
        Brm     = mag_mid*(np.cos(ph_mid) + 1j*np.sin(ph_mid))
        cur2 = xr.DataArray(np.cos(pi*phasing/180.)+conv*1j*np.sin(pi*phasing/180.),
                            {'phasing2':phasing})
        Bres_mid = cur_mid*xr.DataArray(Brm,[('Psi',Psi)])
        Bres = Bres + cur2*Bres_mid
    
    if phase:
        Bres.data = np.angle(Bres.data,deg=True) % 360.
    else:
        Bres.data = np.abs(Bres.data)
    q = xr.DataArray(q,[('Psi',Psi)])
    
    return (Bres, q)
