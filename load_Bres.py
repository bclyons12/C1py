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
               cur_low=1.,ntor=None):
                   
    fup  = folder+'/Bres_upper-'+str(slice)+'.txt'
    flow = folder+'/Bres_lower-'+str(slice)+'.txt'    
    
    if machine is 'diiid':
        conv = 1.0
    elif machine is 'aug':
        conv = -1.0
        
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
    
    
    cur = xr.DataArray(np.cos(pi*phasing/180.)+conv*1j*np.sin(pi*phasing/180.),
                             {'phasing':phasing})
                             
    Bres_up  = cur_up*xr.DataArray(Bru,[('Psi',Psi)])
    Bres_low = cur_low*xr.DataArray(Brl,[('Psi',Psi)])
    Bres = np.abs(cur*Bres_up + Bres_low)
    q = xr.DataArray(q,[('Psi',Psi)])
    
    return (Bres, q)