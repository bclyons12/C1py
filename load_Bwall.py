# -*- coding: utf-8 -*-
"""
load_Bwall

Author:       Brendan Carrick Lyons
Date created: Thu Jan 28 13:54:22 2016
Date edited:  
"""
import numpy as np
import xarray as xr
from numpy import pi
from C1py.geofac import geofac

def load_Bwall(folder='./',slice=0,phasing=0.,machine='diiid',cur_up=1.,
               cur_low=1.,ntor=None):
                   
    file = folder+'/Bwall-'+str(slice)+'.txt'
    
    if machine is 'diiid':
        conv = 1.0
    elif machine is 'aug':
        conv = -1.0
        
    fac = geofac(machine=machine,ntor=ntor)
    
    cur_up  = fac*cur_up
    cur_low = fac*cur_low
        
    A  = np.loadtxt(file)
    Z  = A[:,0]
    Bu = A[:,1] + 1j*A[:,2]
    Bl = A[:,3] + 1j*A[:,4]
    
    cur = xr.DataArray(np.cos(pi*phasing/180.)+conv*1j*np.sin(pi*phasing/180.),
                       coords=[('phasing',phasing)])
                             
    Bwu = cur_up*xr.DataArray(Bu,[('Z',Z)])
    Bwl = cur_low*xr.DataArray(Bl,[('Z',Z)])
    Bwall = cur*Bwu + Bwl
    
    return Bwall
