# -*- coding: utf-8 -*-
"""
calc_puc

Calculate magnetic signal versus phasing

Author:       Brendan Carrick Lyons
Date created: Thu Jan 21 12:53:28 2016
Date edited:  
"""
import numpy as np
import xarray as xr
from numpy import pi
from geofac import geofac

def calc_puc(hfs='hfs.txt', lfs='lfs.txt', phasing=0., cur_up=1., cur_low=1.,
              machine=None, field=None,  ntor=None):
                  
    if machine is 'diiid':
        conv = 1.0
    elif machine is 'aug':
        conv = -1.0
    elif machine is 'kstar':
        conv = 1.0
        
    phasing = np.asarray(phasing)
    
    fac = geofac(machine=machine,ntor=ntor)            
    cur_up  = fac*cur_up
    cur_low = fac*cur_low
    
        
    dBz    = np.loadtxt(hfs)
    dBh    = dBz[:,0] + 1j*dBz[:,1]
    dBh_tu = dBh[0]
    dBh_tl = dBh[1]
    dBh_vu = dBh[2]
    dBh_vl = dBh[3]
    
    dBz    = np.loadtxt(lfs)
    dBl    = dBz[:,0] + 1j*dBz[:,1]
    dBl_tu = dBl[0]
    dBl_tl = dBl[1]
    dBl_vu = dBl[2]
    dBl_vl = dBl[3]
    
    cur = xr.DataArray(np.cos(pi*phasing/180.)+conv*1j*np.sin(pi*phasing/180.),
                             {'phasing':phasing})
                             
    if   field == 0:
        dBl = cur*cur_up*dBl_vu + cur_low*dBl_vl
        dBh = cur*cur_up*dBh_vu + cur_low*dBh_vl
    elif field == 1:
        dBl = cur*cur_up*(dBl_tu - dBl_vu) + cur_low*(dBl_tl - dBl_vl)
        dBh = cur*cur_up*(dBh_tu - dBh_vu) + cur_low*(dBh_tl - dBh_vl)
    elif field == 2:
        dBl = cur*cur_up*dBl_tu + cur_low*dBl_tl
        dBh = cur*cur_up*dBh_tu + cur_low*dBh_tl
    else:
        print("Error: field = " + str(field) + " unacceptable")
        return None
        
    
    return (dBh, dBl)