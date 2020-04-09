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
from C1py.geofac import geofac

def calc_puc_3(hfsa='hfsa.txt', hfsb='hfsb.txt', lfs='lfs.txt', phasing=0., cur_up=1., cur_low=1.,
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
    
        
    dBz    = np.loadtxt(hfsa)
    dBha    = dBz[:,0] + 1j*dBz[:,1]
    dBha_tu = dBha[0]
    dBha_tl = dBha[1]
    dBha_vu = dBha[2]
    dBha_vl = dBha[3]
    
    dBz    = np.loadtxt(hfsb)
    dBhb    = dBz[:,0] + 1j*dBz[:,1]
    dBhb_tu = dBhb[0]
    dBhb_tl = dBhb[1]
    dBhb_vu = dBhb[2]
    dBhb_vl = dBhb[3]
    
    dBz    = np.loadtxt(lfs)
    dBl    = dBz[:,0] + 1j*dBz[:,1]
    dBl_tu = dBl[0]
    dBl_tl = dBl[1]
    dBl_vu = dBl[2]
    dBl_vl = dBl[3]
    
    cur = xr.DataArray(np.cos(pi*phasing/180.)+conv*1j*np.sin(pi*phasing/180.),
                       coords=[('phasing',phasing)])
                             
    if   field == 0:
        dBl = cur*cur_up*dBl_vu   + cur_low*dBl_vl
        dBha = cur*cur_up*dBha_vu + cur_low*dBha_vl
        dBhb = cur*cur_up*dBhb_vu + cur_low*dBhb_vl
    elif field == 1:
        dBl = cur*cur_up*(dBl_tu - dBl_vu)    + cur_low*(dBl_tl - dBl_vl)
        dBha = cur*cur_up*(dBha_tu - dBha_vu) + cur_low*(dBha_tl - dBha_vl)
        dBhb = cur*cur_up*(dBhb_tu - dBhb_vu) + cur_low*(dBhb_tl - dBhb_vl)
    elif field == 2:
        dBl = cur*cur_up*dBl_tu   + cur_low*dBl_tl
        dBha = cur*cur_up*dBha_tu + cur_low*dBha_tl
        dBhb = cur*cur_up*dBhb_tu + cur_low*dBhb_tl
    else:
        print("Error: field = " + str(field) + " unacceptable")
        return None
        
    
    return (dBha, dBhb, dBl)
