#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 16:14:55 2017

@author: blyons
"""

import numpy as np
import xarray as xr
from numpy import pi
from C1py.geofac import geofac
from C1py.calc_puc import calc_puc

def calc_puc_multi(hfs=['hfs.txt'], lfs=['lfs.txt'], phasing=0., cur_up=1., cur_low=1.,
             machine=None, field=None,  ntor=None,
             label='time',values=[0]):
    
    
    dBhs, dBls = calc_puc(hfs=hfs[0], lfs=lfs[0], phasing=phasing,
                          cur_up=cur_up, cur_low=cur_low,
                          machine=machine,field=field,ntor=ntor)
    dBhs[label] = values[0]
    dBls[label] = values[0]


    for i,val in enumerate(values):
    
        if i==0:
            continue
        
        dBh, dBl = calc_puc(hfs=hfs[i], lfs=lfs[i], phasing=phasing,
                            cur_up=cur_up, cur_low=cur_low,
                            machine=machine,field=field,ntor=ntor)
        dBh[label] = val
        dBl[label] = val
        dBhs = xr.concat([dBhs,dBh],label)
        dBls = xr.concat([dBls,dBl],label)
    
    return dBhs, dBls
