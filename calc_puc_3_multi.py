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
from C1py.calc_puc_3 import calc_puc_3

def calc_puc_3_multi(hfsa=['hfsa.txt'], hfsb=['hfsb.txt'], lfs=['lfs.txt'], phasing=0., cur_up=1., cur_low=1.,
                     machine=None, field=None,  ntor=None,
                     label='time',values=[0]):
    
    dBhas, dBhbs, dBls = calc_puc_3(hfsa=hfsa[0], hfsb=hfsb[0], lfs=lfs[0], phasing=phasing,
                                    cur_up=cur_up, cur_low=cur_low,
                                    machine=machine,field=field,ntor=ntor)
    dBhas[label] = values[0]
    dBhbs[label] = values[0]
    dBls[label] = values[0]


    for i,val in enumerate(values):
    
        if i==0:
            continue
        
        dBha, dBhb, dBl = calc_puc_3(hfsa=hfsa[i], hfsb=hfsb[i], lfs=lfs[i], phasing=phasing,
                                   cur_up=cur_up, cur_low=cur_low,
                                   machine=machine,field=field,ntor=ntor)
        dBha[label] = val
        dBhb[label] = val
        dBl[label] = val
        dBhas = xr.concat([dBhas,dBha],label)
        dBhbs = xr.concat([dBhbs,dBhb],label)
        dBls = xr.concat([dBls,dBl],label)
    
    return dBhas, dBhbs, dBls
