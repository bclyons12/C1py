# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 15:02:21 2016

@author: lyonsbc
"""

import xarray as xr
from C1py.load_Bres import load_Bres

def multi_Bres(cases,index,head='./',foot='/',slice=0,phasing=0.,
               machine='diiid',cur_up=1.,cur_low=1.,ntor=None,phase=False,Jres=False):
                  

    key  = index.keys()[0]
    vals = index[key]                  

     
    fout = head+'/'+cases[0]+'/'+foot
    (Brs,q) = load_Bres(folder=fout,slice=slice,phasing=phasing,ntor=ntor,
                        machine=machine,cur_up=cur_up,cur_low=cur_low,phase=phase,Jres=Jres)
    Brs[key]=vals[0]
    
    N = len(cases)
    
    for i in range(1,N):
        case = cases[i]
        fout = head+'/'+case+'/'+foot
        (Br,q) = load_Bres(folder=fout,slice=slice,phasing=phasing,ntor=ntor, 
                           machine=machine,cur_up=cur_up,cur_low=cur_low,phase=phase,Jres=Jres)
        Br.Psi.data = 1.0*Brs.Psi.data
        Br[key]=vals[i]
        Brs = xr.concat([Brs,Br],key)
    
    return (Brs,q)
