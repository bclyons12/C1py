# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 15:02:21 2016

@author: lyonsbc
"""

import xarray as xr
from load_Bres import load_Bres

def multi_Bres(z,c,w,index,head='./',foot='/',slice=0,phasing=0.,
               machine='diiid',cur_up=1.,cur_low=1.,ntor=None,phase=False):
                  

    key  = index.keys()[0]
    vals = index[key]                  

    r = "{:.3f}".format(z[0])+'_'+"{:.3f}".format(c[0])+'_'+"{:.3f}".format(w[0])
    
    fout = head+'/om_'+r+'/'+foot
    (Brs,q) = load_Bres(folder=fout,slice=slice,phasing=phasing,ntor=ntor,
                        machine=machine,cur_up=cur_up,cur_low=cur_low,phase=phase)
    Brs[key]=vals[0]
    
    N = len(z)
    
    for i in range(1,N):
        if i==10:
            r = "{:.4f}".format(z[i])+'_'+"{:.4f}".format(c[i])+'_'+"{:.3f}".format(w[i])
        else:
            r = "{:.3f}".format(z[i])+'_'+"{:.3f}".format(c[i])+'_'+"{:.3f}".format(w[i])
        fout = head+'/om_'+r+'/'+foot
        (Br,q) = load_Bres(folder=fout,slice=slice,phasing=phasing,ntor=ntor, 
                           machine=machine,cur_up=cur_up,cur_low=cur_low,phase=phase)
        Br.Psi.data = 1.0*Brs.Psi.data
        Br[key]=vals[i]
        Brs = xr.concat([Brs,Br],key)
    
    return (Brs,q)