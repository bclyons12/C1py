# -*- coding: utf-8 -*-
"""
load_Bmns

Reads Fourier harmonics of B from M3D-C1 output in netCDF format

Author:       Brendan Carrick Lyons
Date created: Wed Jan 20 17:25:40 2016
Date edited:  
"""

import numpy as np
from numpy import pi
import xarray as xr
import xarray.ufuncs as xru
from scipy.interpolate import interp1d
from scipy.interpolate import RectBivariateSpline as RBS
from geofac import geofac

def load_Bmns(folder='./', phasing=0., slice=0, cur_up=1., cur_low=1.,
              cur_mid=None, machine=None, iplasma=None, code='m3dc1',
              ntor=None,phase=False,Jmn=False):
    
    if machine is 'diiid':
        conv = 1.0
    elif machine is 'aug':
        conv = -1.0
    elif machine is 'kstar':
        conv = 1.0
        
    if phase is not True:
        phase = False
        
    phasing = np.asarray(phasing)
    
    if code is 'all':
        
        Bc1 = load_Bmns(folder=folder, phasing=phasing, slice=slice, 
                        cur_up=cur_up, cur_low=cur_low, machine=machine, 
                        iplasma=iplasma, code='m3dc1', ntor=ntor)
        Bmars = load_Bmns(folder=folder, phasing=phasing, slice=slice, 
                          cur_up=cur_up, cur_low=cur_low, machine=machine, 
                          iplasma=iplasma, code='mars', ntor=ntor)
        Bipec = load_Bmns(folder=folder, phasing=phasing, slice=slice, 
                          cur_up=cur_up, cur_low=cur_low, machine=machine, 
                          iplasma=iplasma, code='ipec', ntor=ntor)
        
        return (Bc1, Bmars, Bipec)
    
    elif code is 'm3dc1':
        
        if Jmn:
            file_up  = folder+'/jmn3_upper-'+str(slice)+'.cdf'
            file_low = folder+'/jmn3_lower-'+str(slice)+'.cdf'
        else:
            file_up  = folder+'/bmn_upper-'+str(slice)+'.cdf'
            file_low = folder+'/bmn_lower-'+str(slice)+'.cdf'
        
        ds_up  = xr.open_dataset(file_up)
        ds_low = xr.open_dataset(file_low)
        
        ntor = ds_up.attrs['ntor']
        
        fac = geofac(machine=machine,ntor=ntor)            
        cur_up  = fac*cur_up
        cur_low = fac*cur_low

        m    = ds_up.m
        Psi  = ds_up.psi        
        
        Bup  = ds_up.bmn_real.data + 1j*ds_up.bmn_imag.data
        Blow = ds_low.bmn_real.data + 1j*ds_low.bmn_imag.data

        if cur_mid is not None:
            if Jmn:
                file_mid = folder+'/jmn_middle-'+str(slice)+'.cdf'
            else:
                file_mid = folder+'/bmn_middle-'+str(slice)+'.cdf'
            ds_mid = xr.open_dataset(file_mid)
            cur_mid = fac*cur_mid
            Bmid = ds_mid.bmn_real.data + 1j*ds_mid.bmn_imag.data
            
        if (iplasma is not None) and (slice > 0):
            
            if Jmn:
                file_up  = folder+'/jmn3_upper-0.cdf'
                file_low = folder+'/jmn3_lower-0.cdf'
            else:
                file_up  = folder+'/bmn_upper-0.cdf'
                file_low = folder+'/bmn_lower-0.cdf'
            
            ds_vu = xr.open_dataset(file_up)
            ds_vl = xr.open_dataset(file_low)
            
            Bup  -= ds_vu.bmn_real.data
            Bup  -= 1j*ds_vu.bmn_imag.data
            Blow -= ds_vl.bmn_real.data
            Blow -= 1j*ds_vl.bmn_imag.data
            
            if cur_mid is not None:
                if Jmn:
                    file_mid= folder+'/jmn3_middle-0.cdf'
                else:
                    file_mid= folder+'/bmn_middle-0.cdf'
                ds_vm = xr.open_dataset(file_mid)
                Bmid -= ds_vm.bmn_real.data
                Bmid -= 1j*ds_vm.bmn_imag.data
                
        
        q = ds_up.q.data        
        
    elif code is 'mars':
        
        file = folder + '/mars_bmn.nc'
        dset  = xr.open_dataset(file)
        
    
        m = dset.m 
        Psi = dset.s**2.0
    
        if slice is 0:
            Bup  = dset.vacuum_upper_real.data + 1j*dset.vacuum_upper_imag.data
            Blow = dset.vacuum_lower_real.data + 1j*dset.vacuum_lower_imag.data
        elif slice is 1:
            Bup  = dset.total_upper_real.data + 1j*dset.total_upper_imag.data
            Blow = dset.total_lower_real.data + 1j*dset.total_lower_imag.data
            
            if (iplasma is not None):
                Bup  -= dset.vacuum_upper_real.data 
                Bup  -= 1j*dset.vacuum_upper_imag.data
                Blow -= dset.vacuum_lower_real.data
                Blow -= 1j*dset.vacuum_lower_imag.data
                
        else:
            return NotImplemented
            
        Psi1 = dset.q_s**2.0
        q1 = dset.q
        f = interp1d(Psi1,q1,bounds_error=False)
        q = f(Psi) 
        
        # flip signs because theta is flipped
        conv = -conv
        m = -m
        q = -q
        
    elif code is 'ipec':

        if (slice == 0) or (iplasma is not None):
            
            # get the vacuum field
            
            base = folder+'ipec_vbnormal'
            sr = 7; # skip this many rows
            cr = 3; # column of real Bmn (zero-based)
            ci = 4; # column of imaginary Bmn (zero-based)
        
            (Pv,mv,Bvu,Bvl,qv) = B_ipec(base,sr,cr,ci)
        
        if slice == 1:
            
            # get the total field
            
            base = folder+'ipec_xbnormal';
            sr = 8; # skip this many rows
            cr = 7; # column of real Bmn
            ci = 8; # column of imaginary Bmn
            
            (Pt,mt,Btu,Btl,qt) = B_ipec(base,sr,cr,ci)
        
        if slice == 0:
            # just the vacuum field
            (Psi,m,Bup,Blow,q) = (Pv,mv,Bvu,Bvl,qv)
        elif (slice == 1) and (iplasma is None):
            # just the total field
            (Psi,m,Bup,Blow,q) = (Pt,mt,Btu,Btl,qt)
            
        else:
            # we need to subtract the vacuum from the total field
            Rv = RBS(Pv,mv,Bvu.real)
            Iv = RBS(Pv,mv,Bvu.imag)
            Bvu = Rv(Pt,mt) + 1j*Iv(Pt,mt)
            
            Rv = RBS(Pv,mv,Bvl.real)
            Iv = RBS(Pv,mv,Bvl.imag)
            Bvl = Rv(Pt,mt) + 1j*Iv(Pt,mt)
                        
            (Psi,m,Bup,Blow,q) = (Pt,mt,Btu-Bvu,Btl-Bvl,qt)
            
            
        conv = -conv
        m = -m
        q = -q
        
    else:
        return NotImplemented
    
    # Quantities used for all code results
    
    
    cur = xr.DataArray(np.cos(pi*phasing/180.) + conv*1j*np.sin(pi*phasing/180.),
                             {'phasing':phasing})
    
    Bmn_up  = cur_up*xr.DataArray(Bup,[('Psi',Psi),('m',m)])
    Bmn_low = cur_low*xr.DataArray(Blow,[('Psi',Psi),('m',m)]) 
    Bmn = cur*Bmn_up + Bmn_low
    
    if cur_mid is not None:
        cur2 = xr.DataArray(np.cos(pi*phasing/180.) + conv*1j*np.sin(pi*phasing/180.),
                             {'phasing2':phasing})
        Bmn_mid = cur_mid*xr.DataArray(Bmid,[('Psi',Psi),('m',m)])
        Bmn = Bmn + cur2*Bmn_mid
    
    if phase:
        Bmn = xru.angle(Bmn,deg=True) % 360.
    else:
        Bmn = np.abs(Bmn)     
    
    q = xr.DataArray(q,[('Psi',Psi)])
        
    if Jmn:
        Bmns = xr.Dataset({'Jmn':Bmn,'q':q})
    else:
        Bmns = xr.Dataset({'Bmn':Bmn,'q':q})
    
    Bmns.attrs['ntor'] = ntor
    Bmns.attrs['phase'] = phase
        
    
    return Bmns
    
    
def B_ipec(base,sr,cr,ci):

    file_up  = base+'_iu.out'
    file_low = base+'_il.out'
    
    Aup = np.loadtxt(file_up,skiprows=sr)
    Psi = Aup[:,0]
    q   = Aup[:,1]
    m   = Aup[:,2].astype(int)
    Bup = Aup[:,cr] + 1j*Aup[:,ci]
    
    Alow = np.loadtxt(file_low,skiprows=sr)
    Blow = Alow[:,cr] + 1j*Alow[:,ci]
    
    Mm = max(m) - min(m) + 1
    
    Psi = Psi.reshape(-1,Mm)[:,0]
    q   = q.reshape(-1,Mm)[:,0]
    m   = m.reshape(-1,Mm)[0,:]
    
    Bup  = 1.e4*Bup.reshape(-1,Mm)
    Blow = 1.e4*Blow.reshape(-1,Mm)
    
    return (Psi,m,Bup,Blow,q)
    