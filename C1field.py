# -*- coding: utf-8 -*-
"""
C1field class

Author:       Brendan Carrick Lyons
Date created: Tue Jan 12 14:04:27 2016
Date edited:  
"""
import numpy as np
import fio_py
import copy
from scipy.interpolate import UnivariateSpline

class C1field:
    
    def __init__(self):
        
        self.R       = None
        self.Z       = None
        self.rank    = None
        self.data    = None
        self.type    = None
        self.handle  = None
        self.species = None
        self.dim     = None
        self.ntor    = None
        
    def eval(self,x):
        
        
        
        if self.rank == 1:
            
            if self.dim is None:
                # Scalar field
                R = fio_py.eval_scalar_field(self.handle, x)
            else:
                # One component of a vector field
                v =  fio_py.eval_vector_field(self.handle, x)
                R = v[self.dim]
        
        elif self.rank == 2:
            
            # Vector field
            v =  fio_py.eval_vector_field(self.handle, x)
            R =  np.asarray(v)
        
        if self.ntor is None:
            # return the real part
            return R
        else:
            # want complex, so need to add in the imaginary part
            y = (x[0], x[1] + 3.0*np.pi/(2.0*self.ntor), x[2])
            
            if self.rank == 1:
                
                if self.dim is None:
                    # Scalar field
                    I = fio_py.eval_scalar_field(self.handle, y)
                else:
                    # One component of a vector field
                    v =  fio_py.eval_vector_field(self.handle, y)
                    I = v[self.dim]
            
            elif self.rank == 2:
                
                # Vector field
                v =  fio_py.eval_vector_field(self.handle, y)
                I =  np.asarray(v)
            
            return R + 1j*I
        
    def component(self,dim):
        
        if self.rank != 2:
            return NotImplemented
        
        new = copy.deepcopy(self)
        
        new.rank = 1
        new.dim = dim
        new.data = self.data[dim,:,:]
        
        return new
        
    def enumerate(self):
        
        if self.rank == 1:
            return np.ndenumerate(self.data)
        elif self.rank == 2:
            return np.ndenumerate(self.data[0,:,:])
            
        return NotImplemented
        
    def dR(self):
        
        if self.rank == 1:
            
            new = copy.deepcopy(self)
            
            for j in range(len(self.Z)):
                
                y = self.data[:,j]
                x = self.R
                
                spl = UnivariateSpline(x,y)
        
        elif self.rank == 2:
            return NotImplemented
        
        return NotImplemented

###############################################################################

    def __str__(self):
        return 'C1field\n' + np.str(self.data)
    
    def __radd__(self,other):
        new = copy.deepcopy(self)
        if isinstance(other,C1field):
            new.data = self.data + other.data
        else:
            new.data = self.data + other
        return new
    
    def __add__(self,other):
        new = copy.deepcopy(self)
        if isinstance(other,C1field):
            new.data = self.data + other.data
        else:
            new.data = self.data + other
        return new
    
    def __rsub__(self,other):
        new = copy.deepcopy(self)
        if isinstance(other,C1field):
            new.data = other.data - self.data
        else:
            new.data = other - self.data
        return new
    
    def __sub__(self,other):
        new = copy.deepcopy(self)
        if isinstance(other,C1field):
            new.data = self.data - other.data
        else:
            new.data = self.data - other
        return new
        
    def __rmul__(self,other):
        new = copy.deepcopy(self)
        if isinstance(other,C1field):
            new.data = np.multiply(self.data, other.data)
        else:
            new.data = np.multiply(self.data, other)
        return new
        
    def __mul__(self,other):
        new = copy.deepcopy(self)
        if isinstance(other,C1field):
            new.data = np.multiply(self.data, other.data)
        else:
            new.data = np.multiply(self.data, other)
        return new
        
    def __rdiv__(self,other):
        new = copy.deepcopy(self)
        if isinstance(other,C1field):
            new.data = np.divide(other.data,self.data)
        else:
            new.data = np.divide(other,self.data)
        return new
        
    def __div__(self,other):
        new = copy.deepcopy(self)
        if isinstance(other,C1field):
            new.data = np.divide(self.data, other.data)
        else:
            new.data = np.divide(self.data, other)
        return new
    
    def __rpow__(self,other):
        new = copy.deepcopy(self)
        if isinstance(other,C1field):
            new.data = np.power(other.data,self.data)
        else:
            new.data = np.power(other,self.data)
        return new
        
    def __pow__(self,other):
        new = copy.deepcopy(self)
        if isinstance(other,C1field):
            new.data = np.power(self.data, other.data)
        else:
            new.data = np.power(self.data, other)
        return new
        
    