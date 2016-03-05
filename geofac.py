# -*- coding: utf-8 -*-
"""
geofac

Author:       Brendan Carrick Lyons
Date created: Thu Jan 28 14:32:30 2016
Date edited:  
"""
import numpy as np
from numpy import pi

def geofac(machine=None,ntor=None):
    
    if machine is 'diiid':
                        
        if ntor==1:
            fac = 3.0/pi
        elif ntor==2:
            fac = 3.0*np.sqrt(3.0)/(2.0*pi)
        elif ntor==3:
            fac = 4.0/pi
        else:
            print('Error: ntor  unacceptable for ' + machine)
            return None
        
    elif machine is 'aug':
        
        if ntor==1:
            fac = 8.0*np.sin(pi/8.0)/pi
        elif ntor==2:
            fac = 2.0*np.sqrt(2.0)/pi
        elif ntor==3:
            fac = 8.0*np.cos(pi/8.0)/(3.0*pi)
        elif ntor==4:
            fac = 4.0/pi
        else:
            print('Error: ntor  unacceptable for ' + machine)
            return None

    elif machine is 'kstar':
        
        if ntor==1:
            fac = 2.0*np.sqrt(2.0)/pi
        elif ntor==2:
            fac = 4.0/pi
        else:
            print('Error: ntor  unacceptable for ' + machine)
            return None
    
    else:
        print('Error: No such machine ' + machine)
        return None
        
    return fac