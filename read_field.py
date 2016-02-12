"""
read_field

Read fields from a C1.h5 file

M3D-C1 suite

Brendan Carrick Lyons
12/23/2015
"""
import numpy as np
import fio_py
import xarray as xr
#from C1field import C1field

# constants
import scipy.constants as spc

def read_field(name, slice=0, filename='C1.h5', points=200, phi=0.,
               rrange=None, zrange=None, linfac=1., iabs=None, iphase=None,
               isum=None, iavg=None, idiff=None,ilinear=None, iequil=None, 
               icomplex=None, ntor=None):
    
    # Redefine zero options to None for convenience
    if iabs==0:
        iabs = None
    if iphase==0:
        iphase = None
    if isum==0:
        isum = None
    if iavg==0:
        iavg = None
    if idiff==0:
        idiff = None
    if ilinear==0:
        ilinear = None
    if iequil==0:
        iequil = None
    if icomplex==0:
        icomplex = None
    
    if (isum is not None) or (iavg is not None) or (idiff is not None):
        
        if isinstance(filename,basestring):
            filename = [filename]
        if isinstance(slice,int):
            slice = [slice]
        
        N_files  = len(filename)
        N_slices = len(slice)
        N = N_files*N_slices
        
        if (N_files > 1) and (N_slices > 1):
            print("Error:  can't combine mutliple filenames and slices")
            return None
        elif N_files > 1:
            slice = slice*N
        elif N_slices > 1:
            filename = filename*N
            
        if not hasattr(linfac,"__len__"):
            # Duplicate constant linear factor                
            linfac = np.ones(N)*linfac
        elif len(linfac) != N:
            print("Error:  incorrect number of linfacs")
            return None
            
        if not hasattr(phi,"__len__"):
            # Duplicate constant linear factor                
            phi = np.ones(N)*phi
        elif len(phi) != N:
            print("Error:  incorrect number of phis")
            return None
            
        field = read_field(name, slice=slice[0], filename=filename[0], 
                           points=points, phi=phi[0], rrange=rrange, 
                           zrange=zrange, linfac=linfac[0], 
                           ilinear=ilinear, iequil=iequil, 
                           icomplex=icomplex, ntor=ntor)
                               
        if N==1:
            return field
                                 
        for i in [x+1 for x in range(N-1)]:
            field2 = read_field(name, slice=slice[i], filename=filename[i], 
                                points=points, phi=phi[i], rrange=rrange, 
                                zrange=zrange, linfac=linfac[i], 
                                ilinear=ilinear, iequil=iequil, 
                                icomplex=icomplex, ntor=ntor)
            if (isum is not None) or (iavg is not None):
                field = field + field2
            elif (idiff is not None):
                field = field + field2*((-1)**i)
                
            
        if iavg is not None:
            return field/N
        else:
            return field
            
    else:
        if not isinstance(filename,basestring):
            print('Warning:  Only considering first filename')
            filename = filename[0]
        if  hasattr(slice,"__len__"):
            print('Warning:  Only considering first slice')
            slice = slice[0]
        if hasattr(linfac,"__len__"):
            print('Warning:  Only considering first linfac')
            linfac = linfac[0]
        if hasattr(phi,"__len__"):
            print('Warning:  Only considering first phi')
            phi = phi[0]
            
        
    
    # initialize the M3D-C1 file
    isrc = fio_py.open_source(fio_py.FIO_M3DC1_SOURCE,filename)
    fio_py.get_options(isrc)
    
    # Set up options given as imput
    fio_py.set_int_option(fio_py.FIO_TIMESLICE,slice)
    
    if ilinear is not None:    
        fio_py.set_int_option(fio_py.FIO_PART,fio_py.FIO_PERTURBED_ONLY)
    elif iequil is not None:
        fio_py.set_int_option(fio_py.FIO_PART,fio_py.FIO_EQUILIBRIUM_ONLY)
    else:
        fio_py.set_int_option(fio_py.FIO_PART,fio_py.FIO_TOTAL)
    
    
    if rrange is None:
        rrange = [1.,2.5]      # I want to read this from the file somehow
    if zrange is None:
        zrange = [-1.25,1.25]  # I want to read this from the file somehow

        
    R = np.linspace(rrange[0],rrange[1],points)
    Z = np.linspace(zrange[0],zrange[1],points)
    
    # Are we reading a real field?
    if icomplex is None:
        ntor = None
    elif ntor is None:
        print('Warning:  ntor must be assigned for complex evaluation')
        print('Reading real field')
        icomplex = None
    elif ilinear is None:
        print('Warning:  Trying to access complex values but ilinear=None')
        print('          This has not been properly implemented yet')
        print('Reading real field')
        icomplex = None
        ntor = None
    
    if icomplex is None:
        dtype = np.float64
    else:
        dtype = np.complex128
    data = np.zeros((points,points),dtype=dtype)   
    
    # Lists for scalar fields read from scalars
    total_pressure    = ['total pressure','p']
    ion_pressure      = ['ion pressure','pi']
    electron_pressure = ['electron pressure','pe']
    ion_density       = ['ion density', 'ni', 'den']
    electron_density  = ['electron density', 'ne']
    
    # Dictionaries for scalar fields read  from vectors
    B_comp =  {'bx':'R', 'by':'phi', 'bt':'phi', 'bphi':'phi', 'bz':'Z'}
    J_comp =  {'jx':'R', 'jy':'phi', 'jt':'phi', 'jphi':'phi', 'jz':'Z'}
    E_comp =  {'ex':'R', 'ey':'phi', 'et':'phi', 'ephi':'phi', 'ez':'Z'}
    v_comp =  {'vx':'R', 'vy':'phi', 'vt':'phi', 'vphi':'phi', 'vz':'Z'}
    A_comp =  {'ax':'R', 'ay':'phi', 'at':'phi', 'aphi':'phi', 'az':'Z'}
    comps  = ['R', 'phi', 'Z']
    
    # Lists for vector fields
    magnetic_field   = ['magnetic field', 'b field', 'bfield', 'b_field']
    current_density  = ['current density', 'j']
    electric_field   = ['electric field', 'e', 'efield', 'e_field']
    fluid_velocity   = ['fluid velocity', 'v']
    vector_potential = ['vector potential', 'a']
    
    # Lists for composite fields
    major_radius         = ['major radius', 'r', 'x']
    height               = ['height', 'z']
    ion_temperature      = ['ion_temperature', 'ti']
    electron_temperature = ['electron_temperature', 'te']

    name_lc = name.lower()

    # Primitive scalar fields
    if name_lc in total_pressure:

        field = xr.DataArray(data,[('R',R),('Z',Z)])
        field.attrs['type']    = fio_py.FIO_TOTAL_PRESSURE
        field.name = total_pressure[0]
    
    elif name_lc in ion_pressure:

        field = xr.DataArray(data,[('R',R),('Z',Z)])
        field.attrs['species'] = fio_py.FIO_MAIN_ION
        field.attrs['type']    = fio_py.FIO_PRESSURE
        field.name = ion_pressure[0]

    elif name_lc in ion_density:

        field = xr.DataArray(data,[('R',R),('Z',Z)])
        field.attrs['species'] = fio_py.FIO_MAIN_ION
        field.attrs['type']    = fio_py.FIO_DENSITY
        field.name = ion_density[0]

    elif name_lc in electron_pressure:

        field = xr.DataArray(data,[('R',R),('Z',Z)])
        field.attrs['species'] = fio_py.FIO_ELECTRON
        field.attrs['type']    = fio_py.FIO_PRESSURE
        field.name = electron_pressure[0]
        
    elif name_lc in electron_density:

        field = xr.DataArray(data,[('R',R),('Z',Z)])
        field.attrs['species'] = fio_py.FIO_ELECTRON
        field.attrs['type']    = fio_py.FIO_DENSITY
        field.name = electron_density[0]
        
    # Primitive vector fields
    elif name_lc in magnetic_field:

        data = np.array([data,data,data]) 
        
        field = xr.DataArray(data,[('component',comps),('R',R),('Z',Z)])    
        field.attrs['type'] = fio_py.FIO_MAGNETIC_FIELD
        field.name = magnetic_field[0]

    elif name_lc in current_density:
        
        data = np.array([data,data,data]) 
        
        field = xr.DataArray(data,[('component',comps),('R',R),('Z',Z)])  
        field.attrs['type'] = fio_py.FIO_CURRENT_DENSITY
        field.name = current_density[0]
    
    elif name_lc in electric_field:
        
        data = np.array([data,data,data])
        
        field = xr.DataArray(data,[('component',comps),('R',R),('Z',Z)])    
        field.attrs['type'] = fio_py.FIO_ELECTRIC_FIELD
        field.name = electric_field[0]
        
    elif name_lc in fluid_velocity:
        
        data = np.array([data,data,data])
        
        field = xr.DataArray(data,[('component',comps),('R',R),('Z',Z)])    
        field.attrs['type'] = fio_py.FIO_FLUID_VELOCITY
        field.name = fluid_velocity[0]
        
    elif name_lc in vector_potential:
        
        data = np.array([data,data,data])
        
        field = xr.DataArray(data,[('component',comps),('R',R),('Z',Z)])    
        field.attrs['type'] = fio_py.FIO_VECTOR_POTENTIAL
        field.name = vector_potential[0]
                
    #  Vector components
    elif name_lc in B_comp:

        B = read_field('bfield', slice=slice, filename=filename, points=points, 
                       phi=phi, rrange=rrange, zrange=zrange, ilinear=ilinear,
                       iequil=iequil, icomplex=icomplex, ntor=ntor)
        
        field = B.sel(component=B_comp[name_lc])
        field.attrs['type'] = 'component'

    elif name_lc in J_comp:

        J = read_field('j', slice=slice, filename=filename, points=points, 
                       phi=phi, rrange=rrange, zrange=zrange, ilinear=ilinear,
                       iequil=iequil, icomplex=icomplex, ntor=ntor)
        
        field = J.sel(component=J_comp[name_lc])
        field.attrs['type'] = 'component'
    
    elif name_lc in E_comp:

        E = read_field('efield', slice=slice, filename=filename, points=points, 
                       phi=phi, rrange=rrange, zrange=zrange, ilinear=ilinear,
                       iequil=iequil, icomplex=icomplex, ntor=ntor)
        
        field = E.sel(component=E_comp[name_lc])
        field.attrs['type'] = 'component'
        
    elif name_lc in v_comp:

        v = read_field('v', slice=slice, filename=filename, points=points, 
                       phi=phi, rrange=rrange, zrange=zrange, ilinear=ilinear,
                       iequil=iequil, icomplex=icomplex, ntor=ntor)
        
        field = v.sel(component=v_comp[name_lc])
        field.attrs['type'] = 'component'
        
    elif name_lc in A_comp:

        A = read_field('a', slice=slice, filename=filename, points=points, 
                       phi=phi, rrange=rrange, zrange=zrange, ilinear=ilinear,
                       iequil=iequil, icomplex=icomplex, ntor=ntor)
        
        field = A.sel(component=A_comp[name_lc])
        field.attrs['type'] = 'component'
                      
    # Composite fields
    elif name_lc == 'zero':
        field = xr.DataArray(data,[('R',R),('Z',Z)])
        field.attrs['type'] = 'composite'
        field.name = name_lc
        
    elif name_lc in major_radius:
        
        field = xr.DataArray(data,[('R',R),('Z',Z)])
        
        for j in range(points):

            field[dict(Z=j)]=field.coords['R']
        
        field.attrs['type'] = 'composite'
        field.name = major_radius[0]
    
    elif name_lc in height:
        
        for i in range(points):
            field[dict(R=i)]=field.coords['Z']
        
        field.attrs['type'] = 'composite'
        field.name = height[0]
        
    elif name_lc == 'beta':
        
        if ilinear is not None:
            
            print(name + ': Cannot calculate linear version')
            print('Running for total instead')   
        
        P = read_field('p', slice=slice, filename=filename, points=points,
                       phi=phi, rrange=rrange, zrange=zrange,iequil=iequil)
        B = read_field('b', slice=slice, filename=filename, points=points,
                       phi=phi, rrange=rrange, zrange=zrange,iequil=iequil)
                       
        field = 2.0*spc.mu_0*P/(B**2)
        
        field.attrs['type'] = 'composite'
        field.name = name_lc
        
    elif name_lc in ion_temperature:
        
        Pi = read_field('pi', slice=slice, filename=filename, points=points,
                        phi=phi, rrange=rrange, zrange=zrange, ilinear=ilinear,
                        iequil=iequil, icomplex=icomplex, ntor=ntor)
                        
        ni = read_field('ni', slice=slice, filename=filename, points=points,
                        phi=phi, rrange=rrange, zrange=zrange, ilinear=ilinear,
                        iequil=iequil, icomplex=icomplex, ntor=ntor)
                        
        if ilinear is None:
            Ti = Pi/ni            
        else:
            Pi0 = read_field('pi', slice=slice, filename=filename, points=points,
                             phi=phi, rrange=rrange, zrange=zrange, iequil=1)
                        
            ni0 = read_field('ni', slice=slice, filename=filename, points=points,
                             phi=phi, rrange=rrange, zrange=zrange, iequil=1)
            
            Ti = Pi/ni0 - (Pi0*ni)/(ni0**2)
            
        J_eV = spc.physical_constants['electron volt'][0]
        
        field = Ti/J_eV
        
        field.attrs['type'] = 'composite'
        field.name = ion_temperature[0]
        
    elif name_lc in electron_temperature:
        
        Pe = read_field('pe', slice=slice, filename=filename, points=points,
                        phi=phi, rrange=rrange, zrange=zrange, ilinear=ilinear,
                        iequil=iequil, icomplex=icomplex, ntor=ntor)
                        
        ne = read_field('ne', slice=slice, filename=filename, points=points,
                        phi=phi, rrange=rrange, zrange=zrange, ilinear=ilinear,
                        iequil=iequil, icomplex=icomplex, ntor=ntor)
                        
        if ilinear is not None:
            Pe0 = read_field('pe', slice=slice, filename=filename, points=points,
                             phi=phi, rrange=rrange, zrange=zrange, iequil=1)
                        
            ne0 = read_field('ne', slice=slice, filename=filename, points=points,
                             phi=phi, rrange=rrange, zrange=zrange, iequil=1)
            Te = Pe/ne0 - (Pe0*ne)/(ne0**2)
        else:
            Te = Pe/ne  
            
            
        J_eV = spc.physical_constants['electron volt'][0]
        
        field = Te/J_eV
        
        field.attrs['type'] = 'composite'
        field.name = electron_temperature[0]
        
    elif name_lc == 'b':
        
        B = read_field('bfield', slice=slice, filename=filename,
                       points=points, phi=phi, rrange=rrange, zrange=zrange, 
                       ilinear=ilinear, iequil=iequil, icomplex=icomplex,
                       ntor=ntor)
                      
        Bx = B.sel(component='R')
        By = B.sel(component='phi')
        Bz = B.sel(component='Z')
        
        if ilinear is not None:
            B0 = read_field('bfield', slice=slice, filename=filename, phi=phi, 
                            points=points, rrange=rrange, zrange=zrange, 
                            iequil=1)
            Bx0 = B0.sel(component='R')
            By0 = B0.sel(component='phi')
            Bz0 = B0.sel(component='Z')
        
            b0 = (Bx0**2 + By0**2 + Bz0**2)**0.5

            b = (Bx*Bx0 + By*By0 + Bz*Bz0)/b0
            
        else:
            b = (Bx**2 + By**2 + Bz**2)**0.5
            
        field = b.drop('component')
        field.attrs['type'] = 'composite'
        field.name = name_lc
        
    elif name_lc == 'b2':
        
        B = read_field('bfield', slice=slice, filename=filename, 
                       points=points, phi=phi, rrange=rrange, zrange=zrange, 
                       ilinear=ilinear, iequil=iequil, icomplex=icomplex,
                       ntor=ntor)
                      
        Bx = B.sel(component='R')
        By = B.sel(component='phi')
        Bz = B.sel(component='Z')
        
        if ilinear is not None:
            B0 = read_field('bfield', slice=slice, filename=filename, phi=phi, 
                            points=points, rrange=rrange, zrange=zrange, 
                            iequil=1)
            Bx0 = B0.sel(component='R')
            By0 = B0.sel(component='phi')
            Bz0 = B0.sel(component='Z')
        
            b2 = 2.0*(Bx*Bx0 + By*By0 + Bz*Bz0)
            
        else:
            b2 = Bx**2 + By**2 + Bz**2
            
        field = b2.drop('component')
        field.attrs['type'] = 'composite'
        field.name = name_lc
    
    elif name_lc == 'omega':

        vy = read_field('vy', slice=slice, filename=filename, phi=phi, 
                        points=points, rrange=rrange, zrange=zrange, 
                        ilinear=ilinear, iequil=iequil, icomplex=icomplex,
                       ntor=ntor)
        R = read_field('r', slice=slice, filename=filename, phi=phi, 
                       points=points, rrange=rrange, zrange=zrange, 
                       ilinear=ilinear, iequil=iequil, icomplex=icomplex,
                       ntor=ntor)
                        
        field = vy/R
        
        field.attrs['type'] = 'composite'
        field.name = name_lc
        
    elif name_lc == 'psi':
        
        Aphi = read_field('ay', slice=slice, filename=filename, phi=phi, 
                        points=points, rrange=rrange, zrange=zrange, 
                        ilinear=ilinear, iequil=iequil, icomplex=icomplex,
                       ntor=ntor)
        R = read_field('r', slice=slice, filename=filename, phi=phi, 
                       points=points, rrange=rrange, zrange=zrange, 
                       ilinear=ilinear, iequil=iequil, icomplex=icomplex,
                       ntor=ntor)

        field = R*Aphi
        field = field.drop('component')
        
        field.attrs['type'] = 'composite'
        field.name = name_lc
        
    else:
        print "Field '" + name_lc + "' is+- not defined"
        return None
        
    field.attrs['ntor'] = ntor
    
    #**************************************************************************
    
    if not isinstance(field.attrs['type'],basestring):
        
        if 'species' in field.attrs:
            fio_py.set_int_option(fio_py.FIO_SPECIES, field.attrs['species'])
    
        handle = fio_py.get_field(isrc,field.attrs['type'])
        
        if len(field.shape) == 2:
            enum = np.ndenumerate(field)
            eval_field = fio_py.eval_scalar_field
        elif len(field.shape) == 3:
            enum = np.ndenumerate(field.data[0,:,:])
            eval_field = fio_py.eval_vector_field
        
        for (i,j), dummy in enum:
    
            # Real part
    
            x = (R[i], phi*np.pi/180., Z[j])
            
            val = np.asarray(eval_field(handle, x),dtype=dtype)
            
            # Imaginary part
            if icomplex is not None:
                y = (x[0], x[1] + 3.0*np.pi/(2.0*ntor), x[2])
                val += 1j*np.asarray(eval_field(handle, y),dtype=dtype)
            
            field[dict(R=i,Z=j)] = val
    
    return linfac*field

