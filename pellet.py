import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, CubicHermiteSpline
from numbers import Number
import sys
try:
    import seaborn as sns
except ModuleNotFoundError:
    pass

def plot_resolution(Rp,dY,Np,dist='blend', cf=0.0, ntor=None, fig=None, ax=None):
    try:
        sns.set_style('white')
        sns.set_palette('colorblind')
    except NameError:
        pass

    if (fig is None) and (ax is None):
        fig,ax = plt.subplots()

    th = np.linspace(-np.pi, np.pi, 100001)
    # plot location of uniformly spaced planes
    ps = np.linspace(-np.pi, np.pi, Np+1)
    ymin = [-1,0][dist!='fourier']
    for p in ps:
        ax.plot([p,p], [ymin,1], 'k--', lw=1)
    ax.plot([-np.pi,np.pi], [0,0], 'k', lw=1)

    # plot the von Mises distribution through the center of the pellet
    if dist == 'vonmises':
        f = np.exp(-(Rp/dY)**2*(1. - np.cos(th)))
        dfdt = -(Rp/dY)**2*np.sin(th)*f
    elif dist == 'cauchy':
        gamma = dY/Rp
        f = (np.cosh(gamma)-1.)/(np.cosh(gamma)-np.cos(th))
        dfdt = -(np.cosh(gamma)-1.)*np.sin(th)/(np.cosh(gamma)-np.cos(th))**2
    elif dist == 'blend':
        fv = np.exp(-(Rp/dY)**2*(1. - np.cos(th)))
        sv = np.trapz(fv,th)
        fv /= sv
        dfvdt = -(Rp/dY)**2*np.sin(th)*fv
        gamma = dY/Rp
        fc = (np.cosh(gamma)-1.)/(np.cosh(gamma)-np.cos(th))
        sc = np.trapz(fc,th)
        fc /= sc
        dfcdt = -fc*np.sin(th)/(np.cosh(gamma)-np.cos(th))

        f = (1.-cf)*fv + cf*fc
        sf = np.trapz(f,th)
        f /= sf
        dfdt = ((1.-cf)*dfvdt + cf*dfcdt)/sf
    elif dist == 'fourier':
        f = np.cos(ntor*th)
        dfdt = -ntor*np.sin(ntor*th)


    ax.plot(th, f, lw=3)

    # plot Cubic Hermite interpolation using plane locations
    fp    = interp1d(th, f)(ps)
    dfpdt = interp1d(th, dfdt)(ps)
    fi = CubicHermiteSpline(ps, fp, dfpdt)(th)

    if dist!='fourier' and np.any(fi<0):
       print("Warning: interpolated pellet distribution goes negative: %.2e"%fi.min())

    ax.plot(th, fi, lw=3)

    ax.set_xlim([-np.pi, np.pi])
    ax.set_xticks([-np.pi, -np.pi/2, 0., np.pi/2, np.pi])
    ax.set_xticklabels([r'$-\pi$',r'$-\frac{\pi}{2}$','0',
                        r'$\frac{\pi}{2}$',r'$-\pi$'])
    fig.tight_layout()
    plt.show()

    return (fig, ax)



def create_plume(D_p=1e-3, L_D = 1.5, v=200., device='diii-d',
                 th_pl = 5., L_pl = 0.2, R_front=None, seed=None,
                 sdist='uniform', tdist='uniform', pdist='sunflower',
                 N_s=None, species='Ar',
                 facct=1., sublimate='total',
                 pellet_rate=0., pellet_var=1., pellet_var_tor=0.,
                 cloud_pel=1., pellet_mix=0., cauchy_fraction=0.,
                 v_A=1.54235e6, debug=False):

    # SPI configurations
    if device == 'diii-d':
        R_torus = 1.75
        Z_torus = .0
        a_torus = .55
        R_SPI = 2.262  # radial location of injector
        PHI_SPI = 0.0  # phi location of injector
        Z_SPI = 0.664  # height location of injector
        azi_SPI = 0.   # angle of toroidal injection velocity
        inc_SPI = 49. # angle of poloidal injection velocity
        bend_SPI = 20.
    elif device == 'jet':
        R_torus = 2.9
        Z_torus = .2
        a_torus = .9
        R_SPI = 3.239  # radial location of injector
        PHI_SPI = 0.0  # phi location of injector
        Z_SPI = 1.885  # height location of injector
        azi_SPI = 0.   # angle of toroidal injection velocity
        inc_SPI = 65.59 # angle of poloidal injection velocity
        bend_SPI = 20.
    elif device == 'kstar':
        R_torus = 1.8
        Z_torus = 0.0
        a_torus = .45
        R_SPI = 2.63818  # radial location of injector
        PHI_SPI = 0.0  # phi location of injector
        Z_SPI = -0.45954  # height location of injector
        azi_SPI = 0.   # angle of toroidal injection velocity
        inc_SPI = -20. # angle of poloidal injection velocity
        bend_SPI = 20.
    elif device == 'iter':
        R_torus = 6.225
        Z_torus = .55
        a_torus = 1.975
        R_SPI = 8.454232668
        PHI_SPI = 0.
        Z_SPI = 0.6855
        azi_SPI = 0.
        inc_SPI = 0.
        bend_SPI = 20.
    else:
        raise ValueError("Device %s not implemented"%device)

    inc_SPI *= np.pi/180.
    azi_SPI *= np.pi/180.
    PHI_SPI *= np.pi/180.
    th_pl *= np.pi/180.

    np.random.seed(seed)
    # shatter pellet into shards of radius r_s
    r_s = shatter(D_p=D_p, L_D=L_D, v=v, bend=bend_SPI, sdist=sdist, N_s=N_s,
                  species=species, facct=facct, sublimate=sublimate)

    # distribute shards into cone
    dist = None
    while dist is None:
        np.random.shuffle(r_s)
        dist = distribute(r_s, v=v, th_pl=th_pl, L_pl=L_pl,
                          pdist=pdist, tdist=tdist, angle='rad',debug=debug)
    t_s, th_s, phi_s = dist

    def rotate(xpp,ypp,zpp):
        xp = -xpp*np.sin(inc_SPI) - zpp*np.cos(inc_SPI)
        yp = ypp
        zp = xpp*np.cos(inc_SPI) - zpp*np.sin(inc_SPI)
        x = xp*np.cos(azi_SPI-PHI_SPI) + yp*np.sin(azi_SPI-PHI_SPI)
        y = -xp*np.sin(azi_SPI-PHI_SPI)  + yp*np.cos(azi_SPI-PHI_SPI)
        z = zp
        return x,y,z

    # get cartesian velocities (v'')
    vx_s = v*np.sin(th_s)*np.cos(phi_s)
    vy_s = v*np.sin(th_s)*np.sin(phi_s)
    vz_s = v*np.cos(th_s)
    vx_s,vy_s,vz_s = rotate(vx_s,vy_s,vz_s)

    # Advance plume to R_front
    if R_front is not None:
        tadv = ((R_SPI-L_pl*np.cos(inc_SPI)*np.cos(azi_SPI))-R_front)/(v*np.cos(inc_SPI)*np.cos(azi_SPI))
        t_s += tadv

    X_s = R_SPI*np.cos(PHI_SPI) + vx_s*t_s
    Y_s = R_SPI*np.sin(PHI_SPI) + vy_s*t_s
    R_s = np.sqrt(X_s**2 + Y_s**2)
    Phi_s = np.arctan2(Y_s,X_s)
    Z_s = Z_SPI + vz_s*t_s

    vr_s = (X_s*vx_s + Y_s*vy_s)/R_s
    vphi_s = (X_s*vy_s - Y_s*vx_s)/R_s
    vz_s = vz_s

    if debug:
        from mpl_toolkits import mplot3d
        plt.close('all')
        f = plt.figure()
        ax = plt.axes(projection='3d')
        phi = np.linspace(0.,2.*np.pi,10001)
        x = (R_torus + a_torus*np.cos(300.*phi))*np.cos(phi)
        y = (R_torus + a_torus*np.cos(300.*phi))*np.sin(phi)
        z = Z_torus + 1.5*a_torus*np.sin(300.*phi)
        ax.plot3D(x,y,z,'k',alpha=.25)

        t = np.linspace(0,1,10001)
        if R_front is not None:
            z = (R_SPI - R_front)*t/(np.cos(inc_SPI)*np.cos(azi_SPI))
        else:
            z = L_pl*t
        x = z*np.tan(th_pl)*np.cos(100*2.*np.pi*t)
        y = z*np.tan(th_pl)*np.sin(100*2.*np.pi*t)
        x,y,z = rotate(x,y,z)
        x += R_SPI*np.cos(PHI_SPI)
        y += R_SPI*np.sin(PHI_SPI)
        z += Z_SPI
        ax.plot3D(x,y,z,'b',alpha=0.2)

        ax.set_xlim(R_SPI*np.cos(PHI_SPI) + np.array([-L_pl,L_pl]))
        ax.set_xlabel('x')
        ax.set_ylim(R_SPI*np.sin(PHI_SPI) + np.array([-L_pl,L_pl]))
        ax.set_ylabel('Y')
        ax.set_zlim(Z_SPI + np.array([-L_pl,L_pl]))
        ax.set_zlabel('Z')

        if sdist=='uniform':
            ax.scatter3D(X_s,Y_s,Z_s,c='red')
        else:
            ax.scatter3D(X_s,Y_s,Z_s,c=r_s,cmap='viridis')


    plume = np.zeros((len(r_s),13))
    plume[:,0] = R_s
    plume[:,1] = Phi_s
    plume[:,2] = Z_s
    plume[:,3] = pellet_rate
    plume[:,4] = pellet_var
    plume[:,5] = pellet_var_tor
    plume[:,6] = vr_s/v_A
    plume[:,7] = vphi_s/v_A
    plume[:,8] = vz_s/v_A
    plume[:,9] = r_s
    plume[:,10] = cloud_pel
    plume[:,11] = pellet_mix
    plume[:,12] = cauchy_fraction

    # sort by R location
    plume = plume[np.argsort(plume[:,0])]

    return plume



def shatter(D_p=1e-3, L_D = 1.5, v=200., bend=20., device='diii-d',
            sdist='uniform', N_s=None, species='Ar', seed=None,
            facct=1., sublimate='total'):
    from scipy.special import kn, erf
    from scipy.interpolate import interp1d

    def uniform(x):
        Vtot =  0.25*L_D*np.pi*D_p**3
        return (0.75*(Vtot/N_s)/np.pi)**(1./3.)

    def fragmentation(x):

        if (sublimate=='small') and (x < (1.-facct)):
            return 0.

        if species == 'Ar':
            v_th = 6. # m/s for argon
            C = 8. # pellet material constant for argon
        elif species == 'Ne':
            v_th = 8.
            C = 5.
        elif species == 'D2':
            v_th = 20.
            C = 2.5
        else:
            #species is mol(D)/[mol(D)+mol(Ne)]
            w_Ne = (1.0-species)*20.1797
            w_Ne /= (1.0-species)*20.1797 + species*2.014
            v_th = 8.*w_Ne + 20.*(1.-w_Ne)
            C = 2.5*(1.0 + w_Ne)

        D_cm = D_p*1e2

        # initial velocity: v_i = v/cos(bend)
        # normal velocity: v_n= v_i*sin(bend)= v*tan(bend)
        v_n = v*np.tan(np.pi*bend/180.)
        X_R = (v_n/v_th)**2
        b = X_R/(L_D*(D_cm)*C)

        d = np.linspace(0.,D_cm,10000)
        K1 = kn(1,b*d)
        K1[0] = 0.
        F = 1. - b*d*K1
        F[0] = 0.

        try:
            dp = interp1d(F,d)(x)/1e2
        except ValueError:
            #print('Skipping x=%f, greater than pellet diameter at x=%f'%(x,F[-1]))
            dp = 0.

        return dp/2 # divide by two to get radius


    # Calcluate the fragment distribution
    Vtot = 0.25*L_D*np.pi*D_p**3
    if sublimate in ['total','small']:
        Vtot *= facct
    V = 0.
    r_s = []
    if seed is not None:
        np.random.seed(seed)
    while V < Vtot:
        r = eval(sdist)(np.random.random())
        if r==0.:
            continue
        if V + (4./3.)*np.pi*r**3 > Vtot:
            if (N_s is not None) and (len(r_s) == N_s):
                break
            else:
                r = (0.75*(Vtot - V)/np.pi)**(1./3.)
        V += (4./3.)*np.pi*r**3
        r_s += [r]
    r_s = np.array(r_s)

    # Sublimate some of the pellet away weighted by surface area
    if sublimate == 'surf_area':
        Vacct = Vtot
        dV = (1.-facct)*Vtot/100.
        while Vacct > facct*Vtot:
            Atot = 4.*np.pi*np.sum(r_s**2)
            Vsub = (4.*np.pi*(r_s**2)/Atot)*dV
            V_s = (4.*np.pi*r_s**3/3. - Vsub)
            V_s = V_s[V_s > 0.]
            Vacct = np.sum(V_s)
            r_s = (0.75*V_s/np.pi)**(1./3.)

    return r_s

def distribute(r_s, v=200., device='diii-d', th_pl = 5., L_pl = 0.2,
               pdist='uniform', tdist='uniform', angle='deg', seed=None,
               debug=False):

    # time length of plume
    t_pl = L_pl/v

    if angle == 'deg':
        th_pl *= np.pi/180.
    sigma = th_pl/np.sqrt(2.)

    # we call it "mass" but it's volume here
    ms = (4./3.)*np.pi*r_s**3
    M = sum(ms)
    N = len(ms)

    if debug:
        print("%d shards totaling %.2e m^-3"%(N,M))

    ths = np.zeros(N)
    phis = np.zeros(N)
    ts = np.zeros(N)
    xc = 0.
    yc = 0.
    zc = 0.
    dc = 0.
    imax = np.argmax(ms)

    if seed is not None:
        np.random.seed(seed)

    if tdist == 'uniform':
        # randomize order
        ir = np.arange(N)
        np.random.shuffle(ir)

    for i,m in enumerate(ms):
        if tdist == 'uniform':
            ts[i] = t_pl*ir[i]/(N-1)
        elif tdist == 'random':
            ts[i] = t_pl*np.random.rand()
        elif tdist == 'normal':
            ts[i] = t_pl*(0.5 + 0.25*np.random.randn())
        else:
            raise ValueError("tdist=%s not valid"%tdist)

        if pdist == 'sunflower':
            phis[i] = 8.*np.pi*i/(np.sqrt(5.)+1.)**2
            ths[i] = th_pl*np.sqrt(i+0.5)/np.sqrt(N+0.5)
        elif pdist == 'normal':
            phis[i] = 2.*np.pi*np.random.rand()
            x = np.random.rand()
            ths[i] = np.sqrt(-2*sigma**2*np.log(1.-x))
        else:
            raise ValueError("pdist=%s not valid"%pdist)
        d = v*ts[i]
        dc += m*d
        xc += m*d*np.sin(ths[i])*np.cos(phis[i])
        yc += m*d*np.sin(ths[i])*np.sin(phis[i])
        zc += m*d*np.cos(ths[i])

    dc /= M
    xc /= M
    yc /= M
    zc /= M
    d = v*t_pl
    rc = np.sqrt(xc**2+yc**2)

    dr = 0.05*(0.5*d*np.sin(th_pl))
    if debug:
        print("Center of mass at th=%.2e, limit %.2e"%(rc,dr))
    if rc > dr:
        if debug:
            print("th too far, try again")
        return

    z0 = 0.5*d*np.cos(th_pl)
    dz = 0.05*d
    #d0 = 0.5*d
    #dd = 0.05*d
    if debug:
        #print("Center of mass at d-d0=%.2e, limit %.2e"%(dc-d0, dd))
        print("Center of mass at z-z0=%.2e, limit %.2e"%(zc-z0, dz))
    if np.abs(zc-z0) > dz:
        if debug:
            print("d too far, try again")
        return
    # Else we've got a good distribution

    if debug:
        print("Good")
        cols = sns.color_palette('viridis',N)
        f,ax = plt.subplots()

        size = 3.*(ms*N/M)**(1./3.)
        for i,m in enumerate(ms):
            col = cols[i]
            ax.plot([ths[i]*np.cos(phis[i])],[ths[i]*np.sin(phis[i])],'o',ms=size[i],color=col)
        ax.plot([xc],[yc],'m*',ms=10)

        print("1 sigma: ",np.sum(np.where(ths<=1*sigma,ms,0.*ms))/M)
        print("th_pl: ",np.sum(np.where(ths<=1*th_pl,ms,0.*ms))/M)
        print("2 sigma: ",np.sum(np.where(ths<=2*sigma,ms,0.*ms))/M)
        print("2*th_pl: ",np.sum(np.where(ths<=2*th_pl,ms,0.*ms))/M)

        plt.gca().set_aspect('equal',adjustable='box')
        th = np.linspace(0,2*np.pi,1001)
        ax.plot(th_pl*np.cos(th),th_pl*np.sin(th),'k--')
        ax.plot(2*th_pl*np.cos(th),2*th_pl*np.sin(th),'k--')

    return ts, ths, phis

if __name__ == "__main__":
    Rp = float(sys.argv[1])
    dY = float(sys.argv[2])
    Np = int(sys.argv[3])
    plot_resolution(Rp, dY, Np)
