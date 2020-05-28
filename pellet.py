import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, CubicHermiteSpline
import sys
try:
    import seaborn as sns
except ModuleNotFoundError:
    pass

def plot_resolution(Rp,dY,Np,dist='blend', cf=0.0, fig=None, ax=None):
    try:
        sns.set_style('white')
        sns.set_palette('colorblind')
    except NameError:
        pass

    if (fig is None) and (ax is None):
        fig,ax = plt.subplots()

    th = np.linspace(-np.pi, np.pi, 1001)

    # plot location of uniformly spaced planes
    ps = np.linspace(-np.pi, np.pi, Np+1)
    for p in ps:
        ax.plot([p,p], [0,1], 'k--', lw=1)
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
        
        
    ax.plot(th, f, lw=3)

    # plot Cubic Hermite interpolation using plane locations
    fp    = interp1d(th, f)(ps)
    dfpdt = interp1d(th, dfdt)(ps)
    fi = CubicHermiteSpline(ps, fp, dfpdt)(th)

    if np.any(fi<0):
       print("Warning: interpolated pellet distribution goes negative: %.2e"%fi.min())

    ax.plot(th, fi, lw=3)

    ax.set_xlim([-np.pi, np.pi])
    ax.set_xticks([-np.pi, -np.pi/2, 0., np.pi/2, np.pi])
    ax.set_xticklabels([r'$-\pi$',r'$-\frac{\pi}{2}$','0',
                        r'$\frac{\pi}{2}$',r'$-\pi$'])
    fig.tight_layout()
    plt.show()
    
    return (fig, ax)

if __name__ == "__main__":
    Rp = float(sys.argv[1])
    dY = float(sys.argv[2])
    Np = int(sys.argv[3])
    plot_resolution(Rp, dY, Np)
