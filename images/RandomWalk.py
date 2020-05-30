from LaunchPhoton import LaunchPhoton
from RandomStep import RandomStep
from BoundaryCheck import BoundaryCheck
from TakeStep import TakeStep
from AbsorbPhoton import AbsorbPhoton
from ScatterPhoton import ScatterPhoton

def RandomWalk(mu_att, absFact, g, n_t, n_m, alpha_critical, d, dz, dDR, dDT, pDW):
    '''
    Take a random walk within tissue untill it crosses the boundary or 
    losses all its energy in absorption.

    Parameters
    ----------
    mu_att : double
        Attenuation coefficient of tissue. 
        Sum of absorption coefficient and scattering coefficient.
    absFact : double
        Energy absorption factor.
    g : double
        Anisotropy of scatter.
    n_t : double
        Refractive index of tissue.
    n_m : double
        Refractive index of surrounding medium.
    alpha_critical : double or NaN
        Critical angle for total internal reflectance.
        Has a double value when n_t > n_m. Otherwise, NaN.
    d : double
        Thickness of tissue.
    dz : double
        Distance between points in pDW array.
    dDR : double
        Current diffuse reflectance.
    dDT : double
        Current diffuse transmittance.
    pDW : array of double
        Current energy absorption at various depth.

    Returns
    -------
    dDR : double
        Updated diffuse reflectance.
    dDT : double
        Updated diffuse transmittance.
    pDW : array of double
        Updated energy absorption at various depth.

    '''
    
    
    w, x, y, z, dpx, dpy, dpz = LaunchPhoton()
    
    while True:
        s = RandomStep(mu_att)
        w, x, y, z, dpx, dpy, dpz, s, dDR, dDT = BoundaryCheck(w, x, y, z, dpx, dpy, dpz, s, dDR, dDT, d, n_t, n_m, alpha_critical)
        
        if w <= 0:
            break
        
        x, y, z = TakeStep(x, y, z, dpx, dpy, dpz, s)
        w, pDW = AbsorbPhoton(w, z, pDW, absFact, dz)
        
        if w <= 0:
            break
        
        dpx, dpy, dpz = ScatterPhoton(dpx, dpy, dpz, g)
    
    return dDR, dDT, pDW



if __name__ == '__main__':
    
    from random import seed
    from math import asin
    import numpy as np
    
    seed(1)
    
    mu_a = 1 # aborption coefficient
    mu_s = 1 # scattering coefficient
    mu_att = mu_a + mu_s
    absFact = mu_a / mu_att
    g = 0.5
    
    n_t = 1.4
    n_m = 1.5
    if n_t > n_m:
        alpha_critical = asin(n_m / n_t)
    else:
        alpha_critical = float('nan')
    
    N_bins = 1000 # number of bins
    d = 1
    dz = d / N_bins
    
    dDR = 0
    dDT = 0
    pDW = np.zeros(N_bins)
    
    dDR, dDT, pDW = RandomWalk(mu_att, absFact, g, n_t, n_m, alpha_critical, d, dz, dDR, dDT, pDW)