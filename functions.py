from random import seed
from random import random
from math import log
from math import acos
from math import asin
from math import cos
from math import sin
from math import tan
from math import floor
from math import pi
from math import copysign
from math import sqrt
from matplotlib.pyplot import figure
from matplotlib.pyplot import plot
from matplotlib.pyplot import semilogy
import numpy as np

seed(1)



def LaunchPhoton():
    '''
    

    Returns
    -------
    w : TYPE
        DESCRIPTION.
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    z : TYPE
        DESCRIPTION.
    dpx : TYPE
        DESCRIPTION.
    dpy : TYPE
        DESCRIPTION.
    dpz : TYPE
        DESCRIPTION.

    '''
    
    w = 1
    
    x = 0
    y = 0
    z = 0
    
    dpx = 0
    dpy = 0
    dpz = 1
    
    return w, x, y, z, dpx, dpy, dpz



def RandomStep(mu_att):
    '''
    

    Parameters
    ----------
    mu_att : TYPE
        DESCRIPTION.

    Returns
    -------
    s : TYPE
        DESCRIPTION.

    '''
    
    s = - log(random()) / mu_att
    
    return s



def BoundaryCheck(x, y, z, dpx, dpy, dpz, s, d, n_t, n_m, w, alpha_critical, dDR, dDT, epsilon = 1e-9):
    '''
    

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    z : TYPE
        DESCRIPTION.
    dpx : TYPE
        DESCRIPTION.
    dpy : TYPE
        DESCRIPTION.
    dpz : TYPE
        DESCRIPTION.
    s : TYPE
        DESCRIPTION.
    d : TYPE
        DESCRIPTION.
    n_t : TYPE
        DESCRIPTION.
    n_m : TYPE
        DESCRIPTION.
    w : TYPE
        DESCRIPTION.
    alpha_critical : TYPE
        DESCRIPTION.
    dDR : TYPE
        DESCRIPTION.
    dDT : TYPE
        DESCRIPTION.
    epsilon : TYPE, optional
        DESCRIPTION. The default is 1e-9.

    Returns
    -------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    z : TYPE
        DESCRIPTION.
    dpx : TYPE
        DESCRIPTION.
    dpy : TYPE
        DESCRIPTION.
    dpz : TYPE
        DESCRIPTION.
    s : TYPE
        DESCRIPTION.
    w : TYPE
        DESCRIPTION.
    dDR : TYPE
        DESCRIPTION.
    dDT : TYPE
        DESCRIPTION.

    '''
    
    while True:
        
        if dpz < 0:
            sbound = (0 - z) / dpz
        elif dpz > 0:
            sbound = (d - z) / dpz
        else:
            break
        
        if sbound > s:
            break
        else:
            x += sbound * dpx
            y += sbound * dpy
            z += sbound * dpz
            
            s -= sbound
            
            alpha_i = acos(abs(dpz))
            
            if n_t > n_m:
                if alpha_i > alpha_critical:
                    dpz = - dpz
                    break
            else:
                alpha_t = asin(n_t * sin(alpha_i) / n_m)
                
                R = 0.5 * ((sin(alpha_i - alpha_t))**2 / ((sin(alpha_i + alpha_t))**2 + epsilon) + (tan(alpha_i - alpha_t))**2 / ((tan(alpha_i + alpha_t))**2 + epsilon))
                
                if R < random():
                    if dpz < 0:
                        dDR += w
                    else:
                        dDT += w
                    
                    w = 0
                
                dpz = - dpz
    
    return x, y, z, dpx, dpy, dpz, s, w, dDR, dDT




def TakeStep(x, y, z, dpx, dpy, dpz, s):
    '''
    

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    z : TYPE
        DESCRIPTION.
    dpx : TYPE
        DESCRIPTION.
    dpy : TYPE
        DESCRIPTION.
    dpz : TYPE
        DESCRIPTION.
    s : TYPE
        DESCRIPTION.

    Returns
    -------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    z : TYPE
        DESCRIPTION.

    '''
    
    x += s * dpx
    y += s * dpy
    z += s * dpz
    
    return x, y, z



def AbsorbPhoton(w, absFact, z, dz, pDW):
    '''
    

    Parameters
    ----------
    w : TYPE
        DESCRIPTION.
    absFact : TYPE
        DESCRIPTION.
    z : TYPE
        DESCRIPTION.
    dz : TYPE
        DESCRIPTION.
    pDW : TYPE
        DESCRIPTION.

    Returns
    -------
    w : TYPE
        DESCRIPTION.
    pDW : TYPE
        DESCRIPTION.

    '''
    
    dw = w * absFact
    
    if (w - dw) < 0:
        dw = w
    w -= dw
    
    nz = floor(z / dz)
    pDW[nz] += dw
    
    return w, pDW





def ScatterPhoton(dpx, dpy, dpz, g):
    '''
    

    Parameters
    ----------
    dpx : TYPE
        DESCRIPTION.
    dpy : TYPE
        DESCRIPTION.
    dpz : TYPE
        DESCRIPTION.
    g : TYPE
        DESCRIPTION.

    Returns
    -------
    dpx : TYPE
        DESCRIPTION.
    dpy : TYPE
        DESCRIPTION.
    dpz : TYPE
        DESCRIPTION.

    '''
    
    
    if g == 0:
        theta = acos(2 * random() - 1)
    else:
        theta = acos((1 + g**2 - ((1 - g**2) / (1 - g + 2 * g * random()))**2) / (2 * g))
    
    psi = 2 * pi * random()
    
    
    if abs(dpz) > 0.9999:
        dpx = sin(theta) * cos(psi)
        dpy = sin(theta) * sin(psi)
        dpz = copysign(cos(theta), dpz)
    else:
        dpx = sin(theta) * (dpx * dpz * cos(psi) - dpy * sin(psi)) / sqrt(1 - dpz**2) + dpx * cos(theta)
        dpy = sin(theta) * (dpy * dpz * cos(psi) - dpx * sin(psi)) / sqrt(1 - dpz**2) + dpy * cos(theta)
        dpz = -sin(theta) * cos(psi) * sqrt(1 - dpz**2) + dpz * cos(theta)
    
    return dpx, dpy, dpz



def RandomWalk(mu_att, absFact, g, n_t, n_m, d, dDR, dDT, pDW, alpha_critical, dz):
    '''
    

    Parameters
    ----------
    mu_att : TYPE
        DESCRIPTION.
    absFact : TYPE
        DESCRIPTION.
    g : TYPE
        DESCRIPTION.
    n_t : TYPE
        DESCRIPTION.
    n_m : TYPE
        DESCRIPTION.
    d : TYPE
        DESCRIPTION.
    dDR : TYPE
        DESCRIPTION.
    dDT : TYPE
        DESCRIPTION.
    pDW : TYPE
        DESCRIPTION.
    alpha_critical : TYPE
        DESCRIPTION.
    dz : TYPE
        DESCRIPTION.

    Returns
    -------
    dDR : TYPE
        DESCRIPTION.
    dDT : TYPE
        DESCRIPTION.
    pDW : TYPE
        DESCRIPTION.

    '''
    
    
    w, x, y, z, dpx, dpy, dpz = LaunchPhoton()
    
    while True:
        s = RandomStep(mu_att)
        x, y, z, dpx, dpy, dpz, s, w, dDR, dDT = BoundaryCheck(x, y, z, dpx, dpy, dpz, s, d, n_t, n_m, w, dDR, dDT, alpha_critical)
        
        if w <= 0:
            break
        
        x, y, z = TakeStep(x, y, z, dpx, dpy, dpz, s)
        w, pDW = AbsorbPhoton(w, absFact, z, dz, pDW)
        
        if w <= 0:
            break
        
        dpx, dpy, dpz = ScatterPhoton(dpx, dpy, dpz, g)
    
    return dDR, dDT, pDW




mu_as = np.array([1.63, 6.85, 16.13])
mu_s = 74.7
g = 0.759
n_t = 1.447
n_m = 1.5
d = 0.1

N_bins = int(1e3)
dz = d / N_bins
dalpha = pi / (2 * N_bins)

N_photons = int(5e6)

mu_atts = mu_as + mu_s
absFacts = mu_as / mu_atts


pDW = np.zeros((3, N_bins))
dDR = np.zeros(3)
dDT = np.zeros(3)


for n, (mu_att, absFact) in enumerate(zip(mu_atts, absFacts)):
    for nph in range(N_photons):
        dDR[n], dDT[n], pDW[n] = RandomWalk(mu_att, absFact, g, n_t, n_m, d, dDR[n], dDT[n], pDW[n], dz)



dDR /= N_photons
dDT /= N_photons
dDW = np.sum(pDW, axis = 1) / N_photons


zs = np.arange(0, d, dz)
pDW = pDW / (dz * N_photons)
figure()
semilogy(zs, pDW[0], 'r', zs, pDW[1], 'g', zs, pDW[2], 'b')




mu_as = np.array([30.16, 13.5, 6.85, 3.68])
mu_ss = np.array([210.4, 121.6, 74.7, 55.48])
gs = np.array([0.702, 0.728, 0.759, 0.787])
n_ts = np.array([1.488, 1.448, 1.447, 1.433])

mu_atts = mu_as + mu_ss
absFacts = mu_as / mu_atts


pDW = np.zeros((4, N_bins))
dDR = np.zeros(4)
dDT = np.zeros(4)


for n, (mu_att, absFact, g, n_t) in enumerate(zip(mu_atts, absFacts, gs, n_ts)):
    for nph in range(N_photons):
        dDR[n], dDT[n], pDW[n] = RandomWalk(mu_att, absFact, g, n_t, n_m, d, dDR[n], dDT[n], pDW[n], dz)



dDR /= N_photons
dDT /= N_photons
dDW = np.sum(pDW, axis = 1) / N_photons


zs = np.arange(0, d, dz)
pDW = pDW / (dz * N_photons)
figure()
semilogy(zs, pDW[0], 'r', zs, pDW[1], 'g', zs, pDW[2], 'b', pDW[3], 'k')


wavelengths = np.array([350, 450, 550, 650])
p = np.polyfit(wavelengths, dDW, 3)
pd = np.poly1d(p)
waves = np.arange(350, 650, 10)
fitted = pd(waves)
figure()
plot(wavelengths, dDW, 'ro', waves, fitted, 'b')