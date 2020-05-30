from random import random
from math import sin
from math import tan
from math import asin
from math import acos

def BoundaryCheck(w, x, y, z, dpx, dpy, dpz, s, dDR, dDT, d, n_t, n_m, alpha_critical, epsilon = 1e-9):
    '''
    Check whether the photon hits the boundary or not, and update 
    photon energy, step size and direction of propagation accordingly. 
    Register diffuse reflectance and diffuse transmittance if the photon 
    crosses the top and bottom boundary, respectively.

    Parameters
    ----------
    w : double
        Current energy of photon.
    x : double
        Current x position of photon.
    y : double
        Current y position of photon.
    z : double
        Current z position of photon.
    dpx : double
        Current direction of propagation along the x direction.
    dpy : double
        Current direction of propagation along the y direction.
    dpz : double
        Current direction of propagation along the z direction.
    s : double
        Current step size.
    dDR : double
        Current diffuse reflectance.
    dDT : double
        Current diffuse transmittance.
    d : double
        Thickness of tissue.
    n_t : double
        Refractive index of tissue.
    n_m : double
        Refractive index of surrounding medium.
    alpha_critical : double or NaN
        Critical angle for total internal reflectance. 
        Has a double value when n_t > n_m. Otherwise, NaN.
    epsilon : double, optional
        A small value for avoiding division by zero. 
        The default is 1e-9.

    Returns
    -------
    w : double
        Updated energy of photon
    x : double
        Updated x position of photon
    y : double
        Updated y position of photon
    z : double
        Updated z position of photon
    dpx : double
        Updated direction of propagation along the x direction
    dpy : double
        Updated direction of propagation along the y direction
    dpz : double
        Updated direction of propagation along the z direction
    s : double
        Updated step size
    dDR : double
        Updated diffuse reflectance
    dDT : double
        Updated diffuse transmittance

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
    
    return w, x, y, z, dpx, dpy, dpz, s, dDR, dDT



if __name__ == "__main__":
    
    from random import seed
    
    seed(1)
    
    d = 1.0
    n_t = 1.4
    n_m = 1.5
    if n_t > n_m:
        alpha_critical = asin(n_m / n_t)
    else:
        alpha_critical = float('nan')
    
    w = 1.0
    x = 0.0
    y = 0.0
    z = 0.0
    dpx = 0.0
    dpy = 0.0
    dpz = 1.0
    s = 0.1
    dDR = 0.0
    dDT = 0.0
    
    
    w, x, y, z, dpx, dpy, dpz, s, dDR, dDT = BoundaryCheck(w, x, y, z, dpx, dpy, dpz, s, dDR, dDT, d, n_t, n_m, alpha_critical)