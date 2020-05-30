from random import random
from math import pi
from math import sin
from math import cos
from math import acos
from math import copysign
from math import sqrt

def ScatterPhoton(dpx, dpy, dpz, g):
    '''
    Scatter photon and update direction of propagation accordingly.

    Parameters
    ----------
    dpx : double
        Current direction of propagation along the x direction.
    dpy : double
        Current direction of propagation along the y direction.
    dpz : double
        Current direction of propagation along the z direction.
    g : double
        Anisotropy of scatter

    Returns
    -------
    dpx : double
        Updated direction of propagation along the x direction.
    dpy : double
        Updated direction of propagation along the y direction.
    dpz : double
        Updated direction of propagation along the z direction.

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



if __name__ == '__main__':
    
    g = 0.5
    
    dpx = 0.0
    dpy = 0.0
    dpz = 1.0
    
    dpx, dpy, dpz = ScatterPhoton(dpx, dpy, dpz, g)