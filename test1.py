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
import numpy as np

seed(1)


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


pDW = np.zeros((3, N_bins + 1))
dDR = np.zeros(3)
dDT = np.zeros(3)



epsilon = 1e-9



for n, (mu_att, absFact) in enumerate(zip(mu_atts, absFacts)):
    for nph in range(N_photons):
        
        w = 1
    
        x = 0
        y = 0
        z = 0
        
        dpx = 0
        dpy = 0
        dpz = 1
        
        
        while True:
            
        
            s = - log(random()) / mu_att
            
            
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
                    
                    alpha_t = asin(n_t * sin(alpha_i) / n_m)
                    
                    R = 0.5 * ((sin(alpha_i - alpha_t))**2 / ((sin(alpha_i + alpha_t))**2 + epsilon) + (tan(alpha_i - alpha_t))**2 / ((tan(alpha_i + alpha_t))**2 + epsilon))
                    
                    if R < random():
                        if dpz < 0:
                            dDR[n] += w
                        else:
                            dDT[n] += w
                        
                        w = 0
                    
                    dpz = - dpz
            
            
            if w <= 0:
                break
            
            
            x += s * dpx
            y += s * dpy
            z += s * dpz
            
            
            dw = w * absFact
            if (w - dw) < 0:
                dw = w
            w -= dw
            
            nz = floor(z / dz)
            pDW[n, nz] += dw
            
            
            if w <= 0:
                break
            
            
            
            
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
            
