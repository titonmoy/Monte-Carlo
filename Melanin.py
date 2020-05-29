from RandomWalk import RandomWalk

from random import seed
import numpy as np
from math import asin

from tabulate import tabulate
from matplotlib import pyplot as plt



seed(1)


def EffectOfMelaninConcentration():
    

    sMelaninLevel = ['Low', 'Medium', 'High']
    mu_as = np.array([1.63, 6.85, 16.13])
    mu_s = 74.7
    g = 0.759
    n_t = 1.447
    n_m = 1.5
    d = 0.1
    
    mu_atts = mu_as + mu_s
    absFacts = mu_as / mu_atts
    
    if n_t > n_m:
            alpha_critical = asin(n_m / n_t)
    else:
        alpha_critical = float('nan')
    
    N_bins = int(1e3)
    dz = d / N_bins
    
    N_photons = int(1e4)
    
    pDW = np.zeros((3, N_bins))
    dDR = np.zeros(3)
    dDT = np.zeros(3)
    
    
    for n, (mu_att, absFact) in enumerate(zip(mu_atts, absFacts)):
        for nph in range(N_photons):
            dDR[n], dDT[n], pDW[n] = RandomWalk(mu_att, absFact, g, n_t, n_m, alpha_critical, d, dz, dDR[n], dDT[n], pDW[n])
    
    
    
    dDR /= N_photons
    dDT /= N_photons
    dDW = np.sum(pDW, axis = 1) / N_photons
    
    print('\n')
    print(tabulate({'Melanin Level': sMelaninLevel, 'Diffuse Transmittance': dDT, 'Diffuse Reflectance': dDR, 'Absorption': dDW}, headers = 'keys'))
    print('\n')
    
    pz = np.arange(0, d, dz)
    pDW = pDW / (dz * N_photons)
    plt.figure()
    plt.semilogy(pz, pDW[0], 'r', label = sMelaninLevel[0])
    plt.semilogy(pz, pDW[1], 'g', label = sMelaninLevel[1])
    plt.semilogy(pz, pDW[2], 'b', label = sMelaninLevel[2])
    plt.xlabel('Depth (cm)')
    plt.ylabel('Absorption (cm$^{-1}$)')
    plt.legend()



def EffectOfWavelength():
    
    
    nWavelengths = np.array([350, 450, 550, 650])
    mu_as = np.array([30.16, 13.5, 6.85, 3.68])
    mu_ss = np.array([210.4, 121.6, 74.7, 55.48])
    gs = np.array([0.702, 0.728, 0.759, 0.787])
    n_ts = np.array([1.488, 1.448, 1.447, 1.433])
    n_m = 1.5
    d = 0.1
    
    mu_atts = mu_as + mu_ss
    absFacts = mu_as / mu_atts
    
    alpha_criticals = np.zeros(4)
    for n, n_t in enumerate(n_ts):
        if n_t > n_m:
                alpha_criticals[n] = asin(n_m / n_t)
        else:
            alpha_criticals[n] = float('nan')
    
    N_bins = int(1e3)
    dz = d / N_bins
    
    N_photons = int(1e4)
    
    pDW = np.zeros((4, N_bins))
    dDR = np.zeros(4)
    dDT = np.zeros(4)
    
    for n, (mu_att, absFact, g, n_t, alpha_critical) in enumerate(zip(mu_atts, absFacts, gs, n_ts, alpha_criticals)):
        for nph in range(N_photons):
            dDR[n], dDT[n], pDW[n] = RandomWalk(mu_att, absFact, g, n_t, n_m, alpha_critical, d, dz, dDR[n], dDT[n], pDW[n])
    
    dDR /= N_photons
    dDT /= N_photons
    dDW = np.sum(pDW, axis = 1) / N_photons
    
    print('\n')
    print(tabulate({'Wavelengths': nWavelengths, 'Diffuse Transmittance': dDT, 'Diffuse Reflectance': dDR, 'Absorption': dDW}, headers = 'keys'))
    print('\n')
    
    pz = np.arange(0, d, dz)
    pDW = pDW / (dz * N_photons)
    plt.figure()
    plt.semilogy(pz, pDW[0], 'r', label = str(nWavelengths[0]) + ' nm')
    plt.semilogy(pz, pDW[1], 'g', label = str(nWavelengths[1]) + ' nm')
    plt.semilogy(pz, pDW[2], 'b', label = str(nWavelengths[2]) + ' nm')
    plt.semilogy(pz, pDW[3], 'k', label = str(nWavelengths[3]) + ' nm')
    plt.xlabel('Depth (cm)')
    plt.ylabel('Absorption (cm$^{-1}$)')
    plt.legend()
    
    
    p = np.polyfit(nWavelengths, dDW, 3)
    nWaves = np.arange(350, 660, 10)
    fitted = np.polyval(p, nWaves)
    plt.figure()
    plt.plot(nWavelengths, dDW, 'ro', label = 'Original')
    plt.plot(nWaves, fitted, 'b', label = 'Fitted')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Absorption')
    plt.legend()




if __name__ == '__main__':
    EffectOfMelaninConcentration()
    EffectOfWavelength()