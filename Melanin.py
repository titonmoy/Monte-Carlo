from functions.RandomWalk import RandomWalk

from random import seed
import numpy as np
from math import asin

from tabulate import tabulate
from matplotlib import pyplot as plt

seed(1)



def EffectOfMelaninConcentration():
    
    # optical properties of epidermis for 550 nm light
    sMelaninLevel = ['Low', 'Medium', 'High'] # melanin levels
    mu_as = np.array([1.63, 6.85, 16.13]) # absorption coefficient for different melanin levels
    mu_s = 74.7 # scattering coefficient
    g = 0.759 # anisotropy of scatter
    n_t = 1.447 # refrctive index of epidermis
    n_m = 1.5 # refractive index of surrounding media (coverslip)
    d = 0.1 # tickness of epidermis
    
    mu_atts = mu_as + mu_s # attenuation coefficient for different melanin levels
    absFacts = mu_as / mu_atts # energy absorption factor for different melanin levels
    
    # critical angle for total internal reflection
    if n_t > n_m:
            alpha_critical = asin(n_m / n_t)
    else:
        alpha_critical = float('nan')    
    
    N_bins = int(1e3) # number of bins in pDW
    dz = d / N_bins # distance between bins in pDW
    
    # variables for recording absorption, diffused reflectance, and diffused transmittance
    pDW = np.zeros((3, N_bins))
    dDR = np.zeros(3)
    dDT = np.zeros(3)
    
    N_photons = int(5e6) # number of photons to be launched    
    
    for n, (mu_att, absFact) in enumerate(zip(mu_atts, absFacts)):
        for nph in range(N_photons):
            dDR[n], dDT[n], pDW[n] = RandomWalk(mu_att, absFact, g, n_t, n_m, alpha_critical, d, dz, dDR[n], dDT[n], pDW[n])
    
    
    # calculating avergae absorption, diffused reflectance, and diffused transmittance
    dDR /= N_photons
    dDT /= N_photons
    dDW = np.sum(pDW, axis = 1) / N_photons
    
    print('\n')
    print(tabulate({'Melanin Level': sMelaninLevel, 'Diffuse Transmittance': dDT, 'Diffuse Reflectance': dDR, 'Absorption': dDW}, headers = 'keys'))
    print('\n')
    
    # plotting absorption depth profile
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
    
    # optical properties of epidermis for medium melanin concentration
    nWavelengths = np.array([350, 450, 550, 650]) # wavelengths of light
    mu_as = np.array([30.16, 13.5, 6.85, 3.68]) # absorption coefficient for different wavelengths
    mu_ss = np.array([210.4, 121.6, 74.7, 55.48]) # scattering coeffcient for different wavelengths
    gs = np.array([0.702, 0.728, 0.759, 0.787]) # anisotropy of scatter for different wavelengths
    n_ts = np.array([1.488, 1.448, 1.447, 1.433]) # refractive index of epidermis for different wavelengths
    n_m = 1.5 # refractive index of surrounding medium (coverslip)
    d = 0.1 # thickness of epidermis
    
    mu_atts = mu_as + mu_ss # attenuation coefficient for different wavelengths
    absFacts = mu_as / mu_atts # energy absorption factor for different wavelengths
    
    # critical angle for total internal reflection for different wavelengths
    alpha_criticals = np.zeros(4)
    for n, n_t in enumerate(n_ts):
        if n_t > n_m:
                alpha_criticals[n] = asin(n_m / n_t)
        else:
            alpha_criticals[n] = float('nan')
    
    N_bins = int(1e3) # number of bins in pDW
    dz = d / N_bins # distance between bins in pDW
    
    # variables for recording absorption, diffused reflectance, and diffused transmittance
    pDW = np.zeros((4, N_bins))
    dDR = np.zeros(4)
    dDT = np.zeros(4)
    
    N_photons = int(5e6) # number of photons to be launched
    
    for n, (mu_att, absFact, g, n_t, alpha_critical) in enumerate(zip(mu_atts, absFacts, gs, n_ts, alpha_criticals)):
        for nph in range(N_photons):
            dDR[n], dDT[n], pDW[n] = RandomWalk(mu_att, absFact, g, n_t, n_m, alpha_critical, d, dz, dDR[n], dDT[n], pDW[n])
    
    
    # calculating avergae absorption, diffused reflectance, and diffused transmittance
    dDR /= N_photons
    dDT /= N_photons
    dDW = np.sum(pDW, axis = 1) / N_photons
    
    print('\n')
    print(tabulate({'Wavelengths': nWavelengths, 'Diffuse Transmittance': dDT, 'Diffuse Reflectance': dDR, 'Absorption': dDW}, headers = 'keys'))
    print('\n')
    
    # plotting absorption depth profile
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
    
    # plotting absorption vs wavelength
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