from math import floor

def AbsorbPhoton(w, z, pDW, absFact, dz):
    '''
    Register absorption of photon energy at current depth within tissue.

    Parameters
    ----------
    w : double
        Current energy of photon.
    z : double
        Current z position of photon.
    pDW : array of doubles
        Current energy absorption at various depth.
    absFact : double
        Energy absorption factor.
    dz : double
        Distance between points in pDW array.

    Returns
    -------
    w : double
        Updated energy of photon.
    pDW : array of double
        Updated energy absorption at various depth.

    '''
    
    dw = w * absFact
    
    if (w - dw) < 0:
        dw = w
    w -= dw
    
    nz = floor(z / dz)
    pDW[nz] += dw
    
    return w, pDW



if __name__ == '__main__':
    
    import numpy as np
    
    mu_a = 1.0 # aborption coefficient
    mu_s = 1.0 # scattering coefficient
    d = 1.0 # tissue thickness
    N_bins = 1000 # number of bins
    
    mu_att = mu_a + mu_s
    absFact = mu_a / mu_att
    dz = d / N_bins
    
    w = 1.0
    z = 0.0
    pDW = np.zeros(N_bins)

    w, pDW = AbsorbPhoton(w, z, pDW, absFact, dz)