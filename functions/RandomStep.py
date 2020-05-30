from random import random
from math import log

def RandomStep(mu_att):
    '''
    Generate a random step size for the photon.

    Parameters
    ----------
    mu_att : double
        Attenuation coefficient of tissue.
        Sum of absorption coefficient and scattering coefficient.

    Returns
    -------
    s : double
        Random step size.

    '''
    
    s = - log(random()) / mu_att
    
    return s



if __name__ == "__main__":
    
    mu_a = 1.0 # aborption coefficient
    mu_s = 1.0 # scattering coefficient
    mu_att = mu_a + mu_s
    s = RandomStep(mu_att);