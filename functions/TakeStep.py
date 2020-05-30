def TakeStep(x, y, z, dpx, dpy, dpz, s):
    '''
    Take random step within tissue.

    Parameters
    ----------
    x : double
        Current x position of photon.
    y : double
        Current y position of photon.
    z : double
        Current z position of photon.
    dpx : double
        Direction of propagation along the x direction.
    dpy : double
        Direction of propagation along the y direction.
    dpz : double
        Direction of propagation along the z direction.
    s : double
        Step size

    Returns
    -------
    x : double
        Updated x position of photon.
    y : double
        Updated y position of photon.
    z : double
        Updated z position of photon.

    '''
    
    x += s * dpx
    y += s * dpy
    z += s * dpz
    
    return x, y, z



if __name__ == '__main__':
    
    x = 0.0
    y = 0.0
    z = 0.0
    dpx = 0.0
    dpy = 0.0
    dpz = 1.0
    s = 0.1
    
    x, y, z = TakeStep(x, y, z, dpx, dpy, dpz, s)