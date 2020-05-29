def LaunchPhoton():
    '''
    Launch a collimated photon from the origin directing downward in the tissue.

    Returns
    -------
    w : double
        Energy of photon.
    x : double
        x position of photon.
    y : double
        y position of photon.
    z : double
        z position of photon.
    dpx : double
        Direction of propagation along the x direction.
    dpy : double
        Direction of propagation along the y direction.
    dpz : double
        Direction of propagation along the z direction.

    '''
    
    w = 1.0
    
    x = 0.0
    y = 0.0
    z = 0.0
    
    dpx = 0.0
    dpy = 0.0
    dpz = 1.0
    
    return w, x, y, z, dpx, dpy, dpz



if __name__ == "__main__":
    
    w, x, y, z, dpx, dpy, dpz = LaunchPhoton()