import numpy as np
import random as rnd
def radial_profile(data, center):
    """
    Calculate the radial profile of a map defined the
    center, results are done in pixel scales
    """
    x, y = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)
    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile

def pop_sources(clid, fov, sources_density):
    """
    Populate the field of view of the map and sample it
    with point sources with a given density in square arcmin
    fov in deg
    clid the seed of the random number generator
    """
    arcsec = fov*3600.0
    arcmin2 = arcsec*arcsec/60.0/60.0
    nsources = sources_density*arcmin2
    nsources = int(nsources)
    t = np.zeros([nsources])
    xs = np.zeros([nsources])
    ys = np.zeros([nsources])
    rnd.seed(clid)
    for i in range(0, nsources):
        a = rnd.uniform(-0.5, 0.5)
        b = rnd.uniform(-0.5, 0.5)
        xs[i] = a
        ys[i] = b
    t = np.arctan(xs/ys)
    print('   - - -  ')
    print('   Sampling the fov with point like sources  ')
    print('   - - -  ')
    print('   - - -  ')
    print(' min and max values of the sources in the fov [-0.5,0,5] ')
    print(np.amin(xs), np.amax(xs))
    print(np.amin(ys), np.amax(ys))
    print('   - - -  ')
    #...center is the center of the map
    rs = np.sqrt(xs*xs+ys*ys)*fov
    print( ' min and max values from the center of the map')
    print(np.amin(rs), np.amax(rs))
    print('   - - -  ')
    return xs, ys, t, rs

def bin_prof(x, y, nbins, ngal, sigmagal):
    """
    Bin the profile, nbins set the number of bins
    ngal number of background galaxies for the error
    sigmagal is the intrinsic ellipticity distribution
    """
    x = np.array(x)
    y = np.array(y)
    y2 = y*y
    if nbins>0:
        a, b = np.histogram(x, bins=nbins)
        c, b = np.histogram(x, bins=nbins, weights=y)
        e, d = np.histogram(x, bins=nbins, weights=y2)
    else:
        xx = np.logspace(np.log10(x[0]), np.log10(x[-1]), -nbins)
        a, b0 = np.histogram(x, bins=xx)
        c, b = np.histogram(x, bins=xx, weights=y)
        e, d = np.histogram(x, bins=xx, weights=y2)
    x = ((b[1:] + b[:-1])/2.0)
    s = sigmagal/np.sqrt(3600*np.pi*ngal*(b0[1:]**2-b0[:-1]**2))
    x = x[a != 0]
    c = c[a != 0]
    s = s[a != 0]
    e = e[a != 0]
    a = a[a != 0]
    y = (c / a)
    e = np.sqrt(e/a - y*y)
    s = np.sqrt(s*s+e*e)
    return x, y, s

def bin_prof2(x, y, nbins, ngal, sigmagal):
    """
    Bin the profile, nbins set the number of bins
    ngal number of background galaxies for the error
    sigmagal is the intrinsic ellipticity distribution
    """
    x = np.array(x)
    y = np.array(y)
    y2 = y*y
    if nbins>0:
        a, b = np.histogram(x, bins=nbins)
        c, b = np.histogram(x, bins=nbins, weights=y)
        e, d = np.histogram(x, bins=nbins, weights=y2)
    else:
        xx = np.logspace(np.log10(x[0]), np.log10(x[-1]), -nbins)
        a, b0 = np.histogram(x, bins=xx)
        c, b = np.histogram(x, bins=xx, weights=y)
        e, d = np.histogram(x, bins=xx, weights=y2)
    x = ((b[1:] + b[:-1])/2.0)
    s = sigmagal/np.sqrt(3600*np.pi*ngal*(b0[1:]**2-b0[:-1]**2))
    x = x[a != 0]
    c = c[a != 0]
    s = s[a != 0]
    e = e[a != 0]
    a = a[a != 0]
    y = (c / a)
    e = np.sqrt(e/a - y*y)
    #s = np.sqrt(s*s+e*e)
    return x, y, s, e
