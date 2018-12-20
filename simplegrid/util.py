
import numpy as np

# some useful constants:
edges = (N,S,E,W) = list(range(4))

def lonlat2cart( lons, lats, rad=1.):
    """Convert longitude/latitude to cartesian coordinates.

    Args:
        lons (numpy 2-d array): longitude ("x") values (decimal degrees).
        lats (numpy 2-d array): latitude ("y") values (decimal degrees).
        rad (float): nominal sphere radius (default=1. returns cartesian
            coordinates on the unit sphere).

    Returns:
        cart (numpy 3-d array): cartesian coordinates on the sphere.

    Raises:
        ValueError: If input lons and lats matrix dimensions are not equal.
    """

    if lons.ndim!=2 or lats.ndim!=2 or lons.shape!=lats.shape:
        raise ValueError('lons and lats must be two-dimensional matrices of equal size.')

    dims = list(lons.shape)
    dims.append(3)
    cart = np.zeros(dims)
    cart[:,:,0] = np.cos(np.radians(lons))*np.cos(np.radians(lats))
    cart[:,:,1] = np.sin(np.radians(lons))*np.cos(np.radians(lats))
    cart[:,:,2] = np.sin(np.radians(lats))

    return cart


def nearest(lon,lat,lons,lats,geod):
    """Determine indices of, and distance to, nearest lon/lat point.

    Args:
        lon (float): Longitude of search origin
        lat (float): Latitude of search origin
        lons (numpy 2-d array): Matrix of longitude values used in the nearest
            neighbour search.
        lats (numpy 2-d array): Matrix of latitude values used in the nearest
            neighbour search.
        geod (pyproj.Geod object): Geod to be used as basis for distance
            calculations.

    Returns:
        i (int): Zeros-based row index of nearest neighbour.
        j (int): Zeros-based colum index of nearest neighbour.
        dist (float): Great circle distance from (lon,lat) to
            (lons(i,j), lats(i,j))

    Raises:
        ValueError: If input lons and lats matrix dimensions are not
            equal.
    """

    if lons.ndim!=2 or lats.ndim!=2 or lons.shape!=lats.shape:
        raise ValueError('lons and lats must be two-dimensional matrices of equal size.')

    i,j,dist = -1,-1,1.e10

    it = np.nditer(lons,flags=['multi_index'])
    while not it.finished:
        (fwd_az,back_az,distance) = geod.inv(
            lon,lat,
            lons[it.multi_index[0],it.multi_index[1]],
            lats[it.multi_index[0],it.multi_index[1]])
        if distance<dist:
            i,j,dist = it.multi_index[0],it.multi_index[1],distance
        it.iternext()

    return i,j,dist


def squad_uarea( cart):
    """Compute areas for cartesian array of corner points on the unit sphere.
    
    Args:
        cart (numpy 3-d array): 2-d array of cartesian x,y,z corner points on the unit sphere.

    Returns:
        areas (numpy 2-d array): n-1 x m-1 array of cell areas.
    """

    area = np.zeros((np.size(cart,0)-1,np.size(cart,1)-1))

    for x_idx in np.arange(np.size(cart,0)-1):
        for y_idx in np.arange(np.size(cart,1)-1):

            # quadrilateral corners in counterclockwise direction:
            ptA = cart[x_idx  ,y_idx  ,:]
            ptB = cart[x_idx+1,y_idx  ,:]
            ptC = cart[x_idx+1,y_idx+1,:]
            ptD = cart[x_idx  ,y_idx+1,:]

            # interior angles of first subtriangle:
            at,bt,ct = ptA, ptB, ptC
            ca,cb,cc = np.dot(bt,ct), np.dot(at,ct), np.dot(at,bt)
            sa,sb,sc = np.sin(np.arccos(ca)), np.sin(np.arccos(cb)), np.sin(np.arccos(cc))
            A1 = np.arccos((ca-cb*cc)/(sb*sc))
            B1 = np.arccos((cb-ca*cc)/(sa*sc))
            C1 = np.arccos((cc-ca*cb)/(sa*sb))

            # interior angles of second subtriangle:
            at,bt,ct = ptA, ptC, ptD
            ca,cb,cc = np.dot(bt,ct), np.dot(at,ct), np.dot(at,bt)
            sa,sb,sc = np.sin(np.arccos(ca)), np.sin(np.arccos(cb)), np.sin(np.arccos(cc))
            A2 = np.arccos((ca-cb*cc)/(sb*sc))
            B2 = np.arccos((cb-ca*cc)/(sa*sc))
            C2 = np.arccos((cc-ca*cb)/(sa*sb))

            # area:
            area[x_idx,y_idx] = A1+B1+C1 + A2+B2+C2 - 2*np.pi

    return area

