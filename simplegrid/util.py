
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
        cart (numpy 3-d array): cartesian coordinates on the sphere.  Rows and
            columns are arranged according to the input lons/lats arrays, while
            along the third axis, the 0, 1, and 2 indices correspond to the x, y
            and z components, respectively.

    Raises:
        ValueError: If input lons and lats matrix dimensions are not equal.

    """

    if lons.ndim!=2 or lats.ndim!=2 or lons.shape!=lats.shape:
        raise ValueError('lons and lats must be two-dimensional matrices of equal size.')

    dims = list(lons.shape)
    dims.append(3)
    cart = np.zeros(dims)
    cart[:,:,0] = np.cos(np.radians(lons))*np.cos(np.radians(lats)) # x
    cart[:,:,1] = np.sin(np.radians(lons))*np.cos(np.radians(lats)) # y
    cart[:,:,2] = np.sin(np.radians(lats))                          # z

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
    """Compute quadrilateral surface areas for cartesian array of corner points
    on the unit sphere.

    Args:
        cart (numpy array): 3-d array (m x n x 3) of cartesian x,y,z corner
            points on the unit sphere (for every (i,j), (x,y,z) = (i,j,0),
            (i,j,1), (i,j,2)).

    Returns:
        areas (numpy array): 2-d array (n-1 x m-1) of cell areas.

    Note:
        One of many possible approaches, the algorithm implemented here is based
        on Girard's spherical excess formula, with direct calculation of angles
        using the spherical law of cosines.  Note that, due to numerical
        round-off, its results can be inaccurate for small angles/areas (edge
        lengths less than roughly 5km on the scaled sphere).  In such cases,
        pquad_uarea is recommended.

    """

    # only because some geometries are poorly conditioned:
    import warnings
    warnings.filterwarnings("ignore")

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


def pquad_uarea( cart):
    """ Compute planar quadrilateral (faceted) areas for cartesian array of
    corner points on the unit sphere.

    Args:
        cart (numpy array): 3-d array (m x n x 3) of cartesian x,y,z corner
            points on the unit sphere (for every (i,j), (x,y,z) = (i,j,0),
            (i,j,1), (i,j,2)).

    Returns:
        areas (numpy array): 2-d array (n-1 x m-1) of cell areas.

    Note:
        The algorithm implemented here computes areas using edge vector cross
        products and is recommented for small angles/areas (edge lengths less
        than roughly 5km on the scaled sphere).  If areas are larger,
        squad_uarea is recommended.

    """

    area = np.zeros((np.size(cart,0)-1,np.size(cart,1)-1))

    for x_idx in np.arange(np.size(cart,0)-1):
        for y_idx in np.arange(np.size(cart,1)-1):

            # quadrilateral corners in counterclockwise direction:
            ptA = cart[x_idx  ,y_idx  ,:]
            ptB = cart[x_idx+1,y_idx  ,:]
            ptC = cart[x_idx+1,y_idx+1,:]
            ptD = cart[x_idx  ,y_idx+1,:]

            # edge vectors:
            ab = ptB-ptA
            ac = ptC-ptA
            ad = ptD-ptA

            area[x_idx,y_idx] = 0.5 * np.linalg.norm(np.cross(ab,ac)+np.cross(ac,ad))

    return area

