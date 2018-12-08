
import numpy as np

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


def matchedges( aXG, aYG, bXG, bYG, geod, verbose=False):
    """Given two tiles, find their matching edges.

    Args:
        aXG, aYG (numpy arrays): Grid point longitudes and latitudes,
            respectively, for the first tile (mitgridfile format).
        bXG, bYG (numpy arrays): Grid point longitudes and latitudes,
            respectively, for the second tile (mitgridfile format).
        geod (pyproj.Geod object): Geod to be used as basis for distance
            calculations.
        verbose (logical): verbose output

    Returns:
        (a_edge,b_edge) tuple of slice object pairs that will yield the input
            tile matching edge grids, e.g, for an east-west match between 'a'
            and 'b', ((numpy.s_[-1],numpy.s_[:]),(numpy.s_[0],numpy.s_[:]))

    Raises:
        ValueErrors if matching edges are not found, are not sufficiently close,
        or are inconsistent (e.g., incorrect alignment)

    Note:
        Input tiles are assumed to have the same orientation.
        
    """

    # some useful parameters:
    edges   = (N,S,E,W)     = (1,0,1,0)                     # j_n,j_s,i_e,i_w
    corners = (SW,SE,NE,NW) = (0,0),(-1,0),(-1,-1),(0,-1)   # full matrix idxs
    CLOSE_ENOUGH = 1.e-6                                    # error checking
    a_rtn  = b_rtn = None                                   # return vals

    # determine nominal diagonal distance for error checking:
    _,_,a_swne_dist = geod.inv(aXG[SW],aYG[SW],aXG[NE],aYG[NE])
    _,_,a_senw_dist = geod.inv(aXG[SE],aYG[SE],aXG[NW],aYG[NW])
    _,_,b_swne_dist = geod.inv(bXG[SW],bYG[SW],bXG[NE],bYG[NE])
    _,_,b_senw_dist = geod.inv(bXG[SE],bYG[SE],bXG[NW],bYG[NW])
    nom_diag = (a_swne_dist+a_senw_dist+b_swne_dist+b_senw_dist)/4.

    # build 2x2 corner matrices of minimum distances between a and b:
    a_to_b_dist = np.full((2,2),np.inf)
    b_to_a_dist = a_to_b_dist.copy()

    for corner_a in corners:
        for corner_b in corners:
            _,_,dist = geod.inv(
                aXG[corner_a], aYG[corner_a], bXG[corner_b], bYG[corner_b])
            if dist < a_to_b_dist[corner_a]:
                # keep track of minimum distances without regard to the actual
                # corner in b that's closest...:
                a_to_b_dist[corner_a] = dist
                # ...and, by reciprocity:
                b_to_a_dist[corner_b] = a_to_b_dist[corner_a]

    # make sure minimums are sufficently close to zero:
    if np.amin([np.amin(a_to_b_dist),np.amin(b_to_a_dist)]) > nom_diag*CLOSE_ENOUGH:
        raise ValueError('Input tiles not sufficiently close to determine common edge')

    # edge detection:

    a_to_b_rowmins = np.argmin(a_to_b_dist,axis=0).tolist()
    a_to_b_colmins = np.argmin(a_to_b_dist,axis=1).tolist()
    if a_to_b_rowmins == [W,W]:
        a_edge, a_rtn = W, (np.s_[0],np.s_[:])
        if verbose:
            msg = "western edge of tile 'a' "
    elif a_to_b_rowmins == [E,E]:
        a_edge, a_rtn = E, (np.s_[-1],np.s_[:])
        if verbose:
            msg = "eastern edge of tile 'a' "
    elif a_to_b_colmins == [S,S]:
        a_edge, a_rtn = S, (np.s_[:],np.s_[0])
        if verbose:
            msg = "southern edge of tile 'a' "
    elif a_to_b_colmins == [N,N]:
        a_edge, a_rtn = N, (np.s_[:],np.s_[-1])
        if verbose:
            msg = "northern edge of tile 'a' "
    else:
        raise ValueError('No matching edge found for first input tile')

    b_to_a_rowmins = np.argmin(b_to_a_dist,axis=0).tolist()
    b_to_a_colmins = np.argmin(b_to_a_dist,axis=1).tolist()
    if b_to_a_rowmins == [W,W]:
        b_edge, b_rtn = W, (np.s_[0],np.s_[:])
        if verbose:
            msg = msg + "matches western edge of tile 'b'"
    elif b_to_a_rowmins == [E,E]:
        b_edge, b_rtn = E, (np.s_[-1],np.s_[:])
        if verbose:
            msg = msg + "matches eastern edge of tile 'b'"
    elif b_to_a_colmins == [S,S]:
        b_edge, b_rtn = S, (np.s_[:],np.s_[0])
        if verbose:
            msg = msg + "matches southern edge of tile 'b'"
    elif b_to_a_colmins == [N,N]:
        b_edge, b_rtn = N, (np.s_[:],np.s_[-1])
        if verbose:
            msg = msg + "matches northern edge of tile 'b'"
    else:
        raise ValueError('No matching edge found for second input tile')

    if a_edge+b_edge != 1:
        raise ValueError('Inconsistent edge match for input grids')

    if verbose:
        print(msg)

    return (a_rtn,b_rtn)


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

