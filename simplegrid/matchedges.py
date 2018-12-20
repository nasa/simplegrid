
import numpy as np
import pyproj
from . import computegrid

from .util import (N,S,E,W)

def matchedges( aXG, aYG, bXG, bYG, geod, verbose=False):

    """Determine whether or not a common edge exists between two tiles.

    Args:
        aXG, aYG (numpy arrays): Grid point longitudes and latitudes,
            respectively, for the first tile (mitgridfile format).
        bXG, bYG (numpy arrays): Grid point longitudes and latitudes,
            respectively, for the second tile (mitgridfile format).
        geod (pyproj.Geod object): Geod to be used as basis for distance
            calculations.
        verbose (logical): verbose output

    Returns:
        (a_edge, a_edge_slice, b_edge, b_edge_slice,
        compute_grid_edge_xg, compute_grid_edge_yg,
        compute_grid_edge_join_slice):
            tuple of tile edge indicators, slice objects, and interpolated
            compute edge values and corresponding slice object, that can
            facilitate subsequent tile join and edge value interpolation
            operations.
            - a_edge: indicator denoting the matching edge of tile 'A' (N,S,E,
                or W)
            - a_edge_slice: slice operator on tile 'A' (aXG, aYG) that can be
                used to produce edge in common with tile 'B'
            - b_edge: indicator denoting the matching edge of tile 'B' (N,S,E,
                or W)
            - b_edge_slice: slice operator on tile 'B' (bXG, bYG) that can be
                used to produce edge in common with tile 'A'. If tiles match
                exactly, note that aXG[a_edge_slice] = bXG[b_edge_slice] and
                aYG[a_edge_slice] = bYG[b_edge_slice]
            - compute_edge_xg: vector of grid X (longitude) values that have
                been interpolated from the midpoints of the tile 'B' cells
                bordering tile 'A', the purpose of which is to provide compute
                grid boundary edge values for updates to tile 'A'.
            - compute_edge_yg - same as compute_edge_xg, but with respect to Y
                (latitude) values
            - compute_grid_edge_join_slice: slice object that can be used to
                place compute_edge_xg and compute_edge_yg into a compute grid
                for tile 'A', e.g.,
                    compute_grid_xg_A[compute_grid_edge_join_slice] = \\
                        compute_edge_xg

    Raises:
        ValueErrors if matching edges are not found, are not sufficiently close,
        or are inconsistent (e.g., incorrect alignment)

    Note:
        Input tiles are assumed to have the same orientation.
        
    """

    # some useful parameters:
    corners = (SW,SE,NE,NW) = (0,0),(-1,0),(-1,-1),(0,-1)   # full matrix idxs
    CLOSE_ENOUGH = 1.e-6                                    # error checking

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

    for corner_b in corners:
        for corner_a in corners:
            _,_,dist = geod.inv(
                aXG[corner_a], aYG[corner_a], bXG[corner_b], bYG[corner_b])
            if dist < b_to_a_dist[corner_b]:
                # keep track of minimum distances without regard to the actual
                # corner in a that's closest...:
                b_to_a_dist[corner_b] = dist

    # make sure minimums are sufficently close to zero:
    if np.amin([np.amin(a_to_b_dist),np.amin(b_to_a_dist)]) > nom_diag*CLOSE_ENOUGH:
        raise ValueError('Input tiles not sufficiently close to determine common edge')

    # edge detection:

    a_to_b_closest = np.less_equal(a_to_b_dist,nom_diag*CLOSE_ENOUGH)
    b_to_a_closest = np.less_equal(b_to_a_dist,nom_diag*CLOSE_ENOUGH)

    if np.all(a_to_b_closest[0,:]):
        a_edge          = W
        a_edge_slice    = (np.s_[0],np.s_[:])
        if verbose:
            msg = "western edge of tile 'a' "
    elif np.all(a_to_b_closest[-1,:]):
        a_edge          = E
        a_edge_slice    = (np.s_[-1],np.s_[:])
        if verbose:
            msg = "eastern edge of tile 'a' "
    elif np.all(a_to_b_closest[:,0]):
        a_edge          = S
        a_edge_slice    = (np.s_[:],np.s_[0])
        if verbose:
            msg = "southern edge of tile 'a' "
    elif np.all(a_to_b_closest[:,-1]):
        a_edge          = N
        a_edge_slice    = (np.s_[:],np.s_[-1])
        if verbose:
            msg = "northern edge of tile 'a' "
    else:
        raise ValueError('No matching edge found for first input tile')

    if np.all(b_to_a_closest[0,:]):
        b_edge              = W
        b_edge_slice        = (np.s_[0],np.s_[:])
        b_edge_next_slice   = (np.s_[1],np.s_[:])
        if verbose:
            msg = msg + "matches western edge of tile 'b'"
    elif np.all(b_to_a_closest[-1,:]):
        b_edge              = E
        b_edge_slice        = (np.s_[-1],np.s_[:])
        b_edge_next_slice   = (np.s_[-2],np.s_[:])
        if verbose:
            msg = msg + "matches eastern edge of tile 'b'"
    elif np.all(b_to_a_closest[:,0]):
        b_edge              = S
        b_edge_slice        = (np.s_[:],np.s_[0])
        b_edge_next_slice   = (np.s_[:],np.s_[1])
        if verbose:
            msg = msg + "matches southern edge of tile 'b'"
    elif np.all(b_to_a_closest[:,-1]):
        b_edge              = N
        b_edge_slice        = (np.s_[:],np.s_[-1])
        b_edge_next_slice   = (np.s_[:],np.s_[-2])
        if verbose:
            msg = msg + "matches northern edge of tile 'b'"
    else:
        raise ValueError('No matching edge found for second input tile')

    if (a_edge+b_edge!=N+S) and (a_edge+b_edge!=E+W):
        raise ValueError('Inconsistent edge match for input grids')

    if verbose:
        print(msg)

    # use edge information to make a single boundary edge only "compute grid",
    # subdivide by two, and extract the middle row or column, as the case may
    # be, along with an assembly slice operator, so that it may be joined to
    # tile A's appropriate edge for subsequent mitgrid operations:

    if a_edge==E or a_edge==W:
        _,xgyg_cols = aXG.shape
        compute_grid_xg = np.zeros((3,2*xgyg_cols-1))
        compute_grid_yg = np.zeros((3,2*xgyg_cols-1))
        ilb,iub = 0,2
        jlb,jub = 0,2*(xgyg_cols-1)
        if a_edge==E:
            for i in range(xgyg_cols):
                compute_grid_xg[0,2*i] = aXG[a_edge_slice][i]
                compute_grid_yg[0,2*i] = aYG[a_edge_slice][i]
                compute_grid_xg[2,2*i] = bXG[b_edge_next_slice][i]
                compute_grid_yg[2,2*i] = bYG[b_edge_next_slice][i]
        elif a_edge==W:
            for i in range(xgyg_cols):
                compute_grid_xg[2,2*i] = aXG[a_edge_slice][i]
                compute_grid_yg[2,2*i] = aYG[a_edge_slice][i]
                compute_grid_xg[0,2*i] = bXG[b_edge_next_slice][i]
                compute_grid_yg[0,2*i] = bYG[b_edge_next_slice][i]

    elif a_edge==N or a_edge==S:
        xgyg_rows,_ = aXG.shape
        compute_grid_xg = np.zeros((2*xgyg_rows-1,3))
        compute_grid_yg = np.zeros((2*xgyg_rows-1,3))
        ilb,iub = 0,2*(xgyg_rows-1)
        jlb,jub = 0,2
        if a_edge==N:
            for i in range(xgyg_rows):
                compute_grid_xg[2*i,0] = aXG[a_edge_slice][i]
                compute_grid_yg[2*i,0] = aYG[a_edge_slice][i]
                compute_grid_xg[2*i,2] = bXG[b_edge_next_slice][i]
                compute_grid_yg[2*i,2] = bYG[b_edge_next_slice][i]
        elif a_edge==S:
            for i in range(xgyg_rows):
                compute_grid_xg[2*i,2] = aXG[a_edge_slice][i]
                compute_grid_yg[2*i,2] = aYG[a_edge_slice][i]
                compute_grid_xg[2*i,0] = bXG[b_edge_next_slice][i]
                compute_grid_yg[2*i,0] = bYG[b_edge_next_slice][i]

    # produce subdivided compute grid:

    (compute_grid_xg,compute_grid_yg) = computegrid.fill(
        compute_grid_xg, compute_grid_yg, ilb, iub, jlb, jub,
        1, 1, geod, verbose)

    # extract middle row or column, as the case may be, along with tile A
    # compute grid assembly slice operator:

    if a_edge==E or a_edge==W:
        compute_grid_edge_xg = compute_grid_xg[1,:]
        compute_grid_edge_yg = compute_grid_yg[1,:]
        if a_edge==E:
            compute_grid_edge_join_slice = (np.s_[-1],np.s_[1:-1]) 
        elif a_edge==W:
            compute_grid_edge_join_slice = (np.s_[0],np.s_[1:-1])

    elif a_edge==N or a_edge==S:
        compute_grid_edge_xg = compute_grid_xg[:,1]
        compute_grid_edge_yg = compute_grid_yg[:,1]
        if a_edge==N:
            compute_grid_edge_join_slice = (np.s_[1:-1],np.s_[-1])
        elif a_edge==S:
            compute_grid_edge_join_slice = (np.s_[1:-1],np.s_[0])

    return (
        a_edge, a_edge_slice,
        b_edge, b_edge_slice,
        compute_grid_edge_xg, compute_grid_edge_yg,
        compute_grid_edge_join_slice)

