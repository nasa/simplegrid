
import numpy as np
import pyproj
from . import util

def cgareas( compute_grid_xg, compute_grid_yg, geod, verbose=False):
    """Using spherical geoid assumption, compute sub-areas for the compute grid.

    Args:
        compute_grid_xg, compute_grid_yg (numpy arrays): Arrays which have been
            dimensioned according to lon_subscale*2 and lat_subscale*2 (possibly
            with "boundary ring" of elements, depending on the application), and
            with (x,y) values filled in according to cgfill().
        geod (pyproj.Geod): current geoid instance (only the semi-major, or
            equatorial axis radius, geod.a is used here).
        verbose (logical): verbose output

    Returns:
        compute_areas: numpy array of sub-areas for the compute grid.

    Note:
        If any of the grid corner points are undefined (NaN), the corresponding
        area will also be undefined.

    """

    compute_areas = util.squad_uarea(
        util.lonlat2cart(compute_grid_xg,compute_grid_yg)) \
        * np.power(geod.a,2)

    if verbose:
        print('compute_areas:')
        print(compute_areas)

    return compute_areas


def cgedges( compute_grid_xg, compute_grid_yg, geod, verbose=False):
    """
    Use great circle arcs to compute distances between all compute_grid points

    Args:
        compute_grid_xg, compute_grid_yg (numpy arrays): Arrays which have been
            dimensioned according to lon_subscale*2 and lat_subscale*2 (possibly
            with "boundary ring" of elements, depending on the application), and
            with (x,y) values filled in according to cgfill().
        geod (pyproj.Geod): current geoid instance (only the semi-major, or
            equatorial axis radius, geod.a is used here).
        verbose (logical): verbose output

    Returns:
        (compute_edges_x, compute_edges_y) tuple of numpy arrays of grid edge lengths.
        compute_edges_x(i,j) contains the edge length from compute_grid_x /
        compute_grid_y (i,j) to (i+1,j) and compute_edges_y(i,j) contains the
        edge length from compute_grid_x / compute_grid_y (i,j) to (i,j+1)

    Note:
        If either of the grid corner points is undefined (NaN), the
        corresponding edge length will also be undefined.

    """


    # initialization can be based on either *_xg or *_yg since they're the same
    # shape:
    (cg_rows,cg_cols) = compute_grid_xg.shape
    compute_edges_x = np.empty((cg_rows-1,cg_cols  ))*np.nan
    compute_edges_y = np.empty((cg_rows  ,cg_cols-1))*np.nan

    # x-direction edge lengths:

    it = np.nditer(compute_edges_x,flags=['multi_index'],op_flags=['readwrite'])
    while not it.finished:
        lon1 = compute_grid_xg[it.multi_index[0]  ,it.multi_index[1]  ]
        lat1 = compute_grid_yg[it.multi_index[0]  ,it.multi_index[1]  ]
        lon2 = compute_grid_xg[it.multi_index[0]+1,it.multi_index[1]  ]
        lat2 = compute_grid_yg[it.multi_index[0]+1,it.multi_index[1]  ]
        if not any(np.isnan((lon1,lat1,lon2,lat2))):
            _,_,compute_edges_x[it.multi_index[0],it.multi_index[1]] = geod.inv(
                lon1,lat1,lon2,lat2)
        it.iternext()

    # y-direction edge lengths:

    it = np.nditer(compute_edges_y,flags=['multi_index'],op_flags=['readwrite'])
    while not it.finished:
        lon1 = compute_grid_xg[it.multi_index[0]  ,it.multi_index[1]  ]
        lat1 = compute_grid_yg[it.multi_index[0]  ,it.multi_index[1]  ]
        lon2 = compute_grid_xg[it.multi_index[0]  ,it.multi_index[1]+1]
        lat2 = compute_grid_yg[it.multi_index[0]  ,it.multi_index[1]+1]
        if not any(np.isnan((lon1,lat1,lon2,lat2))):
            _,_,compute_edges_y[it.multi_index[0],it.multi_index[1]] = geod.inv(
                lon1,lat1,lon2,lat2)
        it.iternext()

    if verbose:
        print('compute_edges_x:')
        print( compute_edges_x)
        print('compute_edges_y:')
        print( compute_edges_y)

    return (compute_edges_x,compute_edges_y)


def cgfill( compute_grid_xg, compute_grid_yg, ilb, iub, jlb, jub,
    lon_subscale, lat_subscale, geod, verbose=False):
    """Use great circle subdivisions to fill in compute grid intermediate points
    according to x/y subdivision levels and index ranges.

    Args:
        compute_grid_xg, compute_grid_yg (numpy arrays): Arrays which have been
            dimensioned according to lon_subscale*2 and lat_subscale*2 (possibly
            with "boundary ring" of elements, depending on the application), and
            with existing XG, YG, grid values already mapped to their
            corresponding positions.
        ilb, iub (int): compute_grid x-direction indices that define range of
            fill-in.
        jlb, jub (int): compute_grid y-direction indices that define range of
            fill-in.
        lon_subscale (int): user-specified x-direction cell subdivision level.
        lat_subscale (int): user-specified y-direction cell subdivision level.
        geod (pyproj.Geod): current geoid instance.
        verbose (logical): verbose output

    Returns:
        (compute_grid_xg_out, compute_grid_yg_out) tuple of input arrays with
        intermediate points filled in.

    Note:
        If either of the endpoints is undefined (NaN), the corresponding
        intermediate grids will also be undefined.

    """

    # explicit output quantities to avoid function "side effects":
    compute_grid_xg_out = np.copy(compute_grid_xg)
    compute_grid_yg_out = np.copy(compute_grid_yg)

    # since matrix slice operations are from start:(stop-1):stride, define a
    # lambda function to make it clear when the stop point should be inclusive,
    # i.e., start:incl(stop):stride
    incl = lambda idx : idx+1

    #
    # Step 1: x-edge fill-in according to user-specified subdivision level:
    #

    # compute_grid range, stride:
    cg_first_i  = ilb
    cg_last_i   = iub               # (not inclusive)
    cg_stride_i = 2*lon_subscale
    cg_first_j  = jlb
    cg_last_j   = incl(jub)         # (inclusive)
    cg_stride_j = 2*lat_subscale

    # transformation from partitioned, strided nditer space (i_n,j_n) to
    # underlying compute_grid space (i_cg,j_cg):
    i_cg = lambda i_n,ilb,lon_subscale : ilb + i_n*2*lon_subscale
    j_cg = lambda j_n,jlb,lat_subscale : jlb + j_n*2*lat_subscale

    it = np.nditer(
        [compute_grid_xg_out[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_j],
         compute_grid_yg_out[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_j]],
        flags=['multi_index'])
    while not it.finished:
        # compute equally-spaced x-edge subdivisions:
        lon1 = it[0]
        lat1 = it[1]
        lon2 = compute_grid_xg_out[
            i_cg(it.multi_index[0]+1,ilb,lon_subscale),
            j_cg(it.multi_index[1]  ,jlb,lat_subscale)]
        lat2 = compute_grid_yg_out[
            i_cg(it.multi_index[0]+1,ilb,lon_subscale),
            j_cg(it.multi_index[1]  ,jlb,lat_subscale)]
        if not any(np.isnan((lon1,lat1,lon2,lat2))):
            x_edge_subdivided_lonlats = geod.npts(
                lon1,lat1,lon2,lat2,
                lon_subscale*2-1)                   # n intermediate points
            # ...and store to updated compute_grid:
            compute_grid_xg_out[
                i_cg(it.multi_index[0],ilb,lon_subscale)+1:
                i_cg(it.multi_index[0]+1,ilb,lon_subscale),
                j_cg(it.multi_index[1],jlb,lat_subscale)] = \
                np.array(x_edge_subdivided_lonlats)[:,0]
            compute_grid_yg_out[
                i_cg(it.multi_index[0],ilb,lon_subscale)+1:
                i_cg(it.multi_index[0]+1,ilb,lon_subscale),
                j_cg(it.multi_index[1],jlb,lat_subscale)] = \
                np.array(x_edge_subdivided_lonlats)[:,1]
        it.iternext()

    if verbose:
        print('compute_grid after x-edge subdivision:')
        print('compute_grid_xg:')
        print(compute_grid_xg_out)
        print('compute_grid_yg:')
        print(compute_grid_yg_out)

    #
    # Step 2: y-direction fill-in for every x-direction subdivision:
    #

    # compute_grid range, stride:
    cg_first_i  = ilb
    cg_last_i   = incl(iub)         # (inclusive)
    cg_stride_i = 1
    cg_first_j  = jlb
    cg_last_j   = jub               # (not inclusive)
    cg_stride_j = 2*lat_subscale

    # transformation from partitioned, strided nditer space (i_n,j_n) to
    # underlying compute_grid space (i_cg,j_cg) (j_cg same as above):
    i_cg = lambda i_n,ilb : ilb + i_n

    it = np.nditer(
        [compute_grid_xg_out[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_j],
         compute_grid_yg_out[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_j]],
        flags=['multi_index'])
    while not it.finished:
        # compute equally-spaced y-direction subdivisions...:
        lon1 = it[0]
        lat1 = it[1]
        lon2 = compute_grid_xg_out[
            i_cg(it.multi_index[0],ilb),
            j_cg(it.multi_index[1]+1,jlb,lat_subscale)]
        lat2 = compute_grid_yg_out[
            i_cg(it.multi_index[0],ilb),
            j_cg(it.multi_index[1]+1,jlb,lat_subscale)]
        if not any(np.isnan((lon1,lat1,lon2,lat2))):
            y_subdivided_lonlats = geod.npts(
                lon1,lat1,lon2,lat2,
                lat_subscale*2-1)                   # n intermediate points
            # ...and store to compute_grid:
            compute_grid_xg_out[
                i_cg(it.multi_index[0],ilb),
                j_cg(it.multi_index[1],jlb,lat_subscale)+1:
                j_cg(it.multi_index[1]+1,jlb,lat_subscale)] = \
                np.array(y_subdivided_lonlats)[:,0]
            compute_grid_yg_out[
                i_cg(it.multi_index[0],ilb),
                j_cg(it.multi_index[1],jlb,lat_subscale)+1:
                j_cg(it.multi_index[1]+1,jlb,lat_subscale)] = \
                np.array(y_subdivided_lonlats)[:,1]
        it.iternext()

    if verbose:
        print('compute_grid after y_direction subdivision fill-in:')
        print('compute_grid_xg:')
        print(compute_grid_xg_out)
        print('compute_grid_yg:')
        print(compute_grid_yg_out)

    return compute_grid_xg_out, compute_grid_yg_out

