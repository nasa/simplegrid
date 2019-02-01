
import numpy as np
import pyproj
from . import config
from . import util
from . import mitgridfilefields

def areas( compute_grid_xg, compute_grid_yg,
    compute_grid_edges_x, compute_grid_edges_y, geod, verbose=False):
    """
    Compute sub-areas for the compute grid.

    Args:
        compute_grid_xg, compute_grid_yg (numpy arrays): Arrays which have been
            dimensioned according to lon_subscale*2 and lat_subscale*2 (possibly
            with "boundary ring" of elements, depending on the application), and
            with (x,y) values filled in according to cgfill().
        compute_grid_edges_x, compute_grid_edges_y (numpy arrays): Arrays of
            grid edge lengths (as computed by edge()), used to determine whether
            a planar or a spherical surface approximation should be used for
            areas.
        geod (pyproj.Geod): current geoid instance (only the semi-major, or
            equatorial axis radius, geod.a is used here).
        verbose (logical): verbose output

    Returns:
        compute_areas: numpy array of sub-areas for the compute grid.

    """

    x_edge_mean = compute_grid_edges_x[compute_grid_edges_x.nonzero()].mean()
    y_edge_mean = compute_grid_edges_y[compute_grid_edges_y.nonzero()].mean()

    if  x_edge_mean<config.PLANAR_SPHER_TRANSITION and \
        y_edge_mean<config.PLANAR_SPHER_TRANSITION:
        # areas are nearly planar; use flat faceted surface approximation:
        compute_areas = util.pquad_uarea(
            util.lonlat2cart(compute_grid_xg,compute_grid_yg)) \
            * np.power(geod.a,2)
    else:
        # areas are surface triangles; use spherical excess formulation:
        compute_areas = util.squad_uarea(
            util.lonlat2cart(compute_grid_xg,compute_grid_yg)) \
            * np.power(geod.a,2)

    if verbose:
        print('compute_areas nonzero subset:')
        print(compute_areas)

    return compute_areas


def edges( compute_grid_xg, compute_grid_yg,
    ilb, iub, jlb, jub, geod, verbose=False):
    """
    Use great circle arcs to compute distances between all compute_grid points

    Args:
        compute_grid_xg, compute_grid_yg (numpy arrays): Arrays which have been
            dimensioned according to lon_subscale*2 and lat_subscale*2 (possibly
            with "boundary ring" of elements, depending on the application), and
            with (x,y) values filled in according to cgfill().
        ilb, iub, jlb, jub (ints): row (i) and column (j) lower and upper index
            bounds that define compute_grid_xg and compute_grid_yg range for
            edge calculations (avoids referenct to zero boundary values during
            length calculations).
        geod (pyproj.Geod): current geoid instance (only the semi-major, or
            equatorial axis radius, geod.a is used here).
        verbose (logical): verbose output

    Returns:
        (compute_edges_x, compute_edges_y) tuple of numpy arrays of grid edge lengths.
        compute_edges_x(i,j) contains the edge length from compute_grid_x /
        compute_grid_y (i,j) to (i+1,j) and compute_edges_y(i,j) contains the
        edge length from compute_grid_x / compute_grid_y (i,j) to (i,j+1)

    """

    # initialization can be based on either *_xg or *_yg since they're the same
    # shape:
    (cg_rows,cg_cols) = compute_grid_xg.shape
    compute_edges_x = np.zeros((cg_rows-1,cg_cols  ))
    compute_edges_y = np.zeros((cg_rows  ,cg_cols-1))

    # since matrix slice operations are from start:(stop-1):stride, define a
    # lambda function to make it clear when the stop point should be inclusive,
    # i.e., start:incl(stop):stride
    incl = lambda idx : idx+1

    # x-direction edge lengths:
    for i in range(ilb,iub):
        for j in range(jlb,incl(jub)):
            lon1 = compute_grid_xg[i  ,j]
            lat1 = compute_grid_yg[i  ,j]
            lon2 = compute_grid_xg[i+1,j]
            lat2 = compute_grid_yg[i+1,j]
            _,_,compute_edges_x[i,j] = geod.inv( lon1,lat1,lon2,lat2)

    # y-direction edge lengths:
    for i in range(ilb,incl(iub)):
        for j in range(jlb,jub):
            lon1 = compute_grid_xg[i,j  ]
            lat1 = compute_grid_yg[i,j  ]
            lon2 = compute_grid_xg[i,j+1]
            lat2 = compute_grid_yg[i,j+1]
            _,_,compute_edges_y[i,j] = geod.inv( lon1,lat1,lon2,lat2)

    if verbose:
        print('compute_edges_x:')
        print( compute_edges_x)
        print('compute_edges_y:')
        print( compute_edges_y)

    return (compute_edges_x,compute_edges_y)


def fill( compute_grid_xg, compute_grid_yg, ilb, iub, jlb, jub,
    lon_subscale, lat_subscale, geod, verbose=False):
    """
    Use great circle subdivisions to fill in compute grid intermediate points
    according to x/y subdivision levels and index ranges.

    Args:
        compute_grid_xg, compute_grid_yg (numpy arrays): Arrays which have been
            dimensioned according to lon_subscale*2 and lat_subscale*2 (possibly
            with "boundary ring" of elements, depending on the application), and
            with existing XG, YG, grid values already mapped to their
            corresponding positions.
        ilb, iub (int): compute_grid indices that define range of x-direction
            fill-in.
        jlb, jub (int): compute_grid indices that define range of y-direction
            fill-in.
        lon_subscale (int): user-specified x-direction cell subdivision level.
        lat_subscale (int): user-specified y-direction cell subdivision level.
        geod (pyproj.Geod): current geoid instance.
        verbose (logical): verbose output

    Returns:
        (compute_grid_xg_out, compute_grid_yg_out) tuple of input arrays with
        intermediate points filled in.

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
        x_edge_subdivided_lonlats = geod.npts(
            lon1,lat1,lon2,lat2,
            lon_subscale*2-1)   # n intermediate points
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
        y_subdivided_lonlats = geod.npts(
            lon1,lat1,lon2,lat2,
            lat_subscale*2-1)   # n intermediate points
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


def tomitgrid(compute_grid_xg,compute_grid_yg,
    iLB,ilb,iub,iUB,jLB,jlb,jub,jUB,
    geod,verbose=False):
    """Generates an mit grid from a compute grid.

    Args:
        compute_grid_xg, compute_grid_yg (numpy arrays): Arrays which have been
            dimensioned according to lon_subscale*2 and lat_subscale*2 (possibly
            with "boundary ring" of elements, depending on the application), and
            with (x,y) values filled in according to cgfill().
        iLB, iUB (int): indices that define absolute lower, upper range of
            x-direction values for compute_grid_xg and compute_grid_yg. Used for
            compute_grid supporting matrix dimensioning purposes.
        jLB, jUB (int): indices that define absolute lower, upper range of
            y-direction values for compute_grid_xg and compute_grid_yg. Used for
            compute_grid supporting matrix dimensioning purposes.
        ilb, iub (int): compute_grid indices that define inner range (i.e.,
            non-boundary grids) of x-direction values for compute_grid_xg and
            compute_grid_yg.
        jlb, jub (int): compute_grid indices that define inner range (i.e.,
            non-boundary grids)of y-direction values for compute_grid_xg and
            compute_grid_yg.
        geod (pyproj.Geod): current geoid instance.
        verbose (logical): verbose output.

    Returns:
        mitgrid (dict): name/value (numpy 2-d array) pairs
            corresponding to matrix name and ordering convention listed in
            mitgridfilefields module.

    """

    outgrid = {key:None for key in mitgridfilefields.names}

    # deduce overall lon/lat subdivisions (equal to tracer cell counts) from
    # input compute_grid ranges:
    lon_subdiv = (iub-ilb)//2   # integer division in python 3x
    lat_subdiv = (jub-jlb)//2   # "                           "

    # since matrix slice operations are from start:(stop-1):stride, define a
    # lambda function to make it clear when the stop point should be inclusive,
    # i.e., start:incl(stop):stride
    incl = lambda idx: idx+1

    # compute x- and y-direction compute_grid edge lengths:
    (compute_grid_edges_x,compute_grid_edges_y) = edges(
        compute_grid_xg,compute_grid_yg,ilb,iub,jlb,jub,geod,verbose)

    # compute subgrid areas:
    compute_grid_areas = np.zeros((iUB-iLB,jUB-jLB))
    compute_grid_areas[ilb:iub,jlb:jub] = areas(
        compute_grid_xg[ilb:incl(iub),jlb:incl(jub)],
        compute_grid_yg[ilb:incl(iub),jlb:incl(jub)],
        compute_grid_edges_x,compute_grid_edges_y,
        geod,verbose)

    # grid (tracer) cell location data:
    #   XC, YC - tracer cell center longitudes and latitudes
    #   XG, YG - tracer cell corner longitudes and latitudes

    # XC, YC directly from compute grid partitions:

    # compute grid partitioning:
    cg_first_i  = ilb+1
    cg_last_i   = incl(iub-1)
    cg_stride_i = 2
    cg_first_j  = jlb+1
    cg_last_j   = incl(jub-1)
    cg_stride_j = 2
    outgrid['XC'] = compute_grid_xg[
        cg_first_i:cg_last_i:cg_stride_i,
        cg_first_j:cg_last_j:cg_stride_j]
    outgrid['YC'] = compute_grid_yg[
        cg_first_i:cg_last_i:cg_stride_i,
        cg_first_j:cg_last_j:cg_stride_j]
    if verbose:
        print("outgrid['XC']:")
        print(outgrid['XC'])
        print("outgrid['YC']:")
        print(outgrid['YC'])

    # XG, YG directly from compute grid partitions:

    # compute grid partitioning:
    cg_first_i  = ilb
    cg_last_i   = incl(iub)
    cg_stride_i = 2
    cg_first_j  = jlb
    cg_last_j   = incl(jub)
    cg_stride_j = 2
    outgrid['XG'] = compute_grid_xg[
        cg_first_i:cg_last_i:cg_stride_i,
        cg_first_j:cg_last_j:cg_stride_j]
    outgrid['YG'] = compute_grid_yg[
        cg_first_i:cg_last_i:cg_stride_i,
        cg_first_j:cg_last_j:cg_stride_j]
    if verbose:
        print("outgrid['XG']:")
        print(outgrid['XG'])
        print("outgrid['YG']:")
        print(outgrid['YG'])

    # tracer cell-related quantities, RAC, DXG, DYG:
    #   DXG - (tracer) cell face separation in X along southern cell wall
    #   DYG - (tracer) cell face separation in Y along western cell wall
    #   RAC - tracer cell area presented in the vertical direction

    # DXG tracer cell southern edge from edge summations:

    outgrid['DXG'] = np.zeros((lon_subdiv,lat_subdiv+1))
    for i in range(0,lon_subdiv):
        for j in range(0,incl(lat_subdiv)):
            if not any(np.PZERO==np.array((
                compute_grid_edges_x[ ilb + 2*i    , jlb + 2*j],
                compute_grid_edges_x[ ilb + 2*i + 1, jlb + 2*j]))):
                outgrid['DXG'][i,j] = \
                    compute_grid_edges_x[ ilb + 2*i    , jlb + 2*j] + \
                    compute_grid_edges_x[ ilb + 2*i + 1, jlb + 2*j]
    if verbose:
        print("outgrid['DXG']:")
        print(outgrid['DXG'])

    # DYG tracer cell western edge from edge summations:

    outgrid['DYG'] = np.zeros((lon_subdiv+1,lat_subdiv))
    for i in range(0,incl(lon_subdiv)):
        for j in range(0,lat_subdiv):
            if not any(np.PZERO==np.array((
                compute_grid_edges_y[ ilb + 2*i, jlb + 2*j    ],
                compute_grid_edges_y[ ilb + 2*i, jlb + 2*j + 1]))):
                outgrid['DYG'][i,j] = \
                    compute_grid_edges_y[ ilb + 2*i, jlb + 2*j    ] + \
                    compute_grid_edges_y[ ilb + 2*i, jlb + 2*j + 1]
    if verbose:
        print("outgrid['DYG']:")
        print(outgrid['DYG'])

    # RAC from subcell area sums:

    outgrid['RAC'] = np.zeros((lon_subdiv,lat_subdiv))
    for i in range(0,lon_subdiv):
        for j in range(0,lat_subdiv):
            if not any(np.PZERO==np.array((
                compute_grid_areas[ilb + 2*i    , jlb + 2*j    ],
                compute_grid_areas[ilb + 2*i + 1, jlb + 2*j    ],
                compute_grid_areas[ilb + 2*i + 1, jlb + 2*j + 1],
                compute_grid_areas[ilb + 2*i    , jlb + 2*j + 1]))):
                outgrid['RAC'][i,j] = \
                    compute_grid_areas[ilb + 2*i    , jlb + 2*j    ] + \
                    compute_grid_areas[ilb + 2*i + 1, jlb + 2*j    ] + \
                    compute_grid_areas[ilb + 2*i + 1, jlb + 2*j + 1] + \
                    compute_grid_areas[ilb + 2*i    , jlb + 2*j + 1]
    if verbose:
        print("outgrid['RAC']:")
        print(outgrid['RAC'])

    # vorticity cell-related quantities, DXC, DYC, RAZ
    #   DXC - vorticity cell southern edge length in X-direction
    #   DYC - vorticity cell western edge length in Y-direction
    #   RAZ - vorticity cell area presented in the vertical direction

    # DXC vorticity cell edge lengths from x-direction edge summations:

    outgrid['DXC'] = np.zeros((lon_subdiv+1,lat_subdiv))
    for i in range(0,incl(lon_subdiv)):
        for j in range(0,lat_subdiv):
            if not any(np.PZERO==np.array((
                compute_grid_edges_x[ 2*i    , jlb + 2*j + 1],
                compute_grid_edges_x[ 2*i + 1, jlb + 2*j + 1]))):
                outgrid['DXC'][i,j] = \
                    compute_grid_edges_x[ 2*i    , jlb + 2*j + 1] + \
                    compute_grid_edges_x[ 2*i + 1, jlb + 2*j + 1]
    if verbose:
        print("outgrid['DXC']:")
        print(outgrid['DXC'])

    # DYC vorticity cell edge lengths from y-direction edge summations:

    outgrid['DYC'] = np.zeros((lon_subdiv,lat_subdiv+1))
    for i in range(0,lon_subdiv):
        for j in range(0,incl(lat_subdiv)):
            if not any(np.PZERO==np.array((
                compute_grid_edges_y[ ilb + 2*i + 1, 2*j    ],
                compute_grid_edges_y[ ilb + 2*i + 1, 2*j + 1]))):
                outgrid['DYC'][i,j] = \
                    compute_grid_edges_y[ ilb + 2*i + 1, 2*j    ] + \
                    compute_grid_edges_y[ ilb + 2*i + 1, 2*j + 1]
    if verbose:
        print("outgrid['DYC']:")
        print(outgrid['DYC'])

    # RAZ vorticity cell areas computed from subcell area sums:

    outgrid['RAZ'] = np.zeros((lon_subdiv+1,lat_subdiv+1))
    for i in range(0,incl(lon_subdiv)):
        for j in range(0,incl(lat_subdiv)):
            if not any(np.PZERO==np.array((
                compute_grid_areas[2*i    , 2*j    ],
                compute_grid_areas[2*i + 1, 2*j    ],
                compute_grid_areas[2*i + 1, 2*j + 1],
                compute_grid_areas[2*i    , 2*j + 1]))):
                outgrid['RAZ'][i,j] = \
                    compute_grid_areas[2*i    , 2*j    ] + \
                    compute_grid_areas[2*i + 1, 2*j    ] + \
                    compute_grid_areas[2*i + 1, 2*j + 1] + \
                    compute_grid_areas[2*i    , 2*j + 1]
    if verbose:
        print("outgrid['RAZ']:")
        print(outgrid['RAZ'])

    # "U" cell-related quantities, DXV, DYF, RAW
    #   DXV - U cell edge length in X-direction between v-points
    #   DYF - U cell edge length in Y-direction between tracer cell faces
    #   RAW - U cell area presented in the vertical direction

    # DXV U cell edge lengths from x-direction edge summations:

    outgrid['DXV'] = np.zeros((lon_subdiv+1,lat_subdiv+1))
    for i in range(0,incl(lon_subdiv)):
        for j in range(0,incl(lat_subdiv)):
            if not any(np.PZERO==np.array((
                compute_grid_edges_x[ 2*i    , jlb + 2*j],
                compute_grid_edges_x[ 2*i + 1, jlb + 2*j]))):
                outgrid['DXV'][i,j] = \
                    compute_grid_edges_x[ 2*i    , jlb + 2*j] + \
                    compute_grid_edges_x[ 2*i + 1, jlb + 2*j]
    if verbose:
        print("outgrid['DXV']:")
        print(outgrid['DXV'])

    # DYF U cell edge lengths from y-direction edge summations:

    outgrid['DYF'] = np.zeros((lon_subdiv,lat_subdiv))
    for i in range(0,lon_subdiv):
        for j in range(0,lat_subdiv):
            if not any(np.PZERO==np.array((
                compute_grid_edges_y[ ilb + 2*i + 1, jlb + 2*j    ],
                compute_grid_edges_y[ ilb + 2*i + 1, jlb + 2*j + 1]))):
                outgrid['DYF'][i,j] = \
                    compute_grid_edges_y[ ilb + 2*i + 1, jlb + 2*j    ] + \
                    compute_grid_edges_y[ ilb + 2*i + 1, jlb + 2*j + 1]
    if verbose:
        print("outgrid['DYF']:")
        print(outgrid['DYF'])

    # RAW U cell areas from subcell area sums:

    outgrid['RAW'] = np.zeros((lon_subdiv+1,lat_subdiv))
    for i in range(0,incl(lon_subdiv)):
        for j in range(0,lat_subdiv):
            if not any(np.PZERO==np.array((
                compute_grid_areas[2*i    , jlb + 2*j    ],
                compute_grid_areas[2*i + 1, jlb + 2*j    ],
                compute_grid_areas[2*i + 1, jlb + 2*j + 1],
                compute_grid_areas[2*i    , jlb + 2*j + 1]))):
                outgrid['RAW'][i,j] = \
                    compute_grid_areas[2*i    , jlb + 2*j    ] + \
                    compute_grid_areas[2*i + 1, jlb + 2*j    ] + \
                    compute_grid_areas[2*i + 1, jlb + 2*j + 1] + \
                    compute_grid_areas[2*i    , jlb + 2*j + 1]
    if verbose:
        print("outgrid['RAW']:")
        print(outgrid['RAW'])

    # "V" cell-related quantities, DXF, DYU, RAS
    #   DXF - V cell northern edge length in X-direction between u-points
    #   DYU - V cell western edge length in Y-direction
    #   RAS - V cell area presented in the vertical direction

    # DXF V cell edge lengths from x-direction edge summations:

    outgrid['DXF'] = np.zeros((lon_subdiv,lat_subdiv))
    for i in range(0,lon_subdiv):
        for j in range(0,lat_subdiv):
            if not any(np.PZERO==np.array((
                compute_grid_edges_x[ ilb + 2*i    , jlb + 2*j + 1],
                compute_grid_edges_x[ ilb + 2*i + 1, jlb + 2*j + 1]))):
                outgrid['DXF'][i,j] = \
                    compute_grid_edges_x[ ilb + 2*i    , jlb + 2*j + 1] + \
                    compute_grid_edges_x[ ilb + 2*i + 1, jlb + 2*j + 1]
    if verbose:
        print("outgrid['DXF']:")
        print(outgrid['DXF'])

    # DYU V cell western edge lengths from y-direction edge summations:

    outgrid['DYU'] = np.zeros((lon_subdiv+1,lat_subdiv+1))
    for i in range(0,incl(lon_subdiv)):
        for j in range(0,incl(lat_subdiv)):
            if not any(np.PZERO==np.array((
                compute_grid_edges_y[ ilb + 2*i, 2*j    ],
                compute_grid_edges_y[ ilb + 2*i, 2*j + 1]))):
                outgrid['DYU'][i,j] = \
                    compute_grid_edges_y[ ilb + 2*i, 2*j    ] + \
                    compute_grid_edges_y[ ilb + 2*i, 2*j + 1]
    if verbose:
        print("outgrid['DYU']:")
        print(outgrid['DYU'])

    # RAS V cell areas from subcell area sums:

    outgrid['RAS'] = np.zeros((lon_subdiv,lat_subdiv+1))
    for i in range(0,lon_subdiv):
        for j in range(0,incl(lat_subdiv)):
            if not any(np.PZERO==np.array((
                compute_grid_areas[ilb + 2*i    , 2*j    ],
                compute_grid_areas[ilb + 2*i + 1, 2*j    ],
                compute_grid_areas[ilb + 2*i + 1, 2*j + 1],
                compute_grid_areas[ilb + 2*i    , 2*j + 1]))):
                outgrid['RAS'][i,j] = \
                    compute_grid_areas[ilb + 2*i    , 2*j    ] + \
                    compute_grid_areas[ilb + 2*i + 1, 2*j    ] + \
                    compute_grid_areas[ilb + 2*i + 1, 2*j + 1] + \
                    compute_grid_areas[ilb + 2*i    , 2*j + 1]
    if verbose:
        print("outgrid['RAS']:")
        print(outgrid['RAS'])

    return outgrid

