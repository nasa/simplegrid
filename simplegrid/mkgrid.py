#!/usr/bin/env python

import argparse
import math
import numpy as np
import pyproj
from . import computegrid
from . import gridio
from . import mitgridfilefields
from . import util


def create_parser():
    """Set up the list of arguments to be provided to mkgrid.
    """
    parser = argparse.ArgumentParser(
        description="""
            Create a simple rectangular grid over a specified domain.""",
        epilog="""
            Corner point information (lon1/lat1, lon2/lat2) can be either
            literal northwest/southeast corners, or diagonal corners provided
            simply for orientation purposes. Either way, the output grid will be
            aligned such that 'northwest' is at min i, max j, and 'southeast' is
            at max i, min j.""")
    parser.add_argument('--lon1', type=float, help="""
        longitude of northwest corner point""")
    parser.add_argument('--lat1', type=float, help="""
        latitude of northwest corner point""")
    parser.add_argument('--lon2', type=float, help="""
        longitude of southeast corner point""")
    parser.add_argument('--lat2', type=float, help="""
        latitude of southeast corner point""")
    parser.add_argument('--lon_subscale', type=int, help="""
        number of longitudinal subdivisions in the resulting grid (number of
        x-direction tracer cells)""")
    parser.add_argument('--lat_subscale', type=int, help="""
        number of latitudinal subdivisions in the resulting grid (number of
        y-direction tracer cells)""")
    parser.add_argument('--outfile', help="""
        file to which grid matrices will be written (mitgridfile format)""")
    parser.add_argument('-v','--verbose',action='store_true',help="""
        verbose output""")
    return parser


def mkgrid(lon1,lat1,lon2,lat2,lon_subscale,lat_subscale,verbose=False):
    """Creates a rectangular lon/lat grid using simple great circle-based
    subdivisions. A normal spherical geoid is currently assumed.

    Args:
        lon1 (float): longitude of northwest corner point.
        lat1 (float): latitude of northwest corner point.
        lon2 (float): longitude of southeast corner point.
        lat2 (float): latitude of southeast corner point.
        lon_subscale (int): desired number of longitudinal subdivisions in the
            resulting grid (number of x-direction tracer cells).
        lat_subscale (int): desired number of latitudinal subdivisions in the
            resulting grid (number of y-direction tracer cells).
        verbose (bool): True for diagnostic output, False otherwise.

    Returns:
        (newgrid,newgrid_ni,newgrid_nj): For consistency with regrid(), tuple
            consisting of dictionary of grid matrices and corresponding tracer
            cell counts ni (=lon_subscale input) and nj (=lat_subscale input).

    Note:
        Corner point information (lon1/lat1, lon2/lat2) can be either literal
        northwest/southeast corners, or diagonal corners provided simply for
        orientation purposes. Either way, the output grid will be aligned such
        that 'northwest' is at min i, max j, and 'southeast' is at max i, min j.

    """

    # for now, assume spherical geoid (perhaps user-specified later):
    geod = pyproj.Geod(ellps='sphere')

    if verbose:
        print('creating simple {0} x {1} rectangular grid between the extents'.format(
            lon_subscale,lat_subscale))
        print('({0:9.4f},{1:9.4f}) and ({2:9.4f},{3:9.4f})'.format(
            lon1,lat1,lon2,lat2))


    # since matrix slice operations are from start:(stop-1):stride, define a
    # lambda function to make it clear when the stop point should be inclusive,
    # i.e., start:incl(stop):stride

    incl = lambda idx: idx+1

    #
    # Step 1:
    #
    # Based on user-selected discretization level, create a "compute grid" that
    # spans the expected range, plus a boundary "ring" one compute cell wide.
    # (Note that this differs slightly from the logic in regrid.py where the
    # boundary ring is the width of one tracer cell at the initial resolution,
    # which translates into 2*lat_subscale x 2*lon_subscale.)
    # Initializing grid values to NaN will allow us to compute values using
    # consistent indexing, while naturally producing undefined values at the
    # boundaries.
    #

    # mimicing the logic in regrid(), index bounds based on lon/lat_subscale
    # "plus one" grid extents:

    # "plus one" extents:
    iLB, iUB = 0, 2*(lon_subscale+1)    # (i.e., iUB=iub+1)
    jLB, jUB = 0, 2*(lat_subscale+1)    # (i.e., jUB=jub+1)
    # "user-provided" extents (relative to "plus one" extents):
    ilb, iub = 1, 1+2*lon_subscale
    jlb, jub = 1, 1+2*lat_subscale

    # "compute grid" dimensions, allocation:
    num_compute_grid_rows = 2*(lon_subscale+1) + 1
    num_compute_grid_cols = 2*(lat_subscale+1) + 1
    compute_grid_xg = np.empty((num_compute_grid_rows,num_compute_grid_cols))
    compute_grid_xg[:,:] = np.nan
    compute_grid_yg = np.empty((num_compute_grid_rows,num_compute_grid_cols))
    compute_grid_yg[:,:] = np.nan

    # "compute areas" dimensions, allocation:
    num_compute_areas_rows = 2*(lon_subscale+1)
    num_compute_areas_cols = 2*(lat_subscale+1)
    compute_areas = np.empty((num_compute_areas_rows,num_compute_areas_cols))
    compute_areas[:,:] = np.nan

    #
    # Step 2: Populate compute grid points that map to user-selected range:
    #

    # 2a: user-specified northwest, southeast corners:

    min_lon, max_lon = lon1, lon2   # nw, se
    min_lat, max_lat = lat2, lat1   # se, nw
    compute_grid_xg[ilb,jlb] = compute_grid_xg[ilb,jub] = min_lon
    compute_grid_xg[iub,jlb] = compute_grid_xg[iub,jub] = max_lon
    compute_grid_yg[ilb,jlb] = compute_grid_yg[iub,jlb] = min_lat
    compute_grid_yg[ilb,jub] = compute_grid_yg[iub,jub] = max_lat

    # 2b: fill in intermediate points per user-specified subdivision level:

    (compute_grid_xg,compute_grid_yg) = computegrid.cgfill(
        compute_grid_xg,compute_grid_yg,ilb,iub,jlb,jub,
        lon_subscale,lat_subscale,geod,verbose)

    #
    # Step 3: Generate areas for sub-quads at the compute_grid array resolution
    #

    # areas, from cartesian coordinates on the unit sphere, scaled according to
    # mean spherical ellipsoid radius:

    compute_areas[1:incl(2*lon_subscale),1:incl(2*lat_subscale)] = util.squad_uarea(
        util.lonlat2cart(
            compute_grid_xg[ilb:incl(iub),jlb:incl(jub)],
            compute_grid_yg[ilb:incl(iub),jlb:incl(jub)])) \
            * np.power(geod.a,2)

    if verbose:
        print('compute_areas:')
        print(compute_areas)

    #
    # Step 4: Create and fill in output quantities based on compute grid data:
    #

    outgrid = {key:None for key in mitgridfilefields.names}

    # grid location data:
    #   XC - longitude east of center of grid (tracer) cell
    #   YC - latitude north of center of grid (tracer) cell
    #   XG - latitude east of southwest corner of grid (tracer) cell
    #   YC - latitude north of southwest corner of grid (tracer) cell

    # XC, YC directly from compute grid partitions:

    # compute grid partitioning:
    cg_first_i  = ilb+1
    cg_last_i   = incl(cg_first_i + 2*(lon_subscale-1))
    cg_stride_i = 2
    cg_first_j  = jlb+1
    cg_last_j   = incl(cg_first_j + 2*(lat_subscale-1))
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
    cg_last_i   = incl(cg_first_i + 2*lon_subscale)
    cg_stride_i = 2
    cg_first_j  = jlb
    cg_last_j   = incl(cg_first_j + 2*lat_subscale)
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

    # DXG tracer cell southern edge computed from compute_grid grid point
    # distances:

    outgrid['DXG'] = np.empty((lon_subscale,lat_subscale+1))
    outgrid['DXG'][:,:] = np.nan    # error check, mostly

    # compute grid partitioning:
    cg_first_i  = ilb
    cg_last_i   = incl(cg_first_i + 2*(lon_subscale-1))
    cg_stride_i = 2
    cg_first_j  = jlb
    cg_last_j   = incl(cg_first_j + 2*lat_subscale)
    cg_stride_j = 2

    # transformation from partitioned, strided nditer space (i_n,j_n) to
    # underlying compute_grid space (i_cg,j_cg):
    i_cg = lambda i_n,ilb : ilb + 2*i_n
    j_cg = lambda j_n,jlb : jlb + 2*j_n

    it = np.nditer(
        [compute_grid_xg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_i],
         compute_grid_yg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_j]],
        flags=['multi_index'])
    while not it.finished:
        _,_,outgrid['DXG'][it.multi_index[0],it.multi_index[1]] = geod.inv(
            it[0],                                  # lon1
            it[1],                                  # lat1
            compute_grid_xg[
                i_cg(it.multi_index[0]+1,ilb),
                j_cg(it.multi_index[1],jlb)],       # lon2
            compute_grid_yg[
                i_cg(it.multi_index[0]+1,ilb),
                j_cg(it.multi_index[1],jlb)])       # lat2
        it.iternext()
    if verbose:
        print("outgrid['DXG']:")
        print(outgrid['DXG'])

    # DYG tracer cell western edge computed from compute_grid grid point
    # distances:

    outgrid['DYG'] = np.empty((lon_subscale+1,lat_subscale))
    outgrid['DYG'][:,:] = np.nan    # error check, mostly

    # compute grid partitioning:
    cg_first_i  = ilb
    cg_last_i   = incl(cg_first_i + 2*(lon_subscale))
    cg_stride_i = 2
    cg_first_j  = jlb
    cg_last_j   = incl(cg_first_j + 2*(lat_subscale-1))
    cg_stride_j = 2

    # (i_n,j_n) to (i_cg,j_cg) transformations are the same as for DXG above.

    it = np.nditer(
        [compute_grid_xg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_i],
         compute_grid_yg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_j]],
        flags=['multi_index'])
    while not it.finished:
        _,_,outgrid['DYG'][it.multi_index[0],it.multi_index[1]] = geod.inv(
            it[0],                                  # lon1
            it[1],                                  # lat1
            compute_grid_xg[
                i_cg(it.multi_index[0],ilb),
                j_cg(it.multi_index[1]+1,jlb)],     # lon2
            compute_grid_yg[
                i_cg(it.multi_index[0],ilb),
                j_cg(it.multi_index[1]+1,jlb)])     # lat2
        it.iternext()
    if verbose:
        print("outgrid['DYG']:")
        print(outgrid['DYG'])

    # RAC computed from subcell area sums in compute_areas:

    outgrid['RAC'] = np.empty((lon_subscale,lat_subscale))
    outgrid['RAC'][:,:] = np.nan    # error check, mostly

    # compute_areas partitioning:
    ca_first_i  = 1
    ca_last_i   = incl(ca_first_i + 2*(lon_subscale-1))
    ca_stride_i = 2
    ca_first_j  = 1
    ca_last_j   = incl(ca_first_j + 2*(lat_subscale-1))
    ca_stride_j = 2

    # transformation from partitioned, strided compute_areas indices (i_n,j_n) to
    # underlying compute_areas indices:
    i_ca = lambda i_n : 2*i_n + 1
    j_ca = lambda j_n : 2*j_n + 1

    it = np.nditer(
        compute_areas[ca_first_i:ca_last_i:ca_stride_i,ca_first_j:ca_last_j:ca_stride_j],
        flags=['multi_index'])
    while not it.finished:
        outgrid['RAC'][it.multi_index[0],it.multi_index[1]] = \
            it[0] + \
            compute_areas[i_ca(it.multi_index[0])+1,j_ca(it.multi_index[1])  ] + \
            compute_areas[i_ca(it.multi_index[0])+1,j_ca(it.multi_index[1])+1] + \
            compute_areas[i_ca(it.multi_index[0])  ,j_ca(it.multi_index[1])+1]
        it.iternext()
    if verbose:
        print("outgrid['RAC']:")
        print(outgrid['RAC'])

    # vorticity cell-related quantities, DXC, DYC, RAZ
    #   DXC - vorticity cell southern edge length in X-direction
    #   DYC - vorticity cell western edge length in Y-direction
    #   RAZ - vorticity cell area presented in the vertical direction

    # DXC vorticity cell edge lengths computed from compute_grid points:

    outgrid['DXC'] = np.empty((lon_subscale+1,lat_subscale))
    outgrid['DXC'][:,:] = np.nan

    # compute grid partitioning:
    cg_first_i  = 0
    cg_last_i   = incl(cg_first_i + 2*lon_subscale)
    cg_stride_i = 2
    cg_first_j  = jlb + 1
    cg_last_j   = incl(cg_first_j + 2*(lat_subscale-1))
    cg_stride_j = 2

    # transformation from partitioned, strided nditer space (i_n,j_n) to
    # underlying compute_grid space (i_cg,j_cg):
    i_cg = lambda i_n : 2*i_n
    j_cg = lambda j_n : 2*j_n + 2

    it = np.nditer(
        [compute_grid_xg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_i],
         compute_grid_yg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_j]],
        flags=['multi_index'])
    while not it.finished:
        lon1 = it[0]
        lat1 = it[1]
        lon2 = compute_grid_xg[
            i_cg(it.multi_index[0]+1),
            j_cg(it.multi_index[1]  )]
        lat2 = compute_grid_yg[
            i_cg(it.multi_index[0]+1),
            j_cg(it.multi_index[1]  )]
        if not any(np.isnan((lon1,lat1,lon2,lat2))):
            _,_,outgrid['DXC'][it.multi_index[0],it.multi_index[1]] = geod.inv(
                lon1,lat1,lon2,lat2)
        it.iternext()
    if verbose:
        print("outgrid['DXC']:")
        print(outgrid['DXC'])

    # DYC vorticity cell edge lengths computed from compute_grid points:

    outgrid['DYC'] = np.empty((lon_subscale,lat_subscale+1))
    outgrid['DYC'][:,:] = np.nan

    # compute grid partitioning:
    cg_first_i  = ilb + 1
    cg_last_i   = incl(cg_first_i + 2*(lon_subscale-1))
    cg_stride_i = 2
    cg_first_j  = 0
    cg_last_j   = incl(cg_first_j + 2*lat_subscale)
    cg_stride_j = 2

    # transformation from partitioned, strided nditer space (i_n,j_n) to
    # underlying compute_grid space (i_cg,j_cg):
    i_cg = lambda i_n : 2*i_n + 2
    j_cg = lambda j_n : 2*j_n

    it = np.nditer(
        [compute_grid_xg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_i],
         compute_grid_yg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_j]],
        flags=['multi_index'])
    while not it.finished:
        lon1 = it[0]
        lat1 = it[1]
        lon2 = compute_grid_xg[
            i_cg(it.multi_index[0]  ),
            j_cg(it.multi_index[1]+1)]
        lat2 = compute_grid_yg[
            i_cg(it.multi_index[0]  ),
            j_cg(it.multi_index[1]+1)]
        if not any(np.isnan((lon1,lat1,lon2,lat2))):
            _,_,outgrid['DYC'][it.multi_index[0],it.multi_index[1]] = geod.inv(
                lon1,lat1,lon2,lat2)
        it.iternext()
    if verbose:
        print("outgrid['DYC']:")
        print(outgrid['DYC'])

    # RAZ vorticity cell areas computed from subcell areas in compute_areas:

    outgrid['RAZ'] = np.empty((lon_subscale+1,lat_subscale+1))
    outgrid['RAZ'][:,:] = np.nan

    # compute_areas partitioning:
    ca_first_i  = 0
    ca_last_i   = incl(2*lon_subscale)
    ca_stride_i = 2
    ca_first_j  = 0
    ca_last_j   = incl(2*lat_subscale)
    ca_stride_j = 2

    # transformation from partitioned, strided compute_areas indices (nditer) to
    # underlying compute_areas indices:
    i_ca = lambda i_n : 2*i_n
    j_ca = lambda j_n : 2*j_n

    it = np.nditer(
        compute_areas[ca_first_i:ca_last_i:ca_stride_i,ca_first_j:ca_last_j:ca_stride_j],
        flags=['multi_index'])
    while not it.finished:
        outgrid['RAZ'][it.multi_index[0],it.multi_index[1]] = \
            it[0] + \
            compute_areas[i_ca(it.multi_index[0])+1,j_ca(it.multi_index[1])  ] + \
            compute_areas[i_ca(it.multi_index[0])+1,j_ca(it.multi_index[1])+1] + \
            compute_areas[i_ca(it.multi_index[0])  ,j_ca(it.multi_index[1])+1]
        it.iternext()
    if verbose:
        print("outgrid['RAZ']:")
        print(outgrid['RAZ'])

    # "U" cell-related quantities, DXC, DYC, RAZ
    #   DXV - U cell edge length in X-direction between v-points
    #   DYF - U cell edge length in Y-direction between tracer cell faces
    #   RAW - U cell area presented in the vertical direction

    # DXV U cell edge lengths computed from compute_grid points: 

    outgrid['DXV'] = np.empty((lon_subscale+1,lat_subscale+1))
    outgrid['DXV'][:,:] = np.nan

    # compute grid partitioning:
    cg_first_i  = 0
    cg_last_i   = incl(lon_subscale*2)
    cg_stride_i = 2
    cg_first_j  = jlb
    cg_last_j   = incl(cg_first_j + 2*lat_subscale)
    cg_stride_j = 2

    # transformation from partitioned, strided nditer space (i_n,j_n) to
    # underlying compute_grid space (i_cg,j_cg):
    i_cg = lambda i_n : 2*i_n
    j_cg = lambda j_n : 2*j_n + 1

    it = np.nditer(
        [compute_grid_xg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_i],
         compute_grid_yg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_j]],
        flags=['multi_index'])
    while not it.finished:
        lon1 = it[0]
        lat1 = it[1]
        lon2 = compute_grid_xg[
            i_cg(it.multi_index[0]+1),
            j_cg(it.multi_index[1]  )]
        lat2 = compute_grid_yg[
            i_cg(it.multi_index[0]+1),
            j_cg(it.multi_index[1]  )]
        if not any(np.isnan((lon1,lat1,lon2,lat2))):
            _,_,outgrid['DXV'][it.multi_index[0],it.multi_index[1]] = geod.inv(
                lon1,lat1,lon2,lat2)
        it.iternext()
    if verbose:
        print("outgrid['DXV']:")
        print(outgrid['DXV'])

    # DYF U cell edge lengths computed from compute_grid points:

    outgrid['DYF'] = np.empty((lon_subscale,lat_subscale))
    outgrid['DYF'][:,:] = np.nan

    # compute grid partitioning:
    cg_first_i  = ilb + 1
    cg_last_i   = incl(cg_first_i + 2*(lon_subscale-1))
    cg_stride_i = 2
    cg_first_j  = jlb
    cg_last_j   = incl(cg_first_j + 2*(lat_subscale-1))
    cg_stride_j = 2

    # transformation from partitioned, strided nditer space (i_n,j_n) to
    # underlying compute_grid space (i_cg,j_cg):
    i_cg = lambda i_n : 2*i_n + 2
    j_cg = lambda j_n : 2*j_n + 1

    it = np.nditer(
        [compute_grid_xg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_i],
         compute_grid_yg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_j]],
        flags=['multi_index'])
    # note isnan checks not required for DYF
    while not it.finished:
        _,_,outgrid['DYF'][it.multi_index[0],it.multi_index[1]] = geod.inv(
            it[0],                              # lon1
            it[1],                              # lat1
            compute_grid_xg[
                i_cg(it.multi_index[0]  ),
                j_cg(it.multi_index[1]+1)],     # lon2
            compute_grid_yg[
                i_cg(it.multi_index[0]  ),
                j_cg(it.multi_index[1]+1)])     # lat2
        it.iternext()
    if verbose:
        print("outgrid['DYF']:")
        print(outgrid['DYF'])

    # RAW vertical face area of U cell computed from subcell areas in
    # compute_areas:

    outgrid['RAW'] = np.empty((lon_subscale+1,lat_subscale))
    outgrid['RAW'][:,:] = np.nan

    # compute areas partitioning:
    ca_first_i  = 0
    ca_last_i   = incl(2*lon_subscale)
    ca_stride_i = 2
    ca_first_j  = 1
    ca_last_j   = incl(ca_first_j + 2*(lat_subscale-1))
    ca_stride_j = 2

    # transformation from partitioned, strided compute_areas indices (nditer) to
    # underlying compute_areas indices:
    i_ca = lambda i_n : 2*i_n
    j_ca = lambda j_n : 2*j_n + 1

    it = np.nditer(
        compute_areas[ca_first_i:ca_last_i:ca_stride_i,ca_first_j:ca_last_j:ca_stride_j],
        flags=['multi_index'])
    while not it.finished:
        outgrid['RAW'][it.multi_index[0],it.multi_index[1]] = \
            it[0] + \
            compute_areas[i_ca(it.multi_index[0])+1,j_ca(it.multi_index[1])  ] + \
            compute_areas[i_ca(it.multi_index[0])+1,j_ca(it.multi_index[1])+1] + \
            compute_areas[i_ca(it.multi_index[0])  ,j_ca(it.multi_index[1])+1]
        it.iternext()
    if verbose:
        print("outgrid['RAW']:")
        print(outgrid['RAW'])

    # "V" cell-related quantities, DXF, DYU, RAS
    #   DXF - V cell northern edge length in X-direction between u-points
    #   DYU - V cell western edge length in Y-direction
    #   RAS - V cell area presented in the vertical direction

    # DXF V cell edge lengths computed using compute_grid point distances:

    outgrid['DXF'] = np.empty((lon_subscale,lat_subscale))
    outgrid['DXF'][:,:] = np.nan

    # compute grid partitioning:
    cg_first_i  = ilb
    cg_last_i   = incl(cg_first_i + 2*(lon_subscale-1))
    cg_stride_i = 2
    cg_first_j  = jlb + 1
    cg_last_j   = incl(cg_first_j + 2*(lat_subscale-1))
    cg_stride_j = 2

    # transformation from partitioned, strided nditer space (i_n,j_n) to
    # underlying compute_grid space (i_cg,j_cg):
    i_cg = lambda i_n : 2*i_n + 1
    j_cg = lambda j_n : 2*j_n + 2

    it = np.nditer(
        [compute_grid_xg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_i],
         compute_grid_yg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_j]],
        flags=['multi_index'])
    # note isnan checks not required for DXF
    while not it.finished:
        _,_,outgrid['DXF'][it.multi_index[0],it.multi_index[1]] = geod.inv(
            it[0],                              # lon1
            it[1],                              # lat1
            compute_grid_xg[
                i_cg(it.multi_index[0]+1),
                j_cg(it.multi_index[1]  )],     # lon2
            compute_grid_yg[
                i_cg(it.multi_index[0]+1),
                j_cg(it.multi_index[1]  )])     # lat2
        it.iternext()
    if verbose:
        print("outgrid['DXF']:")
        print(outgrid['DXF'])

    # DYU V cell western edge length computed from compute_grid point distances:

    outgrid['DYU'] = np.empty((lon_subscale+1,lat_subscale+1))
    outgrid['DYU'][:,:] = np.nan

    # compute grid partitioning:
    cg_first_i  = ilb
    cg_last_i   = incl(cg_first_i + 2*lon_subscale)
    cg_stride_i = 2
    cg_first_j  = 0
    cg_last_j   = incl(2*lat_subscale)
    cg_stride_j = 2

    # transformation from partitioned, strided nditer space (i_n,j_n) to
    # underlying compute_grid space (i_cg,j_cg):
    i_cg = lambda i_n : 2*i_n + 1
    j_cg = lambda j_n : 2*j_n

    it = np.nditer(
        [compute_grid_xg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_i],
         compute_grid_yg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_j]],
        flags=['multi_index'])
    while not it.finished:
        lon1 = it[0]
        lat1 = it[1]
        lon2 = compute_grid_xg[
            i_cg(it.multi_index[0]  ),
            j_cg(it.multi_index[1]+1)]
        lat2 = compute_grid_yg[
            i_cg(it.multi_index[0]  ),
            j_cg(it.multi_index[1]+1)]
        if not any(np.isnan((lon1,lat1,lon2,lat2))):
            _,_,outgrid['DYU'][it.multi_index[0],it.multi_index[1]] = geod.inv(
                lon1,lat1,lon2,lat2)
        it.iternext()
    if verbose:
        print("outgrid['DYU']:")
        print(outgrid['DYU'])

    # RAS vertical face area of V cell computed from subcell areas in
    # compute_areas:

    outgrid['RAS'] = np.empty((lon_subscale,lat_subscale+1))
    outgrid['RAS'][:,:] = np.nan

    # compute areas partitioning:
    ca_first_i  = ilb
    ca_last_i   = incl(ilb + 2*(lon_subscale-1))
    ca_stride_i = 2
    ca_first_j  = 0
    ca_last_j   = incl(2*lat_subscale)
    ca_stride_j = 2

    # transformation from partitioned, strided compute_areas indices (i_n,j_n) to
    # underlying compute_areas indices:
    i_ca = lambda i_n : 2*i_n + 1
    j_ca = lambda j_n : 2*j_n

    it = np.nditer(
        compute_areas[ca_first_i:ca_last_i:ca_stride_i,ca_first_j:ca_last_j:ca_stride_j],
        flags=['multi_index'])
    while not it.finished:
        outgrid['RAS'][it.multi_index[0],it.multi_index[1]] = \
            it[0] + \
            compute_areas[i_ca(it.multi_index[0])+1,j_ca(it.multi_index[1])  ] + \
            compute_areas[i_ca(it.multi_index[0])+1,j_ca(it.multi_index[1])+1] + \
            compute_areas[i_ca(it.multi_index[0])  ,j_ca(it.multi_index[1])+1]
        it.iternext()
    if verbose:
        print("outgrid['RAS']:")
        print(outgrid['RAS'])

    return outgrid, lon_subscale, lat_subscale


def main():
    """Command-line entry point."""

    parser = create_parser()
    args = parser.parse_args()
    (newgrid,newgrid_ni,newgrid_nj) = mkgrid(
        args.lon1,
        args.lat1,
        args.lon2,
        args.lat2,
        args.lon_subscale,
        args.lat_subscale,
        args.verbose)
    if args.verbose:
        print('writing {0:s} with ni={1:d}, nj={2:d}...'.
            format(args.outfile,newgrid_ni,newgrid_nj))
    gridio.write_mitgridfile(args.outfile,newgrid,newgrid_ni,newgrid_nj)
    if args.verbose:
        print('...done.')

if __name__ == '__main__':
    main()

