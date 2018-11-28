#!/usr/bin/env python

import argparse
import math
import numpy as np
import os
import pyproj
from . import computegrid
from . import gridio
from . import mitgridfilefields
from . import util


def create_parser():
    """Set up the list of arguments to be provided to regrid.
    """
    parser = argparse.ArgumentParser(
        description="""
            Regrid subdomain of an existing grid.""",
        epilog="""
            Corner point information (lon1/lat1, lon2/lat2) can be either
            literal northwest/southeast corners, or diagonal corners provided
            simply for orientation purposes. Either way, the output grid will be
            aligned such that 'northwest' is at min i, max j, and 'southeast' is
            at max i, min j.""")
    parser.add_argument('--mitgridfile', help="""
        mitgrid (path and) file name""")
    parser.add_argument('--xg_file', help="""
        XG (path and) file input alternative to mitgridfile (csv if .csv
        extension (one matrix row per line), double-precision column-ordered
        binary otherwise; --yg_file must also be provided)""")
    parser.add_argument('--yg_file', help="""
        YG (path and) file input alternative to mitgridfile (see --xg_file
        comments)""")
    parser.add_argument('--ni', type=int, help="""
        number of tracer points in model grid 'x' direction""")
    parser.add_argument('--nj', type=int, help="""
        number of tracer points in model grid 'y' direction""")
    parser.add_argument('--lon1', type=float, help="""
        longitude of northwest corner point""")
    parser.add_argument('--lat1', type=float, help="""
        latitude of northwest corner point""")
    parser.add_argument('--lon2', type=float, help="""
        longitude of southeast corner point""")
    parser.add_argument('--lat2', type=float, help="""
        latitude of southeast corner point""")
    parser.add_argument('--lon_subscale', type=int, help="""
        subscale factor to be applied to each cell in the model grid 'x'
        direction (e.g., '2' doubles the number of x-direction cells)
        (integer>=1)""")
    parser.add_argument('--lat_subscale', type=int, help="""
        subscale factor to be applied to each cell in the model grid 'y'
        direction (see --lon_subscale comments) (integer>=1)""")
    parser.add_argument('--outfile', help="""
        file to which regridded matrices will be written (mitgridfile
        format)""")
    parser.add_argument('-v','--verbose',action='store_true',help="""
        verbose output""")
    return parser


def regrid( mitgridfile,xg_file,yg_file,ni,nj,
        lon1,lat1,lon2,lat2,lon_subscale,lat_subscale,verbose=False):
    """Regrids a rectangular lon/lat region using simple great circle-based
    subdivision, preserving any corner grids that may already exist within the
    region. A normal spherical geoid is currently assumed.

    Args:
        mitgridfile (str): (path and) filename of data to be regridded.
        xg_file (str): xg (path and) file input alternative to mitgridfile (csv
            if .csv extension (one matrix row per line), double-precision
            column-ordered binary otherwise; --yg_file must also be
            provided).
        yg_file (str): yg (path and) file input alternative to mitgridfile (see
            --xg_file comments)
        ni (int): number of tracer points in the model grid 'x' direction.
        nj (int): number of tracer points in the model grid 'y' direction.
        lon1 (float): longitude of northwest corner point.
        lat1 (float): latitude of northwest corner point.
        lon2 (float): longitude of southeast corner point.
        lat2 (float): latitude of southeast corner point.
        lon_subscale (int): subscale factor to be applied to each cell in the
            model grid 'x' direction (e.g., '2' doubles the number of
            x-direction cells; int>=1).
        lat_subscale (int): subscale factor to be applied to each cell in the
            model grid 'y' direction (see lon_subscale comments; int>=1).
        verbose (bool): True for diagnostic output, False otherwise.

    Returns:
        (newgrid,newgrid_ni,newgrid_nj): Tuple consisting of dictionary of
            newly-regridded subdomain matrices (ref. mitgridfilefields.py for
            names and ordering), and regridded ni, nj subdomain cell counts.

    Note:
        Corner point information (lon1/lat1, lon2/lat2) can be either literal
        northwest/southeast corners, or diagonal corners provided simply for
        orientation purposes. Either way, the output grid will be aligned such
        that 'northwest' is at min i, max j, and 'southeast' is at max i, min j.
    """

    # read XG, YG data from source provided:
    if mitgridfile:
        mitgrid = gridio.read_mitgridfile( mitgridfile, ni, nj, verbose)
    elif xg_file and yg_file:
        mitgrid = {key:None for key in mitgridfilefields.names}
        if os.path.splitext(xg_file)[1] == '.csv' and \
           os.path.splitext(yg_file)[1] == '.csv':
            # read .csv data, store in dictionary fields that read_mitgridfile
            # would have produced:
            mitgrid['XG'] = np.loadtxt(xg_file,delimiter=',')
            mitgrid['YG'] = np.loadtxt(yg_file,delimiter=',')
        else:
            # read binary column-ordered data, store in dictionary fields that
            # read_mitgridfile would have produced:
            XG_raw = np.fromfile(xg_file,mitgridfilefields.datatype)
            mitgrid['XG'] = np.reshape(XG_raw,(ni+1,nj+1),order='F')
            YG_raw = np.fromfile(yg_file,mitgridfilefields.datatype)
            mitgrid['YG'] = np.reshape(YG_raw,(ni+1,nj+1),order='F')
    else:
        raise ValueError(
            "Either an mitgridfile or xg_file/yg_file pair must be provided.")

    # for now, assume spherical geoid (perhaps user-specified later):
    geod = pyproj.Geod(ellps='sphere')

    # determine original XG, YG matrix indices corresponding to user input
    # lat/lon corners:

    i_nw,j_nw,_ = util.nearest(lon1,lat1,mitgrid['XG'],mitgrid['YG'],geod)
    i_se,j_se,_ = util.nearest(lon2,lat2,mitgrid['XG'],mitgrid['YG'],geod)

    # make sure diagonal corners and not degenerate sides have been input:
    if i_nw==i_se or j_nw==j_se:
        raise ValueError("lon/lat input pairs must be diagonally-opposed.")

    # if necessary, rotate xg, yg into user-directed nw/se alignment:
    while not (j_nw>j_se and i_se>i_nw):
        if verbose:
            print('performing 90 degree XG, YG rotation...')
            print('...start: (i,j)_nw = ({0},{1}), (i,j)_se = ({2},{3})'.format(i_nw,j_nw,i_se,j_se))
        mitgrid['XG'] = np.rot90(mitgrid['XG'])
        mitgrid['YG'] = np.rot90(mitgrid['YG'])
        i_nw,j_nw,_ = util.nearest(lon1,lat1,mitgrid['XG'],mitgrid['YG'],geod)
        i_se,j_se,_ = util.nearest(lon2,lat2,mitgrid['XG'],mitgrid['YG'],geod)
        if verbose:
            print('...end:   (i,j)_nw = ({0},{1}), (i,j)_se = ({2},{3})'.format(i_nw,j_nw,i_se,j_se))
    # assign to range-defining indices understood by subsequent code:
    i1,i2 = i_nw, i_se
    j1,j2 = j_nw, j_se

    if verbose:
        print('remeshing {0} x {1} cell grid based on'.format(
            abs(i2-i1),abs(j2-j1)))
        print('located corner points ({0:9.4f},{1:9.4f}) and ({2:9.4f},{3:9.4f})'.format(
            lon1,lat1,lon2,lat2))
        print('resulting grid will be {0} x {1} cells (lon/lat subscale = {2}/{3})'.format(
            abs(i2-i1)*lon_subscale,abs(j2-j1)*lat_subscale,lon_subscale,lat_subscale))

    # since matrix slice operations are from start:(stop-1):stride, define a
    # lambda function to make it clear when the stop point should be inclusive,
    # i.e., start:incl(stop):stride

    incl = lambda idx : idx+1

    #
    # Step 1:
    #
    # Based on user-selected extents, create a "plus one" matrix partition
    # consisting of original cell range plus a boundary "ring" one cell wide
    # to allow for computing boundary grid edge values.
    #
    # In this, and subsequent operations, the following grid index notation will
    # be useful:
    #
    #    (i, j)        - original (user-selected) grid indices
    #    (i_cg, jc_cg) - "compute grid" indices (doubled resolution)
    #    (i_ca, jc_ca) - "compute area" indices (doubled resolution)
    #

    # original grid index bounds based on user lon/lat selections:
    ilb, jlb = min(i1,i2), min(j1,j2)
    iub, jub = max(i1,i2), max(j1,j2)

    # "plus one" grid extents (note that if tile boundary points have been
    # selected, lower bounds may be negative, and/or upper bounds may be out of
    # range; both are ok in the current logic):
    iLB, iUB = ilb-1, iub+1
    jLB, jUB = jlb-1, jub+1

    # "compute grid" dimensions, allocation:
    num_compute_grid_rows = (iUB-iLB)*lon_subscale*2 + 1
    num_compute_grid_cols = (jUB-jLB)*lat_subscale*2 + 1
    compute_grid_xg = np.empty((num_compute_grid_rows,num_compute_grid_cols))
    compute_grid_xg[:,:] = np.nan
    compute_grid_yg = np.empty((num_compute_grid_rows,num_compute_grid_cols))
    compute_grid_yg[:,:] = np.nan

    # map mitgrid values to corresponding locations in compute_grid
    # (note: if user-selected range includes a tile boundary, "plus one" extents
    # may be out of range in XG and YG.  For purposes of compute_grid
    # initialization, make sure this does not happen):

    iLB_bounded = max(iLB,0)
    jLB_bounded = max(jLB,0)
    iUB_bounded = min(iUB,mitgrid['XG'].shape[0])
    jUB_bounded = min(jUB,mitgrid['XG'].shape[1])
    
    # index transformation from partitioned grid space to compute grid space:
    i_cg = lambda i,iLB,iLB_bounded,lon_subscale : \
            (iLB_bounded-iLB+i)*lon_subscale*2
    j_cg = lambda j,jLB,jLB_bounded,lat_subscale : \
            (jLB_bounded-jLB+j)*lat_subscale*2

    it = np.nditer(
        [mitgrid['XG'][iLB_bounded:incl(iUB_bounded),jLB_bounded:incl(jUB_bounded)],
         mitgrid['YG'][iLB_bounded:incl(iUB_bounded),jLB_bounded:incl(jUB_bounded)]],
        flags=['multi_index'])
    while not it.finished:
        compute_grid_xg[
            i_cg(it.multi_index[0],iLB,iLB_bounded,lon_subscale),
            j_cg(it.multi_index[1],jLB,jLB_bounded,lat_subscale)] = it[0]
        compute_grid_yg[
            i_cg(it.multi_index[0],iLB,iLB_bounded,lon_subscale),
            j_cg(it.multi_index[1],jLB,jLB_bounded,lat_subscale)] = it[1]
        it.iternext()

    if verbose:
        print('user-selected range:')
        print("mitgrid['XG']:")
        print(mitgrid['XG'][ilb:incl(iub),jlb:incl(jub)])
        print("mitgrid['YG']:")
        print(mitgrid['YG'][ilb:incl(iub),jlb:incl(jub)])

        print('user-selected range, "plus one":')
        print("mitgrid['XG'] plus one:")
        print(mitgrid['XG'][iLB_bounded:incl(iUB_bounded),jLB_bounded:incl(jUB_bounded)])
        print("mitgrid['YG'] plus one:")
        print(mitgrid['YG'][iLB_bounded:incl(iUB_bounded),jLB_bounded:incl(jUB_bounded)])

        print('user-selected range, plus one, mapped to compute_grid:')
        print('compute_grid_xg:')
        print(compute_grid_xg)
        print('compute_grid_yg:')
        print(compute_grid_yg)

    #
    # Step 2: "fill in" grid cell corner points at resolution given by
    # user-specified subdivision level, times two ("compute grid" resolution):
    #

    (compute_grid_xg,compute_grid_yg) = computegrid.cgfill(
        compute_grid_xg,compute_grid_yg,
        0,num_compute_grid_rows-1,0,num_compute_grid_cols-1,
        lon_subscale,lat_subscale,geod,verbose)

    #
    # Step 3: Generate areas for sub-quads at the compute_grid array resolution
    #

    # areas, from cartesian coordinates on the unit sphere, scaled according to
    # mean spherical ellipsoid radius:

    compute_areas = util.squad_uarea(
        util.lonlat2cart(compute_grid_xg,compute_grid_yg)) \
        * np.power(geod.a,2)

    if verbose:
        print('compute_areas:')
        print(compute_areas)

    #
    # Step 4: Create and fill in output quantities based on compute grid data:
    #

    outgrid = {key:None for key in mitgridfilefields.names}

    # compute regridded grid location quantities:
    #   XC - longitude east of center of grid (tracer) cell
    #   YC - latitude north of center of grid (tracer) cell
    #   XG - latitude east of southwest corner of grid (tracer) cell
    #   YG - latitude north of southwest corner of grid (tracer) cell

    # XC, YC directly from compute grid partitions:

    # compute grid partitioning:
    cg_first_i  = lon_subscale*2 + 1
    cg_last_i   = incl(cg_first_i + 2*((iub-ilb)*lon_subscale-1))
    cg_stride_i = 2
    cg_first_j  = lat_subscale*2 + 1
    cg_last_j   = incl(cg_first_j + 2*((jub-jlb)*lat_subscale-1))
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
    cg_first_i  = lon_subscale*2
    cg_last_i   = incl(cg_first_i + 2*(iub-ilb)*lon_subscale)
    cg_stride_i = 2
    cg_first_j  = lat_subscale*2
    cg_last_j   = incl(cg_first_j + 2*(jub-jlb)*lat_subscale)
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

    outgrid['DXG'] = np.empty(((iub-ilb)*lon_subscale, (jub-jlb)*lat_subscale+1))
    outgrid['DXG'][:,:] = np.nan    # error check, mostly

    # compute grid partitioning:
    cg_first_i  = 2*lon_subscale
    cg_last_i   = incl(cg_first_i + (iub-ilb)*2*(lon_subscale-1))
    cg_stride_i = 2
    cg_first_j  = 2*lat_subscale
    cg_last_j   = incl(cg_first_j + (jub-jlb)*2*lat_subscale)
    cg_stride_j = 2

    # transformation from partitioned, strided nditer space (i_n,j_n) to
    # underlying compute_grid space (i_cg,j_cg):
    i_cg = lambda i_n,lon_subscale : 2*i_n + 2*lon_subscale
    j_cg = lambda j_n,lat_subscale : 2*j_n + 2*lat_subscale

    it = np.nditer(
        [compute_grid_xg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_i],
         compute_grid_yg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_j]],
         flags=['multi_index'])
    while not it.finished:
        # (don't need to check for NaNs since all values are within valid index
        # range)
        _,_,outgrid['DXG'][it.multi_index[0],it.multi_index[1]] = geod.inv(
            it[0],                                          # lon1
            it[1],                                          # lat1
            compute_grid_xg[
                i_cg(it.multi_index[0]+1,lon_subscale),
                j_cg(it.multi_index[1]  ,lat_subscale)],    # lon2
            compute_grid_yg[
                i_cg(it.multi_index[0]+1,lon_subscale),
                j_cg(it.multi_index[1]  ,lat_subscale)])    # lat2
        it.iternext()
    if verbose:
        print("outgrid['DXG']:")
        print(outgrid['DXG'])

    # DYG tracer cell western edge computed from compute_grid grid point
    # distances:

    outgrid['DYG'] = np.empty(((iub-ilb)*lon_subscale+1,(jub-jlb)*lat_subscale))
    outgrid['DYG'][:,:] = np.nan    # error check, mostly

    # compute grid partitioning:
    cg_first_i  = 2*lon_subscale
    cg_last_i   = incl(cg_first_i + (iub-ilb)*2*lon_subscale)
    cg_stride_i = 2
    cg_first_j  = 2*lat_subscale
    cg_last_j   = incl(cg_first_j + ((jub-jlb)*lat_subscale-1)*2)
    cg_stride_j = 2

    # (i_n,j_n) to (i_cg,j_cg) transformations are the same as for DXG above.

    it = np.nditer(
        [compute_grid_xg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_j],
         compute_grid_yg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_j]],
        flags=['multi_index'])

    while not it.finished:
        # (don't need to check for NaNs since all values are within valid index
        # range)
        _,_,outgrid['DYG'][it.multi_index[0],it.multi_index[1]] = geod.inv(
            it[0],                                          # lon1
            it[1],                                          # lat1
            compute_grid_xg[
                i_cg(it.multi_index[0]  ,lon_subscale),
                j_cg(it.multi_index[1]+1,lat_subscale)],    # lon2
            compute_grid_yg[
                i_cg(it.multi_index[0]  ,lon_subscale),
                j_cg(it.multi_index[1]+1,lat_subscale)])    # lat2
        it.iternext()
    if verbose:
        print("outgrid['DYG']:")
        print(outgrid['DYG'])

    # RAC computed from subcell area sums in compute_areas:

    outgrid['RAC'] = np.empty(((iub-ilb)*lon_subscale, (jub-jlb)*lat_subscale))
    outgrid['RAC'][:,:] = np.nan    # error check, mostly

    # compute_areas partitioning:
    ca_first_i  = lon_subscale*2
    ca_last_i   = incl(ca_first_i + ((iub-ilb)*lon_subscale-1)*2)
    ca_stride_i = 2
    ca_first_j  = lat_subscale*2
    ca_last_j   = incl(ca_first_j + ((jub-jlb)*lat_subscale-1)*2)
    ca_stride_j = 2

    # transformation from partitioned, strided compute_areas indices (nditer) to
    # underlying compute_areas indices:
    i_ca = lambda i_n,lon_subscale : 2*i_n + 2*lon_subscale
    j_ca = lambda j_n,lat_subscale : 2*j_n + 2*lat_subscale

    it = np.nditer(
        compute_areas[ca_first_i:ca_last_i:ca_stride_i,ca_first_j:ca_last_j:ca_stride_j],
        flags=['multi_index'])
    while not it.finished:
        # (don't need to check for NaNs since all values are within valid index
        # range)
        outgrid['RAC'][it.multi_index[0],it.multi_index[1]] = \
            it[0] + \
            compute_areas[i_ca(it.multi_index[0],lon_subscale)+1,j_ca(it.multi_index[1],lat_subscale)  ] + \
            compute_areas[i_ca(it.multi_index[0],lon_subscale)+1,j_ca(it.multi_index[1],lat_subscale)+1] + \
            compute_areas[i_ca(it.multi_index[0],lon_subscale)  ,j_ca(it.multi_index[1],lat_subscale)+1]
        it.iternext()
    if verbose:
        print("outgrid['RAC']:")
        print(outgrid['RAC'])

    # vorticity cell-related quantities, DXC, DYC, RAZ
    #   DXC - vorticity cell edge length in X-direction
    #   DYC - vorticity cell edge length in Y-direction
    #   RAZ - vorticity cell area presented in the vertical direction

    # DXC vorticity cell edge lengths computed from compute_grid points since
    # endpoints are centered in tracer cells:

    outgrid['DXC'] = np.empty(((iub-ilb)*lon_subscale+1, (jub-jlb)*lat_subscale))
    outgrid['DXC'][:,:] = np.nan

    # compute grid partitioning:
    cg_first_i  = lon_subscale*2 - 1
    cg_last_i   = incl(cg_first_i + (iub-ilb)*lon_subscale*2)
    cg_stride_i = 2
    cg_first_j  = lat_subscale*2 + 1
    cg_last_j   = incl(cg_first_j + ((jub-jlb)*lat_subscale-1)*2)
    cg_stride_j = 2

    # transformation from partitioned, strided nditer space (i_n,j_n) to
    # underlying compute_grid space (i_cg,j_cg):
    i_cg = lambda i_n,lon_subscale : 2*lon_subscale-1 + 2*i_n
    j_cg = lambda j_n,lat_subscale : 2*lat_subscale+1 + 2*j_n

    it = np.nditer(
        [compute_grid_xg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_i],
         compute_grid_yg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_j]],
        flags=['multi_index'])
    while not it.finished:
        lon1 = it[0]
        lat1 = it[1]
        lon2 = compute_grid_xg[
            i_cg(it.multi_index[0]+1,lon_subscale),
            j_cg(it.multi_index[1]  ,lat_subscale)]
        lat2 = compute_grid_yg[
            i_cg(it.multi_index[0]+1,lon_subscale),
            j_cg(it.multi_index[1]  ,lat_subscale)]
        if not any(np.isnan((lon1,lat1,lon2,lat2))):
            _,_,outgrid['DXC'][it.multi_index[0],it.multi_index[1]] = geod.inv(
                lon1,lat1,lon2,lat2)
        it.iternext()
    if verbose:
        print("outgrid['DXC']:")
        print(outgrid['DXC'])

    # DYC vorticity cell edge lengths computed from compute_grid points:

    outgrid['DYC'] = np.empty(((iub-ilb)*lon_subscale, (jub-jlb)*lat_subscale+1))
    outgrid['DYC'][:,:] = np.nan

    # compute grid partitioning:
    cg_first_i  = lon_subscale*2 + 1
    cg_last_i   = incl(cg_first_i + ((iub-ilb)*lon_subscale-1)*2)
    cg_stride_i = 2
    cg_first_j  = lat_subscale*2 - 1
    cg_last_j   = incl(cg_first_j + (jub-jlb)*lat_subscale*2)
    cg_stride_j = 2

    # transformation from partitioned, strided nditer space (i_n,j_n) to
    # underlying compute_grid space (i_cg,j_cg):
    i_cg = lambda i_n,lon_subscale : 2*lon_subscale+1 + 2*i_n
    j_cg = lambda j_n,lat_subscale : 2*lat_subscale-1 + 2*j_n

    it = np.nditer(
        [compute_grid_xg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_i],
         compute_grid_yg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_j]],
        flags=['multi_index'])
    while not it.finished:
        lon1 = it[0]
        lat1 = it[1]
        lon2 = compute_grid_xg[
            i_cg(it.multi_index[0]  ,lon_subscale),
            j_cg(it.multi_index[1]+1,lat_subscale)]
        lat2 = compute_grid_yg[
            i_cg(it.multi_index[0]  ,lon_subscale),
            j_cg(it.multi_index[1]+1,lat_subscale)]
        if not any(np.isnan((lon1,lat1,lon2,lat2))):
            _,_,outgrid['DYC'][it.multi_index[0],it.multi_index[1]] = geod.inv(
                lon1,lat1,lon2,lat2)
        it.iternext()
    if verbose:
        print("outgrid['DYC']:")
        print(outgrid['DYC'])

    # RAZ vorticity cell areas computed from subcell areas in compute_areas:

    outgrid['RAZ'] = np.empty(((iub-ilb)*lon_subscale+1, (jub-jlb)*lat_subscale+1))
    outgrid['RAZ'][:,:] = np.nan

    # compute_areas partitioning:
    ca_first_i  = lon_subscale*2 - 1
    ca_last_i   = incl(ca_first_i + (iub-ilb)*lon_subscale*2)
    ca_stride_i = 2
    ca_first_j  = lat_subscale*2 - 1
    ca_last_j   = incl(ca_first_j + (jub-jlb)*lat_subscale*2)
    ca_stride_j = 2

    # transformation from partitioned, strided compute_areas indices (nditer) to
    # underlying compute_areas indices:
    i_ca = lambda i_n,lon_subscale : 2*lon_subscale - 1 + 2*i_n
    j_ca = lambda j_n,lat_subscale : 2*lat_subscale - 1 + 2*j_n

    it = np.nditer(
        compute_areas[ca_first_i:ca_last_i:ca_stride_i,ca_first_j:ca_last_j:ca_stride_j],
        flags=['multi_index'])
    while not it.finished:
        outgrid['RAZ'][it.multi_index[0],it.multi_index[1]] = \
            it[0] + \
            compute_areas[i_ca(it.multi_index[0],lon_subscale)+1,j_ca(it.multi_index[1],lat_subscale)  ] + \
            compute_areas[i_ca(it.multi_index[0],lon_subscale)+1,j_ca(it.multi_index[1],lat_subscale)+1] + \
            compute_areas[i_ca(it.multi_index[0],lon_subscale)  ,j_ca(it.multi_index[1],lat_subscale)+1]
        it.iternext()
    if verbose:
        print("outgrid['RAZ']:")
        print(outgrid['RAZ'])

    # "U" cell-related quantities, DXC, DYC, RAZ
    #   DXV - U cell edge length in X-direction between v-points
    #   DYF - U cell edge length in Y-direction between tracer cell faces
    #   RAW - U cell area presented in the vertical direction

    # DXV U cell edge lengths computed from compute_grid points: 

    outgrid['DXV'] = np.empty(((iub-ilb)*lon_subscale+1,(jub-jlb)*lat_subscale+1))
    outgrid['DXV'][:,:] = np.nan

    # compute grid partitioning:
    cg_first_i  = lon_subscale*2 -1
    cg_last_i   = incl(cg_first_i + (iub-ilb)*lon_subscale*2)
    cg_stride_i = 2
    cg_first_j  = lat_subscale*2
    cg_last_j   = incl(cg_first_j + (jub-jlb)*lat_subscale*2)
    cg_stride_j = 2

    # transformation from partitioned, strided nditer space (i_n,j_n) to
    # underlying compute_grid space (i_cg,j_cg):
    i_cg = lambda i_n,lon_subscale : 2*lon_subscale - 1 + 2*i_n
    j_cg = lambda j_n,lat_subscale : 2*lat_subscale     + 2*j_n

    it = np.nditer(
        [compute_grid_xg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_i],
         compute_grid_yg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_j]],
        flags=['multi_index'])
    while not it.finished:
        lon1 = it[0]
        lat1 = it[1]
        lon2 = compute_grid_xg[
            i_cg(it.multi_index[0]+1,lon_subscale),
            j_cg(it.multi_index[1]  ,lat_subscale)]
        lat2 = compute_grid_yg[
            i_cg(it.multi_index[0]+1,lon_subscale),
            j_cg(it.multi_index[1]  ,lat_subscale)]
        if not any(np.isnan((lon1,lat1,lon2,lat2))):
            _,_,outgrid['DXV'][it.multi_index[0],it.multi_index[1]] = geod.inv(
                lon1,lat1,lon2,lat2)
        it.iternext()
    if verbose:
        print("outgrid['DXV']:")
        print(outgrid['DXV'])

    # DYF U cell edge lengths computed from compute_grid points:

    outgrid['DYF'] = np.empty(((iub-ilb)*lon_subscale,(jub-jlb)*lat_subscale))
    outgrid['DYF'][:,:] = np.nan

    # compute grid partitioning:
    cg_first_i  = lon_subscale*2 + 1
    cg_last_i   = incl(cg_first_i + ((iub-ilb)*lon_subscale-1)*2)
    cg_stride_i = 2
    cg_first_j  = lat_subscale*2
    cg_last_j   = incl(cg_first_j + ((jub-jlb)*lat_subscale-1)*2)
    cg_stride_j = 2

    # transformation from partitioned, strided nditer space (i_n,j_n) to
    # underlying compute_grid space (i_cg,j_cg):
    i_cg = lambda i_n,lon_subscale : 2*lon_subscale + 1 + 2*i_n
    j_cg = lambda j_n,lat_subscale : 2*lat_subscale     + 2*j_n

    it = np.nditer(
        [compute_grid_xg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_i],
         compute_grid_yg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_j]],
        flags=['multi_index'])
    # note isnan checks not required for DYF
    while not it.finished:
        _,_,outgrid['DYF'][it.multi_index[0],it.multi_index[1]] = geod.inv(
            it[0],                                          # lon1
            it[1],                                          # lat1
            compute_grid_xg[
                i_cg(it.multi_index[0]  ,lon_subscale),
                j_cg(it.multi_index[1]+1,lat_subscale)],    # lon2
            compute_grid_yg[
                i_cg(it.multi_index[0]  ,lon_subscale),
                j_cg(it.multi_index[1]+1,lat_subscale)])    # lat2
        it.iternext()
    if verbose:
        print("outgrid['DYF']:")
        print(outgrid['DYF'])

    # RAW vertical face area of U cell computed from subcell areas in
    # compute_areas:

    outgrid['RAW'] = np.empty(((iub-ilb)*lon_subscale+1,(jub-jlb)*lat_subscale))
    outgrid['RAW'][:,:] = np.nan

    # compute areas partitioning:
    ca_first_i  = lon_subscale*2 - 1
    ca_last_i   = incl(ca_first_i + (iub-ilb)*lon_subscale*2)
    ca_stride_i = 2
    ca_first_j  = lat_subscale*2
    ca_last_j   = incl(ca_first_j + ((jub-jlb)*lat_subscale-1)*2)
    ca_stride_j = 2

    # transformation from partitioned, strided compute_areas indices (nditer) to
    # underlying compute_areas indices:
    i_ca = lambda i_n,lon_subscale : 2*lon_subscale - 1 + 2*i_n
    j_ca = lambda j_n,lat_subscale : 2*lat_subscale     + 2*j_n

    it = np.nditer(
        compute_areas[ca_first_i:ca_last_i:ca_stride_i,ca_first_j:ca_last_j:ca_stride_j],
        flags=['multi_index'])
    while not it.finished:
        outgrid['RAW'][it.multi_index[0],it.multi_index[1]] = \
            it[0] + \
            compute_areas[i_ca(it.multi_index[0],lon_subscale)+1,j_ca(it.multi_index[1],lat_subscale)  ] + \
            compute_areas[i_ca(it.multi_index[0],lon_subscale)+1,j_ca(it.multi_index[1],lat_subscale)+1] + \
            compute_areas[i_ca(it.multi_index[0],lon_subscale)  ,j_ca(it.multi_index[1],lat_subscale)+1]
        it.iternext()
    if verbose:
        print("outgrid['RAW']:")
        print(outgrid['RAW'])

    # "V" cell-related quantities, DXF, DYU, RAS
    #   DXF - V cell northern edge length in X-direction between u-points
    #   DYU - V cell western edge length in Y-direction
    #   RAS - V cell area presented in the vertical direction

    # DXF V cell edge lengths computed using compute_grid point distances:

    outgrid['DXF'] = np.empty(((iub-ilb)*lon_subscale,(jub-jlb)*lat_subscale))
    outgrid['DXF'][:,:] = np.nan

    # compute grid partitioning:
    cg_first_i  = lon_subscale*2
    cg_last_i   = incl(cg_first_i + ((iub-ilb)*lon_subscale-1)*2)
    cg_stride_i = 2
    cg_first_j  = lat_subscale*2 + 1
    cg_last_j   = incl(cg_first_j + ((jub-jlb)*lat_subscale-1)*2)
    cg_stride_j = 2

    # transformation from partitioned, strided nditer space (i_n,j_n) to
    # underlying compute_grid space (i_cg,j_cg):
    i_cg = lambda i_n,lon_subscale : 2*lon_subscale     + 2*i_n
    j_cg = lambda j_n,lat_subscale : 2*lat_subscale + 1 + 2*j_n

    it = np.nditer(
        [compute_grid_xg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_i],
         compute_grid_yg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_j]],
        flags=['multi_index'])
    # note isnan checks not required for DXF
    while not it.finished:
        _,_,outgrid['DXF'][it.multi_index[0],it.multi_index[1]] = geod.inv(
            it[0],                                          # lon1
            it[1],                                          # lat1
            compute_grid_xg[
                i_cg(it.multi_index[0]+1,lon_subscale),
                j_cg(it.multi_index[1]  ,lat_subscale)],    # lon2
            compute_grid_yg[
                i_cg(it.multi_index[0]+1,lon_subscale),
                j_cg(it.multi_index[1]  ,lat_subscale)])    # lat2
        it.iternext()
    if verbose:
        print("outgrid['DXF']:")
        print(outgrid['DXF'])

    # DYU V cell western edge length computed from compute_grid point distances:

    outgrid['DYU'] = np.empty(((iub-ilb)*lon_subscale+1,(jub-jlb)*lat_subscale+1))
    outgrid['DYU'][:,:] = np.nan

    # compute grid partitioning:
    cg_first_i  = 2*lon_subscale
    cg_last_i   = incl(cg_first_i + (iub-ilb)*2*lon_subscale)
    cg_stride_i = 2
    cg_first_j  = 2*lat_subscale-1
    cg_last_j   = incl(cg_first_j + (jub-jlb)*2*lat_subscale)
    cg_stride_j = 2

    # transformation from partitioned, strided nditer space (i_n,j_n) to
    # underlying compute_grid space (i_cg,j_cg):
    i_cg = lambda i_n,lon_subscale : 2*lon_subscale     + 2*i_n
    j_cg = lambda j_n,lat_subscale : 2*lat_subscale - 1 + 2*j_n

    it = np.nditer(
        [compute_grid_xg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_j],
         compute_grid_yg[cg_first_i:cg_last_i:cg_stride_i,cg_first_j:cg_last_j:cg_stride_j]],
        flags=['multi_index'])
    while not it.finished:
        lon1 = it[0]
        lat1 = it[1]
        lon2 = compute_grid_xg[
            i_cg(it.multi_index[0]  ,lon_subscale),
            j_cg(it.multi_index[1]+1,lat_subscale)]
        lat2 = compute_grid_yg[
            i_cg(it.multi_index[0]  ,lon_subscale),
            j_cg(it.multi_index[1]+1,lat_subscale)]
        if not any(np.isnan((lon1,lat1,lon2,lat2))):
            _,_,outgrid['DYU'][it.multi_index[0],it.multi_index[1]] = geod.inv(
                lon1,lat1,lon2,lat2)
        it.iternext()
    if verbose:
        print("outgrid['DYU']:")
        print(outgrid['DYU'])

    # RAS vertical face area of V cell computed from subcell areas in
    # compute_areas:

    outgrid['RAS'] = np.empty(((iub-ilb)*lon_subscale,(jub-jlb)*lat_subscale+1))
    outgrid['RAS'][:,:] = np.nan

    # compute areas partitioning:
    ca_first_i  = 2*lon_subscale
    ca_last_i   = incl(ca_first_i + ((iub-ilb)*lon_subscale-1)*2)
    ca_stride_i = 2
    ca_first_j  = 2*lat_subscale - 1
    ca_last_j   = incl(ca_first_j + (jub-jlb)*lat_subscale*2)
    ca_stride_j = 2

    # transformation from partitioned, strided compute_areas indices (i_n,j_n) to
    # underlying compute_areas indices:
    i_ca = lambda i_n,lon_subscale : 2*lon_subscale     + 2*i_n
    j_ca = lambda j_n,lat_subscale : 2*lat_subscale - 1 + 2*j_n

    it = np.nditer(
        compute_areas[ca_first_i:ca_last_i:ca_stride_i,ca_first_j:ca_last_j:ca_stride_j],
        flags=['multi_index'])
    while not it.finished:
        outgrid['RAS'][it.multi_index[0],it.multi_index[1]] = \
            it[0] + \
            compute_areas[i_ca(it.multi_index[0],lon_subscale)+1,j_ca(it.multi_index[1],lat_subscale)  ] + \
            compute_areas[i_ca(it.multi_index[0],lon_subscale)+1,j_ca(it.multi_index[1],lat_subscale)+1] + \
            compute_areas[i_ca(it.multi_index[0],lon_subscale)  ,j_ca(it.multi_index[1],lat_subscale)+1]
        it.iternext()
    if verbose:
        print("outgrid['RAS']:")
        print(outgrid['RAS'])

    return outgrid, (iub-ilb)*lon_subscale, (jub-jlb)*lat_subscale


def main():
    """Command-line entry point."""

    parser = create_parser()
    args = parser.parse_args()
    (newgrid,ni_regridded,nj_regridded) = regrid(
        args.mitgridfile,
        args.xg_file,
        args.yg_file,
        args.ni,
        args.nj,
        args.lon1,
        args.lat1,
        args.lon2,
        args.lat2,
        args.lon_subscale,
        args.lat_subscale,
        args.verbose)
    if args.verbose:
        print('writing {0:s} with ni={1:d}, nj={2:d}...'.
            format(args.outfile,ni_regridded,nj_regridded))
    gridio.write_mitgridfile(args.outfile,newgrid,ni_regridded,nj_regridded)
    if args.verbose:
        print('...done.')

if __name__ == '__main__':
    main()

