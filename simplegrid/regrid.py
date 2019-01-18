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


def regrid( verbose=False, **kwargs):
    """Regrids a rectangular lon/lat region using simple great circle-based
    subdivision, preserving any corner grids that may already exist within the
    region. A normal spherical geoid is currently assumed.

    Args:
        verbose (bool, optional): True for diagnostic output, False otherwise.
        **kwargs: Arbitrary keyword arguments.

    Kwargs:
        mitgridfile (str, required if no xg_file, yg_file): (path and) filename
            of data to be regridded.
        xg_file (str, required if no mitgridfile): xg (path and) file input
            alternative to mitgridfile (csv if .csv extension (one matrix row
            per line), double-precision column-ordered binary otherwise; yg_file
            must also be provided).
        yg_file (str, required if no mitgridfile): yg (path and) file input
            alternative to mitgridfile (see xg_file comments)
        ni (int, required): number of tracer points in the model grid 'x'
            direction.
        nj (int, required): number of tracer points in the model grid 'y'
            direction.
        lon1 (float, required): longitude of northwest corner point.
        lat1 (float, required): latitude of northwest corner point.
        lon2 (float, required): longitude of southeast corner point.
        lat2 (float, required): latitude of southeast corner point.
        lon_subscale (int, required): subscale factor to be applied to each cell
            in the model grid 'x' direction (e.g., '2' doubles the number of
            x-direction cells; int>=1).
        lat_subscale (int, required): subscale factor to be applied to each cell
            in the model grid 'y' direction (see lon_subscale comments; int>=1).

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

    # kwarg handling:
    mitgridfile     = kwargs.get('mitgridfile')
    xg_file         = kwargs.get('xg_file')
    yg_file         = kwargs.get('yg_file')
    ni              = kwargs.get('ni')
    nj              = kwargs.get('nj')
    lon1            = kwargs.get('lon1')
    lat1            = kwargs.get('lat1')
    lon2            = kwargs.get('lon2')
    lat2            = kwargs.get('lat2')
    lon_subscale    = kwargs.get('lon_subscale')
    lat_subscale    = kwargs.get('lat_subscale')

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

    if verbose:
        print('remeshing {0} x {1} cell grid based on'.format(
            abs(i_nw-i_se),abs(j_nw-j_se)))
        print('located corner points ({0:9.4f},{1:9.4f}) and ({2:9.4f},{3:9.4f})'.format(
            lon1,lat1,lon2,lat2))
        print('resulting grid will be {0} x {1} cells (lon/lat subscale = {2}/{3})'.format(
            abs(i_nw-i_se)*lon_subscale,abs(j_nw-j_se)*lat_subscale,lon_subscale,lat_subscale))

    # since matrix slice operations are from start:(stop-1):stride, define a
    # lambda function to make it clear when the stop point should be inclusive,
    # i.e., start:incl(stop):stride
    incl = lambda idx : idx+1

    #
    # Step 1:
    #
    # Based on user-selected corner points and discretization level, create a
    # "compute grid" that spans the selected NW/SE range, at double the
    # resolution, plus a boundary "ring" one compute cell wide.  Initializing
    # boundary grid values to NaNs will allow us to compute subsequent mitgrid
    # values using consistent indexing, while naturally producing undefined
    # values at the boundaries. The total compute grid ranges are given by
    # LB/UB, while the user-selected range, given by ilb_mitgrid, iub_mitgrid,
    # jlb_mitgrid, jub_mitgrid, is mapped (along with all intermediate mitgrid
    # points) to the corresponding compute grid range, lb/ub; a picture might
    # help:
    #
    #   '<->' == corner point mapping
    #    ...  == range of intermediate grid point mapping
    #
    #                   y,j ^
    #                       |
    #                   jUB +------------------+
    #                       |                  |
    #   jub_mitgrid <-> jub +   o----------+   |<-- ring of NaN values
    #        .              |   |          |   |
    #        .              |   |          |<--|--- area to be regridded
    #        .              |   |          |   |
    #        .              |   |          |   |
    #   jlb_mitgrid <-> jlb +   +----------o   |    o = user-specified NW/SE corners
    #                       |                  |
    #                   jLB +---+----------+---+--> x,i
    #                      iLB ilb        iub iUB
    #                           ^          ^
    #                           |          |
    #                           v          v
    #                  ilb_mitgrid  ...  iub_mitgrid
    #

    # mitgrid index bounds:
    ilb_mitgrid, jlb_mitgrid = min(i_nw,i_se), min(j_nw,j_se)
    iub_mitgrid, jub_mitgrid = max(i_nw,i_se), max(j_nw,j_se)

    # compute grid bounds:
    iLB = 0
    ilb = 1
    iub = ilb + 2*lon_subscale*(iub_mitgrid-ilb_mitgrid)
    iUB = iub + 1
    jLB = 0
    jlb = 1
    jub = jlb + 2*lat_subscale*(jub_mitgrid-jlb_mitgrid)
    jUB = jub + 1

    # compute grid initialization, allocation:
    num_compute_grid_rows = iUB + 1
    num_compute_grid_cols = jUB + 1
    compute_grid_xg = np.empty((num_compute_grid_rows,num_compute_grid_cols))
    compute_grid_xg[:,:] = np.nan
    compute_grid_yg = np.empty((num_compute_grid_rows,num_compute_grid_cols))
    compute_grid_yg[:,:] = np.nan

    # map mitgrid values to corresponding compute_grid locations:
    
    # index transformation from partitioned grid space (e.g., between lower and
    # upper bounds) to compute grid space:
    i_cg = lambda i_mitgrid,lon_subscale : 1 + 2*lon_subscale*i_mitgrid
    j_cg = lambda j_mitgrid,lat_subscale : 1 + 2*lat_subscale*j_mitgrid

    it = np.nditer(
        [mitgrid['XG'][ilb_mitgrid:incl(iub_mitgrid),jlb_mitgrid:incl(jub_mitgrid)],
         mitgrid['YG'][ilb_mitgrid:incl(iub_mitgrid),jlb_mitgrid:incl(jub_mitgrid)]],
        flags=['multi_index'])
    while not it.finished:
        compute_grid_xg[
            i_cg(it.multi_index[0],lon_subscale),
            j_cg(it.multi_index[1],lat_subscale)] = it[0]
        compute_grid_yg[
            i_cg(it.multi_index[0],lon_subscale),
            j_cg(it.multi_index[1],lat_subscale)] = it[1]
        it.iternext()

    if verbose:
        print('user-selected range:')
        print("mitgrid['XG']:")
        print(mitgrid['XG'][
            ilb_mitgrid:incl(iub_mitgrid),
            jlb_mitgrid:incl(jub_mitgrid)])
        print("mitgrid['YG']:")
        print(mitgrid['YG'][
            ilb_mitgrid:incl(iub_mitgrid),
            jlb_mitgrid:incl(jub_mitgrid)])
        print('user-selected range, mapped to compute_grid:')
        print('compute_grid_xg:')
        print(compute_grid_xg)
        print('compute_grid_yg:')
        print(compute_grid_yg)

    #
    # Step 2: "fill in" grid cell corner points at resolution given by
    # user-specified subdivision level, times two ("compute grid" resolution):
    #

    (compute_grid_xg,compute_grid_yg) = computegrid.fill(
        compute_grid_xg,compute_grid_yg, ilb,iub,jlb,jub,
        lon_subscale,lat_subscale,geod,verbose)

    #
    # Step 3: Use compute grid to generate full set of mitgrid data:
    #

    outgrid = computegrid.tomitgrid( compute_grid_xg, compute_grid_yg,
        ilb, iub, jlb, jub, geod, verbose)

    return (
        outgrid,
        (iub_mitgrid-ilb_mitgrid)*lon_subscale,
        (jub_mitgrid-jlb_mitgrid)*lat_subscale)


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

