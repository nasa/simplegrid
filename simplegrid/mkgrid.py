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

    #
    # Step 1:
    #
    # Based on user-selected discretization level, create a "compute grid" that
    # spans the expected range, plus a boundary "ring" one compute cell wide.
    # Initializing boundary grid values to np.PZERO will allow us to compute
    # subsequent mitgrid values using consistent indexing, while naturally
    # producing undefined values at the boundaries. The total compute grid
    # ranges are given by LB/UB, while the user-selected range is given by
    # lb/ub; a picture might help:
    #
    #
    #       y,j ^
    #           |
    #       jUB +------------------+
    #           |                  |
    #       jub +   o----------+   |<-- ring of PZERO values
    #           |   |          |   |
    #           |   |          |<--|--- area to be regridded
    #           |   |          |   |    (lat_subscale x lon_subscale)
    #           |   |          |   |
    #       jlb +   +----------o   |    o = user-specified NW/SE corners
    #           |                  |
    #       jLB +---+----------+---+--> x,i
    #          iLB ilb        iub iUB

    iLB = 0
    ilb = 1
    iub = ilb + 2*lon_subscale  # compute grids are at 2x resolution
    iUB = iub + 1

    jLB = 0
    jlb = 1
    jub = jlb + 2*lat_subscale  # compute grids are at 2x resolution
    jUB = jub + 1

    # "compute grid" dimensions, allocation:
    num_compute_grid_rows = iUB + 1
    num_compute_grid_cols = jUB + 1
    compute_grid_xg = np.zeros((num_compute_grid_rows,num_compute_grid_cols))
    compute_grid_yg = np.zeros((num_compute_grid_rows,num_compute_grid_cols))

    #
    # Step 2: Populate compute grid points corresponding to user-selected corner
    #         point range:
    #

    # 2a: user-specified northwest, southeast corners:

    min_lon, max_lon = lon1, lon2   # nw, se
    min_lat, max_lat = lat2, lat1   # se, nw
    compute_grid_xg[ilb,jlb] = compute_grid_xg[ilb,jub] = min_lon
    compute_grid_xg[iub,jlb] = compute_grid_xg[iub,jub] = max_lon
    compute_grid_yg[ilb,jlb] = compute_grid_yg[iub,jlb] = min_lat
    compute_grid_yg[ilb,jub] = compute_grid_yg[iub,jub] = max_lat

    # 2b: fill in intermediate points per user-specified subdivision level:

    (compute_grid_xg,compute_grid_yg) = computegrid.fill(
        compute_grid_xg,compute_grid_yg,ilb,iub,jlb,jub,
        lon_subscale,lat_subscale,geod,verbose)

    #
    # Step 3: Use compute grid to generate full set of mitgrid data:
    #

    outgrid = computegrid.tomitgrid( compute_grid_xg, compute_grid_yg,
        iLB, ilb, iub, iUB, jLB, jlb, jub, jUB, geod, verbose)

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

