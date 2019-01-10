#!/usr/bin/env python

import argparse
import numpy as np
import pyproj
from . import computegrid
from . import gridio
from . import matchedges
from . import mitgridfilefields as mgf
from . import util


def create_parser():
    """Set up the list of arguments to be provided to addfringe.
    """
    parser = argparse.ArgumentParser(
        description="""
            Determine whether or not a common edge exists between two tiles and,
            if so, use data from the second ('B') to compute boundary grid data
            for the first ('A').""")
    parser.add_argument('--tilea', help="""
        (path and) file name of tile whose boundary grid information is to be
        computed using tile 'B' data (mitgridfile format)""")
    parser.add_argument('--nia', type=int, help="""
        number of tracer points in model grid 'x' direction for tile A""")
    parser.add_argument('--nja', type=int, help="""
        number of tracer points in model grid 'y' direction for tile A""")
    parser.add_argument('--tileb', help="""
        (path and) file name of tile whose data will be used to compute tile 'A'
        boundary grid information (mitgridfile format)""")
    parser.add_argument('--nib', type=int, help="""
        number of tracer points in model grid 'x' direction for tile B""")
    parser.add_argument('--njb', type=int, help="""
        number of tracer points in model grid 'y' direction for tile B""")
    parser.add_argument('--outfile', help="""
        file to which updated tile 'A' will be written (mitgridfile format)""")
    parser.add_argument('-v','--verbose',action='store_true',help="""
        verbose output""")
    return parser


def addfringe( tilea, nia, nja, tileb, nib, njb, verbose=False):
    """Determine whether or not a common edge exists between two tiles and, if
    so, use data from the second ('B') to compute boundary grid data for the
    first ('A').

    Args:
        tilea (str): mitgrid (path and) filename of tile for which boundary data
            are to be computed.
        nia (int): number of tracer points in the model grid 'x' direction for
            tilea.
        nja (int): number of tracer points in the model grid 'y' direction for
            tilea.
        tileb (str): mitgrid (path and) filename of tile that will be used to
            compute boundary data for tilea.
        nib (int): number of tracer points in the model grid 'x' direction for
            tileb.
        njb (int): number of tracer points in the model grid 'y' direction for
            tileb.
        verbose (bool): True for diagnostic output, False otherwise.

    Returns:
        (tilea_edge,tileb_edge,new_tilea_mitgrid): tuple of tilea and tileb
            matching edge indicators (integer 0 (N), 1 (ref. matchedges.matchedges return values),
            and copy of tilea with updated boundary grid data.

        (tilea_edge,tileb_edge,new_tilea_mitgrid): tuple of tilea and tileb
            matching edge indicators (integer 0(N), 1(S), 2(E) or 3(W)), and
            copy of tilea (dictionary of named numpy arrays) with updated
            boundary grid data.

    """

    tilea_mitgrid = gridio.read_mitgridfile( tilea, nia, nja, verbose)
    tileb_mitgrid = gridio.read_mitgridfile( tileb, nib, njb, verbose)

    geod = pyproj.Geod(ellps='sphere')

    (tilea_edge, tilea_edge_slice,
    tileb_edge, tileb_edge_slice,
    compute_grid_edge_xg,compute_grid_edge_yg,
    compute_grid_edge_join_slice) = matchedges.matchedges(
        tilea_mitgrid['XG'], tilea_mitgrid['YG'],
        tileb_mitgrid['XG'], tileb_mitgrid['YG'],
        geod, verbose)

    # possible future checks for tilea_edge, tileb_edge coincidence/averaging...

    # recompute (regrid) tile 'a', augmenting the compute grid with matching
    # edge data from tile 'b'.  terms that had prevously been NaNs define the
    # set of updates to tile 'a'.

    # compute grid initialization, allocation (ref. regrid(), with
    # lon_subscale=lat_subscale=1):
    iLB = 0
    ilb = 1
    iub = ilb + 2*(tilea_mitgrid['XG'].shape[0]-1)  # could use either XG or YG...
    iUB = iub + 1
    jLB = 0
    jlb = 1
    jub = jlb + 2*(tilea_mitgrid['XG'].shape[1]-1)  # ...for dimensioning
    jUB = jub + 1
    cg_rows = iUB + 1
    cg_cols = jUB + 1
    compute_grid_xg = np.empty((cg_rows,cg_cols))
    compute_grid_xg[:,:] = np.nan
    compute_grid_yg = np.empty((cg_rows,cg_cols))
    compute_grid_yg[:,:] = np.nan
    # index transformation from partitioned grid space to compute grid space:
    i_cg = lambda i_mitgrid : 1 + 2*i_mitgrid
    j_cg = lambda j_mitgrid : 1 + 2*j_mitgrid
    it = np.nditer(
        [tilea_mitgrid['XG'],tilea_mitgrid['YG']],
        flags=['multi_index'])
    while not it.finished:
        compute_grid_xg[i_cg(it.multi_index[0]),j_cg(it.multi_index[1])] = it[0]
        compute_grid_yg[i_cg(it.multi_index[0]),j_cg(it.multi_index[1])] = it[1]
        it.iternext()

    # augment tile 'a' compute grid with interpolated tile 'b' edge cell
    # data...:
    compute_grid_xg[compute_grid_edge_join_slice] = compute_grid_edge_xg
    compute_grid_yg[compute_grid_edge_join_slice] = compute_grid_edge_yg

    # ...and recompute:
    lon_subscale = 1
    lat_subscale = 1
    (compute_grid_xg,compute_grid_yg) = computegrid.fill(
        compute_grid_xg, compute_grid_yg, ilb, iub, jlb, jub,
        lon_subscale, lat_subscale, geod, verbose)

    new_tilea_mitgrid = computegrid.tomitgrid(
        compute_grid_xg, compute_grid_yg, ilb, iub, jlb, jub, geod, verbose)

    # compare tilea with new_tilea, replacing NaNs in the former with updated
    # values from the latter (actually performed somewhat in the reverse, so we
    # can return new_tilea merged with tilea computed values):

    for name in mgf.names:
        it = np.nditer(
            [tilea_mitgrid[name],new_tilea_mitgrid[name]],
             flags=['multi_index'],
             op_flags=['readwrite'])
        while not it.finished:
            if not np.isnan(it[0]):
                new_tilea_mitgrid[name][it.multi_index[0],it.multi_index[1]] = \
                    tilea_mitgrid[name][it.multi_index[0],it.multi_index[1]]
            it.iternext()

    return (tilea_edge,tileb_edge,new_tilea_mitgrid)


def main():
    """Command-line entry point."""

    parser = create_parser()
    args = parser.parse_args()

    if args.verbose:
        print('computing boundary grid information for {0:s},'.format(
            args.tilea))
        print('using data from {0:s}...'.format(args.tileb))

    (_,_,tilea_with_new_boundary) = addfringe(
        args.tilea, args.nia, args.nja,
        args.tileb, args.nib, args.njb,
        args.verbose)

    if args.verbose:
        print('writing {0:s} with ni={1:d}, nj={2:d}...'.format(
            args.outfile,args.nia,args.nja))

    gridio.write_mitgridfile(args.outfile,tilea_with_new_boundary,args.nia,args.nja)

    if args.verbose:
        print('...done.')

if __name__ == '__main__':
    main()

