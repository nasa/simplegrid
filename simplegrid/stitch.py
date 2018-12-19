#!/usr/bin/env/python

import argparse
import numpy as np
from . import gridio
from . import util
from .util import (N,S,E,W)


def create_parser():
    """Set up the list of arguments to be provided to stitch.
    """
    parser = argparse.ArgumentParser(
        description="""
            Join two tiles along a common edge.""")
    parser.add_argument('--tilea', help="""
        (path and) file name of first tile (mitgridfile format)""")
    parser.add_argument('--nia', type=int, help="""
        number of tracer points in model grid 'x' direction for the first
        tile""")
    parser.add_argument('--nja', type=int, help="""
        number of tracer points in model grid 'y' direction for the first
        tile""")
    parser.add_argument('--tileb', help="""
        (path and) file name of second tile (mitgridfile format)""")
    parser.add_argument('--nib', type=int, help="""
        number of tracer points in model grid 'x' direction for the second
        tile""")
    parser.add_argument('--njb', type=int, help="""
        number of tracer points in model grid 'y' direction for the second
        tile""")
    parser.add_argument('--outfile', help="""
        file to which combined tiles will be written (mitgridfile format)""")
    parser.add_argument('-v','--verbose',action='store_true',help="""
        verbose output""")
    return parser


def stitch(tilea_file,nia,nja,tileb_file,nib,njb,verbose):
    """Join two tiles along a common edge

    Args:
        tilea_file (str): (path and) filename of first tile (mitgridfile
            format).
        nia (int): number of tracer points in the model grid 'x' direction for
            tile a.
        nja (int): number of tracer points in the model grid 'y' direction for
            tile a.
        tileb_file (str): (path and) filename of second tile (mitgridfile
            format).
        nib (int): number of tracer points in the model grid 'x' direction for
            tile b.
        njb (int): number of tracer points in the model grid 'y' direction for
            tile b.
        verbose (bool): True for diagnostic output, False otherwise.

    Returns:
        (c,nic,njc): combined tile and corresponding tracer point counts in the
            combined model grid 'x' and 'y' directions.

    """

    # calculate "cross-border" terms using addfringe(), store in tile 'a' with
    # updated 'fringe':
    (tilea_edge,tileb_edge,a) = addfringe(
        tilea_file,nia,nja,tileb_file,nib,njb,verbose)

    # read tile b into mitgrid data structure:
    b = gridio.read_mitgridfile( tileb_file, nib, njb, verbose)

    # initialization:
    c = {key:None for key in mitgridfilefields.names}

    # for convenience:
    edges = (tilea_edge,tileb_edge)

    if edges==(E,W):
        # east edge of 'a' and west edge of 'b':
        # grids:
        c['XG'] = np.concatenate(a['XG'],b['XG'][1:,:],0)
        c['YG'] = np.concatenate(a['YG'],b['YG'][1:,:],0)
        # grid cell centers:
        c['XC'] = np.concatenate(a['XC'],b['XC'],0)
        c['YC'] = np.concatenate(a['YC'],b['YC'],0)
        # tracer cells:
        c['RAC'] = np.concatenate(a['RAC'],b['RAC']      ,0)
        c['DXG'] = np.concatenate(a['DXG'],b['DXG']      ,0)
        c['DYG'] = np.concatenate(a['DYG'],b['DYG'][1:,:],0)
        # vorticity cells:
        c['RAZ'] = np.concatenate(a['RAZ'],b['RAZ'][1:,:],0)
        c['DXC'] = np.concatenate(a['DXC'],b['DXC'][1:,:],0)
        c['DYC'] = np.concatenate(a['DYC'],b['DYC']      ,0)
        # "u" (western) cells:
        c['RAW'] = np.concatenate(a['RAW'],b['RAW'][1:,:],0)
        c['DXV'] = np.concatenate(a['DXV'],b['DXV'][1:,:],0)
        c['DYF'] = np.concatenate(a['DYF'],b['DYF']      ,0)
        # "v" (southern) cells:
        c['RAS'] = np.concatenate(a['RAS'],b['RAS']      ,0)
        c['DXF'] = np.concatenate(a['DXF'],b['DXF']      ,0)
        c['DYU'] = np.concatenate(a['DYU'],b['DYU'][1:,:],0)

    elif edges==(N,S):
        # north edge of 'a' and south edge of 'b':
        # grids:
        c['XG'] = np.concatenate(a['XG'],b['XG'][:,1:],1)
        c['YG'] = np.concatenate(a['YG'],b['YG'][:,1:],1)
        # grid cell centers:
        c['XC'] = np.concatenate(a['XC'],b['XC'],1)
        c['YC'] = np.concatenate(a['YC'],b['YC'],1)
        # tracer cells:
        c['RAC'] = np.concatenate(a['RAC'],b['RAC']      ,1)
        c['DXG'] = np.concatenate(a['DXG'],b['DXG'][:,1:],1)
        c['DYG'] = np.concatenate(a['DYG'],b['DYG']      ,1)
        # vorticity cells:
        c['RAZ'] = np.concatenate(a['RAZ'],b['RAZ'][:,1:],1)
        c['DXC'] = np.concatenate(a['DXC'],b['DXC']      ,1)
        c['DYC'] = np.concatenate(a['DYC'],b['DYC'][:,1:],1)
        # "u" (western) cells:
        c['RAW'] = np.concatenate(a['RAW'],b['RAW']      ,1)
        c['DXV'] = np.concatenate(a['DXV'],b['DXV'][:,1:],1)
        c['DYF'] = np.concatenate(a['DYF'],b['DYF']      ,1)
        # "v" (southern) cells:
        c['RAS'] = np.concatenate(a['RAS'],b['RAS'][:,1:],1)
        c['DXF'] = np.concatenate(a['DXF'],b['DXF']      ,1)
        c['DYU'] = np.concatenate(a['DYU'],b['DYU'][:,1:],1)

    elif edges==(W,E):
        # west edge of 'a' and east edge of 'b':
        # grids:
        c['XG'] = np.concatenate(b['XG'][:-1,:],a['XG'],0)
        c['YG'] = np.concatenate(b['YG'][:-1,:],a['YG'],0)
        # grid cell centers:
        c['XC'] = np.concatenate(b['XC'],a['XC'],0)
        c['YC'] = np.concatenate(b['YC'],a['YC'],0)
        # tracer cells:
        c['RAC'] = np.concatenate(b['RAC']       ,a['RAC'],0)
        c['DXG'] = np.concatenate(b['DXG']       ,a['DXG'],0)
        c['DYG'] = np.concatenate(b['DYG'][:-1,:],a['DYG'],0)
        # vorticity cells:
        c['RAZ'] = np.concatenate(b['RAZ'][:-1,:],a['RAZ'],0)
        c['DXC'] = np.concatenate(b['DXC'][:-1,:],a['DXC'],0)
        c['DYC'] = np.concatenate(b['DYC']       ,a['DYC'],0)
        # "u" (western) cells:
        c['RAW'] = np.concatenate(b['RAW'][:-1,:],a['RAW'],0)
        c['DXV'] = np.concatenate(b['DXV'][:-1,:],a['DXV'],0)
        c['DYF'] = np.concatenate(b['DYF']       ,a['DYF'],0)
        # "v" (southern) cells:
        c['RAS'] = np.concatenate(b['RAS']       ,a['RAS'],0)
        c['DXF'] = np.concatenate(b['DXF']       ,a['DXF'],0)
        c['DYU'] = np.concatenate(b['DYU'][:-1,:],a['DYU'],0)

    elif edges==(S,N):
        # south edge of 'a' and north edge of 'b':
        # grids:
        c['XG'] = np.concatenate(b['XG'][:,:-1],a['XG'],1)
        c['YG'] = np.concatenate(b['YG'][:,:-1],a['YG'],1)
        # grid cell centers:
        c['XC'] = np.concatenate(b['XC'],a['XC'],1)
        c['YC'] = np.concatenate(b['YC'],a['YC'],1)
        # tracer cells:
        c['RAC'] = np.concatenate(b['RAC']       ,a['RAC'],1)
        c['DXG'] = np.concatenate(b['DXG'][:,:-1],a['DXG'],1)
        c['DYG'] = np.concatenate(b['DYG']       ,a['DYG'],1)
        # vorticity cells:
        c['RAZ'] = np.concatenate(b['RAZ'][:,:-1],a['RAZ'],1)
        c['DXC'] = np.concatenate(b['DXC']       ,a['DXC'],1)
        c['DYC'] = np.concatenate(b['DYC'][:,:-1],a['DYC'],1)
        # "u" (western) cells:
        c['RAW'] = np.concatenate(b['RAW']       ,a['RAW'],1)
        c['DXV'] = np.concatenate(b['DXV'][:,:-1],a['DXV'],1)
        c['DYF'] = np.concatenate(b['DYF']       ,a['DYF'],1)
        # "v" (southern) cells:
        c['RAS'] = np.concatenate(b['RAS'][:,:-1],a['RAS'],1)
        c['DXF'] = np.concatenate(b['DXF']       ,a['DXF'],1)
        c['DYU'] = np.concatenate(b['DYU'][:,:-1],a['DYU'],1)

    return (c,)+c['XC'].shape


def main():
    """Command-line entry point."""

    parser = create_parser()
    args = parser.parse_args()

    if args.verbose:
        print('joining {0:s} and {1:s} along a common edge...'.format(
            args.tilea,args.tileb))

    (c,nic,njc) = stitch(
        args.tilea, args.nia, args.nja,
        args.tileb, args.nib, args.njb,
        args.verbose)

    if args.verbose:
        print(
            'writing combined tile to {0:s} with ni={1:d}, nj={2:d}...'.format(
            args.outfile,nic,njc))

    gridio.write_mitgridfile(args.outfile,c,nic,njc)

    if args.verbose:
        print('...done.')

if __name__ == '__main__':
    main()

