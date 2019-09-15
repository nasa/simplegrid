
import numpy as np
import os

from . import mds  # (*)
from . import mitgridfilefields as mgf

# (*) MITgcm utility; ref: MITgcm/utils/python/MITgcmutils/MITgcmutils/mds.py


def mds2mitgrid(rundir='./',verbose=False):
    """Recasts MITgcm meta-data files as mitgrid-formatted arrays

    Args:
        rundir (str): MITgcm run directory path
        verbose (bool): True for diagnostic output, False otherwise

    Returns:
        (mitgrid, ni, nj): Tuple consisting of mitgrid dictionary of named numpy
            arrays (ref. mitgridfilefields.py) and corresponding tracer point
            counts in the model grid 'x' and 'y' directions.

    """

    mitgrid = {key:None for key in mgf.names}
    ni = nj = 0

    for (name,ni_del,nj_del) in zip(mgf.names,mgf.ni_delta_sizes,mgf.nj_delta_sizes):

        # MITgcm may not have generated all standard grid files; just ignore
        # those that aren't present:

        try:
            if verbose:
                print('reading {0:>3s}...'.format(name),end='')
            tmp = mds.rdmds(rundir+os.sep+name).T
        except:
            if verbose:
                print('............no data')
            continue

        if verbose:
            print('reformatting...',end='')
        mitgrid[name] = np.zeros((tmp.shape[0]+ni_del,tmp.shape[1]+nj_del))
        mitgrid[name][0:tmp.shape[0],0:tmp.shape[1]] = tmp
        if not ni and not nj:
            (ni,nj) = tmp.shape
        if verbose:
            print('done')

    return (mitgrid,ni,nj)

