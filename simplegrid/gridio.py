
import numpy as np
import struct
import sys
from . import mitgridfilefields as mgf

def read_mitgridfile(filename,ni,nj,verbose=False):
    """Read a serial (plain format) grid definition file.

    *.mitgrid files consist of contiguous binary segments, each of nominal size
    (ni+1)*(nj+1) terms. Since file data structure and sizes are implied,
    expected grid cell counts, ni and nj, must be provided on input.

    Args:
        filename (str): mitgrid (path and) file name
        ni (int): number of expected nominal "east-west" grid cells
        nj (int): number of expected nominal "north-south" grid cells
        verbose (bool): progress reporting to stdio

    Returns:
        mitgrid_matrices (dict): name/value (numpy 2-d array) pairs
            corresponding to matrix name and ordering convention listed in
            mitgridfilefields module

    Raises:
        RuntimeError: If any matrix resize operation attempts to delete nonzero
            row/column data (indicating improper input assumptions)

    Comments:
        - Nominal "north-south"/"east-west" directions depend on the particular
            gridfile orientation
        - Note that ni and nj map to output matrix rows and columns,
            respectively.
        - The data stream in mitgrid files is assumed to read columwise (e.g. in
            Fortran order) into the respective grid descriptor matrices.

    Examples:
        >>> # llc 90 gridfile:
        >>> mitgrid_matrices = read_mitgridfile( './tests/tile005.mitgrid', 270, 90, True)
        reading 394576 doubles from tile005.mitgrid...  done.
        formatting XC  (270x90)...  done.
        formatting YC  (270x90)...  done.
        formatting DXF (270x90)...  done.
        formatting DYF (270x90)...  done.
        formatting RAC (270x90)...  done.
        formatting XG  (271x91)...  done.
        formatting YG  (271x91)...  done.
        formatting DXV (271x91)...  done.
        formatting DYU (271x91)...  done.
        formatting RAZ (271x91)...  done.
        formatting DXC (271x90)...  done.
        formatting DYC (270x91)...  done.
        formatting RAW (271x90)...  done.
        formatting RAS (270x91)...  done.
        formatting DXG (270x91)...  done.
        formatting DYG (271x90)...  done.
        >>> list(mitgrid_matrices)
        ['XC', 'YC', 'DXF', 'DYF', 'RAC', 'XG', 'YG', 'DXV', 'DYU', 'RAZ', 'DXC', 'DYC', 'RAW', 'RAS', 'DXG', 'DYG']
        >>> mitgrid_matrices['XC'] 
        array([[-127.83792318, -127.44420903, -126.96814311, ...,  -39.03185689,
                 -38.55579097,  -38.16207682],
               [-127.77198934, -127.29194279, -126.76402521, ...,  -39.23597479,
                 -38.70805721,  -38.22801066],
               [-127.73007106, -127.17828952, -126.59756402, ...,  -39.40243598,
                 -38.82171048,  -38.26992894],
               ...,
               [-115.85020453, -115.78025365, -115.71079554, ..., -101.42989724,
                -100.4884439 ,  -99.42047829],
               [-115.50567133, -115.46406357, -115.4227498 , ..., -106.83080944,
                -106.24874177, -105.58464976],
               [-115.16698184, -115.15324165, -115.13959868, ..., -112.28604988,
                -112.0900639 , -111.86579135]])

    """
    # read raw data:
    slice_size = (ni+1)*(nj+1)
    dt = np.dtype([(name,mgf.datatype,slice_size) for name in mgf.names])
    if verbose:
        sys.stdout.write('reading {0:d} doubles from {1:s}... '.format(
            len(mgf.names)*slice_size,filename))
    rawdata = np.fromfile(filename,dt)  # scalar array with named fields
    if verbose:
        sys.stdout.write(' done.\n')

    # and store shaped, resized arrays by name/value:
    mitgrid_matrices = dict()
    for (name,ni_del,nj_del) in zip(mgf.names,mgf.ni_delta_sizes,mgf.nj_delta_sizes):
        # reshape, resize separately to allow intermediate error checking:
        mitgrid_matrices[name] = np.reshape(rawdata[0][name],(ni+1,nj+1),order='F')
        if  (not ni_del and any(mitgrid_matrices[name][ni,:])) \
            or \
            (not nj_del and any(mitgrid_matrices[name][:,nj])):
            raise RuntimeError(
                'trying to trim nonzero rows or columns from {0}'.format(name))
        if verbose:
            sys.stdout.write('formatting {0:3s} ({1:d}x{2:d})... '.format(name,ni+ni_del,nj+nj_del))
        mitgrid_matrices[name] = mitgrid_matrices[name][:ni+ni_del,:nj+nj_del]
        if verbose:
            sys.stdout.write(' done.\n')
    return mitgrid_matrices


def write_mitgridfile(filename,griddata,ni,nj,verbose=False):
    """Write a serial (plain format) grid definition file that can be read with
    read_mitgridfile.

    Args:
        filename (str): mitgrid (path and) file to write
        griddata: dictionary of numpy array data to be written (ref.
            mitgridfilefields.py)
        ni (int): number of nominal "east-west" grid cells in the collection of
            griddata matrices
        nj (int): number of nominal "north-south" grid cells in the collection
            of griddata matrices
    Returns:
        bool: True for success, False otherwise.

    """

    fd = open(filename,'wb')
    slice_size = (ni+1)*(nj+1)

    for (name,ni_del,nj_del) in zip(mgf.names,mgf.ni_delta_sizes,mgf.nj_delta_sizes):
        # copy data to standard (ni+1,nj+1)-sized array:
        tmparray = np.zeros((ni+1,nj+1))
        tmparray[:ni+ni_del,:nj+nj_del] = griddata[name]
        # since np.tofile writes in C order, explicity recast array to vector
        # using fortran ordering:
        tmparray = np.reshape(tmparray,slice_size,order='F')
        # finally, ensure double-precision, big-endian representation:
        outarray = tmparray.astype('>f8')
        # append to output file:
        outarray.tofile(fd)

    fd.close
    return True

