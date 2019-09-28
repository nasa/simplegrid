
import glob
import numpy as np
import os.path
import re
import xesmf as xe

from . import gridio
from . import mds
from . import mds2mitgrid
from . import regrid

INTERP_METH = 'bilinear'
NETCDF_EXT = '.nc'


def weightfilename( interp_meth, grid_in, grid_out, response_type, netcdf_ext):
    """Create a descriptive regridder filename consistent with xesmf
    conventions, augmented by response type.

    Args:
        interp_meth (str): interpolation method, typically one of those
            supported by xESMF (though it can be any string): 'bilinear',
            'conservative', 'nearest_s2d', 'nearest_d2s', or 'patch'.
        grid_in (xESMF grid; numpy 2-d arrays stored in dict with 'lon', 'lat'
            keys): input (global) grid, used simply to retrieve input grid
            sizes.
        grid_out (xESMF grid; see grid_in): output (regional) grid, used simply
            to retrieve output grid sizes.
        response_type (str): type of response quantity for which the resulting
            mapping will be used, e.g., 'tracer', 'U', or 'V'.
        netcdf_ext (str): NetCDF file extension (typically '.nc')

    Returns:
        Regridder filename string, e.g., 'bilinear_20x16_40x32_tracer.nc'

    """

    return '{0}_{1}x{2}_{3}x{4}_{5}.nc'.format(
        interp_meth,
        str(grid_in['lon'].shape[0]), str(grid_in['lon'].shape[1]),   # either
        str(grid_out['lon'].shape[0]), str(grid_out['lon'].shape[1]), # lon or lat
        response_type, netcdf_ext)


def getobcs( strict=False, verbose=False, **kwargs):
    """Compute prescribed regional grid open boundary conditions from global
    (parent) simulation results.

    Args:
        strict (bool, optional): True to raise error if nonzero terms found in
            any of the standard mitgrid matrix row/column padding dimensions
            (e.g., last row and column of XG, YG, etc.), False to ignore.)
        verbose (bool, optional): True for diagnostic output, False otherwise.
        **kwargs: Additional keyword arguments.

    Kwargs:
        parent_mitgridfile (str, required if no parent_mitgrid_matrices):
            (path and) filename of global simulation grid data.
        ni_parent (int, required if parent_mitgridfile): number of tracer points
            in the parent model grid 'x' direction.
        nj_parent (int, required if parent_mitgridfile): number of tracer points
            in the parent model grid 'y' direction.
        parent_mitgrid_matrices (dict, required if no parent_mitgridfile):
            global simulation grid data as name/value (numpy 2-d array) pairs
            corresponding to matrix name and ordering convention listed in
            mitgridfilefields module.
        regional_mitgridfile (str, required if no regional_mitgrid_matrices):
            (path and) filename of regional grid data.
        regional_mitgrid_matrices (dict, required if no regional_mitgridfile):
            regional grid data as name/value (numpy 2-d array) pairs
            corresponding to matrix name and ordering convention listed in
            mitgridfilefields module.
        ni_regional (int, required if regional_mitgridfile): number of tracer
            points in the regional model grid 'x' direction.
        nj_regional (int, required if regional_mitgridfile): number of tracer
            points in the regional model grid 'y' direction.
        parent_resultsdir (str, required): MITgcm parent grid results directory
            (e.g., './run').
        OB_Jnorth, OB_Jsouth, OB_Ieast, OB_Iwest (ones-based int vectors):
            North, south, east, and west Nx and Ny index vectors of i and j
            tracer cell indices, per MITgcm OBCS package conventions. Boundaries
            are ignored for any index vector that is identically zero, and
            defaults to edge tracer cell indices if not provided.
        resultsdir (str): Directory in which to store resulting open boundary
            condition arrays.  Defaults to parent_resultsdir+'_obcs'.

    Returns:

    Raises:
        ValueError: If missing parent or regional grid input,

    Note:
        Resulting data are intended for use with MITgcm data.pkg option,
        'useOBCSprescribe = .TRUE.'

    """

    # kwarg handling:
    parent_mitgridfile      = kwargs.get('parent_mitgridfile')
    parent_mitgrid_matrices = kwargs.get('parent_mitgrid_matrices')
    ni_parent               = kwargs.get('ni_parent')
    nj_parent               = kwargs.get('nj_parent')
    regional_mitgridfile    = kwargs.get('regional_mitgridfile')
    regional_mitgrid_matrices = kwargs.get('regional_mitgrid_matrices')
    ni_regional             = kwargs.get('ni_regional')
    nj_regional             = kwargs.get('nj_regional')
    parent_resultsdir       = kwargs.get('parent_resultsdir')
    ob_jnorth               = np.array(kwargs.get('OB_Jnorth'))
    ob_jsouth               = np.array(kwargs.get('OB_Jsouth'))
    ob_ieast                = np.array(kwargs.get('OB_Ieast'))
    ob_iwest                = np.array(kwargs.get('OB_Iwest'))

    #
    # locate global (parent) grid:
    #

    if not parent_resultsdir:
        raise ValueError(
            "Global solution results directory (parent_resultsdir) must be specified.")
    if verbose:
        readmsg = 'reading parent grid from'
    # get parent grid from any one of several possible input sources:
    if parent_mitgridfile:
        if ni_parent and nj_parent:
            if verbose:
                print('{0} {1}...'.format(readmsg,parent_mitgridfile))
            parent_mitgrid = gridio.read_mitgridfile(
                parent_mitgridfile, ni_parent, nj_parent, strict, verbose)
        else:
            raise ValueError(
                "ni_parent and nj_parent required if parent_mitgridfile specified.")
    elif parent_mitgrid_matrices:
        if verbose:
            print('{0} {1}...'.format(readmsg,'parent_mitgrid_matrices'))
        parent_mitgrid = parent_mitgrid_matrices
        # ni, nj from tracer cell counts:
        (ni_parent, nj_parent) = parent_mitgrid['XC'].shape
    else:
        if verbose:
            print('{0} {1}...'.format(readmsg,parent_resultsdir))
        (parent_mitgrid, ni_parent, nj_parent) = mds2mitgrid.mds2mitgrid(
            parent_resultsdir,verbose)
    if verbose:
        print('...successfully read {0}x{1} grid.'.format(ni_parent,nj_parent))

    #
    # locate regional grid:
    #

    if verbose:
        readmsg = 'reading regional grid from'
    if regional_mitgridfile:
        if ni_regional and nj_regional:
            if verbose:
                print('{0} {1}...'.format(readmsg,regional_mitgridfile))
            regional_mitgrid = gridio.read_mitgridfile(
                regional_mitgridfile, ni_regional, nj_regional, strict, verbose)
        else:
            raise ValueError(
                "ni_regional and nj_regional required if regional_mitgridfile specified.")
    elif regional_mitgrid_matrices:
        if verbose:
            print('{0} {1}...'.format(readmsg,'regional_mitgrid_matrices'))
        regional_mitgrid = regional_mitgrid_matrices
        # ni, nj from tracer cell counts:
        (ni_regional, nj_regional) = regional_mitgrid['XC'].shape
    else:
        raise ValueError(
            "Either regional_mitgridfile or regional_mitgrid_matrices must be provided.")
    if verbose:
        print('...successfully read {0}x{1} grid.'.format(ni_regional,nj_regional))

    #
    # boundary index vectors:
    #

    # accommodate the following user controls:
    # - if an open boundary vector isn't defined (blank), assume edge tracer cells
    # - an explicit zero vector (any length) implies no open boundary on that edge
    # - negative indices count backwards from Nx+1, Ny+1; e.g., as documented in
    #   MITgcm manual section 8.3.1.4 OB_Jnorth(3)=-1 means that the point (3,Ny)
    #   is a northern OB

    # defaults:

    do_ob_jnorth, do_ob_jsouth, do_ob_ieast, do_ob_iwest = [True]*4
    if not ob_jnorth:
        ob_jnorth = np.ones(parent_mitgrid['XC'].shape[0],dtype=int)*nj_regional
    elif not any(ob_jnorth!=0):
        do_ob_jnorth = False
    if not ob_jsouth:
        ob_jsouth = np.ones(parent_mitgrid['XC'].shape[0],dtype=int)
    elif not any(ob_jsouth!=0):
        do_ob_jsouth = False
    if not ob_ieast:
        ob_ieast = np.ones(parent_mitgrid['XC'].shape[1],dtype=int)*ni_regional
    elif not any(ob_ieast!=0):
        do_ob_ieast = False
    if not ob_iwest:
        ob_iwest = np.ones(parent_mitgrid['XC'].shape[1],dtype=int)
    elif not any(ob_iwest!=0):
        do_ob_iwest = False

    # apply logic to negative indices and convert vectors to zeros-based:

    if do_ob_jnorth:
        ob_jnorth = np.array([nj_regional+1+j if j<0 else j for j in ob_jnorth]) - 1
    if do_ob_jsouth:
        ob_jsouth = ob_jsouth - 1
    if do_ob_ieast:
        ob_ieast = np.array([ni_regional+1+i if i<0 else i for i in ob_ieast]) - 1
    if do_ob_iwest:
        ob_iwest = ob_iwest - 1

    #
    # interpolators:
    #

    # tracer cell interpolator:
    ds_tracer_parent = {
        'lon':parent_mitgrid['XC'],
        'lat':parent_mitgrid['YC']}
    ds_tracer_regional = {
        'lon':regional_mitgrid['XC'],
        'lat':regional_mitgrid['YC']}
    wgtfname_tracer = weightfilename(
        INTERP_METH, ds_tracer_parent, ds_tracer_regional, 'tracer', NETCDF_EXT)
    tracer_regridder = xe.Regridder( ds_tracer_parent, ds_tracer_regional,
        INTERP_METH, filename=wgtfname_tracer, reuse_weights=True)

    # u, v point interplators - "regrid" both parent and region at 2x
    # resolution, use to extract tracer cell u,v midside nodes:
    (parent_mitgrid_double, ni_parent_double, nj_parent_double) = \
        regrid.regrid(
            mitgrid_matrices=parent_mitgrid, ni=ni_parent, nj=nj_parent,
            lon1=parent_mitgrid['XG'][ 0,-1],   # NW
            lat1=parent_mitgrid['YG'][ 0,-1],   # ""
            lon2=parent_mitgrid['XG'][-1, 0],   # SE
            lat2=parent_mitgrid['YG'][-1, 0],   # ""
            lon_subscale=2, lat_subscale=2)
    (regional_mitgrid_double, ni_regional_double, ni_regional_double) = \
        regrid.regrid(
            mitgrid_matrices=regional_mitgrid, ni=ni_regional, nj=nj_regional,
            lon1=regional_mitgrid['XG'][ 0,-1], # NW
            lat1=regional_mitgrid['YG'][ 0,-1], # ""
            lon2=regional_mitgrid['XG'][-1, 0], # SE
            lat2=regional_mitgrid['YG'][-1, 0], # ""
            lon_subscale=2, lat_subscale=2)
    # U-point interpolator:
    ds_U_parent = {
        'lon': parent_mitgrid_double['XG'][0::2,1::2],
        'lat': parent_mitgrid_double['YG'][0::2,1::2]}
    ds_U_regional = {
        'lon': regional_mitgrid_double['XG'][0::2,1::2],
        'lat': regional_mitgrid_double['YG'][0::2,1::2]}
    wgtfname_U = weightfilename( INTERP_METH, ds_U_parent,
        ds_U_regional, 'U', NETCDF_EXT)
    U_regridder = xe.Regridder( ds_U_parent, ds_U_regional, INTERP_METH,
        filename=wgtfname_U, reuse_weights=True)
    # V-point interpolator:
    ds_V_parent = {
        'lon': parent_mitgrid_double['XG'][1::2,0::2],
        'lat': parent_mitgrid_double['YG'][1::2,0::2]}
    ds_V_regional = {
        'lon': regional_mitgrid_double['XG'][1::2,0::2],
        'lat': regional_mitgrid_double['YG'][1::2,0::2]}
    wgtfname_V = weightfilename( INTERP_METH, ds_V_parent,
        ds_V_regional, 'V', NETCDF_EXT)
    V_regridder = xe.Regridder( ds_V_parent, ds_V_regional, INTERP_METH,
        filename=wgtfname_V, reuse_weights=True)

    #
    # get time steps and number of depths from results database:
    #

    # get integer time steps from filenames:
    re_time = re.compile('\d+')
    times = sorted( [int(re_time.findall(os.path.basename(f))[0])
        for f in glob.glob(os.path.join(parent_resultsdir,'T.*.data'))])
    print('times = {0}'.format(times))
    ntimes = len(times)
    # open up a temperature file, get depths from matrix size:
    T1 = mds.rdmds(os.path.join(parent_resultsdir,'T'),times[0])
    ndepths = T1.shape[0]

    print('ntimes/ndepths = {0}/{1}'.format(ntimes,ndepths))

    #
    # put resulting controls into dictionary to help organize subsequent loops
    # ('None's are data that will be filled in at loop level; included here only
    # for clarity):
    #

    obs = { 'N': {
                'do'        : do_ob_jnorth,
                'partition' : (np.arange(len(ob_jnorth)),ob_jnorth),# "fancy" indexing
                'shape'     : [ni_regional,None,ntimes],            # Nones will be ndepths
                'ob_array'  : None},
            'S': {
                'do'        : do_ob_jsouth,
                'partition' : (np.arange(len(ob_jsouth)),ob_jsouth),
                'shape'     : [ni_regional,None,ntimes],
                'ob_array'  : None},
            'E': {
                'do'        : do_ob_ieast,
                'partition' : (ob_ieast,np.arange(len(ob_ieast))),
                'shape'     : [nj_regional,None,ntimes],
                'ob_array'  : None},
            'W': {
                'do'        : do_ob_iwest,
                'partition' : (ob_iwest,np.arange(len(ob_iwest))),
                'shape'     : [nj_regional,None,ntimes],
                'ob_array'  : None},}

    #
    # Apply interpolators to map from parent->regional, and partition solutions
    # accordingly:
    #

    # tracer cell-related responses:
    responses = ['T','S','Eta']
    for response in responses:
        if glob.glob(os.path.join(parent_resultsdir,response+'.*.data')):
            if verbose:
                print("computing obcs for '{0}'...".format(response))
            # if results exist in database, then for this reponse, and for
            # selected boundaries, assemble boundary matrices for all depths,
            # times:
            for time_idx,time in enumerate(times):
                if time_idx==0:
                    # allocate storage on first pass:
                    # get ndepths from first instance:
                    try:
                        ndepths = mds.rdmds(os.path.join(
                            parent_resultsdir,response),times[time_idx]).T.shape[2]
                    except IndexError:
                        ndepths = 1
                    for ob in obs:
                        if obs[ob]['do']:
                            obs[ob]['shape'][1] = ndepths
                            obs[ob]['ob_array'] = np.zeros(obs[ob]['shape'])
                # get parent results for this response, for this time...:
                parent_results = mds.rdmds(os.path.join(parent_resultsdir,response),time).T
                # for each depth...:
                for depth_idx in range(ndepths):
                    # ... interpolate to the regional grid...:
                    try:
                        regional_results = tracer_regridder(parent_results[:,:,depth_idx])
                    except IndexError:
                        # will be the case if only one depth:
                        parent_results = np.expand_dims(parent_results,2)
                        regional_results = tracer_regridder(parent_results[:,:,depth_idx])
                    # ... and insert the boundary partitions for the requested
                    # N, S, E, W boundary matrices for this depth, time:
                    for ob in obs:
                        if obs[ob]['do']:
                            obs[ob]['ob_array'][:,depth_idx,time_idx] = \
                                regional_results[obs[ob]['partition']]
                        if depth_idx==0 and time_idx==0:
                            print('{0} compare:'.format(ob))
                            if 'N'==ob:
                                print(regional_results[:,-1])
                            elif 'S'==ob:
                                print(regional_results[:,0])
                            elif 'E'==ob:
                                print(regional_results[-1,:])
                            elif 'W'==ob:
                                print(regional_results[0,:])
                            print(obs[ob]['ob_array'][:,depth_idx,time_idx])


    # failed experiment: keep, for now. attempted to interpolate just to string
    # of points as opposed to a grid.
    #if ob_iwest.any():
    #    # tracer cell interpolator:
    #    tracer_region_west_xc = np.zeros(ob_iwest.shape)
    #    tracer_region_west_yc = np.zeros(ob_iwest.shape)
    #    it = np.nditer(ob_iwest,flags=['multi_index'])
    #    while not it.finished:
    #        row_idx = it[0] # ones-based, per MITgcm conventions
    #        col_idx = it.multi_index[0]
    #        if row_idx>0:
    #            tracer_region_west_xc[col_idx] = \
    #                regional_mitgrid['XC'][row_idx-1,col_idx]
    #            tracer_region_west_yc[col_idx] = \
    #                regional_mitgrid['YC'][row_idx-1,col_idx]
    #        it.iternext()
    #    ds_tracer_region_west = {
    #        'lon': tracer_region_west_xc, 'lat': tracer_region_west_yc}
    #    print('tracer_region_west_xc.shape = {0}'.format(tracer_region_west_xc.shape))
    #    print('tracer_region_west_yc.shape = {0}'.format(tracer_region_west_yc.shape))
    #    tracer_region_west_regridder = xe.Regridder(
    #        ds_tracer_parent, ds_tracer_region_west, 'bilinear')

#T, S, U, V, Eta, W, UICE, VICE
#
# runtime parameters:
#OB[N/S/E/W][t/s/u/v]File
#if present,
#OB[N/S/E/W]etaFile
#OB[N/S/E/W]wFile
#OB[N/S/E/W][a,h,sl,sn,uice,vice]


