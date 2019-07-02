
==========
Simplegrid
==========


Simple local grid creation and refinement for ocean circulation models
----------------------------------------------------------------------


simplegrid is a python package for creating, refining, and joining horizontal
quadrilateral grids that are used in connection with the MIT General Circulation
Model (MITgcm).

simplegrid is based on equal great circle arc subdivision (hence the name) with
geodesic computations provided by pyproj/PROJ.4, and implements both
python-callable and command-line functionality for embedded and scripted
solutions. It also contains several useful utilities for grid manipulation, file
i/o, and coincident edge detection.

Grid point location data are specified in decimal longitude and latitude, with
resulting MITgcm grid data in units of meters.


Installation:
-------------

    Requirements:
    ~~~~~~~~~~~~~

        simplegrid is compatible with python 3+, with numpy and pyproj dependencies.

    Installation from github:
    ~~~~~~~~~~~~~~~~~~~~~~~~~

        simplegrid is under active development. To obtain the latest development
        version you may clone the repository and install it::

        git clone https://github.jpl.nasa.gov/gmoore/simplegrid.git
        cd simplegrid
        pip install .


Quick Start Examples:
---------------------

    Creating a grid:
    ~~~~~~~~~~~~~~~~

        The following creates a simple 10x10 grid over a one-by-one degree
        region with northwest/southeast lon/lat corners at (1.,2.)/(2.,1.):

        from python::

            import simplegrid as sg
            (newgrid,newgrid_ni,newgrid_nj) = sg.mkgrid.mkgrid(
                lon1=1., lat1=2.,
                lon2=2., lat2=1.,
                lon_subscale=10, lat_subscale=10)

        from the command line::

            sgmkgrid                            \
                --lon1 1.                       \
                --lat1 2.                       \
                --lon2 2.                       \
                --lat2 1.                       \
                --lon_subscale 10               \
                --lat_subscale 10               \
                --outfile newgrid.mitgrid

        All sixteen horizontal grid quantities (XC, YC, DXF, DYF, RAC, XG, YG,
        DXV, DYU, RAZ, DXC, DYC, RAW, RAS, DXG, and DYG) are either written to a
        dictionary with corresponding key/value pairs in the case of python or,
        in the command-line case, to a contiguous binary file in mitgridfile
        format.


    Refining a grid:
    ~~~~~~~~~~~~~~~~

        The following performs a 10x10 refinement of the llc 90 model, tile005,
        "southwest" corner cell. In the following, northwest/southeast input
        corner lat/lon values for the region of interest have been determined by
        inspection and need not be as precisely defined as shown; simplegrid
        will perform nearest-neighbor checks to determine the closest existing
        corner grids:

        from python::

            import simplegrid as sg
            (newgrid,newgrid_ni,newgrid_nj) = sg.regrid.regrid(
                strict=True, verbose=False,
                mitgridfile='./data/tile005.mitgrid',
                ni=270, nj=90,
                lon1=-127.73445435, lat1=67.56064719,
                lon2=-128., lat2=67.40168504,
                lon_subscale=10, lat_subscale=10)
                       
        from the command line::

            sgregrid                            \
                --mitgridfile tile005.mitgrid   \
                --ni 270                        \
                --nj  90                        \
                --lon1 -127.73445435            \
                --lat1   67.56064719            \
                --lon2 -128.0                   \
                --lat2   67.40168504            \
                --lon_subscale 10               \
                --lat_subscale 10               \
                --outfile regrid005.mitgrid     \
                --strict

        As in the preceding mkgrid case, all horizontal grid quantities are
        either written to a dictionary of name/value pairs in python or, in the
        command-line case, to a contiguous binary file in mitgrid file format.

        In addition to the mitgrid file input, the python and command line
        interfaces to regrid also support binary and comma-separated input
        options; such files would have been produced had an mitgrid file been
        read into matlab, for example, with XG and YG corner grid matrix output
        (the only mitgrid file quantities, in fact, used by regrid) to
        intermediate files:

        from python::

            import simplegrid as sg
            # both *.bin and *.csv supported:
            (newgrid,newgrid_ni,newgrid_nj) = sg.regrid.regrid(
                xg_file='./data/tile005_XG.bin',
                yg_file='./data/tile005_YG.bin',
                ni=270, nj=90,
                lon1=-127.73445435, lat1=67.56064719,
                lon2=-128., lat2=67.40168504,
                lon_subscale=10, lat_subscale=10)

        and, from the command line::

            # both *.bin and *.csv supported:
            sgregrid                            \
                --xg_file tile005_XG.csv        \
                --yg_file tile005_YG.csv        \
                --ni 270                        \
                --nj  90                        \
                --lon1 -127.73445435            \
                --lat1   67.56064719            \
                --lon2 -128.0                   \
                --lat2   67.40168504            \
                --lon_subscale 10               \
                --lat_subscale 10               \
                --outfile regrid005.mitgrid


    Determining boundary quantities:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        In most cases, mitgrid data that spans tracer cells is undefined along
        boundaries (for example, "U" cell quantities RAW and DXV along a tile's
        western and eastern edges).  "addfringe" functionality can compute this
        boundary, or fringe, data if an adjacent tile is provided.

        The following augments a 2x2 tile with boundary data computed using an
        eastern edge-adjacent 2x2 tile:

        from python::

            import simplegrid as sg
            (tilea_edge,tileb_edge,new_tilea_grid) = sg.addfringe.addfringe(
                strict=True,
                tilea='./data/tile_A_2x2.mitgrid',nia=2,nja=2,
                tileb='./data/tile_B_E_2x2.mitgrid',nib=2,njb=2)

        new_tilea_grid is a dictionary of mitgrid name/value pairs containing
        tile_A input data, augmented with eastern edge data computed using
        tile_B.  tilea_edge and tileb_edge are integer indicators confirming the
        A and B edge matches: 0==N, 1==S, 2==E, and 3==W (in this example,
        tilea_edge will be eqal to 2, and tileb_edge, 3).

        from the command line::

            sgaddfringe                         \
                --tilea tile_A_2x2.mitgrid      \
                --nia 2                         \
                --nja 2                         \
                --tileb tile_B_E_2x2.mitgrid    \
                --nib 2                         \
                --njb 2                         \
                --outfile addfringe_A_EW_2x2.mitgrid \
                --strict

        As in the python example, the output file contains tile_A grid
        quantities, augmented with eastern edge data computed using tile_B.
        Output is to a combined binary file in mitgrid file format.


    Joining grids:
    ~~~~~~~~~~~~~~

        Joining, or "stitching", two tiles together produces a single entity,
        assigning common-edge boundary quantities as appropriate.  The following
        joins two 2x2 tiles that match on their northern and southern edges,
        respectively, resulting in a 2x4 mitgrid:

        from python::

            (newgrid,newgrid_ni,newgrid_nj) = sg.stitch.stitch(
                strict=True, verbose=False,
                tilea='./data/tile_A_2x2.mitgrid',nia=2,nja=2,
                tileb='./data/tile_B_N_2x2.mitgrid',nib=2,njb=2)

        As in the previous examples, newgrid is a dictionary of mitgrid
        name/value pairs, and newgrid_ni and newgrid_nj provide the tracer cell
        row and column counts for the combined grid.

        from the command line::

            sgstitch                            \
                --tilea tile_A_2x2.mitgrid      \
                --nia 2                         \
                --nja 2                         \
                --tileb tile_B_N_2x2.mitgrid    \
                --nib 2                         \
                --njb 2                         \
                --outfile stitch_AB_NS_2x4.mitgrid \
                --strict


