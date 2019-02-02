
import unittest
import simplegrid as sg
import numpy as np
import numpy.testing as nptest

class TestRegrid(unittest.TestCase):

    def test_regrid_1(self):
        """ Tests 1x1 'remeshing' of a 1x1 cell from an llc 90 model tile.
        Results are compared against validated mitgrid file created using
        sgregrid command-line call:

        sgregrid                            \\
            --mitgridfile tile005.mitgrid   \\
            --ni 270                        \\
            --nj  90                        \\
            --lon1 -126.95921008            \\
            --lat1   67.31607921            \\
            --lon2 -126.21154908            \\
            --lat2   67.17736362            \\
            --lon_subscale 1                \\
            --lat_subscale 1                \\
            --outfile regrid_test_1.mitgrid
        """

        (newgrid,newgrid_ni,newgrid_nj) = sg.regrid.regrid(
            mitgridfile='./data/tile005.mitgrid',
            ni=270, nj=90,
            lon1=-126.959, lat1=67.316, lon2=-126.212, lat2=67.177,
            lon_subscale=1, lat_subscale=1)

        validated_grid = sg.gridio.read_mitgridfile(
            './data/regrid_test_1.mitgrid',1,1)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(newgrid,validated_grid):
            nptest.assert_almost_equal(newgrid[a],validated_grid[b])


    def test_regrid_2(self):
        """ Tests 30x20 remeshing of a 1x1 cell from an llc 90 model tile.
        Results are compared against validated mitgrid file created using
        sgregrid command-line call:

        sgregrid                            \\
            --mitgridfile tile005.mitgrid   \\
            --ni 270                        \\
            --nj  90                        \\
            --lon1 -126.95921008            \\
            --lat1   67.31607921            \\
            --lon2 -126.21154908            \\
            --lat2   67.17736362            \\
            --lon_subscale 30               \\
            --lat_subscale 20               \\
            --outfile regrid_test_2.mitgrid
        """

        (newgrid,newgrid_ni,newgrid_nj) = sg.regrid.regrid(
            mitgridfile='./data/tile005.mitgrid',
            ni=270, nj=90,
            lon1=-126.95921008, lat1=67.31607921,
            lon2=-126.21154908, lat2=67.17736362,
            lon_subscale=30, lat_subscale=20)

        validated_grid = sg.gridio.read_mitgridfile(
            './data/regrid_test_2.mitgrid',30,20)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(newgrid,validated_grid):
            nptest.assert_almost_equal(newgrid[a],validated_grid[b])


    def test_regrid_3(self):
        """ Tests 1x1 'remeshing' of a 1x1 corner cell from an llc 90 model
        tile.  Results are compared against validated mitgrid file created using
        sgregrid command-line call:

        sgregrid                            \\
            --mitgridfile tile005.mitgrid   \\
            --ni 270                        \\
            --nj  90                        \\
            --lon1 -127.73445435            \\
            --lat1   67.56064719            \\
            --lon2 -128.0                   \\
            --lat2   67.40168504            \\
            --lon_subscale 1                \\
            --lat_subscale 1                \\
            --outfile regrid_test_3.mitgrid
        """

        (newgrid,newgrid_ni,newgrid_nj) = sg.regrid.regrid(
            mitgridfile='./data/tile005.mitgrid',
            ni=270, nj=90,
            lon1=-127.73445435, lat1=67.56064719,
            lon2=-128., lat2=67.40168504,
            lon_subscale=1, lat_subscale=1)

        validated_grid = sg.gridio.read_mitgridfile(
            './data/regrid_test_3.mitgrid',1,1)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(newgrid,validated_grid):
            nptest.assert_almost_equal(newgrid[a],validated_grid[b])


    def test_regrid_4(self):
        """ Tests 10x10 remeshing of a 1x1 corner cell from an llc 90 model
        tile.  Results are compared against validated mitgrid file created using
        sgregrid command-line call:

        sgregrid                            \\
            --mitgridfile tile005.mitgrid   \\
            --ni 270                        \\
            --nj  90                        \\
            --lon1 -127.73445435            \\
            --lat1   67.56064719            \\
            --lon2 -128.0                   \\
            --lat2   67.40168504            \\
            --lon_subscale 10               \\
            --lat_subscale 10               \\
            --outfile regrid_test_4.mitgrid
        """

        (newgrid,newgrid_ni,newgrid_nj) = sg.regrid.regrid(
            mitgridfile='./data/tile005.mitgrid',
            ni=270, nj=90,
            lon1=-127.73445435, lat1=67.56064719,
            lon2=-128., lat2=67.40168504,
            lon_subscale=10, lat_subscale=10)

        validated_grid = sg.gridio.read_mitgridfile(
            './data/regrid_test_4.mitgrid',10,10)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(newgrid,validated_grid):
            nptest.assert_almost_equal(newgrid[a],validated_grid[b])


    def test_regrid_5(self):
        """ Tests 1x1 'remeshing' of a 1x1 cell from an llc 90 model tile.
        Similar to test_regrid_1, but uses --xg_file, --yg_file option with
        binary input.  Results are compared against validated mitgrid file
        created using sgregrid command-line call:

        sgregrid                            \\
            --xg_file tile005_XG.bin        \\
            --yg_file tile005_YG.bin        \\
            --ni 270                        \\
            --nj  90                        \\
            --lon1 -126.95921008            \\
            --lat1   67.31607921            \\
            --lon2 -126.21154908            \\
            --lat2   67.17736362            \\
            --lon_subscale 1                \\
            --lat_subscale 1                \\
            --outfile regrid_test_5.mitgrid
        """

        (newgrid,newgrid_ni,newgrid_nj) = sg.regrid.regrid(
            xg_file='./data/tile005_XG.bin', yg_file='./data/tile005_YG.bin',
            ni=270, nj=90,
            lon1=-126.959, lat1=67.316, lon2=-126.212, lat2=67.177,
            lon_subscale=1, lat_subscale=1)

        validated_grid = sg.gridio.read_mitgridfile(
            './data/regrid_test_5.mitgrid',1,1)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(newgrid,validated_grid):
            nptest.assert_almost_equal(newgrid[a],validated_grid[b])


    def test_regrid_6(self):
        """ Tests 1x1 'remeshing' of a 1x1 cell from an llc 90 model tile.
        Similar to test_regrid_5, but uses --xg_file, --yg_file option with
        *.csv input.  Results are compared against validated mitgrid file
        created using sgregrid command-line call:

        sgregrid                            \\
            --xg_file tile005_XG.csv        \\
            --yg_file tile005_YG.csv        \\
            --ni 270                        \\
            --nj  90                        \\
            --lon1 -126.95921008            \\
            --lat1   67.31607921            \\
            --lon2 -126.21154908            \\
            --lat2   67.17736362            \\
            --lon_subscale 1                \\
            --lat_subscale 1                \\
            --outfile regrid_test_6.mitgrid
        """

        (newgrid,newgrid_ni,newgrid_nj) = sg.regrid.regrid(
            xg_file='./data/tile005_XG.csv', yg_file='./data/tile005_YG.csv',
            ni=270, nj=90,
            lon1=-126.959, lat1=67.316, lon2=-126.212, lat2=67.177,
            lon_subscale=1, lat_subscale=1)

        validated_grid = sg.gridio.read_mitgridfile(
            './data/regrid_test_6.mitgrid',1,1)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(newgrid,validated_grid):
            nptest.assert_almost_equal(newgrid[a],validated_grid[b])

    def test_regrid_7(self):
        """Tests 1x1 'remeshing' of llc 90 tile 005 and compares subset of differences."""

        tile = './data/tile005.mitgrid'
        # mask off the last rows to get a more valid comparison (by inspection,
        # differences are concentrated in this distorted range), and establish
        # (again, by inspection), an allowable percentage difference:
        filter = 10
        maxdiff = 1.5

        # load tile of interest so we can get corner points:
        mg = sg.gridio.read_mitgridfile(tile,270,90)

        # 1x1 "regrid":
        (mg_new,mg_new_ni,mg_new_nj) = sg.regrid.regrid(
            mitgridfile=tile,
            ni=270, nj=90,
            lon1=mg['XG'][0,-1], lat1=mg['YG'][0,-1],   # "NW"
            lon2=mg['XG'][-1,0], lat2=mg['YG'][-1,0],   # "SE"
            lon_subscale=1, lat_subscale=1)

        # compare:
        for name in sg.mitgridfilefields.names:
            mg_new[name][-filter:,:] = 1.
            mg[name][-filter:,:] = 1.
            diffs = (np.nan_to_num(mg_new[name])-np.nan_to_num(mg[name]))/np.nan_to_num(mg[name])
            nptest.assert_array_less(
                np.amax(diffs), maxdiff/100.,
                "array['{0}'] max diff is not less than {1}%".format(name,maxdiff/100.))


if __name__=='__main__':
    unittest.main()

