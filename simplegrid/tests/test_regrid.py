
import unittest
import simplegrid as sg
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


if __name__=='__main__':
    unittest.main()

