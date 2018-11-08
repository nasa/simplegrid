
import unittest
import simplegrid as sg
import numpy.testing as nptest

class TestRegrid(unittest.TestCase):

    def test_regrid_1(self):
        """ Tests 1x1 'remeshing' of a 1x1 cell from an llc 90 model tile.
        Results are compared against validated mitgrid file created using
        simpleregrid command-line call:

        simpleregrid \
            ./data/tile005.mitgrid \
            270 \
            90 \
            -126.95921008 \
            67.31607921 \
            -126.21154908 \
            67.17736362 \
            1 \ 
            1 \ 
            regrid_test_1.mitgrid
        """

        (newgrid,newgrid_ni,newgrid_nj) = sg.regrid.regrid(
            './data/tile005.mitgrid', 270, 90, 
            -126.959, 67.316, -126.212, 67.177,
            1, 1)

        validated_grid = sg.gridio.read_mitgridfile('./data/regrid_test_1.mitgrid',1,1)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(newgrid,validated_grid):
            nptest.assert_array_equal(newgrid[a],validated_grid[b])


    def test_regrid_2(self):
        """ Tests 30x20 remeshing of a 1x1 cell from an llc 90 model tile.
        Results are compared against validated mitgrid file created using
        simpleregrid command-line call:

        simpleregrid \
            ./data/tile005.mitgrid \
            270 \
            90 \
            -126.95921008 \
            67.31607921 \
            -126.21154908 \
            67.17736362 \
            30 \ 
            20 \ 
            regrid_test_2.mitgrid
        """

        (newgrid,newgrid_ni,newgrid_nj) = sg.regrid.regrid(
            './data/tile005.mitgrid', 270, 90, 
            -126.959, 67.316, -126.212, 67.177,
            30, 20)

        validated_grid = sg.gridio.read_mitgridfile('./data/regrid_test_2.mitgrid',30,20)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(newgrid,validated_grid):
            nptest.assert_array_equal(newgrid[a],validated_grid[b])

    def test_regrid_3(self):
        """ Tests 1x1 'remeshing' of a 1x1 corner cell from an llc 90 model
        tile.  Results are compared against validated mitgrid file created using
        simpleregrid command-line call:

        simpleregrid \
            ./data/tile005.mitgrid \
            270 \
            90 \
            -128.0 \
            67.5 \
            -127.58893247 \
            67.42868828 \
            1 \ 
            1 \ 
            regrid_test_3.mitgrid
        """

        (newgrid,newgrid_ni,newgrid_nj) = sg.regrid.regrid(
            './data/tile005.mitgrid', 270, 90, 
            -128., 67.5, -127.58893247, 67.42868828,
            1, 1)

        validated_grid = sg.gridio.read_mitgridfile('./data/regrid_test_3.mitgrid',1,1)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(newgrid,validated_grid):
            nptest.assert_array_equal(newgrid[a],validated_grid[b])

if __name__=='__main__':
    unittest.main()

