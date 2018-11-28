
import unittest
import simplegrid as sg
import numpy.testing as nptest

class TestMkgrid(unittest.TestCase):

    def test_mkgrid_1(self):
        """Tests creation of simple 1x1 grid on 1deg x 1deg region.
        Results are compared against mitgrid file created using simplegrid
        command-line call:

        simplegrid \
            --lon1 1. \
            --lat1 2. \
            --lon2 2. \
            --lat2 1. \
            --lon_subscale 1 \
            --lat_subscale 1 \
            --outfile mkgrid_test_1.mitgrid

        """

        (newgrid,newgrid_ni,newgrid_nj) = sg.mkgrid.mkgrid(
            lon1=1., lat1=2.,
            lon2=2., lat2=1.,
            lon_subscale=1, lat_subscale=1)

        validated_grid = sg.gridio.read_mitgridfile('./data/mkgrid_test_1.mitgrid',1,1)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(newgrid,validated_grid):
            nptest.assert_array_equal(newgrid[a],validated_grid[b])
        

    def test_mkgrid_2(self):
        """Tests creation of simple 3x3 grid on 1deg x 1deg region.
        Results are compared against mitgrid file created using simplegrid
        command-line call:

        simplegrid \
            --lon1 1. \
            --lat1 2. \
            --lon2 2. \
            --lat2 1. \
            --lon_subscale 3 \
            --lat_subscale 3 \
            --outfile mkgrid_test_2.mitgrid

        """

        (newgrid,newgrid_ni,newgrid_nj) = sg.mkgrid.mkgrid(
            lon1=1., lat1=2.,
            lon2=2., lat2=1.,
            lon_subscale=3, lat_subscale=3)

        validated_grid = sg.gridio.read_mitgridfile('./data/mkgrid_test_2.mitgrid',3,3)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(newgrid,validated_grid):
            nptest.assert_array_equal(newgrid[a],validated_grid[b])

if __name__=='__main__':
    unittest.main()

