
import numpy as np
import numpy.testing as nptest
import pyproj
import unittest


class TestPyproj(unittest.TestCase):

    def test_pyproj_1(self):
        """pyproj.Geod.inv and .npts validation test: N-S segment distance sums.
        """

        # does the segment distance sum equal the total endpoint distance?:

        lon1,lat1,lon2,lat2 = 0.,0.,0.,1.
        npts = 10
        geod = pyproj.Geod(ellps='sphere')
        allpts = np.concatenate((
            np.array([[lon1,lat1]]),
            np.array(geod.npts(lon1,lat1,lon2,lat2,npts)),
            np.array([[lon2,lat2]])))
        _,_,segdist = geod.inv( allpts[0:-1,0], allpts[0:-1,1], allpts[1:,0],
            allpts[1:,1])
        _,_,totdist = geod.inv(lon1,lat1,lon2,lat2)
        nptest.assert_almost_equal(np.sum(segdist),totdist)


    def test_pyproj_2(self):
        """pyproj.Geod.inv and .npts validation test: SW-NE segment distance sums.
        """

        # does the segment distance sum equal the total endpoint distance?:

        lon1,lat1,lon2,lat2 = 0.,0.,1.,1.
        npts = 10
        geod = pyproj.Geod(ellps='sphere')
        allpts = np.concatenate((
            np.array([[lon1,lat1]]),
            np.array(geod.npts(lon1,lat1,lon2,lat2,npts)),
            np.array([[lon2,lat2]])))
        _,_,segdist = geod.inv( allpts[0:-1,0], allpts[0:-1,1], allpts[1:,0],
            allpts[1:,1])
        _,_,totdist = geod.inv(lon1,lat1,lon2,lat2)
        nptest.assert_almost_equal(np.sum(segdist),totdist)


    def test_pyproj_3(self):
        """pyproj.Geod.inv and .npts validation test: E-W segment distance sums.
        """

        # does the segment distance sum equal the total endpoint distance?:

        lon1,lat1,lon2,lat2 = 0.,0.,1.,0.
        npts = 10
        geod = pyproj.Geod(ellps='sphere')
        allpts = np.concatenate((
            np.array([[lon1,lat1]]),
            np.array(geod.npts(lon1,lat1,lon2,lat2,npts)),
            np.array([[lon2,lat2]])))
        _,_,segdist = geod.inv( allpts[0:-1,0], allpts[0:-1,1], allpts[1:,0],
            allpts[1:,1])
        _,_,totdist = geod.inv(lon1,lat1,lon2,lat2)
        nptest.assert_almost_equal(np.sum(segdist),totdist)


if __name__=='__main__':
    unittest.main()

