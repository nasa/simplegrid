
import unittest

import pyproj
import simplegrid as sg


class TestNearest(unittest.TestCase):

    def test_nearest_sw_corner(self):
        geod = pyproj.Geod(ellps='sphere')
        mg = sg.gridio.read_mitgridfile('./data/tile005.mitgrid', 270, 90)
        i,j,dist = sg.util.nearest(-128.,67.5,mg['XG'],mg['YG'],geod)
        self.assertEqual((i,j),(0,0))
        self.assertAlmostEqual(dist,1.20941759e-09)

    def test_nearest_ne_corner(self):
        geod = pyproj.Geod(ellps='sphere')
        mg = sg.gridio.read_mitgridfile('./data/tile005.mitgrid', 270, 90)
        i,j,dist = sg.util.nearest(-115.,-88.17570,mg['XG'],mg['YG'],geod)
        self.assertEqual((i,j),(270,90))
        self.assertAlmostEqual(dist,1.14379740)

    def test_nearest_center(self):
        geod = pyproj.Geod(ellps='sphere')
        mg = sg.gridio.read_mitgridfile('./data/tile005.mitgrid', 270, 90)
        i,j,dist = sg.util.nearest(-83.,-24.310,mg['XG'],mg['YG'],geod)
        self.assertEqual((i,j),(135,45))
        self.assertAlmostEqual(dist,6.2719790)


if __name__=='__main__':
    unittest.main()

