
import unittest
import simplegrid as sg
import numpy.testing as nptest

class TestStitch(unittest.TestCase):

    def test_stitch_1(self):
        """Stitch two 1x1 tiles: E edge of first, W edge of second.
        Results are compared with mitgrid file created using 'stitch'
        command-line call:

        # stitch_AB_EW_2x1.sh:
        sgstitch                                \\
          --tilea tile_A_1x1.mitgrid            \\
          --nia 1                               \\
          --nja 1                               \\
          --tileb tile_B_E_1x1.mitgrid          \\
          --nib 1                               \\
          --njb 1                               \\
          --outfile stitch_AB_EW_2x1.mitgrid
        """

        (newgrid,newgrid_ni,newgrid_nj) = sg.stitch.stitch(
            './data/tile_A_1x1.mitgrid',1,1,
            './data/tile_B_E_1x1.mitgrid',1,1)

        validated_grid = sg.gridio.read_mitgridfile(
            './data/stitch_AB_EW_2x1.mitgrid',2,1)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(newgrid,validated_grid):
            nptest.assert_almost_equal(newgrid[a],validated_grid[b])


    def test_stitch_2(self):
        """Stitch two 1x1 tiles: N edge of first, S edge of second.
        Results are compared with mitgrid file created using 'stitch'
        command-line call:

        # stitch_AB_NS_1x2.sh:
        sgstitch                                \\
          --tilea tile_A_1x1.mitgrid            \\
          --nia 1                               \\
          --nja 1                               \\
          --tileb tile_B_N_1x1.mitgrid          \\
          --nib 1                               \\
          --njb 1                               \\
          --outfile stitch_AB_NS_1x2.mitgrid
        """

        (newgrid,newgrid_ni,newgrid_nj) = sg.stitch.stitch(
            './data/tile_A_1x1.mitgrid',1,1,
            './data/tile_B_N_1x1.mitgrid',1,1)

        validated_grid = sg.gridio.read_mitgridfile(
            './data/stitch_AB_NS_1x2.mitgrid',1,2)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(newgrid,validated_grid):
            nptest.assert_almost_equal(newgrid[a],validated_grid[b])


    def test_stitch_3(self):
        """Stitch two 1x1 tiles: W edge of first, E edge of second.
        Results are compared with mitgrid file created using 'stitch'
        command-line call:

        # stitch_AB_WE_2x1.sh:
        sgstitch                                \\
          --tilea tile_A_1x1.mitgrid            \\
          --nia 1                               \\
          --nja 1                               \\
          --tileb tile_B_W_1x1.mitgrid          \\
          --nib 1                               \\
          --njb 1                               \\
          --outfile stitch_AB_WE_2x1.mitgrid
        """

        (newgrid,newgrid_ni,newgrid_nj) = sg.stitch.stitch(
            './data/tile_A_1x1.mitgrid',1,1,
            './data/tile_B_W_1x1.mitgrid',1,1)

        validated_grid = sg.gridio.read_mitgridfile(
            './data/stitch_AB_WE_2x1.mitgrid',2,1)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(newgrid,validated_grid):
            nptest.assert_almost_equal(newgrid[a],validated_grid[b])


    def test_stitch_4(self):
        """Stitch two 1x1 tiles: S edge of first, N edge of second.
        Results are compared with mitgrid file created using 'stitch'
        command-line call:

        # stitch_AB_SN_1x2.sh:
        sgstitch                                \\
          --tilea tile_A_1x1.mitgrid            \\
          --nia 1                               \\
          --nja 1                               \\
          --tileb tile_B_S_1x1.mitgrid          \\
          --nib 1                               \\
          --njb 1                               \\
          --outfile stitch_AB_SN_1x2.mitgrid
        """

        (newgrid,newgrid_ni,newgrid_nj) = sg.stitch.stitch(
            './data/tile_A_1x1.mitgrid',1,1,
            './data/tile_B_S_1x1.mitgrid',1,1)

        validated_grid = sg.gridio.read_mitgridfile(
            './data/stitch_AB_SN_1x2.mitgrid',1,2)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(newgrid,validated_grid):
            nptest.assert_almost_equal(newgrid[a],validated_grid[b])


    def test_stitch_5(self):
        """Stitch two 2x2 tiles: E edge of first, W edge of second.
        Results are compared with mitgrid file created using 'stitch'
        command-line call:

        # stitch_AB_EW_4x2.sh:
        sgstitch                                \\
          --tilea tile_A_2x2.mitgrid            \\
          --nia 2                               \\
          --nja 2                               \\
          --tileb tile_B_E_2x2.mitgrid          \\
          --nib 2                               \\
          --njb 2                               \\
          --outfile stitch_AB_EW_4x2.mitgrid
        """

        (newgrid,newgrid_ni,newgrid_nj) = sg.stitch.stitch(
            './data/tile_A_2x2.mitgrid',2,2,
            './data/tile_B_E_2x2.mitgrid',2,2)

        validated_grid = sg.gridio.read_mitgridfile(
            './data/stitch_AB_EW_4x2.mitgrid',4,2)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(newgrid,validated_grid):
            nptest.assert_almost_equal(newgrid[a],validated_grid[b])


    def test_stitch_6(self):
        """Stitch two 2x2 tiles: N edge of first, S edge of second.
        Results are compared with mitgrid file created using 'stitch'
        command-line call:

        # stitch_AB_NS_2x4.sh:
        sgstitch                                \\
          --tilea tile_A_2x2.mitgrid            \\
          --nia 2                               \\
          --nja 2                               \\
          --tileb tile_B_N_2x2.mitgrid          \\
          --nib 2                               \\
          --njb 2                               \\
          --outfile stitch_AB_NS_2x4.mitgrid
        """

        (newgrid,newgrid_ni,newgrid_nj) = sg.stitch.stitch(
            './data/tile_A_2x2.mitgrid',2,2,
            './data/tile_B_N_2x2.mitgrid',2,2)

        validated_grid = sg.gridio.read_mitgridfile(
            './data/stitch_AB_NS_2x4.mitgrid',2,4)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(newgrid,validated_grid):
            nptest.assert_almost_equal(newgrid[a],validated_grid[b])


    def test_stitch_7(self):
        """Stitch two 2x2 tiles: W edge of first, E edge of second.
        Results are compared with mitgrid file created using 'stitch'
        command-line call:

        # stitch_AB_WE_4x2.sh:
        sgstitch                                \\
          --tilea tile_A_2x2.mitgrid            \\
          --nia 2                               \\
          --nja 2                               \\
          --tileb tile_B_W_2x2.mitgrid          \\
          --nib 2                               \\
          --njb 2                               \\
          --outfile stitch_AB_WE_4x2.mitgrid
        """

        (newgrid,newgrid_ni,newgrid_nj) = sg.stitch.stitch(
            './data/tile_A_2x2.mitgrid',2,2,
            './data/tile_B_W_2x2.mitgrid',2,2)

        validated_grid = sg.gridio.read_mitgridfile(
            './data/stitch_AB_WE_4x2.mitgrid',4,2)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(newgrid,validated_grid):
            nptest.assert_almost_equal(newgrid[a],validated_grid[b])


    def test_stitch_8(self):
        """Stitch two 2x2 tiles: S edge of first, N edge of second.
        Results are compared with mitgrid file created using 'stitch'
        command-line call:

        # stitch_AB_SN_2x4.sh:
        sgstitch                                \\
          --tilea tile_A_2x2.mitgrid            \\
          --nia 2                               \\
          --nja 2                               \\
          --tileb tile_B_S_2x2.mitgrid          \\
          --nib 2                               \\
          --njb 2                               \\
          --outfile stitch_AB_SN_2x4.mitgrid
        """

        (newgrid,newgrid_ni,newgrid_nj) = sg.stitch.stitch(
            './data/tile_A_2x2.mitgrid',2,2,
            './data/tile_B_S_2x2.mitgrid',2,2)

        validated_grid = sg.gridio.read_mitgridfile(
            './data/stitch_AB_SN_2x4.mitgrid',2,4)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(newgrid,validated_grid):
            nptest.assert_almost_equal(newgrid[a],validated_grid[b])


if __name__=='__main__':
    unittest.main()

