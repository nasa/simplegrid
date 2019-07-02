
import unittest
import simplegrid as sg
import numpy.testing as nptest

class TestAddfringe(unittest.TestCase):

    def test_addfringe_1(self):
        """Add fringe (boundary data) to E edge of first tile, using data from
        the W edge of the second; both tiles are 1x1.  Results are compared with
        mitgrid file created using 'sgaddfringe' command-line call.

        sgaddfringe                             \\
          --tilea tile_A_1x1.mitgrid            \\
          --nia 1                               \\
          --nja 1                               \\
          --tileb tile_B_E_1x1.mitgrid          \\
          --nib 1                               \\
          --njb 1                               \\
          --outfile addfringe_A_EW_1x1.mitgrid  \\
          --strict
        """

        (tilea_edge,tileb_edge,new_tilea_grid) = sg.addfringe.addfringe(
            strict=True,
            tilea='./data/tile_A_1x1.mitgrid',nia=1,nja=1,
            tileb='./data/tile_B_E_1x1.mitgrid',nib=1,njb=1)

        validated_grid = sg.gridio.read_mitgridfile(
            './data/addfringe_A_EW_1x1.mitgrid',1,1,True,False)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(new_tilea_grid,validated_grid):
            nptest.assert_almost_equal(new_tilea_grid[a],validated_grid[b])


    def test_addfringe_2(self):
        """Add fringe (boundary data) to N edge of first tile, using data from
        the S edge of the second; both tiles are 1x1.  Results are compared with
        mitgrid file created using 'sgaddfringe' command-line call.

        sgaddfringe                             \\
          --tilea tile_A_1x1.mitgrid            \\
          --nia 1                               \\
          --nja 1                               \\
          --tileb tile_B_N_1x1.mitgrid          \\
          --nib 1                               \\
          --njb 1                               \\
          --outfile addfringe_A_NS_1x1.mitgrid  \\
          --strict
        """

        (tilea_edge,tileb_edge,new_tilea_grid) = sg.addfringe.addfringe(
            strict=True,
            tilea='./data/tile_A_1x1.mitgrid',nia=1,nja=1,
            tileb='./data/tile_B_N_1x1.mitgrid',nib=1,njb=1)

        validated_grid = sg.gridio.read_mitgridfile(
            './data/addfringe_A_NS_1x1.mitgrid',1,1,True,False)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(new_tilea_grid,validated_grid):
            nptest.assert_almost_equal(new_tilea_grid[a],validated_grid[b])


    def test_addfringe_3(self):
        """Add fringe (boundary data) to W edge of first tile, using data from
        the E edge of the second; both tiles are 1x1.  Results are compared with
        mitgrid file created using 'sgaddfringe' command-line call.

        sgaddfringe                             \\
          --tilea tile_A_1x1.mitgrid            \\
          --nia 1                               \\
          --nja 1                               \\
          --tileb tile_B_W_1x1.mitgrid          \\
          --nib 1                               \\
          --njb 1                               \\
          --outfile addfringe_A_WE_1x1.mitgrid  \\
          --strict
        """

        (tilea_edge,tileb_edge,new_tilea_grid) = sg.addfringe.addfringe(
            strict=True,
            tilea='./data/tile_A_1x1.mitgrid',nia=1,nja=1,
            tileb='./data/tile_B_W_1x1.mitgrid',nib=1,njb=1)

        validated_grid = sg.gridio.read_mitgridfile(
            './data/addfringe_A_WE_1x1.mitgrid',1,1,True,False)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(new_tilea_grid,validated_grid):
            nptest.assert_almost_equal(new_tilea_grid[a],validated_grid[b])


    def test_addfringe_4(self):
        """Add fringe (boundary data) to S edge of first tile, using data from
        the N edge of the second; both tiles are 1x1.  Results are compared with
        mitgrid file created using 'sgaddfringe' command-line call.

        sgaddfringe                             \\
          --tilea tile_A_1x1.mitgrid            \\
          --nia 1                               \\
          --nja 1                               \\
          --tileb tile_B_S_1x1.mitgrid          \\
          --nib 1                               \\
          --njb 1                               \\
          --outfile addfringe_A_SN_1x1.mitgrid  \\
          --strict
        """

        (tilea_edge,tileb_edge,new_tilea_grid) = sg.addfringe.addfringe(
            strict=True,
            tilea='./data/tile_A_1x1.mitgrid',nia=1,nja=1,
            tileb='./data/tile_B_S_1x1.mitgrid',nib=1,njb=1)

        validated_grid = sg.gridio.read_mitgridfile(
            './data/addfringe_A_SN_1x1.mitgrid',1,1,True,False)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(new_tilea_grid,validated_grid):
            nptest.assert_almost_equal(new_tilea_grid[a],validated_grid[b])


    def test_addfringe_5(self):
        """Add fringe (boundary data) to E edge of first tile, using data from
        the W edge of the second; both tiles are 2x2.  Results are compared with
        mitgrid file created using 'sgaddfringe' command-line call.

        sgaddfringe                             \\
          --tilea tile_A_2x2.mitgrid            \\
          --nia 2                               \\
          --nja 2                               \\
          --tileb tile_B_E_2x2.mitgrid          \\
          --nib 2                               \\
          --njb 2                               \\
          --outfile addfringe_A_EW_2x2.mitgrid  \\
          --strict
        """

        (tilea_edge,tileb_edge,new_tilea_grid) = sg.addfringe.addfringe(
            strict=True,
            tilea='./data/tile_A_2x2.mitgrid',nia=2,nja=2,
            tileb='./data/tile_B_E_2x2.mitgrid',nib=2,njb=2)

        validated_grid = sg.gridio.read_mitgridfile(
            './data/addfringe_A_EW_2x2.mitgrid',2,2,True,False)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(new_tilea_grid,validated_grid):
            nptest.assert_almost_equal(new_tilea_grid[a],validated_grid[b])


    def test_addfringe_6(self):
        """Add fringe (boundary data) to N edge of first tile, using data from
        the S edge of the second; both tiles are 2x2.  Results are compared with
        mitgrid file created using 'sgaddfringe' command-line call.

        sgaddfringe                             \\
          --tilea tile_A_2x2.mitgrid            \\
          --nia 2                               \\
          --nja 2                               \\
          --tileb tile_B_N_2x2.mitgrid          \\
          --nib 2                               \\
          --njb 2                               \\
          --outfile addfringe_A_NS_2x2.mitgrid  \\
          --strict
        """

        (tilea_edge,tileb_edge,new_tilea_grid) = sg.addfringe.addfringe(
            strict=True,
            tilea='./data/tile_A_2x2.mitgrid',nia=2,nja=2,
            tileb='./data/tile_B_N_2x2.mitgrid',nib=2,njb=2)

        validated_grid = sg.gridio.read_mitgridfile(
            './data/addfringe_A_NS_2x2.mitgrid',2,2,True,False)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(new_tilea_grid,validated_grid):
            nptest.assert_almost_equal(new_tilea_grid[a],validated_grid[b])


    def test_addfringe_7(self):
        """Add fringe (boundary data) to W edge of first tile, using data from
        the E edge of the second; both tiles are 2x2.  Results are compared with
        mitgrid file created using 'sgaddfringe' command-line call.

        sgaddfringe                             \\
          --tilea tile_A_2x2.mitgrid            \\
          --nia 2                               \\
          --nja 2                               \\
          --tileb tile_B_W_2x2.mitgrid          \\
          --nib 2                               \\
          --njb 2                               \\
          --outfile addfringe_A_WE_2x2.mitgrid  \\
          --strict
        """

        (tilea_edge,tileb_edge,new_tilea_grid) = sg.addfringe.addfringe(
            strict=True,
            tilea='./data/tile_A_2x2.mitgrid',nia=2,nja=2,
            tileb='./data/tile_B_W_2x2.mitgrid',nib=2,njb=2)

        validated_grid = sg.gridio.read_mitgridfile(
            './data/addfringe_A_WE_2x2.mitgrid',2,2,True,False)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(new_tilea_grid,validated_grid):
            nptest.assert_almost_equal(new_tilea_grid[a],validated_grid[b])


    def test_addfringe_8(self):
        """Add fringe (boundary data) to S edge of first tile, using data from
        the N edge of the second; both tiles are 2x2.  Results are compared with
        mitgrid file created using 'sgaddfringe' command-line call.

        sgaddfringe                             \\
          --tilea tile_A_2x2.mitgrid            \\
          --nia 2                               \\
          --nja 2                               \\
          --tileb tile_B_S_2x2.mitgrid          \\
          --nib 2                               \\
          --njb 2                               \\
          --outfile addfringe_A_SN_2x2.mitgrid  \\
          --strict
        """

        (tilea_edge,tileb_edge,new_tilea_grid) = sg.addfringe.addfringe(
            strict=True,
            tilea='./data/tile_A_2x2.mitgrid',nia=2,nja=2,
            tileb='./data/tile_B_S_2x2.mitgrid',nib=2,njb=2)

        validated_grid = sg.gridio.read_mitgridfile(
            './data/addfringe_A_SN_2x2.mitgrid',2,2,True,False)

        # individual comparison of dictionary-stored numpy arrays:
        for a,b in zip(new_tilea_grid,validated_grid):
            nptest.assert_almost_equal(new_tilea_grid[a],validated_grid[b])


if __name__=='__main__':
    unittest.main()

