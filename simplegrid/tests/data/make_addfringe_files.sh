#!/bin/bash

# for reference, command-line invocations used to generate validated grid files
# used in ../test_addfringe.py (note: tiles used in these operations were generated
# by ./make_fringe_and_stitch_infiles.sh)


# 1x1 tile tests:

# test_addfringe_1:
sgaddfringe                             \
  --tilea tile_A_1x1.mitgrid            \
  --nia 1                               \
  --nja 1                               \
  --tileb tile_B_E_1x1.mitgrid          \
  --nib 1                               \
  --njb 1                               \
  --outfile addfringe_A_EW_1x1.mitgrid

# test_addfringe_2:
sgaddfringe                             \
  --tilea tile_A_1x1.mitgrid            \
  --nia 1                               \
  --nja 1                               \
  --tileb tile_B_N_1x1.mitgrid          \
  --nib 1                               \
  --njb 1                               \
  --outfile addfringe_A_NS_1x1.mitgrid

# test_addfringe_3:
sgaddfringe                             \
  --tilea tile_A_1x1.mitgrid            \
  --nia 1                               \
  --nja 1                               \
  --tileb tile_B_W_1x1.mitgrid          \
  --nib 1                               \
  --njb 1                               \
  --outfile addfringe_A_WE_1x1.mitgrid

# test_addfringe_4:
sgaddfringe                             \
  --tilea tile_A_1x1.mitgrid            \
  --nia 1                               \
  --nja 1                               \
  --tileb tile_B_S_1x1.mitgrid          \
  --nib 1                               \
  --njb 1                               \
  --outfile addfringe_A_SN_1x1.mitgrid


# 2x2 tile tests:

# test_addfringe_5:
sgaddfringe                             \
  --tilea tile_A_2x2.mitgrid            \
  --nia 2                               \
  --nja 2                               \
  --tileb tile_B_E_2x2.mitgrid          \
  --nib 2                               \
  --njb 2                               \
  --outfile addfringe_A_EW_2x2.mitgrid

# test_addfringe_6:
sgaddfringe                             \
  --tilea tile_A_2x2.mitgrid            \
  --nia 2                               \
  --nja 2                               \
  --tileb tile_B_N_2x2.mitgrid          \
  --nib 2                               \
  --njb 2                               \
  --outfile addfringe_A_NS_2x2.mitgrid

# test_addfringe_7:
sgaddfringe                             \
  --tilea tile_A_2x2.mitgrid            \
  --nia 2                               \
  --nja 2                               \
  --tileb tile_B_W_2x2.mitgrid          \
  --nib 2                               \
  --njb 2                               \
  --outfile addfringe_A_WE_2x2.mitgrid

# test_addfringe_8:
sgaddfringe                             \
  --tilea tile_A_2x2.mitgrid            \
  --nia 2                               \
  --nja 2                               \
  --tileb tile_B_S_2x2.mitgrid          \
  --nib 2                               \
  --njb 2                               \
  --outfile addfringe_A_SN_2x2.mitgrid

