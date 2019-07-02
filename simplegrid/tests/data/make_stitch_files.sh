#!/bin/bash

# for reference, command-line invocations used to generate validated grid files
# used in ../test_stitch.py (note: tiles used in these operations were generated
# by ./make_fringe_and_stitch_infiles.sh)


# 1x1 tile tests:

# test_stitch_1:
sgstitch                                \
  --tilea tile_A_1x1.mitgrid            \
  --nia 1                               \
  --nja 1                               \
  --tileb tile_B_E_1x1.mitgrid          \
  --nib 1                               \
  --njb 1                               \
  --outfile stitch_AB_EW_2x1.mitgrid    \
  --strict

# test_stitch_2:
sgstitch                                \
  --tilea tile_A_1x1.mitgrid            \
  --nia 1                               \
  --nja 1                               \
  --tileb tile_B_N_1x1.mitgrid          \
  --nib 1                               \
  --njb 1                               \
  --outfile stitch_AB_NS_1x2.mitgrid    \
  --strict

# test_stitch_3:
sgstitch                                \
  --tilea tile_A_1x1.mitgrid            \
  --nia 1                               \
  --nja 1                               \
  --tileb tile_B_W_1x1.mitgrid          \
  --nib 1                               \
  --njb 1                               \
  --outfile stitch_AB_WE_2x1.mitgrid    \
  --strict

# test_stitch_4:
sgstitch                                \
  --tilea tile_A_1x1.mitgrid            \
  --nia 1                               \
  --nja 1                               \
  --tileb tile_B_S_1x1.mitgrid          \
  --nib 1                               \
  --njb 1                               \
  --outfile stitch_AB_SN_1x2.mitgrid    \
  --strict


# 2x2 tile tests:

# test_stitch_5:
sgstitch                                \
  --tilea tile_A_2x2.mitgrid            \
  --nia 2                               \
  --nja 2                               \
  --tileb tile_B_E_2x2.mitgrid          \
  --nib 2                               \
  --njb 2                               \
  --outfile stitch_AB_EW_4x2.mitgrid    \
  --strict

# test_stitch_6:
sgstitch                                \
  --tilea tile_A_2x2.mitgrid            \
  --nia 2                               \
  --nja 2                               \
  --tileb tile_B_N_2x2.mitgrid          \
  --nib 2                               \
  --njb 2                               \
  --outfile stitch_AB_NS_2x4.mitgrid    \
  --strict

# test_stitch_7:
sgstitch                                \
  --tilea tile_A_2x2.mitgrid            \
  --nia 2                               \
  --nja 2                               \
  --tileb tile_B_W_2x2.mitgrid          \
  --nib 2                               \
  --njb 2                               \
  --outfile stitch_AB_WE_4x2.mitgrid    \
  --strict

# test_stitch_8:
sgstitch                                \
  --tilea tile_A_2x2.mitgrid            \
  --nia 2                               \
  --nja 2                               \
  --tileb tile_B_S_2x2.mitgrid          \
  --nib 2                               \
  --njb 2                               \
  --outfile stitch_AB_SN_2x4.mitgrid    \
  --strict

