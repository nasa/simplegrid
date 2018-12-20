#!/bin/bash

# join the eastern edge of tile A with the western edge of tile B:

stitch                                  \
  --tilea tile_A_2x2.mitgrid            \
  --nia 2                               \
  --nja 2                               \
  --tileb tile_B_E_2x2.mitgrid          \
  --nib 2                               \
  --njb 2                               \
  --outfile stitch_AB_EW_4x2.mitgrid
