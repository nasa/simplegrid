#!/bin/bash

# join the eastern edge of tile A with the western edge of tile B:

stitch                                  \
  --tilea tile_A_1x1.mitgrid            \
  --nia 1                               \
  --nja 1                               \
  --tileb tile_B_E_1x1.mitgrid          \
  --nib 1                               \
  --njb 1                               \
  --outfile stitch_AB_EW_2x1.mitgrid
