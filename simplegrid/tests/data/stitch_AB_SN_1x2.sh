#!/bin/bash

# join the southern edge of tile A with the northern edge of tile B:

stitch                                  \
  --tilea tile_A_1x1.mitgrid            \
  --nia 1                               \
  --nja 1                               \
  --tileb tile_B_S_1x1.mitgrid          \
  --nib 1                               \
  --njb 1                               \
  --outfile stitch_AB_SN_1x2.mitgrid
