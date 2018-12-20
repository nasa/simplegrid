#!/bin/bash

# join the western edge of tile A with the eastern edge of tile B:

stitch                                  \
  --tilea tile_A_1x1.mitgrid            \
  --nia 1                               \
  --nja 1                               \
  --tileb tile_B_W_1x1.mitgrid          \
  --nib 1                               \
  --njb 1                               \
  --outfile stitch_AB_WE_2x1.mitgrid
