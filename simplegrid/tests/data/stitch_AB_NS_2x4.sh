#!/bin/bash

# join the northern edge of tile A with the southern edge of tile B:

stitch                                  \
  --tilea tile_A_2x2.mitgrid            \
  --nia 2                               \
  --nja 2                               \
  --tileb tile_B_N_2x2.mitgrid          \
  --nib 2                               \
  --njb 2                               \
  --outfile stitch_AB_NS_2x4.mitgrid
