#!/bin/bash

# make a 1x1 grid that lies to the north of a tile with NW,SE corners at
# (0,1),(1,0):

simplegrid                          \
  --lon1  0.0                       \
  --lat1  2.0                       \
  --lon2  1.0                       \
  --lat2  1.0                       \
  --lon_subscale 1                  \
  --lat_subscale 1                  \
  --outfile tile_B_N_1x1.mitgrid

