#!/bin/bash

# make a 2x2 grid that lies to the west of a tile with NW,SE corners at
# (0,1),(1,0):

simplegrid                          \
  --lon1 -1.0                       \
  --lat1  1.0                       \
  --lon2  0.0                       \
  --lat2  0.0                       \
  --lon_subscale 2                  \
  --lat_subscale 2                  \
  --outfile tile_B_W_2x2.mitgrid

