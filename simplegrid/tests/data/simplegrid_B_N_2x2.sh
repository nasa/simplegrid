#!/bin/bash

# make a 2x2 grid that lies to the north of a tile with NW,SE corners at
# (0,1),(1,0):

simplegrid                          \
  --lon1  0.0                       \
  --lat1  2.0                       \
  --lon2  1.0                       \
  --lat2  1.0                       \
  --lon_subscale 2                  \
  --lat_subscale 2                  \
  --outfile tile_B_N_2x2.mitgrid

