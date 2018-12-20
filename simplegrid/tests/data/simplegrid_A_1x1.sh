#!/bin/bash

# make a 1x1 grid with NW,SE corners at (0,1), (1,0):

simplegrid                      \
  --lon1  0.                    \
  --lat1  1.                    \
  --lon2  1.                    \
  --lat2  0.                    \
  --lon_subscale 1              \
  --lat_subscale 1              \
  --outfile tile_A_1x1.mitgrid

