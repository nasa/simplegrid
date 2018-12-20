#!/bin/bash

# make a 2x2 grid with NW,SE corners at (0,1), (1,0):

simplegrid                      \
  --lon1  0.                    \
  --lat1  1.                    \
  --lon2  1.                    \
  --lat2  0.                    \
  --lon_subscale 2              \
  --lat_subscale 2              \
  --outfile tile_A_2x2.mitgrid

