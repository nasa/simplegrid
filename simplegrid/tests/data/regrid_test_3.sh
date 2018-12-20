#!/bin/bash

simpleregrid \
    --mitgridfile tile005.mitgrid   \
    --ni 270                        \
    --nj  90                        \
    --lon1 -127.73445435            \
    --lat1   67.56064719            \
    --lon2 -128.0                   \
    --lat2   67.40168504            \
    --lon_subscale 1                \
    --lat_subscale 1                \
    --outfile regrid_test_3.mitgrid

