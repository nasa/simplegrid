#!/bin/bash

# command-line invocations to generate grid files used in ../test_mkgrid.py

sgmkgrid                            \
    --lon1 1.                       \
    --lat1 2.                       \
    --lon2 2.                       \
    --lat2 1.                       \
    --lon_subscale 1                \
    --lat_subscale 1                \
    --outfile mkgrid_test_1.mitgrid

sgmkgrid                            \
    --lon1 1.                       \
    --lat1 2.                       \
    --lon2 2.                       \
    --lat2 1.                       \
    --lon_subscale 10               \
    --lat_subscale 10               \
    --outfile mkgrid_test_2.mitgrid

