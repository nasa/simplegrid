#!/bin/bash

# for reference, command-line invocations used to generate validated grid files
# used in ../test_regrid.py:

# test_regrid_1:
sgregrid \
    --mitgridfile tile005.mitgrid   \
    --ni 270                        \
    --nj  90                        \
    --lon1 -126.95921008            \
    --lat1   67.31607921            \
    --lon2 -126.21154908            \
    --lat2   67.17736362            \
    --lon_subscale 1                \
    --lat_subscale 1                \
    --outfile regrid_test_1.mitgrid \
    --strict

# test_regrid_2:
sgregrid \
    --mitgridfile tile005.mitgrid   \
    --ni 270                        \
    --nj  90                        \
    --lon1 -126.95921008            \
    --lat1   67.31607921            \
    --lon2 -126.21154908            \
    --lat2   67.17736362            \
    --lon_subscale 30               \
    --lat_subscale 20               \
    --outfile regrid_test_2.mitgrid \
    --strict

# test_regrid_3:
sgregrid \
    --mitgridfile tile005.mitgrid   \
    --ni 270                        \
    --nj  90                        \
    --lon1 -127.73445435            \
    --lat1   67.56064719            \
    --lon2 -128.0                   \
    --lat2   67.40168504            \
    --lon_subscale 1                \
    --lat_subscale 1                \
    --outfile regrid_test_3.mitgrid \
    --strict

# test_regrid_4:
sgregrid \
    --mitgridfile tile005.mitgrid   \
    --ni 270                        \
    --nj  90                        \
    --lon1 -127.73445435            \
    --lat1   67.56064719            \
    --lon2 -128.0                   \
    --lat2   67.40168504            \
    --lon_subscale 10               \
    --lat_subscale 10               \
    --outfile regrid_test_4.mitgrid \
    --strict

# test_regrid_5:
sgregrid \
    --xg_file tile005_XG.bin        \
    --yg_file tile005_YG.bin        \
    --ni 270                        \
    --nj  90                        \
    --lon1 -126.95921008            \
    --lat1   67.31607921            \
    --lon2 -126.21154908            \
    --lat2   67.17736362            \
    --lon_subscale 1                \
    --lat_subscale 1                \
    --outfile regrid_test_5.mitgrid \
    --strict

# test_regrid_6:
sgregrid \
    --xg_file tile005_XG.csv        \
    --yg_file tile005_YG.csv        \
    --ni 270                        \
    --nj  90                        \
    --lon1 -126.95921008            \
    --lat1   67.31607921            \
    --lon2 -126.21154908            \
    --lat2   67.17736362            \
    --lon_subscale 1                \
    --lat_subscale 1                \
    --outfile regrid_test_6.mitgrid \
    --strict

