#!/bin/bash

simpleregrid \
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
    --outfile regrid_test_5.mitgrid

