#!/bin/bash

# create test grid files in dependent order:

make_mkgrid_files.sh && \
make_regrid_files.sh && \
make_fringe_and_stitch_files.sh && \
make_stitch_files.sh && \
make_addfringe_files.sh
