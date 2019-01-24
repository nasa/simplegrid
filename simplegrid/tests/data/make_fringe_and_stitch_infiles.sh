#!/bin/bash

# for reference, scripts used to generate simple validated 1x1 and 2x2 input
# tiles for basic addfringe and stitch testing:


# group of 1x1 tiles:

# 1x1 grid with NW,SE corners at (0,1), (1,0):
sgmkgrid                            \
  --lon1  0.                        \
  --lat1  1.                        \
  --lon2  1.                        \
  --lat2  0.                        \
  --lon_subscale 1                  \
  --lat_subscale 1                  \
  --outfile tile_A_1x1.mitgrid

# 1x1 grid that lies to the west of a tile with NW,SE corners at
# (0,1),(1,0):
sgmkgrid                            \
  --lon1 -1.0                       \
  --lat1  1.0                       \
  --lon2  0.0                       \
  --lat2  0.0                       \
  --lon_subscale 1                  \
  --lat_subscale 1                  \
  --outfile tile_B_W_1x1.mitgrid

# 1x1 grid that lies to the north of a tile with NW,SE corners at
# (0,1),(1,0):
sgmkgrid                            \
  --lon1  0.0                       \
  --lat1  2.0                       \
  --lon2  1.0                       \
  --lat2  1.0                       \
  --lon_subscale 1                  \
  --lat_subscale 1                  \
  --outfile tile_B_N_1x1.mitgrid

# 1x1 grid that lies to the east of a tile with NW,SE corners at
# (0,1),(1,0):
sgmkgrid                            \
  --lon1  1.0                       \
  --lat1  1.0                       \
  --lon2  2.0                       \
  --lat2  0.0                       \
  --lon_subscale 1                  \
  --lat_subscale 1                  \
  --outfile tile_B_E_1x1.mitgrid

# 1x1 grid that lies to the south of a tile with NW,SE corners at
# (0,1),(1,0):
sgmkgrid                            \
  --lon1  0.0                       \
  --lat1  0.0                       \
  --lon2  1.0                       \
  --lat2 -1.0                       \
  --lon_subscale 1                  \
  --lat_subscale 1                  \
  --outfile tile_B_S_1x1.mitgrid


# group of 2x2 tiles:

# 2x2 grid with NW,SE corners at (0,1), (1,0):
sgmkgrid                            \
  --lon1  0.                        \
  --lat1  1.                        \
  --lon2  1.                        \
  --lat2  0.                        \
  --lon_subscale 2                  \
  --lat_subscale 2                  \
  --outfile tile_A_2x2.mitgrid

# 2x2 grid that lies to the west of a tile with NW,SE corners at
# (0,1),(1,0):
sgmkgrid                            \
  --lon1 -1.0                       \
  --lat1  1.0                       \
  --lon2  0.0                       \
  --lat2  0.0                       \
  --lon_subscale 2                  \
  --lat_subscale 2                  \
  --outfile tile_B_W_2x2.mitgrid

# 2x2 grid that lies to the north of a tile with NW,SE corners at
# (0,1),(1,0):
sgmkgrid                            \
  --lon1  0.0                       \
  --lat1  2.0                       \
  --lon2  1.0                       \
  --lat2  1.0                       \
  --lon_subscale 2                  \
  --lat_subscale 2                  \
  --outfile tile_B_N_2x2.mitgrid

# 2x2 grid that lies to the east of a tile with NW,SE corners at
# (0,1),(1,0):
sgmkgrid                            \
  --lon1  1.0                       \
  --lat1  1.0                       \
  --lon2  2.0                       \
  --lat2  0.0                       \
  --lon_subscale 2                  \
  --lat_subscale 2                  \
  --outfile tile_B_E_2x2.mitgrid

# 2x2 grid that lies to the south of a tile with NW,SE corners at
# (0,1),(1,0):
sgmkgrid                            \
  --lon1  0.0                       \
  --lat1  0.0                       \
  --lon2  1.0                       \
  --lat2 -1.0                       \
  --lon_subscale 2                  \
  --lat_subscale 2                  \
  --outfile tile_B_S_2x2.mitgrid

