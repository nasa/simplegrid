
"""
This module collects values used in various numerical and threshold tests, and
which have proven to be useful defaults for a wide range of ocean regridding
tasks. They can, of course, be modified to better suit other specific
applications.
"""

POINTS_ARE_CLOSE_ENOUGH = 1.e-6
"""Coincident edge detection (relative adjacency of grid points). Used in
matchedges.matchedges().
"""

PLANAR_SPHER_TRANSITION = 10.e3
"""Absolute edge length value (meters) used to determine whether a collection of
surface triangles can be considered flat or spherical. If the mean of both the
x- and y-edge length matrices is less than this value a faceted collection of
surfaces is assumed and, if greater, a spherical surface formulation is used.
Used in computegrid.areas().
"""

