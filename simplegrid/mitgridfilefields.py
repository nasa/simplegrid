
"""Attributes that are useful when reading 'packed' binary mitgrid files.

mitgrid files are simply raw streams of binary data without any data field
separators or record length descriptors. To read such a file, then, one has to
know the order of the data conventionally appearing in the file, the data type,
and amount of padding that may have been used when written. For simplicity, all
matrices are written, and thus read, assuming the same dimensions, even if the
underlying data are actually smaller.

For example, a nominal 90 x 270 grid may include fields that are actually up to
91 x 271; an example of the former would be an array of tracer cell center point
longitude coordinates, the latter, the corresponding cell corner longitudes.
However, both will have been written as 91x271 arrays, the former padded with a
single trailing null row and column. Other matrices may have similar row and/or
column padding as noted by delta_sizes.

Attributes:
    names (str): commonly-used field names
    ni_delta_sizes (int): corresponding field ni "padding" (0 implies no
        padding, 1 implies a trailing padded row)
    nj_delta_sizes (int): corresponding field nj "offsets" (0 implies no
        padding, 1 implies a trailing padded column)
    datatype (str): '>f8' for big-endian, 64-bit, floating-point doubles
"""

names = (
 'XC', 'YC','DXF','DYF','RAC', 'XG', 'YG','DXV','DYU','RAZ','DXC','DYC','RAW','RAS','DXG','DYG')

ni_delta_sizes = (
   0,    0,    0,    0,    0,    1,    1,    1,    1,    1,    1,    0,    1,    0,    0,    1 )

nj_delta_sizes = (
   0,    0,    0,    0,    0,    1,    1,    1,    1,    1,    0,    1,    0,    1,    1,    0 )

endianness = '>'        # big
precision_bytes  = '8'  # double
type = 'f'              # float
order = 'F'             # fortran

datatype = endianness + type + precision_bytes

