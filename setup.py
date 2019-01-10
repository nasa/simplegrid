
from setuptools import setup, find_packages
from os import path
from io import open

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name        = 'simplegrid',
    version     = '1.0',
    description = 'Simple local grid creation and refinement for ocean circulation models',
    long_description                = long_description,
    long_description_content_type   = 'text/x-rst',
    packages    = find_packages(exclude=['*tests*']),
    install_requires = [
        'numpy',
        'pyproj',
    ],
    entry_points = {
        'console_scripts': [
            'sgaddfringe    = simplegrid.addfringe:main',
            'sgmkgrid       = simplegrid.mkgrid:main',
            'sgregrid       = simplegrid.regrid:main',
            'sgstitch       = simplegrid.stitch:main',
        ],
    },
)
