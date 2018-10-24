
class datum(object):
    """ Implements key parameters from the WGS84 datum:

    Class Attributes:
        a (float):  equatorial radius, 6378137.0 meters
        f (float):  flattening, 1/298.257223563 meters
        b (float):  polar semi-minor axis, a*(1-f) = 6356752.3142 meters
    """
    a   = 6378137.0
    f   = 1/298.257223563
    b   = a*(1.-f)

    def __init__(self):
        """ Instantiates class attributes only.

        Args:
            none
        """
