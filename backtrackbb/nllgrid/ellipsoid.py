# -*- coding: utf8 -*-
"""
Python porting of c-code from matrix_statistics.c.

Part of the NonLinLoc package, written by Anthony Lomax.

:copyright:
    2013-2018 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement, Version 2.1
    (http://www.cecill.info/index.en.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from math import pi, cos, sin


class Ellipsoid3D():
    """A 3D ellipsoid."""

    # semi-minor axis km
    az1 = None
    dip1 = None
    len1 = None
    # semi-intermediate axis km
    az2 = None
    dip2 = None
    len2 = None
    # semi-major axis km
    len3 = None

    def __str__(self):
        """String representation."""
        s = 'az1: %f\n' % self.az1
        s += 'dip1: %f\n' % self.dip1
        s += 'len1: %f\n' % self.len1
        s += 'az2: %f\n' % self.az2
        s += 'dip2: %f\n' % self.dip2
        s += 'len2: %f\n' % self.len2
        s += 'len3: %f' % self.len3
        return s


class Vect3D():
    """A 3D vector."""

    x = None
    y = None
    z = None

    def __str__(self):
        """String representation."""
        s = 'x: %f y: %f z: %f' % (self.x, self.y, self.z)
        return s


def cross_product_3d(vect_a, vect_b):
    """Cross product between two 3D vectors."""
    product = Vect3D()
    product.x = vect_a.y * vect_b.z - vect_a.z * vect_b.y
    product.y = vect_a.z * vect_b.x - vect_a.x * vect_b.z
    product.z = vect_a.x * vect_b.y - vect_a.y * vect_b.x
    return product


def ellipsiod2Axes(pellipsoid):
    """Get the three axes of an ellipsoid."""
    DE2RA = pi/180

    # strike angles positive CCW from East = 0
    az1 = 90.0 - pellipsoid.az1
    az2 = 90.0 - pellipsoid.az2
    # dip angles increasing downwards from horiz = 0
    dip1 = -pellipsoid.dip1
    dip2 = -pellipsoid.dip2

    # get 3D vector axes
    cosd1 = cos(DE2RA * pellipsoid.dip1)
    paxis1 = Vect3D()
    paxis1.x = cos(DE2RA * az1) * cosd1
    paxis1.y = sin(DE2RA * az1) * cosd1
    paxis1.z = -sin(DE2RA * dip1)

    cosd2 = cos(DE2RA * pellipsoid.dip2)
    paxis2 = Vect3D()
    paxis2.x = cos(DE2RA * az2) * cosd2
    paxis2.y = sin(DE2RA * az2) * cosd2
    paxis2.z = -sin(DE2RA * dip2)

    paxis3 = cross_product_3d(paxis1, paxis2)

    paxis1.x *= pellipsoid.len1
    paxis1.y *= pellipsoid.len1
    paxis1.z *= pellipsoid.len1
    paxis2.x *= pellipsoid.len2
    paxis2.y *= pellipsoid.len2
    paxis2.z *= pellipsoid.len2
    paxis3.x *= pellipsoid.len3
    paxis3.y *= pellipsoid.len3
    paxis3.z *= pellipsoid.len3

    return paxis1, paxis2, paxis3


def toEllipsoid3D(ax1, ax2, center, npts):
    """Get the coordinates of a 3D ellipsoid defined by 2 axes and a center."""
    d_angle = 2.0 * pi / (npts - 1)
    angle = 0.0

    ellArray = []
    for n in range(0, npts):
        cosang = cos(angle)
        sinang = sin(angle)
        vect = Vect3D()
        vect.x = center.x + ax1.x * cosang + ax2.x * sinang
        vect.y = center.y + ax1.y * cosang + ax2.y * sinang
        vect.z = center.z + ax1.z * cosang + ax2.z * sinang
        ellArray.append(vect)
        angle += d_angle

    return ellArray


def main():
    """Test code."""
    import numpy as np
    import matplotlib.pyplot as plt

    ellips = Ellipsoid3D()
    ellips.az1 = 106.592
    ellips.dip1 = -11.4799
    ellips.len1 = 10.7175
    ellips.az2 = 16.9256
    ellips.dip2 = 1.64329
    ellips.len2 = 13.3939
    ellips.len3 = 2.514136e+01

    expect = Vect3D()
    expect.x = 93.1985
    expect.y = 2.36646
    expect.z = 25.8457

    pax1, pax2, pax3 = ellipsiod2Axes(ellips)

    ellArray12 = toEllipsoid3D(pax1, pax2, expect, 100)
    ellArray13 = toEllipsoid3D(pax1, pax3, expect, 100)
    ellArray23 = toEllipsoid3D(pax2, pax3, expect, 100)
    ell12 = np.array([(vect.x, vect.y) for vect in ellArray12])
    ell13 = np.array([(vect.x, vect.y) for vect in ellArray13])
    ell23 = np.array([(vect.x, vect.y) for vect in ellArray23])

    plt.plot(ell12[:, 0], ell12[:, 1])
    plt.plot(ell13[:, 0], ell13[:, 1])
    plt.plot(ell23[:, 0], ell23[:, 1])

    plt.show()


if __name__ == '__main__':
    main()
