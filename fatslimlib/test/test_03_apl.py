# -*- coding: utf-8; Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
#  Copyright (C) 2013-2016  Sébastien Buchoux <sebastien.buchoux@gmail.com>
#
#    This file is part of FATSLiM.
#
#    FATSLiM is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    FATSLiM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with FATSLiM.  If not, see <http://www.gnu.org/licenses/>.
from __future__ import print_function

# Global imports
import numpy as np
from numpy.testing import assert_almost_equal, assert_allclose

# Local imports
from .. import core_geometry
from .. import core_analysis
from . import frame_model_bilayer, frame_bilayer, frame_model_vesicle, frame_vesicle, frame_model_bilayer_prot,\
    frame_bilayer_prot, frame_bilayer_peptide, frame_bilayer_chol

RTOL = 5e-3

triangle_points = [(200., 50),
                   (150, 100),
                   (250, 150)]

square_points = [[100., 100], [300, 100], [300, 300], [100, 300]]

small_square_points = [[150., 100], [250, 100], [250, 200], [150, 200]]

far_square_points = [[1100., 100], [1300, 1100], [1300, 1300], [1100, 1300]]

concave_points = [[50., 150], [200, 50], [350, 150], [350, 300], [250, 300], [200, 250], [150, 350], [100, 250],
                  [100, 200]]


def test_polygon_perimeter_triangle():
    ref_perimeter = 294.3174758686337
    perimeter = core_geometry.get_polygon_perimeter(np.array(triangle_points))

    assert_almost_equal(perimeter, ref_perimeter)


def test_polygon_area_triangle():
    ref_area = 3750.0
    area = core_geometry.get_polygon_area(np.array(triangle_points))

    assert_almost_equal(area, ref_area)


def test_polygon_perimeter_square():
    ref_perimeter = 800.0
    perimeter = core_geometry.get_polygon_perimeter(np.array(square_points))

    assert_almost_equal(perimeter, ref_perimeter)


def test_polygon_area_square():
    ref_area = 40000.0
    area = core_geometry.get_polygon_area(np.array(square_points))

    assert_almost_equal(area, ref_area)


def test_polygon_is_inside():
    poly = np.array(square_points)
    point_inside = np.array([200., 200])
    point_outside = np.array([2000., 2000])

    assert core_geometry.is_inside_polygon(poly, point_inside)
    assert not core_geometry.is_inside_polygon(poly, point_outside)


def test_polygon_clipping_intersecting():
    concave = np.array(concave_points)
    square = np.array(square_points)

    ref_clipped = [[100, 116.66], [125, 100], [275, 100], [300, 116.66],
                   [300, 300], [250, 300], [200, 250], [175, 300], [125, 300], [100, 250]]
    clipped = core_geometry.get_clipped_polygon(concave, square)
    assert_almost_equal(clipped, ref_clipped, 2)

    ref_area = 37083.333
    area = core_geometry.get_polygon_area(clipped)
    assert_almost_equal(area, ref_area, 3)


def test_polygon_clipping_not_intersecting():
    concave = np.array(concave_points)
    square = np.array(far_square_points)

    clipped = core_geometry.get_clipped_polygon(concave, square)
    assert len(clipped) == 0

    ref_area = -1.0
    area = core_geometry.get_polygon_area(clipped)
    assert_almost_equal(area, ref_area, 3)


def test_polygon_clipping_contained():
    concave = np.array(concave_points)
    square = np.array(small_square_points)

    ref_clipped = [small_square_points[i - 1] for i in range(len(small_square_points))]
    clipped = core_geometry.get_clipped_polygon(concave, square)
    assert_almost_equal(clipped, ref_clipped, 2)

    ref_area = 10000
    area = core_geometry.get_polygon_area(clipped)
    assert_almost_equal(area, ref_area, 3)


def test_put_atoms_on_plane_model_bilayer(frame_model_bilayer):
    frame = frame_model_bilayer
    membrane = frame.get_membranes()[0]

    ref_2dpoints = [(-1.60, -0.80), (-1.60, -0.00), (-1.60, 0.80), (-0.80, -1.60), (-0.80, -0.80), (-0.80, -0.00),
                    (-0.80, 0.80), (-0.80, 1.60), (0.00, -1.60), (0.00, -0.80), (0.00, 0.80),
                    (0.00, 1.60), (0.80, -1.60), (0.80, -0.80), (0.80, 0.00), (0.80, 0.80), (0.80, 1.60), (1.60, -0.80),
                    (1.60, 0.00), (1.60, 0.80)]

    for leaflet in membrane:
        leaflet_2dpoints = core_analysis.test_put_atoms_on_plane(leaflet)

        for points in leaflet_2dpoints:
            points = np.array([tuple(val) for val in np.round(points, decimals=2)],
                              dtype=[(str('x'), float), (str('y'), float)])
            points.sort(order=(str('x'), str('y')))

            assert_almost_equal(points.tolist(), ref_2dpoints)


def test_apl_model_bilayer(frame_model_bilayer):
    frame = frame_model_bilayer
    membrane = frame.get_membranes()[0]

    apl = membrane.get_apl(by_type=False, only_average=False)

    ref_apl = [np.array([0.64] * 36),
               np.array([0.64] * 36)]

    for i, value in enumerate(apl[1:]):
        assert_allclose(value, ref_apl[i], rtol=1e-2)


def test_apl_bilayer(frame_bilayer):
    frame = frame_bilayer
    membrane = frame.get_membranes()[0]

    ref_apl = (0.664, 0.664, 0.663)

    apl = membrane.get_apl(by_type=False)

    print("APL: avg:%.3f - L1 (%i lipids): APL=%.3f (min:%.3f - max:%.3f), area=%.3f - L2 (%i lipids): APL=%.3f (min:%.3f - max:%.3f), area=%.3f " % (
        apl[0],
        len(membrane.leaflet1), apl[1][0], apl[1][1], apl[1][2], apl[1][3],
        len(membrane.leaflet2), apl[2][0], apl[2][1], apl[2][2], apl[2][3]))

    assert_allclose(apl[0], ref_apl[0], rtol=RTOL)
    assert_allclose(apl[1][0], ref_apl[1], rtol=RTOL)
    assert_allclose(apl[2][0], ref_apl[2], rtol=RTOL)


def test_apl_bilayer_chol(frame_bilayer_chol):
    frame = frame_bilayer_chol
    membrane = frame.get_membranes()[0]

    ref_apl = (0.490, 0.489, 0.491)
    ref_area = (477.810, 475.138)

    apl = membrane.get_apl(by_type=False)

    print("APL: avg:%.3f - L1 (%i lipids): APL=%.3f (min:%.3f - max:%.3f), area=%.3f - L2 (%i lipids): APL=%.3f (min:%.3f - max:%.3f), area=%.3f " % (
        apl[0],
        len(membrane.leaflet1), apl[1][0], apl[1][1], apl[1][2], apl[1][3],
        len(membrane.leaflet2), apl[2][0], apl[2][1], apl[2][2], apl[2][3]))

    assert_allclose(apl[0], ref_apl[0], rtol=RTOL)
    assert_allclose(apl[1][0], ref_apl[1], rtol=RTOL)
    assert_allclose(apl[2][0], ref_apl[2], rtol=RTOL)
    assert_allclose(apl[1][3], ref_area[0], rtol=RTOL)
    assert_allclose(apl[2][3], ref_area[1], rtol=RTOL)


def test_apl_bilayer_chol_by_type(frame_bilayer_chol):
    frame = frame_bilayer_chol
    membrane = frame.get_membranes()[0]

    ref_apl_chol = (0.389, 0.387)
    ref_area_chol = (114.109, 109.63)
    ref_apl_dppc = (0.456, 0.471)
    ref_area_dppc = (188.035, 196.033)
    ref_apl_dupc = (0.647, 0.633)
    ref_area_dupc = (176.071, 169.683)

    apl = membrane.get_apl()
    for lid in [1, 2]:
        for key, val in apl[lid].items():
            print("Leaflet #%i - %i %s: APL=%.3f (min:%.3f, max:%.3f) - Area=%.3f" %
                  (lid, val[0], key,
                   val[1], val[2], val[3], val[4]))

    # Test Cholesterol
    assert_allclose(apl[1]["CHOL"][1], ref_apl_chol[0], rtol=1e-2)
    assert_allclose(apl[2]["CHOL"][1], ref_apl_chol[1], rtol=1e-2)
    assert_allclose(apl[1]["CHOL"][4], ref_area_chol[0], rtol=1e-2)
    assert_allclose(apl[2]["CHOL"][4], ref_area_chol[1], rtol=1e-2)

    # Test DPPC
    assert_allclose(apl[1]["DPPC"][1], ref_apl_dppc[0], rtol=1e-2)
    assert_allclose(apl[2]["DPPC"][1], ref_apl_dppc[1], rtol=1e-2)
    assert_allclose(apl[1]["DPPC"][4], ref_area_dppc[0], rtol=1e-2)
    assert_allclose(apl[2]["DPPC"][4], ref_area_dppc[1], rtol=1e-2)

    # Test DUPC
    assert_allclose(apl[1]["DUPC"][1], ref_apl_dupc[0], rtol=1e-2)
    assert_allclose(apl[2]["DUPC"][1], ref_apl_dupc[1], rtol=1e-2)
    assert_allclose(apl[1]["DUPC"][4], ref_area_dupc[0], rtol=1e-2)
    assert_allclose(apl[2]["DUPC"][4], ref_area_dupc[1], rtol=1e-2)


def test_apl_vesicle_model(frame_model_vesicle):
    frame = frame_model_vesicle
    membrane = frame.get_membranes()[0]

    ref_apl = (0.570, 0.639, 0.398)
    ref_area = (1254.087, 312.537)

    apl = membrane.get_apl(by_type=False)

    print("APL: avg:%.3f - L1 (%i lipids): APL=%.3f (min:%.3f - max:%.3f), area=%.3f - L2 (%i lipids): APL=%.3f (min:%.3f - max:%.3f), area=%.3f " % (
        apl[0],
        len(membrane.leaflet1), apl[1][0], apl[1][1], apl[1][2], apl[1][3],
        len(membrane.leaflet2), apl[2][0], apl[2][1], apl[2][2], apl[2][3]))

    assert_allclose(apl[0], ref_apl[0], rtol=RTOL)
    assert_allclose(apl[1][0], ref_apl[1], rtol=RTOL)
    assert_allclose(apl[2][0], ref_apl[2], rtol=RTOL)
    assert_allclose(apl[1][3], ref_area[0], rtol=RTOL)
    assert_allclose(apl[2][3], ref_area[1], rtol=RTOL)


def test_apl_vesicle(frame_vesicle):
    frame = frame_vesicle
    membrane = frame.get_membranes()[0]

    ref_apl = (0.681, 0.797, 0.499)
    ref_area = (1474.557, 588.391)

    apl = membrane.get_apl(by_type=False)

    print("APL: avg:%.3f - L1 (%i lipids): APL=%.3f (min:%.3f - max:%.3f), area=%.3f - L2 (%i lipids): APL=%.3f (min:%.3f - max:%.3f), area=%.3f " % (
        apl[0],
        len(membrane.leaflet1), apl[1][0], apl[1][1], apl[1][2], apl[1][3],
        len(membrane.leaflet2), apl[2][0], apl[2][1], apl[2][2], apl[2][3]))

    assert_allclose(apl[0], ref_apl[0], rtol=RTOL)
    assert_allclose(apl[1][0], ref_apl[1], rtol=RTOL)
    assert_allclose(apl[2][0], ref_apl[2], rtol=RTOL)
    assert_allclose(apl[1][3], ref_area[0], rtol=RTOL)
    assert_allclose(apl[2][3], ref_area[1], rtol=RTOL)


def test_apl_bilayer_peptide(frame_bilayer_peptide):
    frame = frame_bilayer_peptide
    membrane = frame.get_membranes()[0]

    ref_apl = (0.614, 0.609, 0.618)
    ref_area = (38.397, 38.934)

    apl = membrane.get_apl(by_type=False)

    print("APL: avg:%.3f - L1 (%i lipids): APL=%.3f (min:%.3f - max:%.3f), area=%.3f - L2 (%i lipids): APL=%.3f (min:%.3f - max:%.3f), area=%.3f " % (
        apl[0],
        len(membrane.leaflet1), apl[1][0], apl[1][1], apl[1][2], apl[1][3],
        len(membrane.leaflet2), apl[2][0], apl[2][1], apl[2][2], apl[2][3]))

    assert_allclose(apl[0], ref_apl[0], rtol=RTOL)
    assert_allclose(apl[1][0], ref_apl[1], rtol=RTOL)
    assert_allclose(apl[2][0], ref_apl[2], rtol=RTOL)
    assert_allclose(apl[1][3], ref_area[0], rtol=RTOL)
    assert_allclose(apl[2][3], ref_area[1], rtol=RTOL)


def test_apl_model_bilayer_prot(frame_model_bilayer_prot):
    frame = frame_model_bilayer_prot
    membrane = frame.get_membranes()[0]

    ref_apl = (0.653, 0.653, 0.653)
    ref_area = (156.627, 156.707)

    apl = membrane.get_apl(by_type=False)

    print("APL: avg:%.3f - L1 (%i lipids): APL=%.3f (min:%.3f - max:%.3f), area=%.3f - L2 (%i lipids): APL=%.3f (min:%.3f - max:%.3f), area=%.3f " % (
        apl[0],
        len(membrane.leaflet1), apl[1][0], apl[1][1], apl[1][2], apl[1][3],
        len(membrane.leaflet2), apl[2][0], apl[2][1], apl[2][2], apl[2][3]))

    assert_allclose(apl[0], ref_apl[0], rtol=RTOL)
    assert_allclose(apl[1][0], ref_apl[1], rtol=RTOL)
    assert_allclose(apl[2][0], ref_apl[2], rtol=RTOL)
    assert_allclose(apl[1][3], ref_area[0], rtol=RTOL)
    assert_allclose(apl[2][3], ref_area[1], rtol=RTOL)


def test_apl_bilayer_prot(frame_bilayer_prot):
    frame = frame_bilayer_prot
    membrane = frame.get_membranes()[0]

    ref_apl = (0.689, 0.672, 0.707)
    ref_area = (37.640, 39.564)

    apl = membrane.get_apl(by_type=False)

    print("APL: avg:%.3f - L1 (%i lipids): APL=%.3f (min:%.3f - max:%.3f), area=%.3f - L2 (%i lipids): APL=%.3f (min:%.3f - max:%.3f), area=%.3f " % (
        apl[0],
        len(membrane.leaflet1), apl[1][0], apl[1][1], apl[1][2], apl[1][3],
        len(membrane.leaflet2), apl[2][0], apl[2][1], apl[2][2], apl[2][3]))

    assert_allclose(apl[0], ref_apl[0], rtol=RTOL)
    assert_allclose(apl[1][0], ref_apl[1], rtol=RTOL)
    assert_allclose(apl[2][0], ref_apl[2], rtol=RTOL)
    assert_allclose(apl[1][3], ref_area[0], rtol=RTOL)
    assert_allclose(apl[2][3], ref_area[1], rtol=RTOL)
