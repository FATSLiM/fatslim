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

# Global imports
import numpy as np
from numpy.testing import assert_allclose

# Local imports
from . import frame_model_bilayer, frame_bilayer, frame_model_vesicle, frame_vesicle, \
    frame_bilayer_chol, frame_model_bilayer_prot, frame_bilayer_prot, frame_bilayer_peptide, \
    traj_vesicle

RTOL = 1e-2
RTOL_MAX = 2e-2
RTOL_MIN = 1e-3


def test_thickness_model_bilayer(frame_model_bilayer):
    frame = frame_model_bilayer
    membrane = frame.get_membranes()[0]

    thickness = membrane.get_thickness(only_average=False)

    print("Thickness: AVG:%.3f - L1:%.3f (min:%.3f, max:%.3f), L2:%.3f (min:%.3f, max:%.3f)" % (
        thickness[0],
        thickness[1].mean(), thickness[1].min(), thickness[1].max(),
        thickness[2].mean(), thickness[2].min(), thickness[2].max(),))

    assert_allclose(thickness[0], 5.664, rtol=RTOL_MIN)
    assert_allclose(thickness[1], np.array([5.664] * 36), rtol=RTOL_MIN)
    assert_allclose(thickness[1], thickness[2], rtol=RTOL_MIN)


def test_thickness_bilayer(frame_bilayer):
    frame = frame_bilayer
    membrane = frame.get_membranes()[0]

    thickness = membrane.get_thickness()

    print("Thickness: AVG:%.3f - L1:%.3f (min:%.3f, max:%.3f), L2:%.3f (min:%.3f, max:%.3f)" % (
        thickness[0],
        thickness[1][0], thickness[1][1], thickness[1][2],
        thickness[2][0], thickness[2][1], thickness[2][2]))

    assert_allclose(thickness[0], 3.948, rtol=RTOL)
    assert_allclose(thickness[1][0], thickness[2][0], rtol=RTOL)


def test_thickness_bilayer_chol(frame_bilayer_chol):
    frame = frame_bilayer_chol
    membrane = frame.get_membranes()[0]

    thickness = membrane.get_thickness()

    print("Thickness: AVG:%.3f - L1:%.3f (min:%.3f, max:%.3f), L2:%.3f (min:%.3f, max:%.3f)" % (
        thickness[0],
        thickness[1][0], thickness[1][1], thickness[1][2],
        thickness[2][0], thickness[2][1], thickness[2][2]))

    assert_allclose(thickness[0], 4.086, rtol=RTOL)
    assert_allclose(thickness[1][0], thickness[2][0], rtol=RTOL)


def test_thickness_model_bilayer_prot(frame_model_bilayer_prot):
    frame = frame_model_bilayer_prot
    membrane = frame.get_membranes()[0]

    thickness = membrane.get_thickness(only_average=False)

    print("Thickness: AVG:%.3f - L1:%.3f (min:%.3f, max:%.3f), L2:%.3f (min:%.3f, max:%.3f)" % (
        thickness[0],
        thickness[1].mean(), thickness[1].min(), thickness[1].max(),
        thickness[2].mean(), thickness[2].min(), thickness[2].max(),))

    assert_allclose(thickness[0], 5.664, rtol=RTOL)
    assert_allclose(thickness[1], np.array([5.664] * 240), rtol=RTOL)
    assert_allclose(thickness[1], thickness[2], rtol=RTOL)


def test_thickness_bilayer_prot(frame_bilayer_prot):
    frame = frame_bilayer_prot
    membrane = frame.get_membranes()[0]

    thickness = membrane.get_thickness()

    print("Thickness: AVG:%.3f - L1:%.3f (min:%.3f, max:%.3f), L2:%.3f (min:%.3f, max:%.3f)" % (
        thickness[0],
        thickness[1][0], thickness[1][1], thickness[1][2],
        thickness[2][0], thickness[2][1], thickness[2][2]))

    assert_allclose(thickness[0], 3.146, rtol=RTOL)
    assert_allclose(thickness[1][0], thickness[2][0], rtol=RTOL)


def test_thickness_bilayer_peptide(frame_bilayer_peptide):
    frame = frame_bilayer_peptide
    membrane = frame.get_membranes()[0]

    thickness = membrane.get_thickness()

    print("Thickness: AVG:%.3f - L1:%.3f (min:%.3f, max:%.3f), L2:%.3f (min:%.3f, max:%.3f)" % (
        thickness[0],
        thickness[1][0], thickness[1][1], thickness[1][2],
        thickness[2][0], thickness[2][1], thickness[2][2]))

    assert_allclose(thickness[0], 3.811, rtol=RTOL)
    assert_allclose(thickness[1][0], thickness[2][0], rtol=RTOL)


def test_thickness_vesicle_model(frame_model_vesicle):
    frame = frame_model_vesicle
    membrane = frame.get_membranes()[0]

    thickness = membrane.get_thickness()

    print("Thickness: AVG:%.3f - L1:%.3f (min:%.3f, max:%.3f), L2:%.3f (min:%.3f, max:%.3f)" % (
        thickness[0],
        thickness[1][0], thickness[1][1], thickness[1][2],
        thickness[2][0], thickness[2][1], thickness[2][2]))

    assert_allclose(thickness[0], 5.0, rtol=RTOL)
    assert_allclose(thickness[1][0], thickness[2][0], rtol=RTOL)


def test_thickness_vesicle(frame_vesicle):
    frame = frame_vesicle
    membrane = frame.get_membranes()[0]

    thickness = membrane.get_thickness()

    print("Thickness: AVG:%.3f - L1:%.3f (min:%.3f, max:%.3f), L2:%.3f (min:%.3f, max:%.3f)" % (
        thickness[0],
        thickness[1][0], thickness[1][1], thickness[1][2],
        thickness[2][0], thickness[2][1], thickness[2][2]))

    assert_allclose(thickness[0], 3.985, rtol=RTOL)
    assert_allclose(thickness[1][0], thickness[2][0], rtol=RTOL)


def test_thickness_vesicle_traj(traj_vesicle):
    avg_thicknesses = [3.985, 3.986, 3.989, 3.981, 3.985, 3.994, 4.004, 4.003, 4.009, 4.001, 3.996]
    for i, frame in enumerate(traj_vesicle):
        membrane = frame.get_membranes()[0]

        thickness = membrane.get_thickness()

        print("Thickness: AVG:%.3f - L1:%.3f (min:%.3f, max:%.3f), L2:%.3f (min:%.3f, max:%.3f)"
              % (
                  thickness[0],
                  thickness[1][0], thickness[1][1], thickness[1][2],
                  thickness[2][0], thickness[2][1], thickness[2][2]))

        assert_allclose(thickness[0], avg_thicknesses[i], rtol=RTOL)
        assert_allclose(thickness[1][0], thickness[2][0], rtol=RTOL_MAX)
