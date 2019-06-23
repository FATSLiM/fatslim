# -*- coding: utf-8; Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
# 
# This file is part of FATSLiM --- http://fatslim.github.io/
# 
# Copyright (c) 2013-2018, Sébastien Buchoux
# Copyright (c) 2019, by the FATSLiM development team (see AUTHORS file)
#
# FATSLiM is free software and is released under the GNU Public Licence,
# v3 or any higher version
#
# If you use FATSLiM in publised work, please cite:
#
# S. Buchoux.
# FATSLiM: a fast and robust software to analyze MD simulations of membranes.
# Bioinformatics 33(1) (2017), 133--134, doi:10.1093/bioinformatics/btw563
#

import numpy as np
from numpy.testing import assert_allclose, assert_equal, assert_almost_equal
import pytest

from .fixtures import *
from .data import MODELS_METADATA

# MD Analysis specific warning
pytestmark = pytest.mark.filterwarnings("ignore: Failed to guess the mass for the following atom")


def test_membrane_model_flat(system_model_flat):
    system = system_model_flat

    expected_leaflets = [
        MODELS_METADATA["flat"]["upper_leaflet_ids"],
        MODELS_METADATA["flat"]["lower_leaflet_ids"]
    ]

    assert len(system.membranes) == 1

    membrane = system.membranes[0]
    for i, leaflet in enumerate(membrane):
        assert_allclose(leaflet.indices, expected_leaflets[i], err_msg="Bad lipid ids for leaflet #{}".format(i))

        assert leaflet.is_planar


def test_membrane_model_vesicle(system_model_vesicle):
    system = system_model_vesicle

    expected_leaflets = [
        MODELS_METADATA["vesicle"]["outer_leaflet_ids"],
        MODELS_METADATA["vesicle"]["inner_leaflet_ids"]
    ]

    assert len(system.membranes) == 1

    membrane = system.membranes[0]
    for i, leaflet in enumerate(membrane):
        assert_allclose(leaflet.indices, expected_leaflets[i], err_msg="Bad lipid ids for leaflet #{}".format(i))

        assert not leaflet.is_planar


def test_membrane_model_curved(system_model_curved):
    system = system_model_curved

    expected_leaflets = [
        MODELS_METADATA["curved"]["upper_leaflet_ids"],
        MODELS_METADATA["curved"]["lower_leaflet_ids"]
    ]

    assert len(system.membranes) == 1

    membrane = system.membranes[0]
    for i, leaflet in enumerate(membrane):
        assert_allclose(leaflet.indices, expected_leaflets[i], err_msg="Bad lipid ids for leaflet #{}".format(i))

        assert leaflet.is_planar


def test_membrane_model_bulged(system_model_bulged):
    system = system_model_bulged

    expected_leaflets = [
        MODELS_METADATA["bulged"]["upper_leaflet_ids"],
        MODELS_METADATA["bulged"]["lower_leaflet_ids"]
    ]

    assert len(system.membranes) == 1

    membrane = system.membranes[0]
    for i, leaflet in enumerate(membrane):
        assert_allclose(leaflet.indices, expected_leaflets[i], err_msg="Bad lipid ids for leaflet #{}".format(i))

        assert leaflet.is_planar


def test_membrane_model_bicelle(system_model_bicelle):
    system = system_model_bicelle

    expected_leaflets = [
        np.array([val for val in MODELS_METADATA["bicelle"]["upper_leaflet_ids"] if val not in MODELS_METADATA["bicelle"]["leaflet_pivot_ids"]]),
        np.array([val for val in MODELS_METADATA["bicelle"]["lower_leaflet_ids"] if val not in MODELS_METADATA["bicelle"]["leaflet_pivot_ids"]]),
    ]

    assert len(system.membranes) == 1

    membrane = system.membranes[0]

    for i, leaflet in enumerate(membrane):
        indices = np.array([val for val in leaflet.indices if val not in MODELS_METADATA["bicelle"]["leaflet_pivot_ids"]])
        assert_allclose(indices, expected_leaflets[i], err_msg="Bad lipid ids for leaflet #{}".format(i))

        assert leaflet.is_planar


def test_membrane_big_deformed(system_big_deformed):
    system = system_big_deformed
    expected_sizes = [11861, 12152]

    assert len(system.membranes) == 1

    membrane = system.membranes[0]
    for i, leaflet in enumerate(membrane):
        assert leaflet.size == expected_sizes[i]
        assert leaflet.is_planar


def test_membrane_ganglio(system_ganglio):
    system = system_ganglio

    expected_sizes = [1703, 1712]

    assert len(system.membranes) == 1

    membrane = system.membranes[0]
    for i, leaflet in enumerate(membrane):
        assert leaflet.size == expected_sizes[i]
        assert leaflet.is_planar
