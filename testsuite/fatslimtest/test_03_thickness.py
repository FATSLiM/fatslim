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

import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_equal

from .data import MODELS_METADATA

# MD Analysis specific warning
pytestmark = pytest.mark.filterwarnings("ignore: Failed to guess the mass for the following atom")

def test_thickness_model_flat(system_model_flat):
    system = system_model_flat

    metadata = MODELS_METADATA["flat"]

    expected_thickness = metadata["thickness"]

    expected_lipid_thicknesses = [
        np.ones(len(metadata["upper_leaflet_ids"])) * metadata["thickness"],
        np.ones(len(metadata["lower_leaflet_ids"])) * metadata["thickness"],
    ]

    membrane = system.membranes[0]

    for i, leaflet in enumerate(membrane):

        thickness = leaflet.thickness

        lipid_thicknesses = leaflet.lipid_thicknesses

        assert_allclose(thickness, expected_thickness, atol=1e-3,
                        err_msg="Bad thickness value for leaflet #{}".format(i))

        assert_allclose(lipid_thicknesses, expected_lipid_thicknesses[i], atol=1e-3)


def test_thickness_model_vesicle(system_model_vesicle):
    system = system_model_vesicle

    metadata = MODELS_METADATA["vesicle"]

    expected_thickness = metadata["thickness"]

    expected_lipid_thicknesses = [
        np.ones(len(metadata["outer_leaflet_ids"])) * metadata["thickness"],
        np.ones(len(metadata["inner_leaflet_ids"])) * metadata["thickness"],
    ]

    membrane = system.membranes[0]

    for i, leaflet in enumerate(membrane):

        thickness = leaflet.thickness

        lipid_thicknesses = leaflet.lipid_thicknesses

        assert_allclose(thickness, expected_thickness, rtol=2e-2,
                        err_msg="Bad thickness value for leaflet #{}".format(i))

        assert_allclose(lipid_thicknesses, expected_lipid_thicknesses[i], rtol=5e-2,
                        err_msg="Bad thicknesses values for leaflet #{}".format(i)
                        )


def test_thickness_model_curved(system_model_curved):
    system = system_model_curved

    metadata = MODELS_METADATA["curved"]

    expected_thickness = metadata["thickness"]

    expected_lipid_thicknesses = []

    lipid_ids = [
        metadata["upper_leaflet_ids"],
        metadata["lower_leaflet_ids"]
    ]
    positions = metadata["positions"]

    for i, vals in enumerate(lipid_ids):
        thicknesses = []

        for j in vals:
            dmin = 10000

            for k in lipid_ids[(i+1)%2]:
                d = np.linalg.norm(positions[j] - positions[k])

                if d < dmin:
                    dmin = d

            thicknesses.append(dmin)

        expected_lipid_thicknesses.append(np.array(thicknesses))

    membrane = system.membranes[0]

    for i, leaflet in enumerate(membrane):

        thickness = leaflet.thickness

        lipid_thicknesses = leaflet.lipid_thicknesses

        assert_allclose(thickness, expected_lipid_thicknesses[i].mean(), rtol=5e-2,
                        err_msg="Bad thickness value for leaflet #{}".format(i))

        assert_allclose(lipid_thicknesses, expected_lipid_thicknesses[i], rtol=10e-2,
                        err_msg="Bad thicknesses values for leaflet #{}".format(i)
                        )


def test_thickness_model_bulged(system_model_bulged):
    system = system_model_bulged

    metadata = MODELS_METADATA["bulged"]

    expected_lipid_thicknesses = []

    lipid_ids = [
        metadata["upper_leaflet_ids"],
        metadata["lower_leaflet_ids"]
    ]
    positions = metadata["positions"]

    for i, vals in enumerate(lipid_ids):
        thicknesses = []

        for j in vals:
            dmin = 10000

            best_id = -1

            for k in lipid_ids[(i + 1) % 2]:
                d = positions[j, :2] - positions[k, :2]
                d = np.linalg.norm(d)

                if d < dmin:
                    dmin = d
                    best_id = k

            thicknesses.append(np.abs(positions[j, 2] - positions[best_id, 2]))

        expected_lipid_thicknesses.append(np.array(thicknesses))

    membrane = system.membranes[0]

    for i, leaflet in enumerate(membrane):

        thickness = leaflet.thickness

        lipid_thicknesses = leaflet.lipid_thicknesses

        assert_allclose(thickness, expected_lipid_thicknesses[i].mean(), rtol=5e-2,
                        err_msg="Bad thickness value for leaflet #{}".format(i))

        assert_allclose(lipid_thicknesses, expected_lipid_thicknesses[i], rtol=10e-2,
                        err_msg="Bad thicknesses values for leaflet #{}".format(i)
                        )


def test_thickness_model_bicelle(system_model_bicelle):
    system = system_model_bicelle

    metadata = MODELS_METADATA["bicelle"]

    expected_thickness = metadata["thickness"]

    membrane = system.membranes[0]

    for i, leaflet in enumerate(membrane):

        indices = []
        for j, ind in enumerate(leaflet.indices):
            if system.lipids[ind].resname == "DPPC":
                indices.append(j)
        indices = np.array(indices, dtype=int)

        lipid_thicknesses = leaflet.lipid_thicknesses[indices]

        assert_allclose(lipid_thicknesses.mean(), expected_thickness, rtol=2e-2,
                        err_msg="Bad thickness value for leaflet #{}".format(i))

        assert_allclose(lipid_thicknesses, expected_thickness, rtol=5e-2,
                        err_msg="Bad thicknesses values for leaflet #{}".format(i)
                        )