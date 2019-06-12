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
from numpy.testing import assert_allclose, assert_equal

from .fixtures import *
from .data import MODELS_METADATA


def test_thickness_model_flat(system_model_flat):
    system = system_model_flat

    metadata = MODELS_METADATA["flat"]

    expected_apl = metadata["apl"]

    expected_area = expected_apl * metadata["size"] * 0.5

    expected_lipid_apls = [
        np.ones(len(metadata["upper_leaflet_ids"])) * expected_apl,
        np.ones(len(metadata["lower_leaflet_ids"])) * expected_apl,
    ]

    membrane = system.membranes[0]

    for i, leaflet in enumerate(membrane):

        apl = leaflet.apl

        area = leaflet.area

        lipid_apls = leaflet.lipid_apls

        assert_allclose(apl, expected_apl, atol=1e-3,
                        err_msg="Bad APL value for leaflet #{}".format(i))

        assert_allclose(area, expected_area, atol=1e-3,
                        err_msg="Bad area value for leaflet #{}".format(i))

        assert_allclose(lipid_apls, expected_lipid_apls[i], atol=1e-3)


def test_thickness_model_vesicle(system_model_vesicle):
    system = system_model_vesicle
    n_ignored = 6

    metadata = MODELS_METADATA["vesicle"]

    expected_apls = [
        metadata["outer_apl"],
        metadata["inner_apl"]
        ]

    expected_areas = [
        metadata["outer_apl"] * metadata["n_outer"],
        metadata["inner_apl"] * metadata["n_inner"]
    ]

    expected_lipid_apls = [
        np.ones(metadata["n_outer"]) * expected_apls[0],
        np.ones(metadata["n_inner"]) * expected_apls[1],
    ]

    membrane = system.membranes[0]

    for i, leaflet in enumerate(membrane):

        apl = leaflet.apl
        area = leaflet.area
        lipid_apls = leaflet.lipid_apls

        for j, val in enumerate(lipid_apls):
            print("resid {}: APL={:.3f} vs {:.3f} - diff: {:.3f}".format(
                system.lipids[leaflet.indices[j]].resid,
                val,
                expected_lipid_apls[i][j],
                abs(val - expected_lipid_apls[i][j]) / expected_lipid_apls[i][j]
            ))

        #  the first and last points are not as homogeneously located as they should so their APL further from the
        #  average value than tolerated. This is not due to the algorithm failing so they are just ignored

        assert_allclose(apl, expected_apls[i], rtol=5e-2,
                        err_msg="Bad APL value for leaflet #{}".format(i))

        assert_allclose(area, expected_areas[i], rtol=5e-2,
                        err_msg="Bad area value for leaflet #{}".format(i))

        assert_allclose(lipid_apls[n_ignored:-n_ignored], expected_lipid_apls[i][n_ignored:-n_ignored], rtol=5e-2)