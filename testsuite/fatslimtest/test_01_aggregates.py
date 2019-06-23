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


def test_aggregate_model_flat(system_model_flat):
    system = system_model_flat

    aggregates = system.aggregates

    expected_aggregate_ids = [
        MODELS_METADATA["flat"]["upper_leaflet_ids"],
        MODELS_METADATA["flat"]["lower_leaflet_ids"]
    ]

    assert len(aggregates) == len(expected_aggregate_ids)

    for i, aggregate in enumerate(aggregates):
        assert_equal(aggregate.indices, expected_aggregate_ids[i], err_msg="Bad lipid ids for aggregate #{}".format(i))


def test_aggregate_neighbours(system_model_flat):
    system = system_model_flat

    for aggid, aggregate in enumerate(system.aggregates):
        tuples = aggregate.lipid_neighbours.tuples

        for i in range(aggregate.size):
            beadid = aggregate.indices[i]

            tuple = []
            for j, d in tuples[i]:
                tuple.append((j + aggregate.indices[0], d))

            assert_almost_equal(system[beadid].neighbours_distances,
                                tuple, decimal=4,
                                err_msg="Bad neigbours for lipid #{} from aggregate #{} (beadid #{})".format(i,
                                                                                                             aggid,
                                                                                                             beadid
                                                                                                             ))


def test_aggregate_positions_raw(system_model_flat):
    system = system_model_flat

    for aggid, aggregate in enumerate(system.aggregates):

        for i in range(aggregate.size):
            beadid = aggregate.indices[i]

            assert_almost_equal(system.lipid_positions[beadid],
                                aggregate.lipid_positions_raw[i], decimal=4,
                                err_msg="Bad raw positions for lipid #{} from aggregate #{} (beadid #{})".format(i,
                                                                                                                 aggid,
                                                                                                                 beadid
                                                                                                                 ))


def test_aggregate_directions(system_model_flat):
    system = system_model_flat

    for aggid, aggregate in enumerate(system.aggregates):

        for i in range(aggregate.size):
            beadid = aggregate.indices[i]

            assert_almost_equal(system.lipid_directions[beadid],
                                aggregate.lipid_directions[i], decimal=4,
                                err_msg="Bad direction for lipid #{} from aggregate #{} (beadid #{})".format(i,
                                                                                                             aggid,
                                                                                                             beadid
                                                                                                             ))


def test_aggregate_normals(system_model_flat):
    system = system_model_flat

    for aggid, aggregate in enumerate(system.aggregates):

        for i in range(aggregate.size):
            beadid = aggregate.indices[i]

            assert_almost_equal(system.lipid_normals[beadid],
                                aggregate.lipid_normals[i], decimal=4,
                                err_msg="Bad normal for lipid #{} from aggregate #{} (beadid #{})".format(i,
                                                                                                             aggid,
                                                                                                             beadid
                                                                                                             ))


def test_aggregate_model_vesicle(system_model_vesicle):
    system = system_model_vesicle

    aggregates = system.aggregates

    expected_aggregate_ids = [
        MODELS_METADATA["vesicle"]["outer_leaflet_ids"],
        MODELS_METADATA["vesicle"]["inner_leaflet_ids"]
    ]

    assert len(aggregates) == len(expected_aggregate_ids)

    for i, aggregate in enumerate(aggregates):
        assert_equal(aggregate.indices, expected_aggregate_ids[i], err_msg="Bad lipid ids for aggregate #{}".format(i))


def test_aggregate_model_bicelle(system_model_bicelle):
    system = system_model_bicelle

    aggregates = system.aggregates

    assert len(aggregates) == 1
    assert len(aggregates[0]) == len(system)


def test_aggregate_vesicle(system_vesicle):
    system = system_vesicle

    # Make sure we reset the trajectory by loading the first frame
    system.trajectory[0]

    aggregates = system.aggregates

    ref_aggregates = [
        1853,
        1179,
        26,
        5,
        3,
        2,
        1,
        1,
        1,
        1
    ]

    assert len(aggregates) == len(ref_aggregates)
    for i, aggregate in enumerate(aggregates):

        assert len(aggregate) == ref_aggregates[i]


def test_clusterize_vesicle(system_vesicle):
    system = system_vesicle

    # List of centroids for each frame
    # For each frame:
    # - First element should correspond to the centroid after clusterization
    # - Second element should correspond to the centroid based on raw coordinates (ie directly from trajectory)
    centroids = np.array(
        [[[[437.875, 342.175, 165.593],
           [277.994, 341.665, 165.593]],

          [[437.989, 341.735, 165.301],
           [328.129, 341.735, 165.301]]],

         [[[-32.889, 342.24, 165.725],
           [276.342, 341.729, 165.725]],

          [[437.889, 341.911, 165.064],
           [331.661, 341.911, 165.064]]],

         [[[437.88, 342.355, 165.621],
           [275.584, 342.355, 165.621]],

          [[438.132, 341.68, 165.52],
           [328.703, 341.68, 165.52]]],

         [[[-32.827, 342.051, 165.458],
           [277.421, 342.051, 165.458]],

          [[437.867, 342.061, 165.697],
           [329.272, 342.061, 165.697]]],

         [[[-32.631, 341.965, 165.436],
           [276.147, 341.965, 165.436]],

          [[437.887, 342.196, 165.811],
           [324.48, 342.196, 165.811]]],

         [[[-32.975, 342.182, 165.705],
           [279.573, 342.182, 165.705]],

          [[437.932, 341.876, 165.368],
           [324.941, 341.876, 165.368]]],

         [[[-32.657, 342.473, 165.52],
           [275.778, 342.473, 165.52]],

          [[437.788, 341.689, 165.69],
           [328.805, 341.689, 165.69]]],

         [[[-33.011, 342.503, 165.602],
           [275.403, 342.503, 165.602]],

          [[438.127, 341.968, 165.611],
           [322.765, 341.968, 165.611]]],

         [[[-33.034, 342.541, 165.662],
           [277.243, 342.541, 165.662]],

          [[438.179, 341.609, 165.426],
           [331.171, 341.609, 165.426]]],

         [[[437.797, 342.426, 165.718],
           [275.544, 342.426, 165.718]],

          [[437.689, 341.731, 165.558],
           [334.279, 341.731, 165.558]]],

         [[[437.904, 342.148, 165.651],
           [274.397, 342.148, 165.651]],

          [[437.742, 341.994, 165.523],
           [327.955, 341.994, 165.523]]]]
    )

    sizes = np.array(
        [[1853, 1179],
         [1853, 1179],
         [1851, 1179],
         [1851, 1179],
         [1851, 1179],
         [1851, 1179],
         [1851, 1179],
         [1851, 1179],
         [1851, 1179],
         [1851, 1179],
         [1851, 1179]]
    )

    for i, ts in enumerate(system.trajectory):

        aggregates = system.aggregates

        # We only check the first two aggregates (ie the biggest ones)
        for j, aggregate in enumerate(aggregates[:2]):
            assert len(aggregate) == sizes[i, j]

            positions = aggregate.lipid_positions
            assert_almost_equal(centroids[i, j, 0], positions.mean(axis=0), decimal=3)
            assert_almost_equal(centroids[i, j, 0], aggregate.position, decimal=3)

            positions = aggregate.lipid_positions_raw
            assert_almost_equal(centroids[i, j, 1], positions.mean(axis=0), decimal=3)


def test_aggregate_big_deformed(system_big_deformed):
    system = system_big_deformed
    expected_sizes = [12152, 11903, 1]

    aggregates = system.aggregates

    assert len(aggregates) == len(expected_sizes)

    for i, aggregate in enumerate(aggregates):
        assert aggregate.size == expected_sizes[i]
