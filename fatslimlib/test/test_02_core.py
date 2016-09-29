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
from . import frame_vesicle, frame_model_bilayer_prot, frame_model_bilayer, frame_model_vesicle, frame_bilayer, \
    frame_bilayer_prot, frame_big_prot, frame_bilayer_chol, frame_model_bilayer_vesicle, frame_model_multibilayer, \
    traj_vesicle, frame_bilayer_allatom, frame_multibilayer, frame_bilayer_ganglio


def dprod(v1, v2):
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]


def test_openmp():
    from ..core_base import set_num_threads, get_max_threads
    from ..core_analysis import test_parallelism

    nthreads = min(get_max_threads(), 8)

    for i in range(nthreads):
        set_num_threads(i)
        assert test_parallelism()


def test_forcefield_identification_allatom(frame_bilayer_allatom):
    from ..core_base import FORCEFIELDTYPE_ALLATOM
    assert frame_bilayer_allatom.get_forcefield_type() == FORCEFIELDTYPE_ALLATOM


def test_forcefield_identification_unified(frame_bilayer_prot):
    from ..core_base import FORCEFIELDTYPE_UNIFIED
    assert frame_bilayer_prot.get_forcefield_type() == FORCEFIELDTYPE_UNIFIED


def test_forcefield_identification_coarse(frame_bilayer):
    from ..core_base import FORCEFIELDTYPE_COARSE
    assert frame_bilayer.get_forcefield_type() == FORCEFIELDTYPE_COARSE


def test_frame_size(frame_model_bilayer):
    frame = frame_model_bilayer

    assert frame.size == len(frame.bead_coords)


def test_core_base_bbox(frame_vesicle):
    frame = frame_vesicle
    coords = frame.hg_group_coords
    coords_bbox = frame.box.put_atoms_in_bbox(frame.hg_group_coords)

    # All the atoms are already inside the bbox
    assert_almost_equal(coords, coords_bbox)


def test_core_base_interacting_bbox(frame_model_bilayer_prot):
    frame = frame_model_bilayer_prot
    coords = frame.interacting_group_coords_bbox
    coords_bbox = frame.box.put_atoms_in_bbox(frame.interacting_group_coords)

    # All the atoms are already inside the bbox
    assert_almost_equal(coords, coords_bbox)


def test_core_base_directions(frame_model_bilayer):
    frame = frame_model_bilayer

    directions = frame.directions

    ref_direction_first = np.array([0, 0, 1.0])
    ref_direction_last = np.array([0, 0, -1.0])

    assert_almost_equal(directions[0], ref_direction_first, 1)
    assert_almost_equal(directions[-1], ref_direction_last, 1)


def test_core_ns_bilayer(frame_model_bilayer):
    from ..core_ns import neighbor_search

    test_bead = 0
    test_result = sorted([35, 30, 31, 5, 1, 11, 6, 7])

    frame = frame_model_bilayer
    coords = frame.bead_coords

    neighbors = neighbor_search(frame.box, coords, cutoff=1.5)

    assert test_result == sorted(neighbors[test_bead])


def test_core_ns_consistency_bilayer(frame_model_bilayer):
    from ..core_ns import neighbor_search

    frame = frame_model_bilayer
    coords = frame.bead_coords

    neighbors = neighbor_search(frame.box, coords, cutoff=2.0)

    for neighborhood in neighbors:
        assert len(neighborhood) == 20


def test_core_ns_vesicle(frame_vesicle):
    from ..core_ns import neighbor_search

    test_bead = 0
    test_result = sorted([624, 655, 780, 1743, 2104, 1981, 2714, 638, 1352])

    frame = frame_vesicle
    coords = frame.hg_group_coords

    neighbors = neighbor_search(frame.box, coords, cutoff=2)

    assert test_result == sorted(neighbors[test_bead])


def test_core_dx(frame_model_bilayer):
    ref_dx = np.array([0.0, -0.002, 4.336])

    frame = frame_model_bilayer
    coords = frame.bead_coords

    for i in range(36):
        l1_bead = i
        l2_bead = 35 + (i // 6 + 1) * 6 - (i % 6)
        dx = frame.box.dx(coords[l1_bead], coords[l2_bead])

        assert_almost_equal(dx, ref_dx)


def test_directions_bilayer(frame_model_bilayer):
    frame = frame_model_bilayer

    directions = frame.directions

    ref_direction_first = np.array([0, 0, 1.0])
    ref_direction_last = np.array([0, 0, -1.0])

    for i in range(36):
        assert_almost_equal(directions[i], ref_direction_first, 1)
        assert_almost_equal(directions[i + 36], ref_direction_last, 1)


def test_direction_ganglio(frame_bilayer_ganglio):
    frame = frame_bilayer_ganglio

    ref_ganglio = 812
    ref_chol = 27

    assert dprod(frame.directions[ref_chol], [0.0, 0.0, 1.0]) > 0
    assert dprod(frame.directions[ref_ganglio], [0.0, 0.0, 1.0]) > 0
    assert dprod(frame.directions[ref_ganglio], frame.directions[ref_ganglio]) > 0


def test_directions_vesicle(frame_model_vesicle):
    frame = frame_model_vesicle

    directions = frame.directions
    coords = frame.bead_coords
    center = frame.box.center

    for i in range(len(directions)):

        # Get the distance between the lipid and the center
        dx = coords[i] - center
        dx_norm = np.sqrt(np.sum(dx * dx))

        # Change the sign if we are dealing with a lipid from the inner leaflet
        if dx_norm < 7:
            dx *= -1

        dx /= dx_norm

        assert dprod(directions[i], dx) > 0.99


def test_normals_model_bilayer(frame_model_bilayer):
    frame = frame_model_bilayer

    normals = frame.normals

    ref_direction_first = np.array([0, 0, 1.0])
    ref_direction_last = np.array([0, 0, -1.0])

    for i in range(36):
        assert_almost_equal(normals[i], ref_direction_first)
        assert_almost_equal(normals[i + 36], ref_direction_last)


def test_normals_bilayer(frame_bilayer):
    frame = frame_bilayer

    normals = frame.normals

    coords = frame.bead_coords
    center_z = 3.5

    ref_direction_up = np.array([0, 0, 1.0])
    ref_direction_down = np.array([0, 0, -1.0])

    for i, normal in enumerate(normals):
        if coords[i][2] > center_z:
            assert dprod(normal, ref_direction_up) > 0.97
        else:
            assert dprod(normal, ref_direction_down) > 0.97


def test_normals_bilayer_prot(frame_bilayer_prot):
    frame = frame_bilayer_prot

    normals = frame.normals
    coords = frame.bead_coords
    center_z = 3.5

    ref_direction_up = np.array([0, 0, 1.0])
    ref_direction_down = np.array([0, 0, -1.0])

    for i, normal in enumerate(normals):
        if coords[i][2] > center_z:
            assert dprod(normal, ref_direction_up) > 0.95
        else:
            assert dprod(normal, ref_direction_down) > 0.95


def test_normals_vesicle_model(frame_model_vesicle):
    frame = frame_model_vesicle

    normals = frame.normals
    coords = frame.bead_coords
    center = frame.box.center

    for i in range(len(normals)):

        # Get the distance between the lipid and the center
        dx = coords[i] - center
        dx_norm = np.sqrt(np.sum(dx * dx))

        # Change the sign if we are dealing with a lipid from the inner leaflet
        if dx_norm < 7:
            dx *= -1

        dx /= dx_norm

        assert dprod(normals[i], dx) > 0.99


def test_normals_big_prot(frame_big_prot):
    test_bead = 6876
    test_neighbors = [21465, 6975, 8575, 8994, 12676, 8731, 8803, 10220, 10573, 12151, 12363, 12366,
                      22388, 6996, 8673, 8786, 8842, 8990, 12283, 12403, 22284, 10471, 21441]
    test_bad = [21441, 10471, 8731]
    frame = frame_big_prot

    normals = frame.normals

    for nid in test_neighbors:
        dp_val = dprod(normals[test_bead], normals[nid])

        if nid in test_bad:
            assert dp_val < 0
        else:
            assert dp_val > 0.90


def test_core_dx_leaflet(frame_model_bilayer):
    ref_dx = np.array([0.000, -0.002, -5.664])

    frame = frame_model_bilayer
    coords = frame.bead_coords
    normals = frame.normals

    for i in range(36):
        l1_bead = i
        l2_bead = 35 + (i // 6 + 1) * 6 - (i % 6)

        dx = frame.box.dx_leaflet(coords[l1_bead], coords[l2_bead], normals[l1_bead])
        assert_almost_equal(dx, ref_dx)

        dx = frame.box.dx_leaflet(coords[l2_bead], coords[l1_bead], normals[l2_bead])
        assert_almost_equal(dx, -ref_dx)


def test_aggregate_model_bilayer(frame_model_bilayer):
    ref_aggregates = [np.arange(36), np.arange(36, 72)]
    frame = frame_model_bilayer
    aggregates = frame.get_aggregates()

    assert len(aggregates) == len(ref_aggregates)

    for i, aggregate in enumerate(aggregates):
        assert_almost_equal(aggregate.beadids, ref_aggregates[i])
        assert_allclose(aggregate.xcm, aggregate.coords.mean(axis=0), 1e-2)

    assert len(frame.get_aggregates(0.5)) == 72


def test_aggregate_bilayer(frame_bilayer):
    ref_aggregates = [169, 169]

    frame = frame_bilayer
    aggregates = frame.get_aggregates()

    assert len(aggregates) == len(ref_aggregates)

    for i, aggregate in enumerate(aggregates):
        assert len(aggregate) == ref_aggregates[i]


def test_aggregate_big_prot(frame_big_prot):
    ref_aggregates = [12152, 11861, 42, 1]
    ref_xcm = [
        [43.19, 43.06, 9.95],
        [43.24, 43.20, 14.27]
    ]
    expected_normal = [
        [0, 0, -1],
        [0, 0, 1]
    ]

    frame = frame_big_prot
    aggregates = frame.get_aggregates()

    assert len(aggregates) == len(ref_aggregates)

    for i, aggregate in enumerate(aggregates):
        assert len(aggregate) == ref_aggregates[i]

        if i < len(ref_xcm):
            assert_allclose(aggregate.xcm, ref_xcm[i], rtol=0.001)

            dp_val = dprod(aggregate.avg_normal, expected_normal[i])

            assert dp_val > 0.90


def test_aggregate_chol(frame_bilayer_chol):
    ref_aggregates = [977, 967]

    frame = frame_bilayer_chol
    aggregates = frame.get_aggregates()

    assert len(aggregates) == len(ref_aggregates)

    for i, aggregate in enumerate(aggregates):
        assert len(aggregate) == ref_aggregates[i]


def test_aggregate_vesicle_model(frame_model_vesicle):
    frame = frame_model_vesicle
    aggregates = frame.get_aggregates()

    ref_aggregates = [1963, 785]

    assert len(aggregates) == len(ref_aggregates)

    for i, aggregate in enumerate(aggregates):
        assert len(aggregate) == ref_aggregates[i]


def test_aggregate_vesicle(frame_vesicle):
    frame = frame_vesicle
    aggregates = frame.get_aggregates()

    ref_aggregates = [1851, 1179, 26, 3, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1]

    assert len(aggregates) == len(ref_aggregates)

    for i, aggregate in enumerate(aggregates):
        assert len(aggregate) == ref_aggregates[i]
        if len(aggregates) > 100:
            assert_allclose(aggregate.xcm, (26.9,  36.7,  15.4), 1e-2)


def test_aggregate_vesicle_traj(traj_vesicle):
    nlipids = len(traj_vesicle[0].bead_coords)

    for frame in traj_vesicle:
        aggregates = frame.get_aggregates()
        assert sum([len(val) for val in aggregates]) == nlipids


def test_aggregate_atomid(frame_model_bilayer):
    ref_atomids = (4, 1804)

    frame = frame_model_bilayer
    aggregates = frame.get_aggregates()

    assert aggregates[0].get_atomid(0) == ref_atomids[0]
    assert aggregates[1].get_atomid(0) == ref_atomids[1]


def test_aggregate_residue(frame_model_bilayer):
    ref_atomid = 1804

    frame = frame_model_bilayer
    aggregates = frame.get_aggregates()
    ref_residue = frame.trajectory.get_topology().get_residue_from_atomid(ref_atomid)

    for key, val in aggregates[1].get_residue(0).items():
        try:
            assert ref_residue[key] == val
        except ValueError:
            assert_almost_equal(ref_residue[key], val)


def test_aggregate_resname(frame_model_bilayer):
    ref_atomid = 1804
    ref_resname = "DPPC"

    frame = frame_model_bilayer
    aggregates = frame.get_aggregates()

    assert aggregates[1].get_resname(0) == ref_resname
    assert aggregates[1].get_resname(0) == frame.trajectory.\
        get_topology().get_residue_from_atomid(ref_atomid)["resname"]


def test_aggregate_same_restype(frame_bilayer_chol):
    frame = frame_bilayer_chol
    aggregates = frame.get_aggregates()

    assert aggregates[0].same_restype(0, 1)
    assert not aggregates[0].same_restype(0, len(aggregates[0])-1)


def test_membranes_model_bilayer(frame_model_bilayer):
    frame = frame_model_bilayer
    nlipids = frame.size

    ref_membrane_beadids = [np.arange(36, 72),
                            np.arange(36)]
    membranes = frame.get_membranes()

    assert len(membranes) == 1
    membrane = membranes[0]

    assert membrane.is_planar

    assert sum([len(val) for val in membrane]) == nlipids

    for i, beadids in enumerate(membrane.beadids):
        assert_almost_equal(beadids, ref_membrane_beadids[i])


def test_membranes_bilayer(frame_bilayer):
    ref_leaflets = [169, 169]

    frame = frame_bilayer
    membranes = frame.get_membranes()

    assert len(membranes) == 1

    membrane = membranes[0]

    assert membrane.is_planar

    for i, leaflet in enumerate(membrane):
        assert len(leaflet) == ref_leaflets[i]


def test_membranes_bilayer_normals(frame_bilayer):
    frame = frame_bilayer

    membrane = frame.get_membranes()[0]

    center_z = 3.5

    ref_direction_up = np.array([0, 0, 1.0])
    ref_direction_down = np.array([0, 0, -1.0])

    for leaflet in membrane:
        coords = leaflet.coords
        normals = leaflet.normals

        for i, normal in enumerate(normals):
            if coords[i][2] > center_z:
                assert dprod(normal, ref_direction_up) > 0.985  # less than 10° off
            else:
                assert dprod(normal, ref_direction_down) > 0.985  # less than 10° off


def test_membranes_big_prot(frame_big_prot):
    ref_leaflets = [12152, 11861]

    frame = frame_big_prot
    membranes = frame.get_membranes()

    assert len(membranes) == 1

    membrane = membranes[0]

    assert membrane.is_planar

    for i, leaflet in enumerate(membrane):
        assert len(leaflet) == ref_leaflets[i]


def test_membranes_chol(frame_bilayer_chol):

    ref_leaflets = [977, 967]

    frame = frame_bilayer_chol
    membranes = frame.get_membranes()

    assert len(membranes) == 1

    membrane = membranes[0]

    assert membrane.is_planar

    for i, leaflet in enumerate(membrane):
        assert len(leaflet) == ref_leaflets[i]


def test_membranes_vesicle_model(frame_model_vesicle):
    ref_leaflets = [1963, 785]

    frame = frame_model_vesicle
    membranes = frame.get_membranes()

    assert len(membranes) == 1

    membrane = membranes[0]

    assert membrane.is_vesicle

    for i, leaflet in enumerate(membrane):
        assert len(leaflet) == ref_leaflets[i]


def test_membranes_vesicle_model_normals(frame_model_vesicle):
    frame = frame_model_vesicle
    membrane = frame.get_membranes()[0]

    for lid, leaflet in enumerate(membrane):
        coords = leaflet.coords
        normals = leaflet.normals
        xcm = leaflet.xcm

        for i, normal in enumerate(normals):
            ref_direction = coords[i] - xcm
            ref_direction /= np.linalg.norm(ref_direction)
            if lid == 1:
                ref_direction *= -1.0
            assert dprod(normal, ref_direction) > 0.998  # less than 4° off


def test_membranes_vesicle(frame_vesicle):
    ref_leaflets = [1851, 1179]

    frame = frame_vesicle
    membranes = frame.get_membranes()

    assert len(membranes) == 1

    membrane = membranes[0]

    assert membrane.is_vesicle

    for i, leaflet in enumerate(membrane):
        assert len(leaflet) == ref_leaflets[i]


def test_membranes_vesicle_traj(traj_vesicle):
    for frame in traj_vesicle:
        membranes = frame.get_membranes()
        assert len(membranes) == 1
        assert membranes[0].is_vesicle


def test_membranes_multibilayer_model(frame_model_multibilayer):

    frame = frame_model_multibilayer
    nlipids = 72

    ref_membrane_beadids = [[np.arange(36, 72),
                            np.arange(36)],
                            [np.arange(72, 108),
                                np.arange(108, 144)]]
    membranes = frame.get_membranes()

    assert len(membranes) == 2

    for mid, membrane in enumerate(membranes):
        assert membrane.is_planar

        assert sum([len(val) for val in membrane]) == nlipids

        for i, beadids in enumerate(membrane.beadids):
            assert_almost_equal(beadids, ref_membrane_beadids[mid][i])


def test_membranes_multibilayer(frame_multibilayer):

    frame = frame_multibilayer
    nlipids = 256

    ref_membrane_beadids_firsts = [
        [
            [198, 240, 254, 366, 506],
            [30, 44, 58, 86, 142]
        ],
        [
            [2, 16, 72, 100, 114],
            [688, 1682, 2022, 2036, 2050]
        ]
    ]

    membranes = frame.get_membranes()

    assert len(membranes) == 2

    for mid, membrane in enumerate(membranes):
        assert membrane.is_planar

        assert sum([len(val) for val in membrane]) == nlipids

        for i, leaflet in enumerate(membrane):
            assert_almost_equal(leaflet.hg_atomids[:5], ref_membrane_beadids_firsts[mid][i])


def test_membranes_bilayer_vesicle(frame_model_bilayer_vesicle):
    frame = frame_model_bilayer_vesicle

    ref_leaflets = [[1963, 785], [1296, 1296]]
    membrane_planar = [False, True]

    membranes = frame.get_membranes()

    assert len(membranes) == 2

    for mid, membrane in enumerate(membranes):

        assert membrane.is_planar == membrane_planar[mid]

        assert membrane.is_vesicle != membrane_planar[mid]

        for i, leaflet in enumerate(membrane):
            assert len(leaflet) == ref_leaflets[mid][i]


def test_membranes_bilayer_ganglio(frame_bilayer_ganglio):
    frame = frame_bilayer_ganglio

    membrane = frame.get_membranes()[0]

    assert len(membrane) == len(frame.bead_coords)
