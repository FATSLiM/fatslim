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

# Global imports
import numpy as np
from numpy.testing import assert_almost_equal, assert_allclose
import pytest
from unittest import TestCase

# Local imports
from ..coreobjects import LipidSystem, Lipid
from . import universe_model_flat, system_model_flat, system_model_vesicle
from .data import MODEL_FLAT_BBOX_GRO

# TODO: Add test for neighbors when position is (0,0)


def test_lipid_bad_init(universe_model_flat):
    atoms = universe_model_flat.select_atoms("resid 1")

    with pytest.raises(TypeError) as excinfo:  # Should trigger a mda.groups.
        Lipid(None, None)
    assert str(excinfo.value).startswith("atoms must be a MDAnalysis.AtomGroup.")

    with pytest.raises(TypeError) as excinfo:
        Lipid(atoms, None)
    assert str(excinfo.value).startswith("hg_atoms must be a MDAnalysis.AtomGroup")

    with pytest.raises(ValueError) as excinfo:
        Lipid(universe_model_flat.select_atoms("resid 1 2"), atoms.select_atoms("name PO4"))
    assert str(excinfo.value).startswith("Only lipids belonging to one single residue are supported")

    with pytest.raises(ValueError) as excinfo:
        Lipid(atoms, atoms.select_atoms("name BAD"))
    assert str(excinfo.value).startswith("'hg_atoms' group is empty")

    with pytest.raises(ValueError) as excinfo:
        Lipid(atoms, universe_model_flat.select_atoms("resid 2 and name PO4"))
    assert str(excinfo.value).startswith("'hg_atoms' group is not consistent with 'atoms' group")


@pytest.fixture(scope="module")
def single_lipid(universe_model_flat):
    atoms = universe_model_flat.select_atoms("resid 7")
    hg_atoms = universe_model_flat.select_atoms("resid 7 and name PO4")

    return Lipid(atoms, hg_atoms)


def test_lipid_size(single_lipid):
    assert len(single_lipid) == 12


def test_lipid_position(single_lipid):
    with pytest.warns(UserWarning) as record:
        assert_almost_equal(single_lipid.position, np.array([0, 6*8., 92.5]), decimal=3)
    assert len(record) == 1
    assert str(record[0].message) == "Lipid does not belong to any registry. No fast calculation " \
                                     "nor PBC-awareness available"


def test_lipid_direction(single_lipid):
    with pytest.warns(UserWarning) as record:
        assert_almost_equal(single_lipid.direction, np.array([0.0, 0.0, 1.0]), decimal=1)

    assert len(record) == 1
    assert str(record[0].message) == "Lipid does not belong to any registry. No fast calculation " \
                                     "nor PBC-awareness available"


def test_lipid_neighbours(single_lipid):
    with pytest.raises(ValueError) as excinfo:
        single_lipid.neighbours_ids
    assert str(excinfo.value) == "Neighbours are not available if lipid does not belong to LipidSystem"


def test_lipid_normal(single_lipid):
    with pytest.raises(ValueError) as excinfo:
        single_lipid.normal
    assert str(excinfo.value) == "Normal is not available if lipid does not belong to LipidSystem"


def test_no_universe():
    with pytest.raises(TypeError) as excinfo:
        LipidSystem(None, None)
    assert "universe must be an instance of MDAnalysis.Universe" in str(excinfo.value)


def test_no_headgroup(universe_model_flat):
    with pytest.raises(TypeError) as excinfo:
        LipidSystem(universe_model_flat, None)
    assert "headgroup_atoms argument must be a string" in str(excinfo.value)


def test_headgroup_selection_bad(universe_model_flat):
    with pytest.raises(ValueError) as excinfo:
        LipidSystem(universe_model_flat, "azerty")
    assert "Bad headgroup selection string: 'azerty'" in str(excinfo.value)


def test_headgroup_selection_empty(universe_model_flat):
    with pytest.raises(ValueError) as excinfo:
        with pytest.warns(UserWarning) as record:
            LipidSystem(universe_model_flat, "")
        assert len(record) == 1
        assert str(record[0].message) == "Empty string to select atoms, empty group returned."
    assert "headgroup selection" in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        LipidSystem(universe_model_flat, "name P")
    assert "Empty headgroup selection" in str(excinfo.value)


def test_headgroup_selection_all(universe_model_flat):
    with pytest.raises(ValueError) as excinfo:
        LipidSystem(universe_model_flat, "all")
    assert "Headgroup selection is whole universe" in str(excinfo.value)


@pytest.mark.filterwarnings("ignore: Failed to guess the mass for the following atom")
def test_headgroup_selection_whole_residue(universe_model_flat):
    with pytest.raises(ValueError) as excinfo:
            LipidSystem(universe_model_flat, "resname DPPC and resid 1")
    assert "Headgroup selection corresponds to whole residue" in str(excinfo.value)


def test_lipid_system_size(system_model_flat):
    assert len(system_model_flat) == 128, "Bad number of lipids"


def test_lipid_system_positions_bbox(system_model_flat):
    import MDAnalysis as mda
    bbox_coords = mda.Universe(MODEL_FLAT_BBOX_GRO).atoms.positions
    assert_almost_equal(system_model_flat.positions_bbox, bbox_coords, decimal=3)


@pytest.mark.filterwarnings("ignore: Lipid does not belong")
def test_lipid_system_positions(system_model_flat, single_lipid):
    assert_almost_equal(system_model_flat[6].position, single_lipid.position, decimal=3)

    base_x = 0
    nlipids = 8
    unit_distance = 8
    base_y = [0, 0]
    base_z = [92.5, 57.5]
    lipid_per_leaflet = 64

    for i in range(2 * lipid_per_leaflet):
        offset_z = int(i//lipid_per_leaflet)
        offset_x = int((i - lipid_per_leaflet * offset_z) // nlipids)
        offset_y = int(i - lipid_per_leaflet * offset_z - nlipids * offset_x)
        if offset_x % 2 == 1:
            offset_y += 0.5

        expected_position = np.array([
            base_x + offset_x * unit_distance,
            base_y[offset_z] + offset_y * unit_distance,
            base_z[offset_z]
        ])

        assert_almost_equal(system_model_flat[i].position, expected_position, decimal=3,
                            err_msg="Bad position for lipid #{} (offsets: {}, {}, {})".format(
                                i,
                                offset_x, offset_y, offset_z
                            ))


def test_lipid_system_centers_naive(system_model_flat):
    # Centroids naively calculated from BBox positions should match MDAnalysis centroid computed using PBC
    u = system_model_flat.universe

    for i, residue in enumerate(u.residues):
        mda_centroid_pbc = residue.atoms.center_of_geometry(pbc=True)

        assert_almost_equal(mda_centroid_pbc,
                            system_model_flat.positions_bbox[system_model_flat[i]._ix].mean(axis=0),
                            decimal=3,
                            err_msg="Difference between MDAnalysis and FATSLiM centroids for lipid #{}".format(i)
                            )


def test_lipid_system_centers_from_headgroups(system_model_flat):
    # Internally lipid centroids are calculated using PBC distance from the headgroup position.
    # This can differ from naive calculation where the headgroup position is close to the box limits.
    # Yet the centroid is ensured to be inside the BBox

    u = system_model_flat.universe

    for i, residue in enumerate(u.residues):
        mda_centroid_nopbc = residue.atoms.center_of_geometry(pbc=False)
        # Make sure the MDA-calculated centroid is inside the bbox
        for dim in range(2):
            if mda_centroid_nopbc[dim] > 64:
                mda_centroid_nopbc[dim] -= 64
            if mda_centroid_nopbc[dim] < 0:
                mda_centroid_nopbc[dim] += 64

        assert_almost_equal(system_model_flat.lipid_centroids[i], mda_centroid_nopbc, decimal=3,
                            err_msg="Bad center for lipid #{}".format(i))


@pytest.mark.filterwarnings("ignore: Lipid does not belong")
def test_lipid_system_directions(system_model_flat, single_lipid):
    assert_almost_equal(system_model_flat[6].direction, single_lipid.direction, decimal=3)

    for i in range(128):
        direction = np.array([0.0, 0.0, 1.0])
        if i > 63:
            direction *= -1

        assert_almost_equal(system_model_flat[i].direction, direction, decimal=1,
                            err_msg="Bad direction for lipid #{}".format(i))


def test_lipid_system_neighbours(system_model_flat):
    from MDAnalysis.lib.distances import self_capped_distance
    cutoff = 20
    nlipids = 128

    pairs, distances = self_capped_distance(system_model_flat.lipid_positions, cutoff,
                                            box=system_model_flat.universe.dimensions,
                                            method="bruteforce")

    mda_neighbours = []
    for i in range(nlipids):
        mda_neighbours.append([])

    for k, [i, j] in enumerate(pairs):
        distance = distances[k]

        mda_neighbours[i].append((j, distance))
        mda_neighbours[j].append((i, distance))

    for beadid in range(nlipids):
        ref_result = sorted(mda_neighbours[beadid])
        print(system_model_flat[beadid].neighbours_distances)
        print(ref_result)
        assert_almost_equal(system_model_flat[beadid].neighbours_distances, ref_result, decimal=4,
                            err_msg="Bad neigbours for lipid #{}".format(beadid))


def test_lipid_system_normals(system_model_flat):

    for i in range(128):
        normal = np.array([0.0, 0.0, 1.0])
        if i > 63:
            normal *= -1

        assert_almost_equal(system_model_flat[i].normal, normal, decimal=4,
                            err_msg="Bad normal for lipid #{}".format(i))


@pytest.fixture(scope="module")
def system_model_flat_multiple_hg(universe_model_flat):
    return LipidSystem(universe_model_flat, "resname DPPC and name PO4 NC3")


def test_lipid_system_multiple_headgroup_atoms_size(system_model_flat_multiple_hg):
    assert len(system_model_flat_multiple_hg) == 128, "Bad number of lipids"


def test_lipid_system_multiple_headgroup_positions(system_model_flat_multiple_hg, system_model_flat):
    mda_hg_atoms = system_model_flat.universe.select_atoms("resname DPPC and name PO4 NC3")

    for resid, atoms in mda_hg_atoms.groupby("resids").items():

        mda_centroid_nopbc = atoms.center_of_geometry(pbc=False)
        # Make sure the MDA-calculated centroid is inside the bbox
        for dim in range(2):
            if mda_centroid_nopbc[dim] > 64:
                mda_centroid_nopbc[dim] -= 64
            if mda_centroid_nopbc[dim] < 0:
                mda_centroid_nopbc[dim] += 64

        assert_almost_equal(system_model_flat_multiple_hg.lipid_positions[resid - 1],
                            mda_centroid_nopbc,
                            decimal=4)


def test_lipid_system_multiple_headgroup_directions(system_model_flat_multiple_hg, system_model_flat):
    assert_almost_equal(system_model_flat_multiple_hg.lipid_directions,
                        system_model_flat.lipid_directions,
                        decimal=1)


def test_lipid_system_multiple_headgroup_normals(system_model_flat_multiple_hg, system_model_flat):
    assert_almost_equal(system_model_flat_multiple_hg.lipid_normals,
                        system_model_flat.lipid_normals,
                        decimal=4)


def test_vesicle_normals(system_model_vesicle):
    system = system_model_vesicle

    center = system.lipid_positions.mean(axis=0)

    threshold_radius = 60 - 35/2

    for i, lipid in enumerate(system):
        v = lipid.position - center
        normv = np.linalg.norm(v)

        normal = v/normv

        if normv < threshold_radius:
            normal *= -1

        assert_almost_equal(lipid.normal, normal, decimal=1,
                            err_msg="Bad normal for lipid #{}".format(i))

