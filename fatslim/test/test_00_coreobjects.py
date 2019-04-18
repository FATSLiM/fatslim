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
from . import universe_model_bilayer, system_model_bilayer, universe_simple_bilayer


def test_lipid_bad_init(universe_model_bilayer):
    atoms = universe_model_bilayer.select_atoms("resid 1")

    with pytest.raises(TypeError) as excinfo:  # Should trigger a mda.groups.
        Lipid(None, None)
    assert str(excinfo.value).startswith("atoms must be a MDAnalysis.AtomGroup.")

    with pytest.raises(TypeError) as excinfo:
        Lipid(atoms, None)
    assert str(excinfo.value).startswith("hg_atoms must be a MDAnalysis.AtomGroup")


@pytest.fixture(scope="module")
def single_lipid(universe_model_bilayer):
    atoms = universe_model_bilayer.select_atoms("resid 7")
    hg_atoms = universe_model_bilayer.select_atoms("resid 7 and name P")

    return Lipid(atoms, hg_atoms)


def test_lipid_size(single_lipid):
    assert len(single_lipid) == 50


def test_lipid_position(single_lipid):
    with pytest.warns(UserWarning) as record:
        assert_almost_equal(single_lipid.position, np.array([11.85, 4.46, 76.25]), decimal=3)
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


def test_no_headgroup(universe_model_bilayer):
    with pytest.raises(TypeError) as excinfo:
        LipidSystem(universe_model_bilayer, None)
    assert "headgroup_atoms argument must be a string" in str(excinfo.value)


def test_headgroup_selection_bad(universe_model_bilayer):
    with pytest.raises(ValueError) as excinfo:
        LipidSystem(universe_model_bilayer, "azerty")
    assert "Bad headgroup selection string: 'azerty'" in str(excinfo.value)


def test_headgroup_selection_empty(universe_model_bilayer):
    with pytest.raises(ValueError) as excinfo:
        with pytest.warns(UserWarning) as record:
            LipidSystem(universe_model_bilayer, "")
        assert len(record) == 1
        assert str(record[0].message) == "Empty string to select atoms, empty group returned."
    assert "headgroup selection" in str(excinfo.value)


def test_headgroup_selection_all(universe_model_bilayer):
    with pytest.raises(ValueError) as excinfo:
        LipidSystem(universe_model_bilayer, "all")
    assert "Headgroup selection is whole universe" in str(excinfo.value)


@pytest.mark.filterwarnings("ignore: Failed to guess the mass for the following atom")
def test_headgroup_selection_whole_residue(universe_simple_bilayer):
    with pytest.raises(ValueError) as excinfo:
            LipidSystem(universe_simple_bilayer, "resname DPC")
    assert "Headgroup selection corresponds to whole residue" in str(excinfo.value)


@pytest.fixture(scope="module")
def lipid_system(universe_model_bilayer):
    return LipidSystem(universe_model_bilayer, "resname DPPC and name P")


def test_lipid_system_size(lipid_system):
    assert len(lipid_system) == 72, "Bad number of lipids"


def test_lipid_system_positions_bbox(lipid_system):
    coords = lipid_system.universe.atoms.positions
    assert_almost_equal(lipid_system.positions_bbox, coords, decimal=3)


@pytest.mark.filterwarnings("ignore: Lipid does not belong")
def test_lipid_system_positions(lipid_system, single_lipid):
    assert_almost_equal(lipid_system[6].position, single_lipid.position, decimal=3)

    base_x = 3.85
    nlipids_x = 6
    unit_distance = 8
    base_y = [4.46, 3.54]
    base_z = [76.25, 23.75]
    lipid_per_leaflet = int(72 / 2)

    for i in range(2 * lipid_per_leaflet):
        offset_z = int(i//lipid_per_leaflet)
        if offset_z == 0:
            offset_x = int((i - lipid_per_leaflet * offset_z) // nlipids_x)
            offset_y = int(i - lipid_per_leaflet * offset_z - nlipids_x * offset_x)
        else:
            offset_x = int((i - lipid_per_leaflet * offset_z) // nlipids_x)
            offset_y = nlipids_x - 1 - int(i - lipid_per_leaflet * offset_z - nlipids_x * offset_x)

        expected_position = np.array([
            base_x + offset_x * unit_distance,
            base_y[offset_z] + offset_y * unit_distance,
            base_z[offset_z]
        ])

        assert_almost_equal(lipid_system[i].position, expected_position, decimal=3,
                            err_msg="Bad position for lipid #{} (offsets: {}, {}, {})".format(
                                i,
                                offset_x, offset_y, offset_z
                            ))


def test_lipid_system_centers(lipid_system):
    u = lipid_system.universe

    #print(u.residues[3].atoms._ix)
    #print(lipid_system[3]._ix)
    assert_almost_equal(u.residues[3].atoms.center_of_geometry(),
                        lipid_system.positions_bbox[lipid_system[3]._ix].mean(axis=0),
                        decimal=3
                        )

    for i, residue in enumerate(u.residues):
        center = residue.atoms.center_of_geometry()

        assert_almost_equal(lipid_system.lipid_centers[i], center, decimal=3,
                            err_msg="Bad center for lipid #{}".format(i))


@pytest.mark.filterwarnings("ignore: Lipid does not belong")
def test_lipid_system_directions(lipid_system, single_lipid):
    assert_almost_equal(lipid_system[6].direction, single_lipid.direction, decimal=3)

    for i in range(72):
        direction = np.array([0.0, 0.0, 1.0])
        if i > 35:
            direction *= -1

        assert_almost_equal(lipid_system[i].direction, direction, decimal=1,
                            err_msg="Bad direction for lipid #{}".format(i))


def test_lipid_system_neighbours(lipid_system):
    nlipids_x = 6
    nlipids_y = 6
    unit_distance = 8.0
    cutoff = 20.0

    def get_results_from_bid(bid):
        if bid >= nlipids_x * nlipids_y:
            zoffset = 1
            bid -= nlipids_x * nlipids_y
        else:
            zoffset = 0

        neighbours = []

        bead_y = bid // nlipids_y
        bead_x = bid - bead_y * nlipids_y

        extent_x = int(cutoff // unit_distance)
        extent_y = int(cutoff // unit_distance)

        for offset_x in range(-extent_x, extent_x + 1):
            x = bead_x + offset_x
            if x < 0:
                x += nlipids_x
            if x >= nlipids_x:
                x -= nlipids_x

            for offset_y in range(-extent_y, extent_y + 1):
                y = bead_y + offset_y
                if y < 0:
                    y += nlipids_y
                if y >= nlipids_y:
                    y -= nlipids_y

                nid = y * nlipids_y + x + zoffset * nlipids_y * nlipids_x
                d = np.sqrt((offset_x * unit_distance)**2 + (offset_y * unit_distance)**2)

                if 0.5 < d < cutoff:
                    neighbours.append((nid, d))

        return sorted(neighbours)

    for bid in range(72):
        assert_almost_equal(lipid_system[bid].neighbours_distances, get_results_from_bid(bid), decimal=4,
                            err_msg="Bad neigbours for lipid #{}".format(bid))


def test_lipid_system_normals(lipid_system):

    for i in range(72):
        normal = np.array([0.0, 0.0, 1.0])
        if i > 35:
            normal *= -1

        assert_almost_equal(lipid_system[i].normal, normal, decimal=4,
                            err_msg="Bad normal for lipid #{}".format(i))


@pytest.fixture(scope="module")
def lipid_system_multiple_hg(universe_model_bilayer):
    return LipidSystem(universe_model_bilayer, "resname DPPC and name P O3*")


def test_lipid_system_multiple_headgroup_atoms_size(lipid_system_multiple_hg):
    assert len(lipid_system_multiple_hg) == 72, "Bad number of lipids"

