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
from numpy.testing import assert_almost_equal
import pytest

# Local imports
from ..coreobjects import LipidSystem, Lipid
from . import universe_model_bilayer, system_model_bilayer, universe_simple_bilayer


def test_lipid_bad_init(universe_model_bilayer):
    atoms = universe_model_bilayer.select_atoms("resid 1")

    with pytest.raises(TypeError):
        Lipid(None, None)
        Lipid(atoms, None)

@pytest.mark.xfail(reason="Directions not implemented")
def test_lipid_init_ok(universe_model_bilayer):
    atoms = universe_model_bilayer.select_atoms("resid 7")
    hg_atoms = universe_model_bilayer.select_atoms("resid 7 and name P")

    lipid = Lipid(atoms, hg_atoms)

    assert len(lipid) == 50

    assert_almost_equal(lipid.position, np.array([11.85, 4.46, 76.25]), decimal=3)

    assert_almost_equal(lipid.direction, np.array([0.0, 0.0, 1.0]))


def test_no_universe():
    with pytest.raises(TypeError) as excinfo:
        LipidSystem(None, None)
    assert "First argument must be an instance of MDAnalysis.Universe" in str(excinfo.value)


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
        LipidSystem(universe_model_bilayer, "")
    assert "Empty headgroup selection" in str(excinfo.value)


def test_headgroup_selection_all(universe_model_bilayer):
    with pytest.raises(ValueError) as excinfo:
        LipidSystem(universe_model_bilayer, "all")
    assert "Headgroup selection is whole universe" in str(excinfo.value)


def test_headgroup_selection_whole_residue(universe_simple_bilayer):
    with pytest.raises(ValueError) as excinfo:
        LipidSystem(universe_simple_bilayer, "resname DPC")
    assert "Headgroup selection corresponds to whole residue" in str(excinfo.value)


def test_headgroup_selection_string_ok(universe_model_bilayer):
    system = LipidSystem(universe_model_bilayer, "resname DPPC and name P")
    assert len(system) == 72, "Bad number of lipids"


def test_headgroup_selection_string_multiple_atoms(universe_model_bilayer):
    system = LipidSystem(universe_model_bilayer, "resname DPPC and name P O3*")
    assert len(system) == 72, "Bad number of lipids"

