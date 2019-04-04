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
import pytest
import MDAnalysis as mda

# Local imports
from .data import BILAYER_GRO, MODEL_BILAYER_GRO, MODEL_BILAYER_NDX, BILAYER_GRO_ALLATOM
from ..coreobjects import LipidSystem


@pytest.fixture(scope="session")
def universe_model_bilayer() -> mda.Universe:
    u = mda.Universe(MODEL_BILAYER_GRO)
    return u


@pytest.fixture(scope="session")
def system_model_bilayer() -> LipidSystem:
    u = mda.Universe(MODEL_BILAYER_GRO)
    system = LipidSystem(u, "name P")
    return system


@pytest.fixture(scope="session")
def universe_simple_bilayer() -> mda.Universe:
    u = mda.Universe(BILAYER_GRO)
    return u
