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
import pytest

# Local imports
from .data import VESICLE_GRO, VESICLE_NDX, VESICLE_XTC, MODEL_BILAYER_NDX, MODEL_BILAYER_GRO, \
    BILAYER_CHOL_GRO, BILAYER_CHOL_NDX, MODEL_VESICLE_GRO, MODEL_VESICLE_NDX, BIG_PROT_GRO, \
    BIG_PROT_NDX, MODEL_MULTIBILAYER_GRO, MODEL_MULTIBILAYER_NDX, MODEL_BILAYER_VESICLE_GRO, \
    MODEL_BILAYER_VESICLE_NDX, BILAYER_PROT_GRO, BILAYER_PROT_NDX, MODEL_BILAYER_PROT_GRO, \
    MODEL_BILAYER_PROT_NDX, BILAYER_GRO, BILAYER_NDX, BILAYER_PEPTIDE_GRO, BILAYER_PEPTIDE_NDX, \
    BILAYER_ALLATOM_GRO, BILAYER_ALLATOM_NDX
from .data import MODEL_BILAYER_THICKNESS_CSV, BILAYER_CHOL_APL_GROUPED, \
    BILAYER_CHOL_APL_NOGROUP, MODEL_BILAYER_APL_CSV, VESICLE_APL_XVG, \
    MODEL_BILAYER_THICKNESS_XVG, VESICLE_THICKNESS_XVG, VESICLE_AREA_XVG, MULTIBILAYER_GRO, \
    MULTIBILAYER_NDX
from ..datareading import load_trajectory


@pytest.fixture(scope="session")
def frame_model_bilayer_vesicle():
    traj = load_trajectory(MODEL_BILAYER_VESICLE_GRO, MODEL_BILAYER_VESICLE_NDX)
    traj.initialize()
    frame = traj[0]
    return frame


@pytest.fixture(scope="session")
def frame_model_multibilayer():
    traj = load_trajectory(MODEL_MULTIBILAYER_GRO, MODEL_MULTIBILAYER_NDX)
    traj.initialize()
    frame = traj[0]
    return frame

@pytest.fixture(scope="session")
def frame_multibilayer():
    traj = load_trajectory(MULTIBILAYER_GRO, MULTIBILAYER_NDX)
    traj.initialize()
    frame = traj[0]
    return frame

@pytest.fixture(scope="session")
def frame_big_prot():
    traj = load_trajectory(BIG_PROT_GRO, BIG_PROT_NDX)
    traj.initialize()
    frame = traj[0]
    return frame


@pytest.fixture(scope="session")
def traj_vesicle():
    traj = load_trajectory(VESICLE_GRO, VESICLE_NDX, VESICLE_XTC)
    traj.initialize()
    return traj


@pytest.fixture(scope="session")
def frame_vesicle():
    traj = load_trajectory(VESICLE_GRO, VESICLE_NDX, VESICLE_XTC)
    traj.initialize()
    frame = traj[0]
    return frame


@pytest.fixture(scope="session")
def frame_model_vesicle():
    traj = load_trajectory(MODEL_VESICLE_GRO, MODEL_VESICLE_NDX)
    traj.initialize()
    frame = traj[0]
    return frame


@pytest.fixture(scope="session")
def frame_model_bilayer():
    traj = load_trajectory(MODEL_BILAYER_GRO, MODEL_BILAYER_NDX)
    traj.initialize()
    frame = traj[0]
    return frame


@pytest.fixture(scope="session")
def frame_bilayer():
    traj = load_trajectory(BILAYER_GRO, BILAYER_NDX)
    traj.initialize()
    frame = traj[0]
    return frame


@pytest.fixture(scope="session")
def frame_bilayer_allatom():
    traj = load_trajectory(BILAYER_ALLATOM_GRO, BILAYER_ALLATOM_NDX)
    traj.initialize()
    frame = traj[0]
    return frame


@pytest.fixture(scope="session")
def frame_model_bilayer_prot():
    traj = load_trajectory(MODEL_BILAYER_PROT_GRO, MODEL_BILAYER_PROT_NDX)
    traj.initialize()
    frame = traj[0]
    return frame


@pytest.fixture(scope="session")
def frame_bilayer_prot():
    traj = load_trajectory(BILAYER_PROT_GRO, BILAYER_PROT_NDX)
    traj.initialize()
    frame = traj[0]
    return frame


@pytest.fixture(scope="session")
def frame_bilayer_peptide():
    traj = load_trajectory(BILAYER_PEPTIDE_GRO, BILAYER_PEPTIDE_NDX)
    traj.initialize()
    frame = traj[0]
    return frame


@pytest.fixture(scope="session")
def frame_bilayer_chol():
    traj = load_trajectory(BILAYER_CHOL_GRO, BILAYER_CHOL_NDX)
    traj.initialize()
    frame = traj[0]
    return frame
