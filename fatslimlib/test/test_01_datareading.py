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
from unittest import TestCase

import pytest
import numpy
from numpy.testing import assert_almost_equal


# Local imports
from .data import MODEL_BILAYER_GRO, MODEL_BILAYER_NDX, VESICLE_GRO, VESICLE_NDX, VESICLE_XTC, \
    VESICLE_TRR, VESICLE_HG_XTC, \
    MODEL_BIG_GRO, MODEL_BIG_NDX, MODEL_BILAYER_PROT_GRO, MODEL_BILAYER_PROT_NDX, \
    BILAYER_ALLATOM_GRO, BILAYER_ALLATOM_NDX, BILAYER_PEPTIDE_GRO, BILAYER_PEPTIDE_H_NDX, \
    BILAYER_CHOL_GRO, BILAYER_CHOL_NDX
from ..datareading import UnknownFileType, get_readers, load_trajectory
from ..core_datareading import NdxReader


def test_index_filetype():
    with pytest.raises(ValueError):
        NdxReader("fake.index")


class TestIndexReader(TestCase):
    def setUp(self):
        self.index = NdxReader(MODEL_BILAYER_NDX)
        self.ref_size = 4
        self.ref_group = "headgroups"
        self.ref_group_size = 144
        self.ref_group_firsts = numpy.array([4, 8, 54, 58, 104], dtype=int)
        self.ref_group_lasts = numpy.array([3458, 3504, 3508, 3554, 3558], dtype=int)

    def test_size(self):
        assert len(self.index) == self.ref_size

    def test_group_size(self):
        assert len(self.index[self.ref_group]) == self.ref_group_size

    def test_group_firsts(self):
        assert_almost_equal(self.index[self.ref_group][:5], self.ref_group_firsts)

    def test_group_lasts(self):
        assert_almost_equal(self.index[self.ref_group][-5:], self.ref_group_lasts)

    def test_fake_group(self):
        with pytest.raises(KeyError):
            assert self.index["fake_group"]


def test_unk_reader():
    with pytest.raises(UnknownFileType):
        get_readers("test.unk", MODEL_BILAYER_NDX)

    with pytest.raises(UnknownFileType):
        get_readers(MODEL_BILAYER_GRO, "test.unk")

    with pytest.raises(UnknownFileType):
        get_readers(MODEL_BILAYER_GRO, MODEL_BILAYER_NDX, MODEL_BILAYER_NDX)

    # This should not raise anything
    get_readers(MODEL_BILAYER_GRO, MODEL_BILAYER_NDX, "test.unk")


def test_gro_reader():
    assert get_readers(MODEL_BILAYER_GRO, MODEL_BILAYER_NDX)


def test_trajectory_initialization():
    traj = load_trajectory(MODEL_BILAYER_GRO, MODEL_BILAYER_NDX)

    with pytest.raises(ValueError):
        assert traj[0]


def test_trajectory_initialization_hydrogen():
    traj = load_trajectory(BILAYER_PEPTIDE_GRO, BILAYER_PEPTIDE_H_NDX)
    with pytest.raises(KeyError):
        traj.initialize()


def test_trajectory_initialization_incoherent_traj():
    with pytest.raises(IndexError):
        traj = load_trajectory(VESICLE_GRO, VESICLE_NDX, VESICLE_HG_XTC)


class TestTopolReading(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.traj = load_trajectory(BILAYER_ALLATOM_GRO, BILAYER_ALLATOM_NDX)
        cls.traj.initialize()
        cls.topology = cls.traj[0].topology

    def test_hydrogens_skipped(self):
        assert self.topology.get_residue(1)["size"] == 46

    def test_solvent_skipped(self):
        assert self.topology.natoms == 128 * 46
        assert self.topology.nresidues == 128


class TestGroReading(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.traj = load_trajectory(MODEL_BILAYER_PROT_GRO, MODEL_BILAYER_PROT_NDX)
        cls.traj.initialize()

    def test_hg_atomids(self):
        frame = self.traj[0]
        atomids = frame.hg_atomids

        ref_group_firsts = numpy.array([2227, 2231, 2277, 2281, 2327], dtype=int)
        ref_group_lasts = numpy.array([26081, 26127, 26131, 26177, 26181], dtype=int)

        assert_almost_equal(atomids[:5], ref_group_firsts)
        assert_almost_equal(atomids[-5:], ref_group_lasts)

    def test_lipid_atomids(self):
        atomids = self.traj[0].lipid_atomids
        ref_group_firsts = numpy.array([2224, 2225, 2226, 2227, 2228], dtype=int)
        ref_group_lasts = numpy.array([26219, 26220, 26221, 26222, 26223], dtype=int)

        assert_almost_equal(atomids[0][:5], ref_group_firsts)
        assert_almost_equal(atomids[-1][-5:], ref_group_lasts)

    def test_timesteps(self):
        assert self.traj[0].timestep == 0.0

    def test_traj_size(self):
        assert len(self.traj) == 1

    def test_pbcbox(self):
        assert_almost_equal(self.traj[0].box.asarray(), numpy.array([[12.8, 0, 0], [0, 12.8, 0], [0, 0, 10.0]]))

    def test_hg_coords(self):
        coords = self.traj[0].hg_group_coords

        ref_coords_first = numpy.array((0.362, 0.356, 8.039))
        ref_coords_last = numpy.array((12.385, 11.554, 2.375))

        assert_almost_equal(coords[0], ref_coords_first)
        assert_almost_equal(coords[-1], ref_coords_last)

    def test_lipid_coords(self):
        coords = self.traj[0].lipid_coords

        ref_coords_first = numpy.array((0.459, 0.281, 8.119))
        ref_coords_last = numpy.array((12.372, 11.480, 4.924))

        assert_almost_equal(coords[0][0], ref_coords_first)
        assert_almost_equal(coords[-1][-1], ref_coords_last)

    def test_bead_coords(self):
        coords = self.traj[0].bead_coords

        ref_coords_first = 0.5 * (numpy.array((0.362, 0.356, 8.039)) + numpy.array((0.385, 0.446, 7.625)))
        ref_coords_last = 0.5 * (numpy.array((12.362, 11.644, 1.961)) + numpy.array((12.385, 11.554, 2.375)))

        assert_almost_equal(coords[0], ref_coords_first)
        assert_almost_equal(coords[-1], ref_coords_last)

    def test_interacting_coords(self):
        coords = self.traj[0].interacting_group_coords

        assert len(coords) > 0

        ref_coords_first = numpy.array((7.258, 6.204, 7.468))
        ref_coords_last = numpy.array((7.312, 5.275, 6.634))

        assert_almost_equal(coords[0], ref_coords_first)
        assert_almost_equal(coords[-1], ref_coords_last)

    @classmethod
    def tearDownClass(cls):
        del cls.traj


class TestBigGroReading(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.traj = load_trajectory(MODEL_BIG_GRO, MODEL_BIG_NDX)
        cls.traj.initialize()

    def test_hg_atomids(self):
        atomids = self.traj[0].hg_atomids

        ref_group_firsts = numpy.array([8, 58, 108, 158, 208], dtype=int)
        ref_group_lasts = numpy.array([129358, 129408, 129458, 129508, 129558], dtype=int)

        assert_almost_equal(atomids[:5], ref_group_firsts)
        assert_almost_equal(atomids[-5:], ref_group_lasts)

    def test_lipid_atomids(self):
        atomids = self.traj[0].lipid_atomids
        ref_group_firsts = numpy.array([1, 2, 3, 4, 5], dtype=int)
        ref_group_lasts = numpy.array([129596, 129597, 129598, 129599, 129600], dtype=int)

        assert_almost_equal(atomids[0][:5], ref_group_firsts)
        assert_almost_equal(atomids[-1][-5:], ref_group_lasts)

    def test_timesteps(self):
        assert self.traj[0].timestep == 0.0

    def test_traj_size(self):
        assert len(self.traj) == 1

    def test_pbcbox(self):
        assert_almost_equal(self.traj[0].box.asarray(), numpy.array([[28.8, 0, 0], [0, 28.8, 0], [0, 0, 10.0]]))

    def test_hg_coords(self):
        coords = self.traj[0].hg_group_coords

        ref_coords_first = numpy.array((0.385, 0.446, 7.625))
        ref_coords_last = numpy.array((28.385, 24.354, 2.375))

        assert_almost_equal(coords[0], ref_coords_first)
        assert_almost_equal(coords[-1], ref_coords_last)

    def test_lipid_coords(self):
        coords = self.traj[0].lipid_coords

        ref_coords_first = numpy.array((0.459, 0.281, 8.119))
        ref_coords_last = numpy.array((28.372, 24.280, 4.924))

        assert_almost_equal(coords[0][0], ref_coords_first)
        assert_almost_equal(coords[-1][-1], ref_coords_last)

    def test_bead_coords(self):
        coords = self.traj[0].bead_coords

        ref_coords_first = numpy.array((0.385, 0.446, 7.625))
        ref_coords_last = numpy.array((28.385, 24.354, 2.375))

        assert_almost_equal(coords[0], ref_coords_first)
        assert_almost_equal(coords[-1], ref_coords_last)

    @classmethod
    def tearDownClass(cls):
        del cls.traj


class TestCholReading(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.traj = load_trajectory(BILAYER_CHOL_GRO, BILAYER_CHOL_NDX)
        cls.traj.initialize()

    def test_hg_atomids(self):
        atomids = self.traj[0].hg_atomids

        ref_group_firsts = numpy.array([2, 14, 26, 38, 50], dtype=int)
        ref_group_lasts = numpy.array([33585, 33593, 33601, 33609, 33617], dtype=int)

        assert_almost_equal(atomids[:5], ref_group_firsts)
        assert_almost_equal(atomids[-5:], ref_group_lasts)

    def test_lipid_atomids(self):
        atomids = self.traj[0].lipid_atomids
        ref_group_firsts = numpy.array([1, 2, 3, 4, 5], dtype=int)
        ref_group_lasts = numpy.array([33620, 33621, 33622, 33623, 33624], dtype=int)

        assert_almost_equal(atomids[0][:5], ref_group_firsts)
        assert_almost_equal(atomids[-1][-5:], ref_group_lasts)

    def test_timesteps(self):
        assert self.traj[0].timestep == 0.0

    def test_traj_size(self):
        assert len(self.traj) == 1

    def test_pbcbox(self):
        assert_almost_equal(self.traj[0].box.asarray(), numpy.array([[21.6588, 0, 0],
                                                                     [0, 21.82714, 0],
                                                                     [0, 0, 7.48249]]))

    def test_hg_coords(self):
        coords = self.traj[0].hg_group_coords

        ref_coords_first = numpy.array((11.950, 9.160, 1.680))
        ref_coords_last = numpy.array((16.720, 0.350, 4.450))

        assert_almost_equal(coords[0], ref_coords_first)
        assert_almost_equal(coords[-1], ref_coords_last)

    def test_lipid_coords(self):
        coords = self.traj[0].lipid_coords

        ref_coords_first = numpy.array((12.300, 9.160, 1.530))
        ref_coords_last = numpy.array((16.590, 0.120, 3.030))

        assert_almost_equal(coords[0][0], ref_coords_first)
        assert_almost_equal(coords[-1][-1], ref_coords_last)

    def test_bead_coords(self):
        coords = self.traj[0].bead_coords

        ref_coords_first = numpy.array((11.950, 9.160, 1.680))
        ref_coords_last = numpy.array((16.720, 0.350, 4.450))

        assert_almost_equal(coords[0], ref_coords_first)
        assert_almost_equal(coords[-1], ref_coords_last)

    @classmethod
    def tearDownClass(cls):
        del cls.traj


class TestXtcReading(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.traj = load_trajectory(VESICLE_GRO, VESICLE_NDX, VESICLE_XTC)
        cls.traj.initialize()

    @classmethod
    def tearDownClass(cls):
        del cls.traj

    def test_hg_atomids(self):
        atomids = self.traj[0].hg_atomids

        ref_group_firsts = numpy.array([2, 14, 26, 38, 50], dtype=int)
        ref_group_lasts = numpy.array([36806, 36818, 36830, 36842, 36854], dtype=int)

        assert_almost_equal(atomids[:5], ref_group_firsts)
        assert_almost_equal(atomids[-5:], ref_group_lasts)

    def test_lipid_atomids(self):
        atomids = self.traj[0].lipid_atomids

        ref_group_firsts = numpy.array([1, 2, 3, 4, 5], dtype=int)
        ref_group_lasts = numpy.array([36860, 36861, 36862, 36863, 36864], dtype=int)

        assert_almost_equal(atomids[0][:5], ref_group_firsts)
        assert_almost_equal(atomids[-1][-5:], ref_group_lasts)

    def test_traj_size(self):
        assert len(self.traj) == 11

    def test_timesteps(self):
        ref_timesteps = [500.00 * i for i in range(len(self.traj))]

        for i, frame in enumerate(self.traj):
            assert frame.timestep == ref_timesteps[i]

    def test_pbcbox(self):
        ref_boxes = [numpy.array([[47.100067138671875, 0.0, 0.0], [0.0, 47.29520797729492, 0.0],
            [23.550033569335938, 23.64760398864746, 31.985841751098633]]),
            numpy.array([[47.08340835571289, 0.0, 0.0], [0.0, 47.30361557006836, 0.0],
                [23.541704177856445, 23.65180778503418, 31.991649627685547]]),
            numpy.array([[47.08655548095703, 0.0, 0.0], [0.0, 47.298133850097656, 0.0],
                [23.543277740478516, 23.649066925048828, 31.99555778503418]]),
            numpy.array([[47.07117462158203, 0.0, 0.0], [0.0, 47.30364227294922, 0.0],
                [23.535587310791016, 23.65182113647461, 31.99390983581543]]),
            numpy.array([[47.079769134521484, 0.0, 0.0], [0.0, 47.30522918701172, 0.0],
                [23.539884567260742, 23.65261459350586, 31.995153427124023]]),
            numpy.array([[47.07292175292969, 0.0, 0.0], [0.0, 47.30583953857422, 0.0],
                [23.536460876464844, 23.65291976928711, 32.00019073486328]]),
            numpy.array([[47.066261291503906, 0.0, 0.0], [0.0, 47.3122444152832, 0.0],
                [23.533130645751953, 23.6561222076416, 31.99447250366211]]),
            numpy.array([[47.06310272216797, 0.0, 0.0], [0.0, 47.31338119506836, 0.0],
                [23.531551361083984, 23.65669059753418, 31.999902725219727]]),
            numpy.array([[47.075565338134766, 0.0, 0.0], [0.0, 47.29836654663086, 0.0],
                [23.537782669067383, 23.64918327331543, 32.0003662109375]]),
            numpy.array([[47.0737419128418, 0.0, 0.0], [0.0, 47.295833587646484, 0.0],
                [23.5368709564209, 23.647916793823242, 32.00419235229492]]),
            numpy.array([[47.06859588623047, 0.0, 0.0], [0.0, 47.29644775390625, 0.0],
                [23.534297943115234, 23.648223876953125, 32.00706100463867]]),
        ]

        for i, frame in enumerate(self.traj):
            assert_almost_equal(frame.box.asarray(), ref_boxes[i], 5)

    def test_hg_coords(self):
        ref_coords = [(numpy.array([46.704, 26.634, 25.838]), numpy.array([0.828, 25.352, 23.336])),
            (numpy.array([0.245, 25.697, 24.964]), numpy.array([1.076, 25.148, 23.05])),
            (numpy.array([46.654, 25.595, 24.466]), numpy.array([0.925, 24.454, 21.854])),
            (numpy.array([0.452, 26.92, 24.795]), numpy.array([1.798, 25.228, 21.268])),
            (numpy.array([0.153, 25.762, 24.212]), numpy.array([1.903, 24.785, 21.654])),
            (numpy.array([0.705, 26.689, 24.573]), numpy.array([1.880, 24.884, 21.318])),
            (numpy.array([0.278, 27.591, 25.560]), numpy.array([1.44, 25.025, 20.538])),
            (numpy.array([0.244, 27.5, 24.755]), numpy.array([1.751, 25.319, 20.476])),
            (numpy.array([0.604, 27.412, 24.827]), numpy.array([1.290, 24.705, 19.278])),
            (numpy.array([46.885, 27.279, 24.879]), numpy.array([1.415, 24.76, 19.318])),
            (numpy.array([47.049, 27.8, 25.764]), numpy.array([1.561, 25.331, 19.281])), ]

        for i, frame in enumerate(self.traj):
            ref_coords_first, ref_coords_last = ref_coords[i]
            coords = frame.hg_group_coords

            assert_almost_equal(coords[0], ref_coords_first, 3)
            assert_almost_equal(coords[-1], ref_coords_last, 3)

    def test_lipid_coords(self):
        ref_coords = [(numpy.array([47.089, 26.709, 25.909]), numpy.array([46.848, 26.583, 21.697])),
            (numpy.array([46.983, 25.571, 24.778]), numpy.array([46.677, 26.159, 21.724])),
            (numpy.array([46.917, 25.729, 24.715]), numpy.array([46.616, 26.104, 21.236])),
            (numpy.array([0.911, 26.731, 24.829]), numpy.array([0.292, 26.295, 21.759])),
            (numpy.array([46.834, 25.704, 24.316]), numpy.array([1.705, 26.92, 20.982])),
            (numpy.array([0.976, 26.966, 24.42]), numpy.array([0.986, 26.455, 20.507])),
            (numpy.array([47.02, 27.213, 25.534]), numpy.array([0.909, 26.701, 19.69])),
            (numpy.array([0.622, 27.231, 24.827]), numpy.array([1.296, 27.084, 19.463])),
            (numpy.array([0.792, 27.017, 24.68]), numpy.array([0.342, 26.490, 19.179])),
            (numpy.array([47.022, 27.047, 25.295]), numpy.array([2.073, 26.168, 19.01])),
            (numpy.array([0.3, 27.619, 25.796]), numpy.array([1.17, 27.371, 19.19])), ]

        for i, frame in enumerate(self.traj):
            ref_coords_first, ref_coords_last = ref_coords[i]
            coords = frame.lipid_coords

            assert_almost_equal(coords[0][0], ref_coords_first, 3)
            assert_almost_equal(coords[-1][-1], ref_coords_last, 3)

    def test_bead_coords(self):
        ref_coords = [(numpy.array([46.704, 26.634, 25.838]), numpy.array([0.828, 25.352, 23.336])),
            (numpy.array([0.245, 25.697, 24.964]), numpy.array([1.076, 25.148, 23.05])),
            (numpy.array([46.654, 25.595, 24.466]), numpy.array([0.925, 24.454, 21.854])),
            (numpy.array([0.452, 26.92, 24.795]), numpy.array([1.798, 25.228, 21.268])),
            (numpy.array([0.153, 25.762, 24.212]), numpy.array([1.903, 24.785, 21.654])),
            (numpy.array([0.705, 26.689, 24.573]), numpy.array([1.880, 24.884, 21.318])),
            (numpy.array([0.278, 27.591, 25.560]), numpy.array([1.44, 25.025, 20.538])),
            (numpy.array([0.244, 27.5, 24.755]), numpy.array([1.751, 25.319, 20.476])),
            (numpy.array([0.604, 27.412, 24.827]), numpy.array([1.290, 24.705, 19.278])),
            (numpy.array([46.885, 27.279, 24.879]), numpy.array([1.415, 24.76, 19.318])),
            (numpy.array([47.049, 27.8, 25.764]), numpy.array([1.561, 25.331, 19.281])), ]

        for i, frame in enumerate(self.traj):
            ref_coords_first, ref_coords_last = ref_coords[i]
            coords = frame.bead_coords

            assert_almost_equal(coords[0], ref_coords_first, 3)
            assert_almost_equal(coords[-1], ref_coords_last, 3)


class TestTrrReading(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.traj = load_trajectory(VESICLE_GRO, VESICLE_NDX, VESICLE_TRR)
        cls.traj.initialize()

    @classmethod
    def tearDownClass(cls):
        del cls.traj

    def test_hg_atomids(self):
        atomids = self.traj[0].hg_atomids

        ref_group_firsts = numpy.array([2, 14, 26, 38, 50], dtype=int)
        ref_group_lasts = numpy.array([36806, 36818, 36830, 36842, 36854], dtype=int)

        assert_almost_equal(atomids[:5], ref_group_firsts)
        assert_almost_equal(atomids[-5:], ref_group_lasts)

    def test_lipid_atomids(self):
        atomids = self.traj[0].lipid_atomids

        ref_group_firsts = numpy.array([1, 2, 3, 4, 5], dtype=int)
        ref_group_lasts = numpy.array([36860, 36861, 36862, 36863, 36864], dtype=int)

        assert_almost_equal(atomids[0][:5], ref_group_firsts)
        assert_almost_equal(atomids[-1][-5:], ref_group_lasts)

    def test_traj_size(self):
        assert len(self.traj) == 11

    def test_timesteps(self):
        ref_timesteps = [500.00 * i for i in range(len(self.traj))]

        for i, frame in enumerate(self.traj):
            assert frame.timestep == ref_timesteps[i]

    def test_pbcbox(self):
        ref_boxes = [numpy.array([[47.100067138671875, 0.0, 0.0], [0.0, 47.29520797729492, 0.0],
            [23.550033569335938, 23.64760398864746, 31.985841751098633]]),
            numpy.array([[47.08340835571289, 0.0, 0.0], [0.0, 47.30361557006836, 0.0],
                [23.541704177856445, 23.65180778503418, 31.991649627685547]]),
            numpy.array([[47.08655548095703, 0.0, 0.0], [0.0, 47.298133850097656, 0.0],
                [23.543277740478516, 23.649066925048828, 31.99555778503418]]),
            numpy.array([[47.07117462158203, 0.0, 0.0], [0.0, 47.30364227294922, 0.0],
                [23.535587310791016, 23.65182113647461, 31.99390983581543]]),
            numpy.array([[47.079769134521484, 0.0, 0.0], [0.0, 47.30522918701172, 0.0],
                [23.539884567260742, 23.65261459350586, 31.995153427124023]]),
            numpy.array([[47.07292175292969, 0.0, 0.0], [0.0, 47.30583953857422, 0.0],
                [23.536460876464844, 23.65291976928711, 32.00019073486328]]),
            numpy.array([[47.066261291503906, 0.0, 0.0], [0.0, 47.3122444152832, 0.0],
                [23.533130645751953, 23.6561222076416, 31.99447250366211]]),
            numpy.array([[47.06310272216797, 0.0, 0.0], [0.0, 47.31338119506836, 0.0],
                [23.531551361083984, 23.65669059753418, 31.999902725219727]]),
            numpy.array([[47.075565338134766, 0.0, 0.0], [0.0, 47.29836654663086, 0.0],
                [23.537782669067383, 23.64918327331543, 32.0003662109375]]),
            numpy.array([[47.0737419128418, 0.0, 0.0], [0.0, 47.295833587646484, 0.0],
                [23.5368709564209, 23.647916793823242, 32.00419235229492]]),
            numpy.array([[47.06859588623047, 0.0, 0.0], [0.0, 47.29644775390625, 0.0],
                [23.534297943115234, 23.648223876953125, 32.00706100463867]]),
        ]

        for i, frame in enumerate(self.traj):
            assert_almost_equal(frame.box.asarray(), ref_boxes[i], 5)

    def test_hg_coords(self):
        ref_coords = [(numpy.array([36.704, 26.634, 25.838]), numpy.array([37.928, 25.352, 23.336])),
            (numpy.array([37.328, 25.697, 24.964]), numpy.array([38.159, 25.148, 23.05])),
            (numpy.array([36.654, 25.595, 24.466]), numpy.array([38.012, 24.454, 21.854])),
            (numpy.array([37.523, 26.92, 24.795]), numpy.array([38.869, 25.228, 21.268])),
            (numpy.array([37.233, 25.762, 24.212]), numpy.array([38.983, 24.785, 21.654])),
            (numpy.array([37.778, 26.689, 24.573]), numpy.array([38.953, 24.884, 21.318])),
            (numpy.array([37.344, 27.591, 25.56]), numpy.array([38.506, 25.025, 20.538])),
            (numpy.array([37.307, 27.5, 24.755]), numpy.array([38.814, 25.319, 20.476])),
            (numpy.array([37.68, 27.412, 24.827]), numpy.array([38.366, 24.705, 19.278])),
            (numpy.array([36.885, 27.279, 24.879]), numpy.array([38.489, 24.76, 19.318])),
            (numpy.array([37.049, 27.8, 25.764]), numpy.array([38.63, 25.331, 19.281])), ]

        for i, frame in enumerate(self.traj):
            ref_coords_first, ref_coords_last = ref_coords[i]
            coords = frame.hg_group_coords

            assert_almost_equal(coords[0], ref_coords_first, 3)
            assert_almost_equal(coords[-1], ref_coords_last, 3)

    def test_lipid_coords(self):
        ref_coords = [(numpy.array([37.089, 26.709, 25.909]), numpy.array([36.848, 26.583, 21.697])),
            (numpy.array([36.983, 25.571, 24.778]), numpy.array([36.677, 26.159, 21.724])),
            (numpy.array([36.917, 25.729, 24.715]), numpy.array([36.616, 26.104, 21.236])),
            (numpy.array([37.982, 26.731, 24.829]), numpy.array([37.363, 26.295, 21.759])),
            (numpy.array([36.834, 25.704, 24.316]), numpy.array([38.785, 26.92, 20.982])),
            (numpy.array([38.049, 26.966, 24.42]), numpy.array([38.059, 26.455, 20.507])),
            (numpy.array([37.02, 27.213, 25.534]), numpy.array([37.975, 26.701, 19.69])),
            (numpy.array([37.685, 27.231, 24.827]), numpy.array([38.359, 27.084, 19.463])),
            (numpy.array([37.868, 27.017, 24.68]), numpy.array([37.418, 26.49, 19.179])),
            (numpy.array([37.022, 27.047, 25.295]), numpy.array([39.147, 26.168, 19.01])),
            (numpy.array([37.369, 27.619, 25.796]), numpy.array([38.239, 27.371, 19.19])), ]

        for i, frame in enumerate(self.traj):
            ref_coords_first, ref_coords_last = ref_coords[i]
            coords = frame.lipid_coords

            assert_almost_equal(coords[0][0], ref_coords_first, 3)
            assert_almost_equal(coords[-1][-1], ref_coords_last, 3)

    def test_bead_coords(self):
        ref_coords = [(numpy.array([36.704, 26.634, 25.838]), numpy.array([37.928, 25.352, 23.336])),
            (numpy.array([37.328, 25.697, 24.964]), numpy.array([38.159, 25.148, 23.05])),
            (numpy.array([36.654, 25.595, 24.466]), numpy.array([38.012, 24.454, 21.854])),
            (numpy.array([37.523, 26.92, 24.795]), numpy.array([38.869, 25.228, 21.268])),
            (numpy.array([37.233, 25.762, 24.212]), numpy.array([38.983, 24.785, 21.654])),
            (numpy.array([37.778, 26.689, 24.573]), numpy.array([38.953, 24.884, 21.318])),
            (numpy.array([37.344, 27.591, 25.56]), numpy.array([38.506, 25.025, 20.538])),
            (numpy.array([37.307, 27.5, 24.755]), numpy.array([38.814, 25.319, 20.476])),
            (numpy.array([37.68, 27.412, 24.827]), numpy.array([38.366, 24.705, 19.278])),
            (numpy.array([36.885, 27.279, 24.879]), numpy.array([38.489, 24.76, 19.318])),
            (numpy.array([37.049, 27.8, 25.764]), numpy.array([38.63, 25.331, 19.281])), ]

        for i, frame in enumerate(self.traj):
            ref_coords_first, ref_coords_last = ref_coords[i]
            coords = frame.bead_coords

            assert_almost_equal(coords[0], ref_coords_first, 3)
            assert_almost_equal(coords[-1], ref_coords_last, 3)

