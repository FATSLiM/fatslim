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
from ..core_base import Atom, AtomGroup, PBCBox, Topology


class TestAtom(TestCase):
    def setUp(self):
        self.ref_x = 1.1
        self.ref_y = 2.2
        self.ref_z = 3.3
        self.ref_atomid = 1
        self.ref_name = "H"
        self.atom = Atom(atomid=self.ref_atomid,
                         name=self.ref_name,
                         x=self.ref_x,
                         y=self.ref_y,
                         z=self.ref_z)

    def test_read_atomid(self):
        assert self.atom.atomid == self.ref_atomid

    def test_ro_atomid(self):
        with pytest.raises(AttributeError):
            self.atom.atomid = 0

    def test_read_name(self):
        assert self.atom.name == self.ref_name

    def test_ro_name(self):
        with pytest.raises(AttributeError):
            self.atom.name = "C"

    def test_read_coords(self):
        assert self.atom.x == self.ref_x
        assert self.atom.y == self.ref_y
        assert self.atom.z == self.ref_z

    def test_ro_coords(self):
        with pytest.raises(AttributeError):
            self.atom.x = 0
        with pytest.raises(AttributeError):
            self.atom.y = 0
        with pytest.raises(AttributeError):
            self.atom.z = 0


class TestAtomGroup(TestCase):
    def setUp(self):
        self.atom1 = Atom(x=1)
        self.atom2 = Atom(x=2)
        self.atom3 = Atom(x=3)
        self.atoms = [self.atom1, self.atom2, self.atom3]

    def test_atomgroup_from_atoms(self):
        group = AtomGroup(atoms=self.atoms)
        assert len(group) == len(self.atoms)

        group = AtomGroup(atoms=set(self.atoms))
        assert len(group) == len(self.atoms)

        group = AtomGroup(atoms=tuple(self.atoms))
        assert len(group) == len(self.atoms)

    def test_atomgroup_creation_error(self):
        with pytest.raises(TypeError):
            AtomGroup(atoms="test")

        # iterable but not af atoms, bad!
        with pytest.raises(TypeError):
            AtomGroup(atoms=(self.atom1, None,))

    def test_atomgroup_append(self):
        group = AtomGroup()

        for atom in self.atoms:
            group.append(atom)

        assert len(group) == len(self.atoms)

        with pytest.raises(TypeError):
            group.append(None)


class TestPBCBox(TestCase):
    def setUp(self):
        self.ref_box = numpy.array([[22.25, 0., 0.], [0., 22.25, 0.], [0., 0., 6.36]])
        self.pbcbox = PBCBox(self.ref_box)

        self.coords_raw = numpy.array(
            [[0.47, 12.28, 5.69], [0.11, 12.61, 5.56], [0.2, 12.39, 5.17], [0.33, 12.43, 4.78], [0.47, 12.4, 4.41],
                [0.37, 12.4, 3.99], [0.47, 12.54, 3.6], [-0.02, 12.21, 5.09], [-0.08, 12.21, 4.65], [0.02, 12.19, 4.23],
                [-0.1, 11.95, 3.89], [-0.35, 11.95, 3.48]])

        self.coords_bbox = numpy.array(
            [[0.47, 12.28, 5.69], [0.11, 12.61, 5.56], [0.2, 12.39, 5.17], [0.33, 12.43, 4.78], [0.47, 12.4, 4.41],
                [0.37, 12.4, 3.99], [0.47, 12.54, 3.6], [22.23, 12.21, 5.09], [22.17, 12.21, 4.65], [0.02, 12.19, 4.23],
                [22.15, 11.95, 3.89], [21.90, 11.95, 3.48]])

    def test_pbc_box(self):
        assert_almost_equal(self.pbcbox.box, self.ref_box)

    def test_pbc_asarray(self):
        assert_almost_equal(self.ref_box, self.pbcbox.asarray())

    def test_pbc_put_in_bbox(self):
        coords_bbox = self.pbcbox.put_atoms_in_bbox(self.coords_raw)
        assert_almost_equal(coords_bbox, self.coords_bbox, 2)

    def test_pbc_dx(self):
        ref_dx = self.pbcbox.dx(self.coords_raw[0], self.coords_raw[1])
        dx = self.pbcbox.dx(self.coords_bbox[0], self.coords_bbox[1])

        assert_almost_equal(dx, ref_dx, 2)

    def test_pbc_distance(self):
        ref_dist = self.pbcbox.distance(self.coords_raw[0], self.coords_raw[1])
        dist = self.pbcbox.distance(self.coords_bbox[0], self.coords_bbox[1])

        assert_almost_equal(dist, ref_dist, 2)

    def test_pbc_xcm(self):
        ref_xcm = self.coords_raw.mean(axis=0)
        xcm = self.pbcbox.get_xcm(self.coords_raw)

        assert_almost_equal(xcm, ref_xcm, 2)

    def test_pbc_xcm_consistency(self):
        ref_xcm = self.pbcbox.get_xcm(self.coords_raw)
        xcm = self.pbcbox.get_xcm(self.coords_bbox)

        assert_almost_equal(xcm, ref_xcm, 2)

class TestTopology(TestCase):
    def setUp(self):
        self.topol = Topology()

        self.atoms = [(1, "DPPC", "NC3", 1),
                      (1, "DPPC", "PO4", 2),
                      (1, "DPPC", "GL0", 3),
                      (1, "DPPC", "CA1", 4),
                      (2, "DPPG", "NC3", 5),
                      (2, "DPPG", "PO4", 6),
                      (3, "DPPE ", "NC3", 7)]
        self.ref_atom = self.atoms[1]
        self.residue = ("DPPC1", [1, 2, 3, 4])
        self.resname = "DPPE3"
        self.size = len(self.atoms)

        for atom in self.atoms:
            self.topol.append(*atom)

    def test_size(self):
        assert len(self.topol) == self.size

    def test_atom_access(self):
        assert self.topol.get_atom(2) == self.ref_atom
        assert self.topol.get_atom(2.2) == self.ref_atom

        with pytest.raises(TypeError):
            self.topol.get_atom("1")

        with pytest.raises(KeyError):
            self.topol.get_atom(0)

        with pytest.raises(KeyError):
            self.topol.get_atom(9999)

    def test_residue_access(self):
        residue = self.topol.get_residue_from_atomid(2)
        assert residue["name"] == self.residue[0]
        assert list(residue["atomids"]) == self.residue[1]

        residue = self.topol.get_residue_from_atomid(2.2)
        assert residue["name"] == self.residue[0]
        assert list(residue["atomids"]) == self.residue[1]

        with pytest.raises(TypeError):
            self.topol.get_residue_from_atomid("1")

        with pytest.raises(KeyError):
            self.topol.get_residue_from_atomid(0)

        with pytest.raises(KeyError):
            self.topol.get_residue_from_atomid(9999)

    def test_name_stripping(self):
        assert self.topol.get_residue_from_atomid(7)["name"] == self.resname
