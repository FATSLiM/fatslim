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

import os
import MDAnalysis as mda
import numpy as np

from ._core import SimplifiedLipid, LipidRegistry


class Lipid(mda.core.groups.ComponentBase, SimplifiedLipid):
    def __new__(cls, atoms: mda.AtomGroup, hg_atoms: mda.AtomGroup):
        # this method is inherited from MDAnalysis.core.groups._MutableBase which updates classes accessible to
        # MDAnalysis.Universe objects.
        # To avoid polluting pure MDAnalysis object with FATSLiM objects, this method is overridden and stripped
        # down to the bare minimum

        return SimplifiedLipid.__new__(cls)

    def __init__(self, atoms: mda.AtomGroup, hg_atoms: mda.AtomGroup):
        try:
            assert isinstance(atoms, mda.AtomGroup)
        except AssertionError:
            raise TypeError("atoms must be a MDAnalysis.AtomGroup. (Actual type: {})".format(
                type(atoms)
            ))
        mda.core.groups.ComponentBase.__init__(self, atoms.ix, atoms.universe)

        SimplifiedLipid.__init__(self, atoms, hg_atoms)

    def __len__(self):
        return self.n_atoms

    @property
    def n_atoms(self):
        return len(self.atoms)

    @property
    def atoms(self) -> mda.AtomGroup:
        """An :class:`MDAnalysis.AtomGroup` of :class:`Atoms<MDAnalysis.core.groups.Atom>` present in this
        :class:`Lipid`.
        """
        return self._atoms

    @property
    def resid(self) -> int:
        return self.residue.resid

    @property
    def resindex(self) -> int:
        return self.residue.resindex

    @property
    def residue(self) -> int:
        return self._atoms.residues[0]


class Leaflet(mda.core.groups.GroupBase):
    pass


class LipidSystem(LipidRegistry):
    def __init__(self, universe: mda.Universe, headgroup_atoms: str, headgroup_index_group: str = "headgroups"):
        super().__init__(universe)

        try:
            assert isinstance(headgroup_atoms, str)
        except AssertionError:
            raise TypeError("headgroup_atoms argument must be a string. (Actual type: {})".format(
                str(type(headgroup_atoms))
            ))

        if os.path.isfile(headgroup_atoms):  # needs to read selection from file
            if os.path.splitext(headgroup_atoms)[1] != ".ndx":
                raise ValueError("Only Gromacs .ndx files are supported!")
            indices = []
            section = None
            with open(headgroup_atoms) as fp:
                for line in fp:
                    line = line.strip()
                    if line.startswith("["):
                        section = line.strip(" []")
                    elif section == headgroup_index_group:
                        indices.extend([int(val) for val in line.split()])
            if len(indices) == 0:
                raise ValueError(".ndx does not contain group named '{}'!".format(headgroup_index_group))
            indices = np.array(indices) - 1  # Substract 1 to the atom indices as gromacs starts counting at 1
            self.hg_selection = self.universe.atoms[indices]

        else:  # headgroup_atoms is a selection string
            try:
                self.hg_selection = self.universe.select_atoms(headgroup_atoms)
            except mda.exceptions.SelectionError as excinfo:
                raise ValueError("Bad headgroup selection string: '{}' (Actual error: '{}')".format(
                    headgroup_atoms, excinfo.args[0]
                ))

        # We need to check if the headgroup selection is coherent
        if len(self.hg_selection) == 0:  # Empty selection
            raise ValueError("Empty headgroup selection")
        elif len(self.hg_selection) == len(self.universe.atoms):  # No selection
            raise ValueError("Headgroup selection is whole universe")

        # We check if the selection correspond to just a few atoms from distinct residues
        for resindex, group in self.hg_selection.groupby("resindices").items():
            lipid = self.universe.residues[resindex].atoms
            try:
                assert len(group) < len(lipid)
            except AssertionError:
                raise ValueError("Headgroup selection corresponds to whole residue")

            self.add_lipid(Lipid(lipid, group))

    def __len__(self):
        return self.nlipids

    def __getitem__(self, item):
        return self.lipids[item]
