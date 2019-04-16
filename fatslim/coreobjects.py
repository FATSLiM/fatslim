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

from ._core import SimplifiedLipid, LipidRegistry


class Lipid(mda.core.groups.ComponentBase, SimplifiedLipid):
    def __new__(cls, *args, **kwargs):
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


class Leaflet(mda.core.groups.GroupBase):
    pass


class LipidSystem(LipidRegistry):
    def __init__(self, universe: mda.Universe, headgroup_atoms: str):
        super().__init__(universe)

        try:
            assert isinstance(headgroup_atoms, str)
        except AssertionError:
            raise TypeError("headgroup_atoms argument must be a string. (Actual type: {})".format(
                str(type(headgroup_atoms))
            ))

        if os.path.isfile(headgroup_atoms):  # needs to read selection from file
            raise NotImplementedError
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
        for resid, group in self.hg_selection.groupby("resids").items():
            lipid = self.universe.residues[resid - 1].atoms  # resids start from 1!
            try:
                assert len(group) < len(lipid)
            except AssertionError:
                raise ValueError("Headgroup selection corresponds to whole residue")

            self.add_lipid(Lipid(lipid, group))

    def __len__(self):
        return self.nlipids

    def __getitem__(self, item):
        return self.lipids[item]
