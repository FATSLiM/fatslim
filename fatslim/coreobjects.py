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


class Lipid(mda.core.groups.ComponentBase):

    def __init__(self, atoms: mda.AtomGroup, hg_atoms: mda.AtomGroup):
        super().__init__(atoms.ix, atoms.universe)

        try:
            assert isinstance(hg_atoms, mda.AtomGroup)
        except AssertionError:
            raise TypeError("hg_atoms must be a MDAnalysis.AtomGroup. (Actual type: {})".format(
                type(atoms)
            ))
        self._hg_atoms = hg_atoms

    def __len__(self):
        return self.n_atoms

    @property
    def n_atoms(self):
        return len(self.atoms)

    @property
    def atoms(self):
        """An :class:`MDAnalysis.AtomGroup` of :class:`Atoms<MDAnalysis.core.groups.Atom>` present in this
        :class:`Lipid`.
        """
        ag = self.universe.atoms[self.ix]
        ag._cache['isunique'] = True
        ag._cache['unique'] = ag
        return ag

    @property
    def position(self):
        return self._hg_atoms.positions.mean(axis=0)

    @property
    def direction(self):
        raise NotImplementedError

    @property
    def normal(self):
        raise NotImplementedError


class Leaflet(mda.core.groups.GroupBase):
    pass


class LipidSystem(object):
    def __init__(self, universe: mda.Universe, headgroup_atoms: str, verbose: bool = False):
        self.verbose = verbose

        try:
            assert isinstance(universe, mda.Universe)
        except AssertionError:
            raise TypeError("First argument must be an instance of MDAnalysis.Universe. (Actual type: {})".format(
                str(type(universe))
            ))
        self.universe = universe

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

        # We need to update universe so it can handle the Lipid custom class
        # See MDAnalysis.core.groups._MutableBase class
        # and _generate_from_topology method from MDAnalysis.core.universe.Universe for details

        self.lipids = []

        # We need to check if the headgroup selection is coherent
        if len(self.hg_selection) == 0:  # Empty selection
            raise ValueError("Empty headgroup selection")
        elif len(self.hg_selection) == len(self.universe.atoms):  # No selection
            raise ValueError("Headgroup selection is whole universe")
        else:
            # We check if the selection correspond to just a few atoms from distinct residues
            for resid, group in self.hg_selection.groupby("resids").items():
                lipid = self.universe.residues[resid - 1].atoms  # resids start from 1!
                try:
                    assert len(group) < len(lipid)
                except AssertionError:
                    raise ValueError("Headgroup selection corresponds to whole residue")

                self.lipids.append(Lipid(lipid, group))

        print(self.hg_selection, len(self.hg_selection))

    def __len__(self):
        return len(self.lipids)

    def __getitem__(self, item):
        return self.lipids[item]
