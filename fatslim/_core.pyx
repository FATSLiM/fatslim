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

# Cython directives
# cython: language_level=3

cimport cython

DEF NOTSET = -12345
DEF DIM = 3
DEF XX = 0
DEF YY = 1
DEF ZZ = 2
DEF EPSILON = 1e-6


import numpy as np
cimport numpy as np

import MDAnalysis as mda
from MDAnalysis.lib.mdamath import triclinic_vectors
import warnings

from ._typedefs cimport matrix, real, rvec, dreal, fsl_int, ivec
from ._typedefs cimport rvec_norm2, rvec_smul, rvec_copy, rvec_inc, rvec_normalize


###############################
# Utility class to handle PBC #
###############################
cdef struct cPBCBox_t:
    matrix     box
    rvec       fbox_diag
    rvec       hbox_diag
    rvec       mhbox_diag
    dreal      max_cutoff2

# Class to handle PBC calculations
@cython.initializedcheck(False)
@cython.boundscheck(False)
cdef class PBCBox(object):
    """
    Cython implementation of
    `PBC-related <https://en.wikipedia.org/wiki/Periodic_boundary_conditions>`_
    operations. This class is used by classes :class:`FastNS`
    and :class:`_NSGrid` to put all particles inside a brick-shaped box
    and to compute PBC-aware distance. The class can also handle
    non-PBC aware distance evaluations through ``periodic`` argument.

    .. warning::
        This class is not meant to be used by end users.

    .. warning::
        Even if MD triclinic boxes can be handled by this class,
        internal optimization is made based on the assumption that
        particles are inside a brick-shaped box. When this is not
        the case, calculated distances are not
        warranted to be exact.
    """

    # Cdef attributes
    cdef cPBCBox_t c_pbcbox
    cdef bint is_triclinic

    def __init__(self, real[:, ::1] box):
        """
        Parameters
        ----------
        box : numpy.ndarray
            box vectors of shape ``(3, 3)`` or
            as returned by ``MDAnalysis.lib.mdamath.triclinic_vectors``
            ``dtype`` must be ``numpy.float32``
        periodic : boolean
            ``True`` for PBC-aware calculations
            ``False`` for non PBC aware calculations
        """

        self.is_triclinic = False
        self.update(box)

    cdef void fast_update(self, real[:, ::1] box) nogil:
        """
        Updates the internal box parameters for
        PBC-aware distance calculations. The internal
        box parameters are used to define the brick-shaped
        box which is eventually used for distance calculations.

        """
        cdef fsl_int i, j
        cdef dreal min_hv2, min_ss, tmp

        # Update matrix
        self.is_triclinic = False
        for i in range(DIM):
            for j in range(DIM):
                self.c_pbcbox.box[i][j] = box[i, j]

                if i != j:
                    # mdamath.triclinic_vectors explicitly sets the off-diagonal
                    # elements to zero if the box is orthogonal, so we can
                    # safely check floating point values for equality here
                    if box[i, j] != 0.0:
                        self.is_triclinic = True

        # Update diagonals
        for i in range(DIM):
            self.c_pbcbox.fbox_diag[i] = box[i, i]
            self.c_pbcbox.hbox_diag[i] = self.c_pbcbox.fbox_diag[i] * 0.5
            self.c_pbcbox.mhbox_diag[i] = - self.c_pbcbox.hbox_diag[i]

        # Update maximum cutoff

        # Physical limitation of the cut-off
        # by half the length of the shortest box vector.
        min_hv2 = min(0.25 * rvec_norm2(&box[XX, XX]), 0.25 * rvec_norm2(&box[YY, XX]))
        min_hv2 = min(min_hv2, 0.25 * rvec_norm2(&box[ZZ, XX]))

        # Limitation to the smallest diagonal element due to optimizations:
        # checking only linear combinations of single box-vectors (2 in x)
        # in the grid search and pbc_dx is a lot faster
        # than checking all possible combinations.
        tmp = box[YY, YY]
        if box[ZZ, YY] < 0:
            tmp -= box[ZZ, YY]
        else:
            tmp += box[ZZ, YY]

        min_ss = min(box[XX, XX], min(tmp, box[ZZ, ZZ]))
        self.c_pbcbox.max_cutoff2 = min(min_hv2, min_ss * min_ss)

    def update(self, real[:, ::1] box):
        """
        Updates internal MD box representation and parameters used for calculations.

        Parameters
        ----------
        box : numpy.ndarray
            Describes the MD box vectors as returned by
            :func:`MDAnalysis.lib.mdamath.triclinic_vectors`.
            `dtype` must be :class:`numpy.float32`

        Note
        ----
        Call to this method is only needed when the MD box is changed
        as it always called when class is instantiated.

        """

        if box.shape[0] != DIM or box.shape[1] != DIM:
            raise ValueError("Box must be a {} x {} matrix. Got: {} x {})".format(
                DIM, DIM, box.shape[0], box.shape[1]))
        if (box[XX, XX] < EPSILON) or (box[YY, YY] < EPSILON) or (box[ZZ, ZZ] < EPSILON):
            raise ValueError("Box does not correspond to PBC=xyz")
        self.fast_update(box)

    cdef void fast_pbc_dx(self, rvec ref, rvec other, rvec dx) nogil:
        """Dislacement between two points for both
        PBC and non-PBC conditions

        Modifies the displacement vector between two points based
        on the minimum image convention for PBC aware calculations.

        For non-PBC aware distance evaluations, calculates the
        displacement vector without any modifications
        """

        cdef fsl_int i, j

        for i in range(DIM):
            dx[i] = other[i] - ref[i]

        for i in range(DIM-1, -1, -1):
            while dx[i] > self.c_pbcbox.hbox_diag[i]:
                for j in range(i, -1, -1):
                    dx[j] -= self.c_pbcbox.box[i][j]

            while dx[i] <= self.c_pbcbox.mhbox_diag[i]:
                for j in range(i, -1, -1):
                    dx[j] += self.c_pbcbox.box[i][j]

    cdef dreal fast_distance2(self, rvec a, rvec b) nogil:
        """Distance calculation between two points
        for both PBC and non-PBC aware calculations

        Returns the distance obeying minimum
        image convention if periodic is set to ``True`` while
        instantiating the :class:`_PBCBox` object.
        """

        cdef rvec dx
        self.fast_pbc_dx(a, b, dx)
        return rvec_norm2(dx)

    cdef void fast_put_atoms_in_bbox(self, real[:, ::1] coords, real[:, ::1] bbox_coords) nogil:
        """Shifts all ``coords`` to an orthogonal brick shaped box

        All the coordinates are brought into an orthogonal
        box. The box vectors for the brick-shaped box
        are defined in ``fast_update`` method.

        """

        cdef fsl_int i, m, d, natoms

        natoms = coords.shape[0]

        if self.is_triclinic:
            for i in range(natoms):
                for m in range(DIM - 1, -1, -1):
                    while bbox_coords[i, m] < 0:
                        for d in range(m+1):
                            bbox_coords[i, d] += self.c_pbcbox.box[m][d]
                    while bbox_coords[i, m] >= self.c_pbcbox.box[m][m]:
                        for d in range(m+1):
                            bbox_coords[i, d] -= self.c_pbcbox.box[m][d]
        else:
            for i in range(natoms):
                for m in range(DIM):
                    while bbox_coords[i, m] < 0:
                        bbox_coords[i, m] += self.c_pbcbox.box[m][m]
                    while bbox_coords[i, m] >= self.c_pbcbox.box[m][m]:
                        bbox_coords[i, m] -= self.c_pbcbox.box[m][m]


    cdef void fast_pbc_xcm(self, real[:, ::1] coords, rvec xcm, fsl_int[:] indices) nogil:
        self.fast_pbc_xcm_from_ref(coords, &coords[0, XX], xcm, indices)

    @cython.cdivision(True)
    cdef void fast_pbc_xcm_from_ref(self, real[:, ::1] coords, rvec ref,
                                    rvec xcm, fsl_int[:] indices) nogil:
        cdef fsl_int i
        cdef fsl_int size
        cdef rvec dx
        cdef real[:,::1] bbox_coords
        cdef fsl_int actual_index

        size = indices.shape[0]

        if size < 1:
            return

        rvec_copy(ref, xcm)
        rvec_smul(size, xcm, xcm)

        # Step 1: Get XCM
        for i in range(size):
            actual_index = indices[i]
            self.fast_pbc_dx(ref, &coords[actual_index, XX], dx)
            rvec_inc(xcm, dx)
        rvec_smul(1.0/size, xcm, xcm)

        # Step 2: Make sure it is inside the brick-shaped box
        for i in range(DIM - 1, -1, -1):
            while xcm[i] < 0:
                for d in range(i+1):
                    xcm[d] += self.c_pbcbox.box[i][d]
            while xcm[i] >= self.c_pbcbox.box[i][i]:
                for d in range(i+1):
                    xcm[d] -= self.c_pbcbox.box[i][d]


cdef class _NSGrid(object):
    # Cdef attributes
    cdef PBCBox box
    cdef dreal cutoff  # cutoff
    cdef dreal optimized_cutoff
    cdef fsl_int size  # total cells
    cdef fsl_int ncoords  # number of coordinates
    cdef fsl_int max_size
    # noinspection PyUnresolvedReferences
    cdef fsl_int[DIM] ncells  # individual cells in every dimension
    # noinspection PyUnresolvedReferences
    cdef fsl_int[DIM] cell_offsets  # Cell Multipliers
    # cellsize MUST be double precision, otherwise coord2cellid() may fail for
    # coordinates very close to the upper box boundaries! See MDAnalysis issue #2132
    # noinspection PyUnresolvedReferences
    cdef dreal[DIM] cellsize  # cell size in every dimension
    cdef fsl_int[:, ::1] beads_in_cell  # size (list of bead in every cell)
    cdef fsl_int max_nbeads
    cdef fsl_int[:] nbeads_in_cell
    cdef fsl_int[:] cell_lastcheckid
    cdef fsl_int[:] cellids  # ncoords (Cell occupation id for every atom)


    def __init__(self, fsl_int ncoords, real cutoff, PBCBox box, fsl_int max_size):

        cdef fsl_int i
        cdef fsl_int ncellx, ncelly, ncellz
        cdef fsl_int xi, yi, zi
        cdef real bbox_vol
        cdef dreal relative_cutoff_margin

        self.box = box

        self.ncoords = ncoords

        self.max_size = max_size

        self.cutoff = cutoff
        self.optimized_cutoff = cutoff

        self.size = -1

        # Allocate memory
        self.cellids = np.empty(self.ncoords, dtype=np.int)

        self.max_nbeads = 50

        self.nbeads_in_cell = None
        self.beads_in_cell = None
        self.cell_lastcheckid = None

        #print("DEBUG: Grid initialized for a maximum of {} beads and a cutoff of {:.3f}".format(self.ncoords, self.cutoff))

    cdef update(self):
        cdef fsl_int i

        # First, we add a small margin to the cell size so that we can safely
        # use the condition d <= cutoff (instead of d < cutoff) for neighbor
        # search.
        relative_cutoff_margin = 1.0e-8
        while self.optimized_cutoff == self.cutoff:
            self.optimized_cutoff = self.cutoff * (1.0 + relative_cutoff_margin)
            relative_cutoff_margin *= 10.0
        bbox_vol = self.box.c_pbcbox.box[XX][XX] * self.box.c_pbcbox.box[YY][YY] * self.box.c_pbcbox.box[YY][YY]
        while bbox_vol/self.optimized_cutoff**3 > self.max_size:
            self.optimized_cutoff *= 1.2

        #if self.optimized_cutoff > self.cutoff:
        #    print("DEBUG: optimized cutoff for grid : {:.3f} instead of requested {:.3f}".format(self.optimized_cutoff,
        #                                                                                        self.cutoff))

        for i in range(DIM):
            self.ncells[i] = <fsl_int> (self.box.c_pbcbox.box[i][i] / self.optimized_cutoff)
            self.cellsize[i] = self.box.c_pbcbox.box[i][i] / self.ncells[i]

        new_size = self.ncells[XX] * self.ncells[YY] * self.ncells[ZZ]

        self.cell_offsets[XX] = 0
        self.cell_offsets[YY] = self.ncells[XX]
        self.cell_offsets[ZZ] = self.ncells[XX] * self.ncells[YY]

        if new_size > self.size:
            self.nbeads_in_cell = np.zeros(new_size, dtype=np.int)
            self.beads_in_cell = np.empty((new_size, self.max_nbeads), dtype=np.int)
            self.cell_lastcheckid = np.empty(new_size, dtype=np.int)

        self.size = new_size

        #print("DEBUG: Grid updated to new size of {}x{}x{}={} cells".format(self.ncells[XX],
        #                                                                    self.ncells[YY],
        #                                                                    self.ncells[ZZ],
        #                                                                    self.size))



    cdef fsl_int coord2cellid(self, rvec coord) nogil:
        """Finds the cell-id for the given coordinate inside the brick shaped box

        Note
        ----
        Assumes the coordinate is already inside the brick shaped box.
        Return wrong cell-id if this is not the case
        """
        return <fsl_int> (coord[ZZ] / self.cellsize[ZZ]) * (self.cell_offsets[ZZ]) +\
               <fsl_int> (coord[YY] / self.cellsize[YY]) * self.cell_offsets[YY] + \
               <fsl_int> (coord[XX] / self.cellsize[XX])

    cdef bint cellid2cellxyz(self, fsl_int cellid, ivec cellxyz) nogil:
        """Finds actual cell position `(x, y, z)` from a cell-id
        """

        if cellid < 0:
            return False
        if cellid >= self.size:
            return False

        cellxyz[ZZ] = <fsl_int> (cellid / self.cell_offsets[ZZ])
        cellid -= cellxyz[ZZ] * self.cell_offsets[ZZ]

        cellxyz[YY] = <fsl_int> (cellid / self.cell_offsets[YY])
        cellxyz[XX] = cellid - cellxyz[YY] * self.cell_offsets[YY]

        return True


    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    cdef bint fill_grid(self, real[:, ::1] coords):
        """Sorts atoms into cells based on their position in the brick shaped box

        Every atom inside the brick shaped box is assigned a
        cell-id based on its position. Another list ``beadids``
        sort the atom-ids in each cell.

        Note
        ----
        The method fails if any coordinate is outside the brick shaped box.

        """

        cdef fsl_int i, cellindex = -1
        cdef fsl_int ncoords = coords.shape[0]
        cdef bint max_reached = False

        #print("DEBUG: fill_grid called with {} coords".format(ncoords))
        if self.ncoords < ncoords:
            return False

        # Update grid according to current box size
        self.update()

        # (Re)initialize buffers
        for i in range(self.size):
            self.nbeads_in_cell[i] = 0
            self.cell_lastcheckid[i] = -1

        # Find cellindex for each bead
        for i in range(ncoords):
            cellindex = self.coord2cellid(&coords[i, XX])

            self.cellids[i] = cellindex

            if self.nbeads_in_cell[cellindex] == self.max_nbeads:
                self.resize_beadlist()

            self.beads_in_cell[cellindex, self.nbeads_in_cell[cellindex]] = i
            self.nbeads_in_cell[cellindex] += 1

        return True

    cdef resize_beadlist(self):
        cdef fsl_int[:, ::1] tmp_memview = self.beads_in_cell
        cdef fsl_int increment = 50

        self.max_nbeads = tmp_memview.shape[1] + increment
        self.beads_in_cell = np.empty((tmp_memview.shape[0], self.max_nbeads), dtype=np.int)
        self.beads_in_cell[:,:-increment] = tmp_memview


cdef class _NSResults:
    cdef fsl_int size
    cdef fsl_int max_nneighbours
    cdef fsl_int[:] nneighbours
    cdef fsl_int[:, ::1] neighbours
    cdef real[:, ::1] distances

    def __init__(self, fsl_int size):

        self.size = size
        self.max_nneighbours = 50

        self.nneighbours = np.zeros(self.size, dtype=np.int)
        self.neighbours = np.empty((self.size, self.max_nneighbours), dtype=np.int)
        self.distances = np.empty((self.size, self.max_nneighbours), dtype=np.float32)

    cdef add_neighbour(self, fsl_int i, fsl_int j, real d2):
        if self.nneighbours[i] == self.max_nneighbours:
            self.resize_neighbourlist()

        self.neighbours[i, self.nneighbours[i]] = j
        self.distances[i, self.nneighbours[i]] = d2
        self.nneighbours[i] += 1

        #if i == 0:
        #    print("DEBUG: Bead#{} has neighbor bead#{} (d2={}, total neighbors={})".format(i,j,d2, self.nneighbours[i]))


    cdef reset(self):
        for i in range(self.size):
            self.nneighbours[i] = 0

    cdef resize_neighbourlist(self):
        cdef fsl_int[:, ::1] tmp_memview = self.neighbours
        cdef real[:, ::1] tmp_memview2 = self.distances
        cdef fsl_int increment = 50

        self.max_nneighbours = tmp_memview.shape[1] + increment

        self.neighbours = np.empty((tmp_memview.shape[0], self.max_nneighbours), dtype=np.int)
        self.neighbours[:, :-increment] = tmp_memview

        self.distances = np.empty((tmp_memview2.shape[0], self.max_nneighbours), dtype=np.float32)
        self.distances[:, :-increment] = tmp_memview2

    def __getitem__(self, fsl_int item):
        return np.asarray(self.neighbours[item, :self.nneighbours[item]]).copy()

    def get_nth_distances(self, fsl_int item):
        return np.sqrt(np.asarray(self.distances[item, :self.nneighbours[item]]))

    def get_tuples(self):
        cdef fsl_int i, j

        tuples = []
        for i in range(self.size):
            tuples.append([])

            for j in range(self.nneighbours[i]):
                tuples[-1].append((self.neighbours[i, j], self.distances[i, j]))


cdef bint self_search(_NSGrid grid, _NSResults results, real[:, ::1] positions):
    cdef fsl_int i, j, m, d
    cdef fsl_int current_beadid, bid, cellindex, cellindex_probe
    cdef fsl_int xi, yi, zi
    cdef rvec probe
    cdef dreal cutoff2 = grid.cutoff * grid.cutoff

    # Empty results
    results.reset()

    # Perform the neighbor search
    for i in range(positions.shape[0]):
        # Start with first search coordinate
        current_beadid = i

        # find the cellindex of the coordinate
        cellindex = grid.cellids[current_beadid]
        for xi in range(DIM):
            for yi in range(DIM):
                for zi in range(DIM):

                    # Calculate and/or reinitialize shifted coordinates
                    # Probe the search coordinates in a brick shaped box
                    probe[XX] = positions[current_beadid, XX] + (xi - 1) * grid.cellsize[XX]
                    probe[YY] = positions[current_beadid, YY] + (yi - 1) * grid.cellsize[YY]
                    probe[ZZ] = positions[current_beadid, ZZ] + (zi - 1) * grid.cellsize[ZZ]

                    # Make sure the shifted coordinates is inside the brick-shaped box
                    for m in range(DIM - 1, -1, -1):
                        while probe[m] < 0:
                            for d in range(m+1):
                                probe[d] += grid.box.c_pbcbox.box[m][d]
                        while probe[m] >= grid.box.c_pbcbox.box[m][m]:
                            for d in range(m+1):
                                probe[d] -= grid.box.c_pbcbox.box[m][d]

                    # Get the cell index corresponding to the probe
                    cellindex_probe = grid.coord2cellid(probe)

                    #if i == 0:
                    #    print("\nDEBUG: Probing cellindex #{} (Probe coordinates: [{:.3f}, {:.3f}, {:.3f}], offsets:[{},{},{}], actual cellindex:{}):".format(
                    #        cellindex_probe,
                    #        probe[XX], probe[YY], probe[ZZ],
                    #        xi-1, yi-1, zi-1,
                    #        cellindex
                    #    ))

                    if grid.cell_lastcheckid[cellindex_probe] == current_beadid:
                        #if i == 0:
                        #    print("DEBUG: No need -> This cell was already checked!")
                        continue
                    grid.cell_lastcheckid[cellindex_probe] = current_beadid

                    # for this cellindex search in grid
                    for j in range(grid.nbeads_in_cell[cellindex_probe]):
                        bid = grid.beads_in_cell[cellindex_probe, j]
                        if bid < current_beadid:
                            continue

                        # find distance between search coords[i] and coords[bid]
                        d2 = grid.box.fast_distance2(&positions[current_beadid, XX], &positions[bid, XX])
                        if EPSILON < d2 <= cutoff2:
                            results.add_neighbour(current_beadid, bid, d2)
                            results.add_neighbour(bid, current_beadid, d2)
    return True


cdef class SimplifiedLipid:
    # Cdef attributes
    cdef fsl_int[:] _ix
    cdef fsl_int _regid
    cdef LipidRegistry _registry

    def __init__(self, atoms: mda.AtomGroup, hg_atoms: mda.AtomGroup):
        self._ix = atoms.ix
        self._atoms = atoms.universe.atoms[self._ix]
        self._atoms._cache['isunique'] = True
        self._atoms._cache['unique'] = self._atoms

        try:
            assert isinstance(hg_atoms, mda.AtomGroup)
        except AssertionError:
            raise TypeError("hg_atoms must be a MDAnalysis.AtomGroup. (Actual type: {})".format(
                type(hg_atoms)
            ))
        self._hg_atoms = hg_atoms

        self._registry = None
        self._regid = -1

    @property
    def hg_atoms(self):
        return self._hg_atoms

    @property
    def position(self):
        if self._registry:
            return self._registry.lipid_positions[self._regid]
        else:
            warnings.warn("Lipid does not belong to any registry. No fast calculation nor PBC-awareness available")
            return self._hg_atoms.positions.mean(axis=0)

    @property
    def direction(self):
        if self._registry:
            return self._registry.lipid_directions[self._regid]
        else:
            warnings.warn("Lipid does not belong to any registry. No fast calculation nor PBC-awareness available")
            direction = self._hg_atoms.positions.mean(axis=0) - self._atoms.positions.mean(axis=0)
            direction /= np.linalg.norm(direction)
            return direction

    @property
    def neighbours_ids(self):
        if self._registry:
            return sorted(self._registry.lipid_neighbours[self._regid])
        else:
            raise ValueError("Neighbours are not available if lipid does not belong to LipidSystem")

    @property
    def neighbours_distances(self):
        if self._registry:
            return sorted(zip(
                self._registry.lipid_neighbours[self._regid],
                self._registry.lipid_neighbours.get_nth_distances(self._regid)
            ))
        else:
            raise ValueError("Neighbours are not available if lipid does not belong to LipidSystem")

    @property
    def normal(self):
        if self._registry:
            return self._registry.lipid_normals[self._regid]
        else:
            raise ValueError("Normal is not available if lipid does not belong to LipidSystem")


cdef class LipidRegistry:
    # Cdef attributes
    cdef fsl_int _nlipids
    cdef fsl_int _lastupdate
    cdef bint _locked
    cdef PBCBox box

    cdef fsl_int max_gridsize
    cdef real ns_cutoff

    cdef real[:, ::1] universe_coords_bbox

    cdef fsl_int[:] lipid_offsets
    cdef fsl_int[:] lipid_indices

    cdef fsl_int[:] hg_offsets
    cdef fsl_int[:] hg_indices

    cdef real[:, ::1] _lipid_positions
    cdef real[:, ::1] _lipid_directions
    cdef real[:, ::1] _lipid_normals
    cdef _NSGrid _lipid_grid
    cdef _NSResults _lipid_neighbours

    def __init__(self, universe: mda.Universe, ns_cutoff=20.0, max_gridsize=5000):
        try:
            assert isinstance(universe, mda.Universe)
        except AssertionError:
            raise TypeError("universe must be an instance of MDAnalysis.Universe. (Actual type: {})".format(
                str(type(universe))
            ))
        self.universe = universe
        self.universe_coords_bbox = None

        self.lipids = []
        self._nlipids = 0

        # Headgroups
        self.hg_coords = None
        self.hg_coords_bbox = None
        self.hg_indices = np.empty(0, dtype=np.int)
        self.hg_offsets = np.empty(0, dtype=np.int)

        # Lipids
        self.lipid_indices = np.empty(0, dtype=np.int)
        self.lipid_offsets = np.empty(0, dtype=np.int)

        # Simplified lipids
        self._lipid_positions = None
        self._lipid_directions = None
        self._lipid_grid = None
        self._lipid_normals = None
        self._lipid_neighbours = None

        # Interacting groups
        self.interacting_group_coords = None
        self.interacting_group_coords_bbox = None

        # NS-related
        self.ns_cutoff = ns_cutoff
        self.max_gridsize = max_gridsize


        self._lastupdate = -1
        self._locked = False

        # PBCBox and NS preparation
        box = triclinic_vectors(self.universe.dimensions)
        self.box = PBCBox(box)

    def add_lipid(self, SimplifiedLipid lipid):
        try:
            assert self.universe == lipid._atoms.universe
        except AssertionError:
            raise ValueError("Lipid does not belong to the same universe")

        if self._locked:
            raise ValueError("Lipid must be added before accessing coordinates!")


        self.lipids.append(lipid)

        lipid._registry = self
        lipid._regid = self._nlipids
        self._nlipids += 1

        # Update indices for fast access
        self.lipid_offsets = np.concatenate((self.lipid_offsets, [self.lipid_indices.shape[0], ]))
        self.lipid_indices = np.concatenate((self.lipid_indices, lipid._ix))

        self.hg_offsets = np.concatenate((self.hg_offsets, [self.hg_indices.shape[0], ]))
        self.hg_indices = np.concatenate((self.hg_indices, lipid._hg_atoms._ix))

        # Force update
        self._lastupdate = -1

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    def update(self, force_update=False):
        cdef fsl_int i, offset
        cdef real[:, ::1] u_pos, positions, hg_positions_bbox
        cdef fsl_int[:] indices
        cdef fsl_int next_offset
        cdef bint bad_normals = False

        if self._lastupdate == self.universe.trajectory.frame and not force_update:
            #warnings.warn("Already uptodate.. No need for update")
            return
        #warnings.warn("Updating lipid registry")

        if not self._locked:
            self._lipid_positions = np.empty((self._nlipids, DIM), dtype=np.float32)
            self._lipid_directions = np.empty((self._nlipids, DIM), dtype=np.float32)
            self._lipid_normals = np.empty((self._nlipids, DIM), dtype=np.float32)
            self._lipid_grid = _NSGrid(self._nlipids, self.ns_cutoff, self.box, self.max_gridsize)
            self._lipid_neighbours = _NSResults(self._nlipids)
        self._locked = True

        # Retrieve coords from MDA Universe
        u_pos = np.ascontiguousarray(self.universe.trajectory.ts.positions[:])
        self.universe_coords_bbox = u_pos.copy()

        # Update PBC box
        box = triclinic_vectors(self.universe.dimensions)
        self.box.fast_update(box)

        self.box.fast_put_atoms_in_bbox(u_pos, self.universe_coords_bbox)

        for i in range(self._nlipids):

            # First: bead position
            if i == self._nlipids-1:
                next_offset = self.hg_indices.shape[0]
            else:
                next_offset = self.hg_offsets[i+1]
            indices = self.hg_indices[self.hg_offsets[i]: next_offset]
            self.box.fast_pbc_xcm(self.universe_coords_bbox, &self._lipid_positions[i, XX], indices)

            # Second: lipid directions
            if i == self._nlipids-1:
                next_offset = self.lipid_indices.shape[0]
            else:
                next_offset = self.lipid_offsets[i+1]
            indices = self.lipid_indices[self.lipid_offsets[i]: next_offset]
            self.box.fast_pbc_xcm(self.universe_coords_bbox, &self._lipid_directions[i, XX], indices)
            self.box.fast_pbc_dx(&self._lipid_directions[i, XX],
                                 &self._lipid_positions[i, XX],
                                 &self._lipid_directions[i, XX])
            rvec_normalize(&self._lipid_directions[i, XX])

        # Third: get neighbours for each lipid
        self._lipid_grid.fill_grid(self._lipid_positions)
        self_search(self._lipid_grid, self._lipid_neighbours, self._lipid_positions)

        # Fourth: get lipid normals

        self._lastupdate = self.universe.trajectory.frame

        return True

    cdef self_search(self, _NSGrid grid, _NSResults results):
        pass


    @property
    def lipid_positions(self) -> np.ndarray:
        self.update()
        return np.asarray(self._lipid_positions).copy()

    @property
    def lipid_directions(self) -> np.ndarray:
        self.update()
        return np.asarray(self._lipid_directions).copy()

    @property
    def lipid_normals(self) -> np.ndarray:
        self.update()
        return np.asarray(self._lipid_normals).copy()

    @property
    def lipid_neighbours(self) -> _NSResults:
        self.update()
        return self._lipid_neighbours

    @property
    def nlipids(self) -> int:
        return self._nlipids
