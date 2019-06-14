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
DEF PI = 3.141592653589793

import numpy as np
cimport numpy as np

import MDAnalysis as mda
from MDAnalysis.lib.mdamath import triclinic_vectors
import warnings

from ._typedefs cimport real, rvec, dreal, fsl_int, ivec, matrix
from ._typedefs cimport rvec_normalize, rvec_clear, mat_clear, rvec_dprod, rvec_inc, rvec_smul, rvec_copy

from ._typedefs cimport rvec_norm
from ._typedefs cimport mat_from_rvec, invert_mat

from ._geometry cimport PBCBox, normal_from_neighbours, curvature_from_neighbours
from ._geometry cimport complete_basis, rvec_to_basis
from ._geometry cimport real_point, Polygon, polygon_new, polygon_get_area, fast_clip_zoi, polygon_destroy, polygon_append, polygon_empty

from libc.math cimport acos


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
        self.cellids = np.empty(self.ncoords, dtype=int)

        self.max_nbeads = 50

        self.nbeads_in_cell = None
        self.beads_in_cell = None
        self.cell_lastcheckid = None

        #print("DEBUG: Grid initialized for a maximum of {} beads and a cutoff of {:.3f}".format(self.ncoords, self.cutoff))

    @cython.cdivision(True)
    cdef int update(self) nogil except -1:
        cdef fsl_int i
        cdef dreal bbox_vol

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
        #    with gil:
        #        print("DEBUG: optimized cutoff for grid : {:.3f} instead of requested {:.3f}".format(self.optimized_cutoff,
        #                                                                                        self.cutoff))

        for i in range(DIM):
            self.ncells[i] = <fsl_int> (self.box.c_pbcbox.box[i][i] / self.optimized_cutoff)
            self.cellsize[i] = self.box.c_pbcbox.box[i][i] / self.ncells[i]

        new_size = self.ncells[XX] * self.ncells[YY] * self.ncells[ZZ]

        self.cell_offsets[XX] = 0
        self.cell_offsets[YY] = self.ncells[XX]
        self.cell_offsets[ZZ] = self.ncells[XX] * self.ncells[YY]

        if new_size > self.size:
            with gil:
                self.nbeads_in_cell = np.zeros(new_size, dtype=int)
                self.beads_in_cell = np.empty((new_size, self.max_nbeads), dtype=int)
                self.cell_lastcheckid = np.empty(new_size, dtype=int)

        self.size = new_size

        #with gil:
        #    print("DEBUG: Grid updated to new size of {}x{}x{}={} cells".format(self.ncells[XX],
        #                                                                    self.ncells[YY],
        #                                                                    self.ncells[ZZ],
        #                                                                    self.size))
        #    print("DEBUG: Grid cell sizes: {}x{}x{}".format(self.cellsize[XX],
        #                                                    self.cellsize[YY],
        #                                                    self.cellsize[ZZ]))


    @cython.cdivision(True)
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

    @cython.cdivision(True)
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
    cdef int fill_grid(self, real[:, ::1] coords) nogil except -1:
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

            #with gil:
            #    print("DEBUG: Filling grid: Bead#{} is put inside cell#{}".format(i, cellindex))

            self.cellids[i] = cellindex

            if self.nbeads_in_cell[cellindex] == self.max_nbeads:
                with gil:
                    self.resize_beadlist()

            self.beads_in_cell[cellindex, self.nbeads_in_cell[cellindex]] = i
            self.nbeads_in_cell[cellindex] += 1

    cdef void resize_beadlist(self):
        cdef fsl_int[:, ::1] tmp_memview = self.beads_in_cell
        cdef fsl_int increment = 50

        self.max_nbeads = tmp_memview.shape[1] + increment
        self.beads_in_cell = np.empty((tmp_memview.shape[0], self.max_nbeads), dtype=int)
        self.beads_in_cell[:,:-increment] = tmp_memview

@cython.initializedcheck(False)
@cython.boundscheck(False)
cdef class _NSResults:
    cdef fsl_int size
    cdef fsl_int max_nneighbours
    cdef fsl_int[:] nneighbours
    cdef fsl_int[:, ::1] neighbours
    cdef real[:, ::1] distances

    def __init__(self, fsl_int size):

        self.size = size
        self.max_nneighbours = 50

        self.nneighbours = np.zeros(self.size, dtype=int)
        self.neighbours = np.empty((self.size, self.max_nneighbours), dtype=int)
        self.distances = np.empty((self.size, self.max_nneighbours), dtype=np.float32)

    cdef int add_neighbour(self, fsl_int i, fsl_int j, real d2) nogil except -1:
        if self.nneighbours[i] == self.max_nneighbours:
            with gil:
                self.resize_neighbourlist()

        self.neighbours[i, self.nneighbours[i]] = j
        self.distances[i, self.nneighbours[i]] = d2
        self.nneighbours[i] += 1

        #if i == 0:
        #    print("DEBUG: Bead#{} has neighbor bead#{} (d2={}, total neighbors={})".format(i,j,d2, self.nneighbours[i]))


    cdef void reset(self) nogil:
        cdef fsl_int i
        for i in range(self.size):
            self.nneighbours[i] = 0

    cdef resize_neighbourlist(self):
        cdef fsl_int[:, ::1] tmp_memview = self.neighbours
        cdef real[:, ::1] tmp_memview2 = self.distances
        cdef fsl_int increment = 50

        self.max_nneighbours = tmp_memview.shape[1] + increment

        self.neighbours = np.empty((tmp_memview.shape[0], self.max_nneighbours), dtype=int)
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

    @property
    def shape(self):
        return (self.size, self.max_nneighbours)

@cython.initializedcheck(False)
@cython.boundscheck(False)
cdef int fast_self_search(_NSGrid grid, _NSResults results, real[:, ::1] positions) nogil except -1:
    cdef fsl_int i, j, m, d
    cdef fsl_int current_beadid, bid, cellindex, cellindex_probe
    cdef ivec cellxyz
    cdef fsl_int xi, yi, zi
    cdef rvec probe
    cdef dreal cutoff2 = grid.cutoff * grid.cutoff
    cdef dreal d2

    # Empty results
    results.reset()

    # Perform the neighbor search
    for i in range(positions.shape[0]):
        # Start with first search coordinate
        current_beadid = i

        # find the cellindex of the coordinate
        cellindex = grid.cellids[current_beadid]
        grid.cellid2cellxyz(cellindex, cellxyz)
        for xi in range(DIM):
            for yi in range(DIM):
                for zi in range(DIM):

                    # Calculate and/or reinitialize shifted coordinates
                    # Probe the search coordinates in a brick shaped box
                    probe[XX] = (cellxyz[XX] + 0.5 + xi - 1) * grid.cellsize[XX]
                    probe[YY] = (cellxyz[YY] + 0.5 + yi - 1) * grid.cellsize[YY]
                    probe[ZZ] = (cellxyz[ZZ] + 0.5 + zi - 1) * grid.cellsize[ZZ]

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

                    #if current_beadid == 0:
                    #    with gil:
                    #        print("\nDEBUG: Probing cellindex #{} (Probe coordinates: [{:.3f}, {:.3f}, {:.3f}], offsets:[{},{},{}], actual cellindex:{}):".format(
                    #            cellindex_probe,
                    #            probe[XX], probe[YY], probe[ZZ],
                    #            xi-1, yi-1, zi-1,
                    #            cellindex
                    #        ))

                    if grid.cell_lastcheckid[cellindex_probe] == current_beadid:
                        #if current_beadid == 0:
                        #    with gil:
                        #        print("DEBUG: No need -> This cell was already checked!")
                        continue
                    grid.cell_lastcheckid[cellindex_probe] = current_beadid

                    # for this cellindex search in grid
                    for j in range(grid.nbeads_in_cell[cellindex_probe]):
                        bid = grid.beads_in_cell[cellindex_probe, j]
                        if bid < current_beadid:
                            #with gil:
                            #    if current_beadid == 0:
                            #        print("DEBUG: bead#{} is ignored (should be already tested and added if necessary".format(bid))
                            continue

                        # find distance between search coords[i] and coords[bid]
                        d2 = grid.box.fast_distance2(&positions[current_beadid, XX], &positions[bid, XX])
                        if EPSILON < d2 <= cutoff2:
                            results.add_neighbour(current_beadid, bid, d2)
                            results.add_neighbour(bid, current_beadid, d2)

                            #if current_beadid==0:
                            #    with gil:
                            #        print("DEBUG: Adding bead#{} as neighbour of bead#0".format(bid))
                        #else:
                        #    if current_beadid==0:
                        #        with gil:
                        #            print("DEBUG: bead#{} is too far from bead#0: d2={:.3f}, cutoff2={:.3f}".format(
                        #                bid,
                        #                d2,
                        #                cutoff2
                        #            ))
                    #if grid.nbeads_in_cell[cellindex_probe] == 0 and current_beadid==0:
                    #    with gil:
                    #        print("DEBUG: Cell is empty!")


cdef class SimplifiedLipid:
    # Cdef attributes
    cdef fsl_int[:] _ix
    cdef fsl_int _regid
    cdef LipidRegistry _registry
    cdef Leaflet _leaflet
    cdef fsl_int _regid_leaflet

    def __init__(self, atoms: mda.AtomGroup, hg_atoms: mda.AtomGroup):
        self._ix = atoms.ix
        self._atoms = atoms.universe.atoms[self._ix]
        self._atoms._cache['isunique'] = True
        self._atoms._cache['unique'] = self._atoms

        if len(atoms) == 0:
            raise ValueError("'atoms' group is empty")

        try:
            assert isinstance(hg_atoms, mda.AtomGroup)
        except AssertionError:
            raise TypeError("hg_atoms must be a MDAnalysis.AtomGroup. (Actual type: {})".format(
                type(hg_atoms)
            ))

        if len(hg_atoms) == 0:
            raise ValueError("'hg_atoms' group is empty")

        residues = atoms.residues
        if len(residues) > 1:
            raise ValueError("Only lipids belonging to one single residue are supported")

        resindex = residues[0].resindex

        hg_residues = hg_atoms.residues
        hg_resindex = hg_residues[0].resindex
        if len(hg_residues) != 1 or hg_resindex != resindex:
            raise ValueError("'hg_atoms' group is not consistent with 'atoms' group")

        self._hg_atoms = hg_atoms

        self._registry = None
        self._regid = -1

        self._leaflet = None
        self._regid_leaflet = -1

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
            direction = self._hg_atoms.positions.mean(axis=0) - \
                        self._atoms[self._atoms.ix >= self._hg_atoms.ix[0]].positions.mean(axis=0)
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

    @property
    def thickness(self):
        if self._leaflet is None:
            raise ("Thickness value is only available when lipid belong to a leaflet")
        else:
            return self._leaflet.lipid_thicknesses[self._regid_leaflet]

    @property
    def leaflet(self):
        return self._leaflet

    @property
    def index(self):
        return self._regid

    @property
    def leaflet_index(self):
        return self._regid_leaflet



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
    cdef readonly fsl_int[:] hg_indices

    cdef real[:, ::1] _lipid_positions
    cdef real[:, ::1] _lipid_centroids
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
        self.hg_indices = np.empty(0, dtype=int)
        self.hg_offsets = np.empty(0, dtype=int)

        # Lipids
        self.lipid_indices = np.empty(0, dtype=int)
        self.lipid_offsets = np.empty(0, dtype=int)

        # Simplified lipids
        self._lipid_positions = None
        self._lipid_centroids = None
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

        self._membranes = None
        self._lastupdate_membrane = -1

        self._aggregates = None
        self._lastupdate_aggregates = -1

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
        cdef fsl_int[:] indices, indices_directions
        cdef rvec cog_directions
        cdef fsl_int next_offset
        cdef bint should_raise = False

        if self._lastupdate == self.universe.trajectory.frame and not force_update:
            #warnings.warn("Already uptodate.. No need for update")
            return
        #warnings.warn("Updating lipid registry")

        if not self._locked:
            self._lipid_positions = np.empty((self._nlipids, DIM), dtype=np.float32)
            self._lipid_centroids = np.empty((self._nlipids, DIM), dtype=np.float32)
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

        with nogil:
            for i in range(self._nlipids):

                # First: bead position
                if i == self._nlipids-1:
                    next_offset = self.hg_indices.shape[0]
                else:
                    next_offset = self.hg_offsets[i+1]
                indices = self.hg_indices[self.hg_offsets[i]: next_offset]
                self.box.fast_pbc_centroid(self.universe_coords_bbox, &self._lipid_positions[i, XX], indices)

                # Second: lipid directions
                if i == self._nlipids-1:
                    next_offset = self.lipid_indices.shape[0]
                else:
                    next_offset = self.lipid_offsets[i+1]
                indices = self.lipid_indices[self.lipid_offsets[i]: next_offset]
                indices_directions = self.lipid_indices[self.hg_indices[self.hg_offsets[i]]: next_offset]
                self.box.fast_pbc_centroid_from_ref(self.universe_coords_bbox,
                                                    &self._lipid_positions[i, XX],
                                                    &self._lipid_centroids[i, XX],
                                                    indices)
                self.box.fast_pbc_centroid_from_ref(self.universe_coords_bbox,
                                                    &self._lipid_positions[i, XX],
                                                    cog_directions,
                                                    indices_directions)
                self.box.fast_pbc_dx(cog_directions,
                                     &self._lipid_positions[i, XX],
                                     &self._lipid_directions[i, XX])
                rvec_normalize(&self._lipid_directions[i, XX])

            # Third: get neighbours for each lipid
            self._lipid_grid.fill_grid(self._lipid_positions)
            fast_self_search(self._lipid_grid, self._lipid_neighbours, self._lipid_positions)

            # Fourth: get lipid normals
            for i in range(self._nlipids):
                normal_from_neighbours(self._lipid_positions,
                                   self._lipid_directions,
                                   i,
                                   self._lipid_neighbours.neighbours[i, :self._lipid_neighbours.nneighbours[i]],
                                   self.box,
                                   &self._lipid_normals[i, XX])
                if self._lipid_neighbours.nneighbours[i] < 3:
                    for j in range(DIM):
                        self._lipid_normals[i, j] = self._lipid_directions[i, j]


        self._lastupdate = self.universe.trajectory.frame

    # TODO: Temp stuff
    @property
    def curvatures(self):
        cdef fsl_int i, dimi, dimj
        cdef matrix eigvecs
        cdef rvec eigvals
        self.update()

        py_eigval = np.empty((self._nlipids, 3), dtype=np.float32)
        py_eigvecs = np.empty((self._nlipids, 3, 3))

        for i in range(self._nlipids):
            rvec_clear(eigvals)
            mat_clear(eigvecs)
            if not curvature_from_neighbours(self._lipid_positions,
                                      self._lipid_directions,
                                      i,
                                      self._lipid_neighbours.neighbours[i, :self._lipid_neighbours.nneighbours[i]],
                                      self.box,
                                      eigvals,
                                      eigvecs):
                print("WARNING: not able to get curvature for lipids #{}".format(i))

            for dimi in range(DIM):
                for dimj in range(DIM):
                    py_eigvecs[i][dimi][dimj] = eigvecs[dimi][dimj]
                py_eigval[i][dimi] = eigvals[dimi]


        return py_eigval, py_eigvecs

    #@cython.initializedcheck(False)
    #@cython.boundscheck(False)
    def find_membranes(self, force_update=False):
        cdef fsl_int bead_id, seed_id, n_edgers, n_leftovers, nid, ref_nid
        cdef real angle, min_angle, max_angle, angle_range, dprod_value
        cdef dreal d_val, best_d
        cdef fsl_int i, j
        cdef fsl_int stack_size
        cdef fsl_int[:] aggregate_lipid_ids
        cdef fsl_int maybe_leaflet_current_id, current_leaflet_size
        cdef LipidAggregate aggregate

        cdef fsl_int[:] edgers, maybe_leaflet_ids, stack, current_leaflet_ids, leftovers

        cdef real[:] normal

        cdef bint in_edger
        # No need to do anything if the membranes are already identified for the current frame
        if self._lastupdate_membrane == self.universe.trajectory.frame and not force_update:
            return

        # Make sure everything is updated
        self.update()

        # Initialize memory
        edgers = np.empty(self._nlipids, dtype=int)
        maybe_leaflet_ids = np.empty(self._nlipids, dtype=int)
        stack = np.empty(self._nlipids, dtype=int)
        current_leaflet_ids = np.empty(self._nlipids, dtype=int)
        leftovers = np.empty(self._nlipids, dtype=int)
        n_leftovers = 0


        # Step 1: Get "naive" aggregates (i.e. aggregates based only on distance)
        self.find_aggregates(force_update=force_update)

        # Step 2: Iterate over the known aggregates to find if they can be splitted into potential leaflets
        maybe_leaflets = []
        for aggregate in self._aggregates:
            aggregate_lipid_ids = aggregate.indices
            # Step 2.1: Find if some lipids are located on edges
            n_edgers = 0

            for i in range(aggregate_lipid_ids.shape[0]):
                bead_id = aggregate_lipid_ids[i]
                min_angle = PI + EPSILON
                max_angle = -EPSILON
                for j in range(self._lipid_neighbours.nneighbours[bead_id]):
                    nid = self._lipid_neighbours.neighbours[bead_id][j]
                    # Check the angle between the local normal and the Z axis
                    # The axis does not really matter as it is just the variation of that angle that matters
                    # But using the Z axis simplify things since dot product is trivial

                    angle = acos(self._lipid_normals[nid, ZZ])

                    if angle > max_angle:
                        max_angle = angle

                    if angle < min_angle:
                        min_angle = angle

                angle_range = max_angle - min_angle

                # if bead_id == DEBUG_BEAD:
                    # print("Checking if resid {} is an edger (angle range: {:.3f}°)".format(
                    #     self.lipids[bead_id].resid,
                    # angle_range))

                if angle_range > PI * 0.5:  #
                    edgers[n_edgers] = bead_id
                    n_edgers += 1

                    # if bead_id == DEBUG_BEAD:
                    #     print("Resid {} is an edger".format(self.lipids[bead_id].resid))
                # else:
                #     if bead_id == DEBUG_BEAD:
                #         print("Resid {} is NOT an edger".format(self.lipids[bead_id].resid))


            # Step 2.2: Loop over the lipids in the aggregate to check if they can fit if a potential leaflet

            maybe_leaflet_current_id = -1
            maybe_leaflet_ids[:] = maybe_leaflet_current_id   # (Re)initialize the found leaflets

            maybe_leaflets_from_aggregates = []
            for i in range(aggregate_lipid_ids.shape[0]):
                seed_id = aggregate_lipid_ids[i]

                if maybe_leaflet_ids[seed_id] > -1:  # This lipid already belong to an aggregate, we can skip it
                    continue

                # If this lipid is on an edge, we will deal with it later
                in_edger = False
                for j in range(n_edgers):
                    if edgers[j] == seed_id:
                        in_edger = True
                        break
                if in_edger:
                    continue

                # Increment the leaflet id as we are dealing with another aggregate
                maybe_leaflet_current_id += 1

                # print("\n\nCreating new aggregate #{} from resid {}".format(maybe_leaflet_current_id,
                #                                                             self.lipids[seed_id].resid))

                # Add the seed bead to the new leaflet
                maybe_leaflet_ids[seed_id] = maybe_leaflet_current_id

                current_leaflet_size = 1
                current_leaflet_ids[0] = seed_id

                # Reset the stack
                stack_size = 1
                stack[0] = seed_id

                # Search for all the lipids that belong to the current leaflet
                # This is done by exhausting the stack that contains all the candidates
                while stack_size > 0:

                    # Pop the stack
                    ref_nid = stack[stack_size - 1]
                    stack_size -= 1

                    for j in range(self._lipid_neighbours.nneighbours[ref_nid]):
                        nid = self._lipid_neighbours.neighbours[ref_nid][j]

                        if maybe_leaflet_ids[nid] > -1:  # This lipid already belong to an aggregate, we can skip it
                            continue

                        # if nid == DEBUG_BEAD:
                        #     print("Dealing with resid {}:".format(self.lipids[nid].resid))

                        # If this lipid is on an edge, we will deal with it later
                        in_edger = False
                        for j in range(n_edgers):
                            if edgers[j] == nid:
                                in_edger = True
                                break
                        if in_edger:
                            # if nid == DEBUG_BEAD:
                            #     print("resid {} is in edger".format(self.lipids[nid].resid))
                            continue

                        dprod_value = rvec_dprod(&self._lipid_normals[ref_nid, XX],
                                                 &self._lipid_normals[nid, XX])

                        # if nid == DEBUG_BEAD:
                        #     print("dot product for resid {}: {:.3f} (ref resid: {} from aggregate {})".format(
                        #         self.lipids[nid].resid,
                        #         dprod_value, self.lipids[ref_nid].resid, maybe_leaflet_current_id
                        #     ))

                        # Angle between the two normals must be smaller than 45°.
                        # Otherwise, it would correspond to a rather small curvature radius which is barely compatible
                        # with a lipid bilayer
                        if dprod_value <= 0.7071067811865476:
                            # if nid == DEBUG_BEAD:
                            #     print("Resid {} rejected!".format(self.lipids[nid].resid))
                            continue

                        # If still here, add bead to current leaflet
                        maybe_leaflet_ids[nid] = maybe_leaflet_current_id
                        current_leaflet_ids[current_leaflet_size] = nid
                        current_leaflet_size += 1

                        # Append bead to stack
                        stack[stack_size] = nid
                        stack_size += 1

                        # if nid == DEBUG_BEAD:
                        #     print("Resid {} accepted in aggregate {}".format(self.lipids[nid].resid, maybe_leaflet_current_id))

                # Store the incomplete leaflet
                maybe_leaflets_from_aggregates.append(current_leaflet_ids[:current_leaflet_size].copy())


            # Step 2.3: Handle the edgers to check if they can be added to current leaflets

            n_added = n_edgers
            n_pass = 0

            leaflet_neighbour = np.empty_like(edgers)
            leaflet_neighbour[:] = -1

            new_members = []
            for l in maybe_leaflets_from_aggregates:
                new_members.append([])

            n_leftovers = 0

            # print("{} edgers to merge into {} leaflet(s):".format(n_edgers, len(maybe_leaflets_from_aggregates)))

            for i in range(n_edgers):
                edger_id = edgers[i]

                if edger_id in range(870, 881):
                    print("\n")

                if True:
                    best_value = 0
                    best_d2min = 1000000

                    best_leaflet = -1

                    for lid, indices in enumerate(maybe_leaflets_from_aggregates):
                        current_stack = np.empty(self._nlipids, dtype=int)
                        current_stack[0] = edger_id
                        current_stack_size = 1

                        next_stack = np.empty(self._nlipids, dtype=int)

                        used_nodes = np.zeros(self._nlipids, dtype=bool)

                        closest_bead_ids = []
                        generation = 0

                        n_checks = 0

                        while len(closest_bead_ids) == 0 and generation < 10:
                            next_stack_size = 0

                            for j in range(current_stack_size):
                                ref_nid = current_stack[j]

                                used_nodes[ref_nid] = True

                                for k in range(self._lipid_neighbours.nneighbours[ref_nid]):
                                    n_checks += 1

                                    nid = self._lipid_neighbours.neighbours[ref_nid][k]

                                    if maybe_leaflet_ids[nid] == lid:
                                        closest_bead_ids.append(nid)
                                    elif maybe_leaflet_ids[nid] == -1 and not used_nodes[nid]:
                                        next_stack[next_stack_size] = nid
                                        next_stack_size += 1

                            generation += 1

                            current_stack = next_stack.copy()
                            current_stack_size = next_stack_size

                        # print("Took {} checks ({} generations) to find {} lipids from leaflet #{}".format(
                        #     n_checks,
                        #     generation,
                        #     len(closest_bead_ids),
                        #     lid
                        # ))

                        if len(closest_bead_ids) == 0:
                            continue


                        d2 = 0
                        dprod_value = 0
                        position = np.zeros(3, dtype=np.float32)
                        normal = np.zeros(3, dtype=np.float32)

                        for bead_id in closest_bead_ids:
                            d2 += self.box.fast_distance2(&self._lipid_positions[edger_id, XX],
                                                         &self._lipid_positions[bead_id, XX])

                            dprod_value += rvec_dprod(&self._lipid_normals[edger_id, XX],
                                                      &self._lipid_normals[bead_id, XX])

                            position += self._lipid_positions[bead_id]

                            for j in range(DIM):
                                normal[j] += self._lipid_normals[bead_id, j]

                        d2 /= len(closest_bead_ids)
                        dprod_value /= len(closest_bead_ids)
                        position /= len(closest_bead_ids)

                        for j in range(DIM):
                            normal[j] /= len(closest_bead_ids)


                        d_leaflet = np.abs(np.dot(normal, position - self._lipid_positions[edger_id]))

                        d2 = d_leaflet
                        if dprod_value > best_value * 0.9 and d2 < best_d2min:
                            best_value = dprod_value
                            best_d2min = d2
                            best_leaflet = lid

                    if best_leaflet > -1:
                        # print("Edger resid {} should be added to leaflet #{}".format(
                        #     self.lipids[edger_id].resid,
                        #     best_leaflet
                        # ))
                        new_members[best_leaflet].append(edger_id)
                    else:
                        # print("Edger resid {} should be added to leftovers".format(
                        #     self.lipids[edger_id].resid,
                        # ))
                        leftovers[n_leftovers] = edger_id
                        n_leftovers += 1

                    # print("")


            for leaflet_id in range(len(maybe_leaflets_from_aggregates)):
                if len(new_members[leaflet_id]) > 0:

                    temp = np.concatenate(
                        (maybe_leaflets_from_aggregates[leaflet_id], new_members[leaflet_id])
                    )
                    maybe_leaflets_from_aggregates[leaflet_id] = temp



            #leftovers = np.unique(leftovers[:n_leftovers])
            # print("{} leftovers: resid {}".format(
            #     n_leftovers,
            #     " ".join([str(self.lipids[val].resid) for val in np.sort(leftovers[:n_leftovers])])
            # ))

            # Step 2.4: Store the completed leaflets
            for ids in maybe_leaflets_from_aggregates:
                maybe_leaflets.append(LipidAggregate(ids, self))

        # Step 3: Inspect potential leaflet to check if they are compatible

        maybe_leaflets.sort(key = len, reverse=True) # Sort leaflets according to their populations
        maybe_leaflets_clean = []
        for maybe_leaflet in maybe_leaflets:
            if maybe_leaflet.size < 50:
                continue

            maybe_leaflets_clean.append(maybe_leaflet)

        # print("{} potential leaflets:".format(len(maybe_leaflets_clean)))
        # for lid, leaflet_ids in enumerate(maybe_leaflets_clean):
        #     print("- leaflet #{}: {} lipids".format(lid, len(leaflet_ids)))

        membranes = []
        while len(maybe_leaflets_clean) > 0:
            ref_leaflet = maybe_leaflets_clean.pop(0)

            compatibles = []
            for leaflet in maybe_leaflets_clean:
                # Todo: check compatibility
                compatibles.append(leaflet)

            companion = None
            for leaflet in compatibles:
                # Todo: check compatibility
                companion = leaflet

            if companion is not None:
                membranes.append(Membrane(ref_leaflet, companion))
                maybe_leaflets_clean.remove(companion)

        self._membranes = membranes
        self._lastupdate_membrane = self.universe.trajectory.frame
        return self._membranes



    def find_aggregates(self, force_update=False) -> [LipidAggregate]:
        # No need to do anything if the membranes are already identified for the current frame
        if self._lastupdate_aggregates == self.universe.trajectory.frame and not force_update:
            return self._aggregates

        # Make sure everything is updated
        self.update()

        current_aggregate_id = -1

        aggregate_ids = np.ones(self._nlipids, dtype=int) * -1

        stack = np.zeros(self._nlipids, dtype=int)
        stack_size = 0

        current_aggregate_lipid_ids = np.zeros(self._nlipids, dtype=int)
        current_aggregate_size = 0


        aggregates = []

        # First step find potential leaflets
        for seed_id in range(self._nlipids):

            if aggregate_ids[seed_id] > -1:  # This lipid already belong to an aggregate, we can skip it
                continue

            # Increment the aggregate id as we are dealing with another aggregate
            current_aggregate_id += 1

            # Add this seed bead to the new aggregate
            aggregate_ids[seed_id] = current_aggregate_id
            current_aggregate_size = 1
            current_aggregate_lipid_ids[0] = seed_id

            # Reset the stack
            stack_size = 1
            stack[0] = seed_id

            # Search for all the lipids that belong to the current aggregate
            # This is done by exhausting the stack that contains all the candidates
            while stack_size > 0:

                # Pop the stack
                ref_nid = stack[stack_size - 1]
                stack_size -= 1

                for nid in self.lipid_neighbours[ref_nid]:

                    if aggregate_ids[nid] > -1:  # This lipid already belong to an aggregate, we can skip it
                        continue

                    # If still here, add bead to current aggregate
                    aggregate_ids[nid] = current_aggregate_id
                    current_aggregate_lipid_ids[current_aggregate_size] = nid
                    current_aggregate_size += 1

                    # Append bead to stack
                    stack[stack_size] = nid
                    stack_size += 1

            # Store the aggregate
            ids = np.sort(current_aggregate_lipid_ids[:current_aggregate_size])
            aggregates.append(ids)

        self._aggregates = []
        aggregates.sort(key = len, reverse=True) # Sort aggregates according to their populations
        for ids in aggregates:
            self._aggregates.append(LipidAggregate(ids, self))

        self._lastupdate_aggregates = self.universe.trajectory.frame
        return self._aggregates

    cdef void compute_weighted_average(self, fsl_int ref_beadid, rvec weighted_position, rvec weighted_normal):
        cdef fsl_int i, j
        cdef rvec dx
        cdef real weight, total_weight, norm


        rvec_clear(weighted_position)
        rvec_clear(weighted_normal)

        total_weight = 1

        for j in range(DIM):
            weighted_normal[j] = self._lipid_normals[ref_beadid, j]
            weighted_position[j] = self._lipid_positions[ref_beadid, j]

        for i in range(self._lipid_neighbours.nneighbours[ref_beadid]):
            beadid = self._lipid_neighbours.neighbours[ref_beadid][i]

            self.box.fast_pbc_dx(
                &self._lipid_positions[ref_beadid, XX],
                &self._lipid_positions[beadid, XX],
                dx
            )

            norm = rvec_norm(dx)


            weight = (norm - self.ns_cutoff)**2 / self.ns_cutoff**2
            total_weight += weight

            for j in range(DIM):
                weighted_normal[j] += weight * self._lipid_normals[beadid, j]
                weighted_position[j] += weight * (self._lipid_positions[ref_beadid, j] + dx[j])

        for j in range(DIM):
            weighted_normal[j] /= total_weight
            weighted_position[j] /= total_weight


    @property
    def positions_bbox(self) -> np.ndarray:
        self.update()
        return np.asarray(self.universe_coords_bbox).copy()

    @property
    def lipid_positions(self) -> np.ndarray:
        self.update()
        return np.asarray(self._lipid_positions).copy()

    @property
    def lipid_directions(self) -> np.ndarray:
        self.update()
        return np.asarray(self._lipid_directions).copy()

    @property
    def lipid_centroids(self) -> np.ndarray:
        self.update()
        return np.asarray(self._lipid_centroids).copy()

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

    @property
    def pbcbox(self) -> PBCBox:
        return self.box

    @property
    def membranes(self) -> [Membrane]:
        self.find_membranes()
        return self._membranes


    @property
    def aggregates(self) -> [LipidAggregate]:
        self.find_aggregates()
        return self._aggregates

    def __eq__(self, other):
        if not isinstance(other, LipidRegistry):
            return False

        return (self.universe == other.universe) and np.all(np.equal(self.hg_indices, other.hg_indices))



cdef class Membrane(object):
    cdef list _leaflets
    cdef readonly LipidRegistry system

    def __init__(self, leaflet1, leaflet2):
        self.system = leaflet1.system

        if leaflet1.system != leaflet2.system:
            raise ValueError("Leaflets do not belong to the same system")

        self._leaflets = list()

        if (leaflet1.is_planar and leaflet1.position[2] > leaflet2.position[2]) or not leaflet1.is_planar:
            self._leaflets.append(Leaflet(leaflet1, self, 0))
            self._leaflets.append(Leaflet(leaflet2, self, 1))
        else:
            self._leaflets.append(Leaflet(leaflet2, self, 0))
            self._leaflets.append(Leaflet(leaflet1, self, 1))

    def __getitem__(self, item):
        if item in (-2, 0):
            return self._leaflets[0]
        elif item in (-1, 1):
            return self._leaflets[1]
        else:
            raise IndexError("A membrane contains two leaflets")


cdef class LipidAggregate:
    cdef readonly LipidRegistry system
    cdef fsl_int[:] _lipid_ids
    cdef fsl_int[:] _is_lipid_id_used
    cdef fsl_int _size

    cdef fsl_int _lastupdate

    cdef real[:] _normal
    cdef real[:] _position
    cdef bint _isplanar

    # cluster-related
    cdef fsl_int[:] _clustered
    cdef real[:, ::1] _positions_clustered_buffer
    cdef real[:, ::1] _positions_clustered
    cdef fsl_int [:, ::1] _cluster_stack


    def __init__(self, fsl_int[:] lipid_ids, LipidRegistry system):
        self.system = system

        self._lipid_ids = np.sort(lipid_ids)
        self._is_lipid_id_used = np.zeros(self.system._nlipids, dtype=int)
        for beadid in self._lipid_ids:
            self._is_lipid_id_used[beadid] = 1

        self._size = self._lipid_ids.shape[0]

        self._lastupdate = -1

        # planarity and average normal and position
        self._isplanar = False
        self._normal = np.empty(DIM, dtype=np.float32)
        self._position = np.empty(DIM, dtype=np.float32)

        # clusterization
        self._clustered = np.empty(self.system._nlipids, dtype=int)
        self._positions_clustered_buffer = np.empty((self.system._nlipids, DIM), dtype=np.float32)
        self._positions_clustered = np.empty((self._size, DIM), dtype=np.float32)
        self._cluster_stack = np.empty((self.system._nlipids, 2), dtype=int)

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.cdivision(True)
    cdef update(self, force_update=False):
        cdef fsl_int i, j, bead_id
        cdef real norm

        # No need to do anything if the membranes are already identified for the current frame
        if self._lastupdate == self.system._lastupdate and not force_update:
            return

        self.system.update(force_update)

        with nogil:
            self._normal[:] = 0

            for i in range(self._size):
                bead_id = self._lipid_ids[i]

                # print("Bead resid {}: position: {}, normal: {}".format(
                #     self.system.lipids[bead_id].resid,
                #     np.asarray(self.system._lipid_positions[bead_id]),
                #     np.asarray(self.system._lipid_normals[bead_id])
                # ))

                for j in range(DIM):
                    self._normal[j] += self.system._lipid_normals[bead_id, j]
                    self._position[j] += self.system._lipid_positions[bead_id, j]


            self.fast_clusterize()


            norm = rvec_norm(&self._normal[XX])

            # print("\nNormal: {}, norm: {}, size: {} - position: {}".format(
            #     np.asarray(self._normal),
            #     norm,
            #     self._size,
            #     np.asarray(self._position)
            # ))

            if norm > 0.75 * self._size:
                self._isplanar = True
            else:
                self._isplanar = False

            for j in range(DIM):
                self._normal[j] /= norm

            self._lastupdate = self.system._lastupdate

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.cdivision(True)
    cdef void fast_clusterize(self, bint force_update=False) nogil:
        cdef fsl_int i, j
        cdef fsl_int beadid, ref_beadid, nid
        cdef fsl_int stack_size = 0, stack_index = 0
        cdef rvec ref_position
        cdef fsl_int counter
        cdef rvec dx, cluster_centroid, raw_dx

        # No need to do anything if self.update was already called
        if self._lastupdate == self.system._lastupdate and not force_update:
            return

        # Reinitialize clutered-related variables
        for i in range(self._size):
            beadid = self._lipid_ids[i]
            self._clustered[beadid] = 0

            for j in range(DIM):
                self._positions_clustered_buffer[beadid, j] = self.system._lipid_positions[beadid, j]


        self._cluster_stack[stack_size][0] = -1
        self._cluster_stack[stack_size][1] = self._lipid_ids[0]
        stack_size = 1
        stack_index = 0

        counter = 0
        while stack_size > stack_index:

            counter += 1

            ref_beadid = self._cluster_stack[stack_index][0]
            beadid  = self._cluster_stack[stack_index][1]
            stack_index += 1

            if self._clustered[beadid] == 2:
                continue

            if ref_beadid < 0:
                for i in range(DIM):
                    ref_position[i] = self._positions_clustered_buffer[beadid, i]
            else:
                for i in range(DIM):
                    ref_position[i] = self._positions_clustered_buffer[ref_beadid, i]

            self.system.box.fast_pbc_dx(ref_position,
                                        &self._positions_clustered_buffer[beadid, XX],
                                        dx)

            for j in range(DIM):
                self._positions_clustered_buffer[beadid, j] = ref_position[j] + dx[j]

            self._clustered[beadid] = 2

            for j in range(self.system._lipid_neighbours.nneighbours[beadid]):
                nid = self.system._lipid_neighbours.neighbours[beadid][j]

                if self._is_lipid_id_used[nid] == 1:
                    if self._clustered[nid] == 0:
                        self._cluster_stack[stack_size][0] = beadid
                        self._cluster_stack[stack_size][1] = nid
                        stack_size += 1

                        self._clustered[nid] = 1

        # Get cluster centroid
        rvec_clear(cluster_centroid)
        for i in range(self._size):
            beadid = self._lipid_ids[i]

            for j in range(DIM):
                cluster_centroid[j] += self._positions_clustered_buffer[beadid, j]
                self._positions_clustered[i][j] = self._positions_clustered_buffer[beadid][j]
        for j in range(DIM):
            cluster_centroid[j] /= self._size

        for i in range(DIM):
            self._position[i] = cluster_centroid[i]


    def __len__(self):
        return self._size

    @property
    def indices(self):
        return np.asarray(self._lipid_ids, dtype=int)

    @property
    def lipids(self):
        lipids = []
        for i in self.indices:
            lipids.append(self.system.lipids[i])
        return lipids

    @property
    def size(self):
        return self._size

    def __getitem__(self, item):
        return self.indices[item]

    @property
    def is_planar(self):
        self.update()
        return self._isplanar

    @property
    def normal(self):
        self.update()
        return np.asarray(self._normal)

    @property
    def centroid(self):
        self.update()
        return np.asarray(self._position)

    @property
    def position(self):
        return self.centroid

    @property
    def positions(self):
        self.update()
        return np.array(self._positions_clustered)

    @property
    def positions_raw(self):
        self.update()
        return np.asarray(self.system.lipid_positions[self._lipid_ids])

    def __str__(self):
        return "Lipid aggregate made of {} lipids".format(self._size)

    def __repr__(self):
        return "<LipidAggregate with {} lipids>".format(self._size)


cdef class Leaflet(LipidAggregate):
    cdef Membrane _membrane
    cdef fsl_int _leaflet_id

    cdef fsl_int _lastupdate_thickness
    cdef real _thickness
    cdef real[:] _lipid_thicknesses
    cdef real[:] _lipid_interleaflet_gaps

    cdef fsl_int _lastupdate_apl
    cdef real _apl
    cdef real[:] _lipid_apls
    cdef real _area

    def __init__(self, ids_or_aggregate, Membrane membrane, fsl_int leaflet_id):
        cdef SimplifiedLipid lipid

        if isinstance(ids_or_aggregate, LipidAggregate):
            ids_or_aggregate = ids_or_aggregate.indices

        system = membrane.system

        super().__init__(ids_or_aggregate, system)

        for i, val in enumerate(self.indices):
            lipid = self.system.lipids[val]

            lipid._regid_leaflet = i
            lipid._leaflet = self

        self._membrane = membrane
        self._leaflet_id = leaflet_id

        # Thickness
        self._lastupdate_thickness = -1
        self._thickness = 0
        self._lipid_thicknesses = np.empty(self._size, dtype=np.float32)
        self._lipid_interleaflet_gaps = np.empty(self._size, dtype=np.float32)

        # APL
        self._lastupdate_apl = -1
        self._apl = 0
        self._lipid_apls = np.empty(self._size, dtype=np.float32)
        self._area = 0


    cdef compute_thickness(self, bint force_update=False):
        cdef fsl_int i, j
        cdef fsl_int ref_beadid, beadid
        cdef Leaflet other_leaflet
        cdef PBCBox box
        cdef rvec dx
        cdef rvec ref_normal, ref_position

        cdef fsl_int locnorm_beadid
        cdef rvec locnorm_position, locnorm_normal, locnorm_dx

        cdef fsl_int closest_beadid
        cdef rvec closest_position, closest_normal, closest_dx

        cdef fsl_int leafnorm_beadid
        cdef rvec leafnorm_position, leafnorm_normal, leafnorm_dx

        cdef fsl_int other_beadid
        cdef rvec other_position, other_normal

        cdef real[:,::] universe_positions


        if self._membrane is None:
            raise ValueError("Can not calculate thickness if leaflet does not belong to a membrane")

        # No need to do anything if the membranes are already identified for the current frame
        if self._lastupdate_thickness == self.system._lastupdate and not force_update:
            return self._thickness, self._lipid_thicknesses


        self.update(force_update)
        box = self.system.box
        universe_positions = self.system.universe_coords_bbox

        if self._leaflet_id == 0:
            other_leaflet = self._membrane._leaflets[1]
        else:
            other_leaflet = self._membrane._leaflets[0]

        for i in range(self._size):
            ref_beadid = self._lipid_ids[i]

            # Get average normal and position for the reference bead
            self.system.compute_weighted_average(ref_beadid, ref_position, ref_normal)

            # print("Weighted position for resid {}: [{:.3f}, {:.3f}, {:.3f}], weighted normal: [{:.3f}, {:.3f}, {:.3f}] - actual position: {}".format(
            #     self.system.lipids[ref_beadid].resid,
            #     ref_position[XX], ref_position[YY], ref_position[ZZ],
            #     ref_normal[XX], ref_normal[YY], ref_normal[ZZ],
            #     np.asarray(self.system._lipid_positions[ref_beadid])
            # ))

            # Get the bead from the other leaflet which can be used to compute thickness
            # The candidate is selected by checking the dot product of dx and reference normal
            beadid_best = -1
            dprod_best = 0

            for j in range(other_leaflet._size):
                beadid = other_leaflet._lipid_ids[j]


                box.fast_pbc_dx_leaflet(
                    &self.system._lipid_positions[ref_beadid, XX],
                    &self.system._lipid_positions[beadid, XX],
                    dx,
                    &ref_normal[XX])

                dprod_val = rvec_dprod(&ref_normal[XX],
                                       dx)

                norm = rvec_norm(dx)

                dprod_val2 = rvec_dprod(&ref_normal[XX],
                                       &self.system._lipid_normals[beadid, XX])

                if dprod_val/norm < dprod_best and dprod_val2 < -0.707:
                    dprod_best = dprod_val/norm
                    beadid_best = beadid

            box.fast_pbc_dx_leaflet(
                    &self.system._lipid_positions[ref_beadid, XX],
                    &self.system._lipid_positions[beadid_best, XX],
                    dx,
                    &ref_normal[XX])

            # print("\nRef bead resid {}, Best bead resid {}: dx:[{}, {}, {}], ref normal:{}, dprod:{}".format(
            #     self.system.lipids[ref_beadid].resid,
            #     self.system.lipids[beadid_best].resid,
            #     dx[XX], dx[YY], dx[ZZ],
            #     np.asarray(ref_normal),
            #     dprod_best
            # ))

            box.fast_pbc_dx_leaflet(
                    &self.system._lipid_positions[ref_beadid, XX],
                    &self.system._lipid_positions[beadid_best, XX],
                    locnorm_dx,
                    &ref_normal[XX])

            self.system.compute_weighted_average(beadid_best, locnorm_position, locnorm_normal)
            locnorm_beadid = beadid_best


            dprod_val = rvec_dprod(&self._normal[XX],
                                   &ref_normal[XX])

            rvec_copy(locnorm_normal, other_normal)
            rvec_copy(locnorm_position, other_position)
            other_beadid = beadid_best

            if self._isplanar and acos(dprod_val) > 10 / 180 * np.pi:
                    # print("lipid resid {} normal is off (by {:.1f}°): {} vs avg normal: {}".format(
                    #     self.system.lipids[ref_beadid].resid,
                    #     acos(dprod_val) / np.pi * 180,
                    #     np.asarray(ref_normal),
                    #     np.asarray(self._normal)
                    # ))


                    # Get closest bead
                    d_min = 100000
                    beadid_best = -1

                    for j in range(other_leaflet._size):
                        beadid = other_leaflet._lipid_ids[j]


                        box.fast_pbc_dx_leaflet(
                            &self.system._lipid_positions[ref_beadid, XX],
                            &self.system._lipid_positions[beadid, XX],
                            dx,
                            &ref_normal[XX])

                        norm = rvec_norm(dx)

                        if norm < d_min:
                            beadid_best = beadid
                            d_min = norm

                    box.fast_pbc_dx_leaflet(
                        &self.system._lipid_positions[ref_beadid, XX],
                        &self.system._lipid_positions[beadid_best, XX],
                        closest_dx,
                        &ref_normal[XX])

                    self.system.compute_weighted_average(beadid_best, closest_position, closest_normal)
                    closest_beadid = beadid_best


                    # Get best using leaflet normal
                    beadid_best = -1
                    dprod_best = 0

                    for j in range(other_leaflet._size):
                        beadid = other_leaflet._lipid_ids[j]


                        box.fast_pbc_dx_leaflet(
                            &self.system._lipid_positions[ref_beadid, XX],
                            &self.system._lipid_positions[beadid, XX],
                            dx,
                            &self._normal[XX])

                        dprod_val = rvec_dprod(&self._normal[XX],
                                               dx)

                        norm = rvec_norm(dx)

                        if dprod_val/norm < dprod_best:
                            dprod_best = dprod_val/norm
                            beadid_best = beadid

                    box.fast_pbc_dx_leaflet(
                    &self.system._lipid_positions[ref_beadid, XX],
                    &self.system._lipid_positions[beadid_best, XX],
                    leafnorm_dx,
                    &ref_normal[XX])


                    self.system.compute_weighted_average(beadid_best, leafnorm_position, leafnorm_normal)
                    leafnorm_beadid = beadid_best

                    delta = abs(rvec_norm(closest_dx) - rvec_norm(locnorm_dx))/rvec_norm(closest_dx)


                    # print("-> Candidates for resid {} (p:[{:.1f},{:.1f},{:.1f}], n: [{:.1f},{:.1f},{:.1f}]):\n  ->locnorm ({}): p:[{:.1f},{:.1f},{:.1f}], n:[{:.1f},{:.1f},{:.1f}], d:{:.1f}\n  ->"
                    #       "leafnorm ({}): p:[{:.1f},{:.1f},{:.1f}], n:[{:.1f},{:.1f},{:.1f}], d:{:.1f}\n  ->"
                    #       "closest ({}): p:[{:.1f},{:.1f},{:.1f}], n:[{:.1f},{:.1f},{:.1f}], d:{:.1f} => delta: {:.3f}\n".format(
                    #     self.system.lipids[ref_beadid].resid,
                    #     ref_position[XX], ref_position[YY], ref_position[ZZ],
                    #     ref_normal[XX], ref_normal[YY], ref_normal[ZZ],
                    #     self.system.lipids[locnorm_beadid].resid,
                    #     locnorm_position[XX], locnorm_position[YY], locnorm_position[ZZ],
                    #     locnorm_normal[XX], locnorm_normal[YY], locnorm_normal[ZZ],
                    #     rvec_norm(locnorm_dx),
                    #     self.system.lipids[leafnorm_beadid].resid,
                    #     leafnorm_position[XX], leafnorm_position[YY], leafnorm_position[ZZ],
                    #     leafnorm_normal[XX], leafnorm_normal[YY], leafnorm_normal[ZZ],
                    #     rvec_norm(leafnorm_dx),
                    #     self.system.lipids[closest_beadid].resid,
                    #     closest_position[XX], closest_position[YY], closest_position[ZZ],
                    #     closest_normal[XX], closest_normal[YY], closest_normal[ZZ],
                    #     rvec_norm(closest_dx),
                    #     delta
                    # ))

                    if delta > 0.2:
                        if delta > 5:
                            self._lipid_thicknesses[i] = np.nan
                            continue
                        else:
                            rvec_copy(leafnorm_normal, other_normal)
                            rvec_copy(leafnorm_position, other_position)
                            rvec_copy(&self._normal[XX], ref_normal)
                            other_beadid = leafnorm_beadid

            box.fast_pbc_dx_leaflet(
            ref_position,
            other_position,
            dx,
            ref_normal)

            norm = rvec_norm(dx)

            dprod_val = abs(rvec_dprod(dx, ref_normal))
            other_dprod_val = abs(rvec_dprod(dx, other_normal))


            self._lipid_thicknesses[i] = 0.5 * (dprod_val + other_dprod_val)

            # print("Thickness for resid {}: ref {:.3f}, dx norm: {:.3f}, other: {:.3f}, avg:{:.3f}".format(
            #     self.system.lipids[ref_beadid].resid,
            #     dprod_val,
            #     norm,
            #     abs(rvec_dprod(dx, other_normal)),
            #     0.5 * (dprod_val + abs(rvec_dprod(dx, other_normal)))
            # ))

            # Compute interleaflet gap
            ref_d = 0

            for ix in self.system.lipids[ref_beadid]._ix:
                box.fast_pbc_dx(
                    &self.system._lipid_positions[ref_beadid, XX],
                    &universe_positions[ix, XX],
                    dx)

                dprod_val = rvec_dprod(dx, ref_normal)

                # print("d ref atom-atom {}: {:.3f} (dz={:.3f})".format(
                #     self.system.universe.atoms[ix].name,
                #     dprod_val,
                #     self.system._lipid_positions[ref_beadid, ZZ] - self.system.universe.atoms[ix].position[ZZ]
                # ))

                if dprod_val < ref_d:
                    ref_d = dprod_val
            # print("")

            other_d = 0
            for ix in self.system.lipids[other_beadid]._ix:
                box.fast_pbc_dx_leaflet(
                    &self.system._lipid_positions[ref_beadid, XX],
                    &universe_positions[ix, XX],
                    dx,
                ref_normal)

                dprod_val = rvec_dprod(dx, ref_normal)

                if dprod_val < other_d:
                    other_d = dprod_val

            self._lipid_interleaflet_gaps[i] = 0.5 * self._lipid_thicknesses[i] - np.abs(other_d)


        self._thickness = np.nanmean(self._lipid_thicknesses)
        self._lastupdate_thickness = self.system._lastupdate

        return self._thickness, self._lipid_thicknesses

    @property
    def thickness(self):
        self.compute_thickness()
        return self._thickness

    @property
    def lipid_thicknesses(self):
        self.compute_thickness()
        return np.asarray(self._lipid_thicknesses).copy()

    @property
    def lipid_interleaflet_gaps(self):
        self.compute_thickness()
        return np.asarray(self._lipid_interleaflet_gaps).copy()


    cdef compute_apl(self, bint force_update=False):
        cdef fsl_int i, j
        cdef fsl_int ref_beadid, beadid
        cdef rvec plane_x, plane_y, plane_z
        cdef real[:] ref_normal, ref_position
        cdef real[:] current_position
        cdef matrix conv_mat_revert, conv_mat
        cdef rvec dx, proj

        cdef Polygon *cell
        cdef real_point tmp_pt
        cdef Polygon *buffer = polygon_new()
        cdef real_point ref_point_2d
        cdef real_point current_position_2d

        if self._membrane is None:
            raise ValueError("Can not calculate APL if leaflet does not belong to a membrane")

        # No need to do anything if the membranes are already identified for the current frame
        if self._lastupdate_apl == self.system._lastupdate and not force_update:
            return self._apl, self._lipid_apls


        self.update(force_update)
        box = self.system.box

        cell = polygon_new()

        ref_point_2d[XX] = 0
        ref_point_2d[YY] = 0

        for i in range(self._size):
            ref_beadid = self._lipid_ids[i]

            ref_normal = self.system._lipid_normals[ref_beadid]
            ref_position = self.system._lipid_positions[ref_beadid]




            # Compute projection matrix
            # Build plane basis
            rvec_copy(&ref_normal[XX], plane_z)
            complete_basis(plane_z, plane_x, plane_y)

            # Get conversion matrix
            mat_from_rvec(plane_x, plane_y, plane_z, conv_mat_revert)
            invert_mat(conv_mat_revert, conv_mat)

            # print("Resid {} -> Convertion matrix:\n{:.3f}, {:.3f}, {:.3f}\n{:.3f}, {:.3f}, {:.3f}\n{:.3f}, {:.3f}, {:.3f}\n".format(
            #     self.system.lipids[ref_beadid].resid,
            #     conv_mat_revert[XX][XX], conv_mat_revert[XX][YY], conv_mat_revert[XX][ZZ],
            #     conv_mat_revert[YY][XX], conv_mat_revert[YY][YY], conv_mat_revert[YY][ZZ],
            #     conv_mat_revert[ZZ][XX], conv_mat_revert[ZZ][YY], conv_mat_revert[ZZ][ZZ],
            # ))

            # (Re)initialize cell
            polygon_empty(cell)

            # Start with a big square
            tmp_pt[XX] = -2 * self.system.ns_cutoff
            tmp_pt[YY] = -2 * self.system.ns_cutoff
            polygon_append(cell, tmp_pt)
            tmp_pt[XX] = 2 * self.system.ns_cutoff
            polygon_append(cell, tmp_pt)
            tmp_pt[YY] = 2 * self.system.ns_cutoff
            polygon_append(cell, tmp_pt)
            tmp_pt[XX] = -2 * self.system.ns_cutoff
            polygon_append(cell, tmp_pt)


            for j in range(self.system._lipid_neighbours.nneighbours[ref_beadid]):
                beadid = self.system._lipid_neighbours.neighbours[ref_beadid][j]

                current_position = self.system._lipid_positions[beadid]

                # Get the closest image of the neighbor
                box.fast_pbc_dx(&ref_position[XX], &current_position[XX], dx)

                # Get the coordinates into the new basis (plane_x, plane_y, plane_z)
                rvec_to_basis(dx, conv_mat, proj)

                # Store the projection onto the (plane_x, plane_y) plane
                current_position_2d[XX] = proj[XX]
                current_position_2d[YY] = proj[YY]

                # print("Resid {} neighbour {}: [{:.3f}, {:.3f}]".format(
                #     self.system.lipids[ref_beadid].resid,
                #     self.system.lipids[beadid].resid,
                #     current_position_2d[XX], current_position_2d[YY]
                # ))

                fast_clip_zoi(cell, ref_point_2d, current_position_2d, buffer)

            self._lipid_apls[i] = polygon_get_area(cell)


        polygon_destroy(cell)
        polygon_destroy(buffer)

        self._lastupdate_apl = self.system._lastupdate
        self._apl = np.mean(self._lipid_apls)
        self._area = np.sum(self._lipid_apls)

        return self._apl, self._lipid_apls, self._area

    @property
    def apl(self):
        self.compute_apl()
        return self._apl

    @property
    def lipid_apls(self):
        self.compute_apl()
        return np.asarray(self._lipid_apls).copy()

    @property
    def area(self):
        self.compute_apl()
        return self._area



