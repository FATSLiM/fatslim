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

DEF NOTSET = -12345
DEF DIM = 3
DEF XX = 0
DEF YY = 1
DEF ZZ = 2
DEF EPSILON = 1e-6
DEF PI = 3.141592653589793

from ._typedefs cimport real, rvec, dreal, fsl_int, ivec

from ._geometry cimport PBCBox

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

    # Cdef methods
    cdef int update(self) nogil except -1
    cdef fsl_int coord2cellid(self, rvec coord) nogil
    cdef bint cellid2cellxyz(self, fsl_int cellid, ivec cellxyz) nogil
    cdef int fill_grid(self, real[:, ::1] coords) nogil except -1
    cdef void resize_beadlist(self)


cdef class _NSResults:
    cdef fsl_int size
    cdef fsl_int max_nneighbours
    cdef fsl_int[:] nneighbours
    cdef fsl_int[:, ::1] neighbours
    cdef real[:, ::1] distances

    # Cdef methods
    cdef int add_neighbour(self, fsl_int i, fsl_int j, real d2) nogil except -1
    cdef void reset(self) nogil
    cdef resize_neighbourlist(self)


cdef int fast_self_search(_NSGrid grid, _NSResults results, real[:, ::1] positions) nogil except -1


cdef class SimplifiedLipid:
    # Cdef attributes
    cdef fsl_int[:] _ix
    cdef fsl_int _regid
    cdef LipidRegistry _registry


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

    # Cdef methods
    cdef void compute_weighted_average(self, fsl_int ref_beadid, rvec weighted_position, rvec weighted_normal)