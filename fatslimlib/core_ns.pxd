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

DEF DIM=3

from .typedefs cimport real, fsl_int

from .core_base cimport PBCBox

# Neighbor search related structures and functions
cdef struct ns_grid:
    fsl_int size
    fsl_int[DIM] ncells
    real[DIM] cellsize
    fsl_int *nbeads
    fsl_int **beadids

cdef struct ns_neighborhood:
    real cutoff
    fsl_int allocated_size
    fsl_int size
    fsl_int *beadids

cdef struct ns_neighborhood_holder:
    fsl_int size
    ns_neighborhood **neighborhoods

cdef ns_neighborhood_holder *fast_neighbor_search(real[:, ::1] ref_coords,
                                                  real[:, ::1] neighbor_coords,
                                                  PBCBox box,
                                                  real cutoff) nogil

cdef void free_neighborhood_holder(ns_neighborhood_holder *holder) nogil
