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

from .typedefs cimport real, rvec, fsl_int

from .core_base cimport Frame
from .core_ns cimport ns_neighborhood_holder

# Classes
cdef class Aggregate:
    # Attributes
    cdef Frame frame
    cdef fsl_int[:] beadids
    cdef fsl_int[:] resids
    cdef fsl_int[:] hg_atomids
    cdef fsl_int[:] lipid_atomids
    cdef real[:, ::1] coords
    cdef real[:, ::1] directions
    cdef real[:, ::1] normals
    cdef real[:, ::1] neighborhood_normals
    cdef real neighborhood_cutoff
    cdef ns_neighborhood_holder *neighborhoods
    cdef rvec xcm
    cdef rvec avg_normal
    cdef bint is_planar
    cdef list lipid_types
    cdef list beadids_by_type
    cdef dict nlipids_by_type

    # Utility
    cdef fsl_int fast_get_atomid(self, fsl_int internalid) nogil
    cpdef fsl_int get_atomid(self, fsl_int internalid)
    cdef bint fast_same_restype(self, fsl_int internalid1, fsl_int internalid2) nogil
    cpdef bint same_restype(self, fsl_int internalid1, fsl_int internalid2)
    cpdef str get_resname(self, fsl_int internalid)
    cpdef fsl_int get_resid(self, fsl_int internalid)

    # APL related
    cdef real apl_cutoff
    cdef real[:] apl_values
    cdef readonly real apl_min
    cdef readonly real apl_max
    cdef readonly real area
    cdef list apl_by_types
    cdef real[:] area_by_types
    cdef readonly real apl_avg

    # Thickness related
    cdef real thickness_cutoff
    cdef real[:] thickness_values
    cdef readonly real thickness_avg
    cdef readonly real thickness_min
    cdef readonly real thickness_max

    # Analysis methods
    cdef void refresh_cache(self) nogil
    cdef fsl_int fast_size(self) nogil
    cdef void compute_neighborhood_normal(self, fsl_int index) nogil
    cdef void set_neighborhood_cutoff(self, real cutoff, bint force=*) nogil
    cdef void compute_apl(self, real cutoff, real area_limit=*, bint force=*) nogil except *
    cdef real fix_apl(self, fsl_int index, real *tmp_apl, bint onlyfix=*) nogil
    cdef void compute_thickness(self, Aggregate other, real interleaflet_cutoff=*, bint force=*) nogil except*
    cdef real fix_thickness(self, fsl_int index, real *tmp_thickness, bint onlyfix=*) nogil


cdef class Membrane:
    # Attributes
    cdef Frame frame
    cdef readonly Aggregate leaflet1
    cdef readonly Aggregate leaflet2
    cdef rvec xcm
    cdef fsl_int type

    # Thickness related
    cdef real thickness

    # APL related
    cdef real apl

    # Analysis methods
    cdef void fast_compute_thickness(self, real interleaflet_cutoff=*, bint force=*) nogil except*
    cdef void fast_compute_apl(self, real cutoff=*, real area_limit=*, bint force=*) nogil except*


# Analysis stuff
cdef list retrieve_aggregates(Frame frame, real cutoff)
cdef list retrieve_membranes(Frame frame, real cutoff)
