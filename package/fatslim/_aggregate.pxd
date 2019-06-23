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

from ._typedefs cimport real, fsl_int

from ._core cimport LipidRegistry, _NSResults


cdef class LipidAggregate:
    cdef readonly LipidRegistry system
    cdef fsl_int[:] _lipid_ids
    cdef fsl_int[:] _system2aggregate_ids
    cdef fsl_int[:] _is_lipid_id_used
    cdef fsl_int _size

    cdef fsl_int _lastupdate

    cdef real[:] _normal
    cdef real[:] _position
    cdef bint _isplanar

    cdef real[:, ::1] _lipid_positions
    cdef real[:, ::1] _lipid_directions
    cdef real[:, ::1] _lipid_normals

    cdef _NSResults _lipid_neighbours

    # cluster-related
    cdef fsl_int[:] _clustered
    cdef real[:, ::1] _positions_clustered_buffer
    cdef real[:, ::1] _lipid_positions_clustered
    cdef fsl_int [:, ::1] _cluster_stack

    # Cdef methods
    cdef update(self, force_update=*)
    cdef void fast_clusterize(self, bint force_update=*) nogil