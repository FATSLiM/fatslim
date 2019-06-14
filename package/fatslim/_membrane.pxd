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


from ._core cimport LipidRegistry

from ._aggregate cimport LipidAggregate


cdef class Membrane(object):
    # Cdef attributes
    cdef list _leaflets
    cdef readonly LipidRegistry system


cdef class Leaflet(LipidAggregate):
    # Cdef attributes
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

    # Cdef methods
    cdef compute_thickness(self, bint force_update=*)
    cdef compute_apl(self, bint force_update=*)
