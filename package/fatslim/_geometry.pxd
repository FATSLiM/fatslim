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

from ._typedefs cimport matrix, rvec, dreal, real, fsl_int

ctypedef real[2] real_point

cdef struct cPBCBox_t:
    matrix     box
    rvec       fbox_diag
    rvec       hbox_diag
    rvec       mhbox_diag
    dreal      max_cutoff2

cdef class PBCBox(object):

    # Cdef attributes
    cdef cPBCBox_t c_pbcbox
    cdef bint is_triclinic

    cdef void fast_update(self, real[:, ::1] box) nogil
    cdef void fast_pbc_dx(self, rvec ref, rvec other, rvec dx) nogil
    cdef void fast_pbc_dx_leaflet(self, rvec ref, rvec other, rvec dx, rvec ref_normal) nogil
    cdef dreal fast_distance2(self, rvec a, rvec b) nogil
    cdef void fast_put_atoms_in_bbox(self, real[:, ::1] coords, real[:, ::1] bbox_coords) nogil
    cdef void fast_pbc_centroid(self, real[:, ::1] coords, rvec xcm, fsl_int[:] indices) nogil
    cdef void fast_pbc_centroid_from_ref(self, real[:, ::1] coords, rvec ref, rvec xcm, fsl_int[:] indices) nogil


cdef bint normal_from_neighbours(real[:, ::1] positions,
                                 real[:, ::1] directions,
                                 fsl_int refid,
                                 fsl_int[:] neighbours_ids,
                                 PBCBox box,
                                 rvec normal) nogil

cdef bint curvature_from_neighbours(real[:, ::1] positions,
                                 real[:, ::1] directions,
                                 fsl_int refid,
                                 fsl_int[:] neighbours_ids,
                                 PBCBox box,
                                 rvec eig_vals,
                                 matrix eig_vecs) nogil

cdef void eigen_33_sym(matrix a, matrix eig_vec, rvec eig_val) nogil
cdef void complete_basis(rvec z_axis, rvec x_axis, rvec y_axis) nogil
cdef void compute_roots(matrix mat, rvec roots) nogil
cdef int compute_rank(matrix m) nogil


# Polygon pseudo class
cdef struct Polygon:
    real_point *points
    int size
    int allocated_size
cdef Polygon *polygon_new() nogil
cdef Polygon *polygon_new_from_memview(real[:, ::1] coords) nogil
cdef Polygon *polygon_new_from_polygon( Polygon *other) nogil
cdef void polygon_copy(Polygon *self, Polygon *other) nogil
cdef void polygon_get_cog(Polygon *self, real_point cog) nogil
cdef void polygon_destroy(Polygon *self) nogil
cdef void polygon_empty(Polygon *self) nogil
cdef real polygon_get_area(Polygon *self) nogil
cdef real polygon_get_perimeter(Polygon *self) nogil
cdef bint polygon_is_inside(Polygon *self, real_point point) nogil
cdef object polygon_as_array(Polygon *self)
cdef void polygon_append(Polygon *self, real_point point) nogil

cdef Polygon *fast_get_clipped_polygon(Polygon *subject, Polygon *clipper) nogil

cdef Polygon *fast_get_zoi(real_point ref_pt, real_point *pts, int size) nogil
cdef void fast_clip_zoi(Polygon *zoi, real_point ref_pt, real_point clipping_pt, Polygon *buffer) nogil

cdef void rvec_to_basis(rvec v, matrix conv_mat, rvec result) nogil
