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

from .typedefs cimport real

ctypedef real[2] real_point

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

