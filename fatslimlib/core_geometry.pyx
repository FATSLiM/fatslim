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
#cython: cdivision=True
#cython: boundscheck=False
#cython: profile=False

from libc.stdlib cimport malloc, realloc, free, abort
from libc.math cimport sqrt, fabs
#cimport numpy as np
import numpy as np

DEF POLYGON_ALLOCATION_INCREMENT = 10
DEF ZOI_DEFAULT_SIZE = 1e6
DEF XX = 0
DEF YY = 1
DEF EPSILON = 1e-6

# Useful functions
cdef bint fast_inside_clipped_polygon(real_point p, real_point centerclip, real_point cp1, real_point cp2) nogil:
    cdef real winding_ok, winding_test
    winding_ok = (cp2[XX]-cp1[XX])*(centerclip[YY]-cp1[YY]) - (cp2[YY]-cp1[YY])*(centerclip[XX]-cp1[XX])
    winding_test = (cp2[XX]-cp1[XX])*(p[YY]-cp1[YY]) - (cp2[YY]-cp1[YY])*(p[XX]-cp1[XX])

    if fabs(winding_test) < EPSILON:
        return False

    return winding_ok / winding_test > 0


cdef void fast_compute_intersection(real_point s, real_point e, real_point cp1, real_point cp2, real_point intersection) nogil:
    cdef real n1, n2, n3
    cdef real_point dc, dp
    dc[XX] = cp1[0] - cp2[0]
    dc[YY] = cp1[1] - cp2[1]
    dp[XX] = s[0] - e[0]
    dp[YY] = s[1] - e[1]

    n1 = cp1[0] * cp2[1] - cp1[1] * cp2[0]
    n2 = s[0] * e[1] - s[1] * e[0]
    n3 = (dc[0] * dp[1] - dc[1] * dp[0])

    if fabs(n3) < EPSILON:
        intersection[XX] = s[XX]
        intersection[YY] = s[YY]
    else:
        n3 = 1.0 / n3
        intersection[XX] = (n1 * dp[0] - n2 * dc[0]) * n3
        intersection[YY] = (n1 * dp[1] - n2 * dc[1]) * n3


# See: http://geomalgorithms.com/a01-_area.html
cdef inline int point_is_left(real_point l0, real_point l1, real_point p) nogil:
    cdef real is_left = (l1[XX] - l0[XX]) * (p[YY] - l0[YY]) - (p[XX] -  l0[XX]) * (l1[YY] - l0[YY])

    if fabs(is_left) < EPSILON:
        return 0 # Point p is on the line

    if is_left > 0:
        return 1 # Point p is on the left

    return -1 # Point p is on the right
########################################################################################################################
# Polygon pseudo class methods
########################################################################################################################
cdef Polygon *polygon_new() nogil:
    cdef Polygon *self = NULL
    self = <Polygon *> malloc(sizeof(Polygon))
    if self == NULL:
        abort()
    self.size = 0
    self.allocated_size = 0
    self.points = NULL
    return self

cdef Polygon *polygon_new_from_memview(real[:, ::1] coords) nogil:
    cdef Polygon *self = polygon_new()
    cdef int size = coords.shape[0]
    cdef int i

    self.points = <real_point *> malloc(size * sizeof(real_point))
    if self.points == NULL:
        abort()

    self.allocated_size = size

    for i in range(size):
        self.points[i][XX] = coords[i, XX]
        self.points[i][YY] = coords[i, YY]
    self.size = size

    return self

cdef Polygon *polygon_new_from_polygon( Polygon *other) nogil:
    cdef Polygon *self = polygon_new()

    polygon_copy(self, other)

    return self

cdef void polygon_copy(Polygon *self, Polygon *other) nogil:
    cdef int size = other.size
    cdef int i

    if self.allocated_size < size:
        if self.allocated_size == 0:
            self.points = <real_point *> malloc(size * sizeof(real_point))
        else:
            self.points = <real_point *> realloc(<void **> self.points, size * sizeof(real_point))

        self.allocated_size = size

    if self.points == NULL:
        abort()

    for i in range(size):
        self.points[i][XX] = other.points[i][XX]
        self.points[i][YY] = other.points[i][YY]
    self.size = size

cdef void polygon_get_cog(Polygon *self, real_point cog) nogil:
    cdef int i

    cog[XX] = 0
    cog[YY] = 0

    for i in range(self.size):
        cog[XX] += self.points[i][XX]
        cog[YY] += self.points[i][YY]

    if self.size == 0:
        return

    cog[XX] /= self.size
    cog[YY] /= self.size


cdef void polygon_destroy(Polygon *self) nogil:
    if self == NULL:
        return
    free(self.points)
    free(self)


cdef void polygon_empty(Polygon *self) nogil:
    self.size = 0


cdef void polygon_append(Polygon *self, real_point point) nogil:
    if self.size + 1 > self.allocated_size:
        self.allocated_size += POLYGON_ALLOCATION_INCREMENT
        if self.points == NULL:
            self.points = <real_point *> malloc(sizeof(real_point) * self.allocated_size)
        else:
            self.points = <real_point *> realloc(<void **> self.points, sizeof(real_point) * self.allocated_size)

        if self.points == NULL:
            abort()

    self.points[self.size][XX] = point[XX]
    self.points[self.size][YY] = point[YY]

    self.size += 1


cdef real polygon_get_area(Polygon *self) nogil:
    cdef int i, previous_i
    cdef real area = 0

    if self == NULL:
        return -1.0

    if self.size < 3:
        return -1.0

    previous_i = self.size - 1
    for i in range(0, self.size):
        area += self.points[previous_i][XX] * self.points[i][YY] - self.points[i][XX] * self.points[previous_i][YY]
        previous_i = i

    if area < 0:
        return area * -0.5
    else:
        return area * 0.5

cdef real polygon_get_perimeter(Polygon *self) nogil:
    cdef int i, previous_i
    cdef real perimeter = 0, delta_x, delta_y

    if self.size < 3:
        return 0.0

    previous_i = self.size - 1
    for i in range(0, self.size):
        delta_x = self.points[i][XX] - self.points[previous_i][XX]
        delta_y = self.points[i][YY] - self.points[previous_i][YY]
        perimeter += sqrt(delta_x * delta_x + delta_y * delta_y)
        previous_i = i

    return perimeter


cdef object polygon_as_array(Polygon *self):
    cdef int i
    cdef real[:, ::1] py_self = np.empty([self.size, 2])

    for i in range(self.size):
        py_self[i, XX] = self.points[i][XX]
        py_self[i, YY] = self.points[i][YY]

    return np.asarray(py_self)



cdef Polygon *fast_get_clipped_polygon(Polygon *subject, Polygon *clipper) nogil:
    cdef real_point clipper_cog, cp1, cp2, s, e, intersect
    cdef int i, j
    cdef Polygon *clipped = NULL
    cdef Polygon *toclip = polygon_new()

    if subject.size < 3 or clipper.size < 3:
        return toclip

    clipped = polygon_new_from_polygon(subject)

    # Compute the center of geometry of the clipping polygon
    polygon_get_cog(clipper, clipper_cog)

    cp1[XX] = clipper.points[clipper.size -1][XX]
    cp1[YY] = clipper.points[clipper.size -1][YY]
    for i in range(clipper.size):
        cp2[XX] = clipper.points[i][XX]
        cp2[YY] = clipper.points[i][YY]

        # Copy clipped polygon to polygon to be clipped
        polygon_copy(toclip, clipped)

        # Check if we need to break
        if toclip.size == 0:
            break

        s[XX] = toclip.points[toclip.size - 1][XX]
        s[YY] = toclip.points[toclip.size - 1][YY]

        # Empty clipped polygon
        polygon_empty(clipped)

        for j in range(toclip.size):
            e[XX] = toclip.points[j][XX]
            e[YY] = toclip.points[j][YY]

            if fast_inside_clipped_polygon(e, clipper_cog, cp1, cp2):
                if not fast_inside_clipped_polygon(s, clipper_cog, cp1, cp2):
                    # Append intersection to clipped polygon
                    fast_compute_intersection(s, e, cp1, cp2, intersect)
                    polygon_append(clipped, intersect)
                # Append current point to clipped polygon
                polygon_append(clipped, e)
            else:
                if fast_inside_clipped_polygon(s, clipper_cog, cp1, cp2):
                    # Append intersection to clipped polygon
                    fast_compute_intersection(s, e, cp1, cp2, intersect)
                    polygon_append(clipped, intersect)

            s[XX] = e[XX]
            s[YY] = e[YY]

        cp1[XX] = cp2[XX]
        cp1[YY] = cp2[YY]

    # Free memory
    polygon_destroy(toclip)

    return clipped

# See: http://geomalgorithms.com/a03-_inclusion.html#wn_PnPoly%28%29
cdef bint polygon_is_inside(Polygon *self, real_point p) nogil:
    cdef int i
    cdef int wn
    cdef real_point next_vertex, current_vertex

    current_vertex[XX] = self.points[self.size - 1][XX] # Last vertex
    current_vertex[YY] = self.points[self.size - 1][YY] # Last vertex

    wn = 0
    for i in range(self.size):
        next_vertex[XX] = self.points[i][XX]
        next_vertex[YY] = self.points[i][YY]

        if current_vertex[YY] <= p[YY]:
            if next_vertex[YY] > p[YY]: # upward crossing
                if point_is_left(current_vertex, next_vertex, p) > 0: # p is on the left side of the edge
                    wn += 1
        else:
            if next_vertex[YY] <= p[YY]: # downward crossing
                if point_is_left(current_vertex, next_vertex, p) < 0: # p is on the right side of the edge
                    wn -= 1

        current_vertex[XX] = next_vertex[XX]
        current_vertex[YY] = next_vertex[YY]

    return wn != 0


########################################################################################################################
#
# Python API
#
########################################################################################################################

def get_polygon_area(real[:,::1] points not None):
    cdef Polygon *polygon = NULL
    cdef real area

    with nogil:
        polygon = polygon_new_from_memview(points)
        area = polygon_get_area(polygon)

        # Free memory
        polygon_destroy(polygon)

    return area

def get_polygon_perimeter(real[:,::1] points not None):
    cdef Polygon *polygon = NULL
    cdef real perimeter

    with nogil:
        polygon = polygon_new_from_memview(points)
        perimeter = polygon_get_perimeter(polygon)

        # Free memory
        polygon_destroy(polygon)

    return perimeter

def get_clipped_polygon(real[:, ::1] subject_points not None, real[:, ::1] clipper_points not None):
    cdef Polygon *subject = polygon_new_from_memview(subject_points)
    cdef Polygon *clipper = polygon_new_from_memview(clipper_points)
    cdef Polygon *clipped = fast_get_clipped_polygon(subject, clipper)

    clipped_py = polygon_as_array(clipped)

    # Free memory
    polygon_destroy(clipper)
    polygon_destroy(subject)
    polygon_destroy(clipped)

    return clipped_py


def is_inside_polygon(real[:, ::1] points, real[:] point):
    cdef Polygon *polygon = polygon_new_from_memview(points)
    cdef real_point p
    cdef bint is_inside
    p[XX] = point[XX]
    p[YY] = point[YY]

    is_inside = polygon_is_inside(polygon, p)

    polygon_destroy(polygon)

    return is_inside

cdef bint same_side(real_point a, real_point b, real_point x_line, real_point y_line) nogil:
    return ((x_line[YY] - y_line[YY])*(a[XX]-x_line[XX]) + (y_line[XX] - x_line[XX])*(a[YY] - x_line[YY])) * \
    ((x_line[YY] - y_line[YY])*(b[XX]-x_line[XX]) + (y_line[XX] - x_line[XX])*(b[YY] - x_line[YY])) >= 0

cdef bint are_parallel(real_point a1, real_point a2, real_point b1, real_point b2) nogil:
    cdef real det = (a1[XX] - a2[XX]) * (b1[YY] - b2[YY]) - (a1[YY] - a2[YY])* (b1[XX] - b2[XX])
    return fabs(det) < EPSILON

cdef void get_intersection(real_point a1, real_point a2, real_point b1, real_point b2,
                           real_point *inter_sect) nogil:
    cdef real det = (a1[XX] - a2[XX]) * (b1[YY] - b2[YY]) - (a1[YY] - a2[YY])* (b1[XX] - b2[XX])
    cdef real val1, val2, x, y

    if fabs(det) < EPSILON: # Lines are coincident
        inter_sect[0][XX] = a1[XX]
        inter_sect[0][YY] = a1[YY]
    else:
        val1 = (a1[0] * a2[1] - a1[1] * a2[0])
        val2 = (b1[0] * b2[1] - b1[1] * b2[0])

        x = val1 * (b1[0] - b2[0]) - (a1[0] - a2[0]) * val2
        x /= det

        y = val1 * (b1[1] - b2[1]) - (a1[1] - a2[1]) * val2
        y /= det

        inter_sect[0][XX] = x
        inter_sect[0][YY] = y

cdef void fast_clip_zoi(Polygon *zoi, real_point ref_pt, real_point clipping_pt, Polygon *buffer) nogil:
    cdef real_point middle_pt, other_pt, delta, line_dir, inter_pt
    cdef int first_vid, second_vid
    cdef bint second_same_side
    cdef bint need_cleanup = False

    if buffer == NULL:
        need_cleanup = True
        buffer = polygon_new()

    # Calculate stuff related to clipping line
    middle_pt[XX] = 0.5 * (clipping_pt[XX] + ref_pt[XX])
    middle_pt[YY] = 0.5 * (clipping_pt[YY] + ref_pt[YY])

    delta[XX] = ref_pt[XX] - clipping_pt[XX]
    delta[YY] = ref_pt[YY] - clipping_pt[YY]

    line_dir[XX] = delta[YY]
    line_dir[YY] = -delta[XX]

    other_pt[XX] = middle_pt[XX] + line_dir[XX]
    other_pt[YY] = middle_pt[YY] + line_dir[YY]

    polygon_empty(buffer)

    # Check if line intersects with the polygon
    for second_vid in range(zoi.size):
        if second_vid == 0:
            first_vid = zoi.size - 1
        else:
            first_vid = second_vid - 1


        second_same_side = same_side(ref_pt, zoi.points[second_vid], middle_pt, other_pt)

        if same_side(zoi.points[first_vid], zoi.points[second_vid], middle_pt, other_pt):
            if not second_same_side:
                continue
        else:

            if not are_parallel(zoi.points[first_vid], zoi.points[second_vid], middle_pt, other_pt):
                # Get the intersection point
                get_intersection(zoi.points[first_vid], zoi.points[second_vid], middle_pt, other_pt, &inter_pt)

                polygon_append(buffer, inter_pt)

        if second_same_side:
            polygon_append(buffer, zoi.points[second_vid])


    # Update zoi
    polygon_copy(zoi, buffer)


    # Destroy buffer polygon if needed
    if need_cleanup:
        polygon_destroy(buffer)


cdef Polygon *fast_get_zoi(real_point ref_pt, real_point *pts, int size) nogil:
    cdef Polygon *zoi = polygon_new()
    cdef Polygon *buffer = polygon_new()
    cdef int i
    cdef real_point tmp_pt

    # Start with a big square
    tmp_pt[XX] = -ZOI_DEFAULT_SIZE
    tmp_pt[YY] = -ZOI_DEFAULT_SIZE
    polygon_append(zoi, tmp_pt)
    tmp_pt[XX] = ZOI_DEFAULT_SIZE
    polygon_append(zoi, tmp_pt)
    tmp_pt[YY] = ZOI_DEFAULT_SIZE
    polygon_append(zoi, tmp_pt)
    tmp_pt[XX] = -ZOI_DEFAULT_SIZE
    polygon_append(zoi, tmp_pt)

    # Clip the ZoI using all the points
    for i in range(size):
        fast_clip_zoi(zoi, ref_pt, pts[i], buffer)

    # Release memory
    polygon_destroy(buffer)
    return zoi



###################################################################################################
#
# Python API
#
###################################################################################################

def get_zoi(real[:] ref_pt, real[:,::1] pts):
    cdef real_point ref_pt_c
    cdef real_point *pts_c = NULL
    cdef int size = pts.shape[0]
    cdef int i
    cdef Polygon *zoi = NULL

    with nogil:
        # Convert Python to C
        ref_pt_c[XX] = ref_pt[XX]
        ref_pt_c[YY] = ref_pt[YY]

        pts_c = <real_point *> malloc(size * sizeof(real_point))
        if pts_c == NULL:
            abort()

        for i in range(size):
            pts_c[i][XX] = pts[i, XX]
            pts_c[i][YY] = pts[i, YY]

        # Get ZoI
        zoi = fast_get_zoi(ref_pt_c, pts_c, size)

        # Free memory
        free(pts_c)

    zoi_py = polygon_as_array(zoi)

    polygon_destroy(zoi)

    return zoi_py

def clip_zoi(real[:, ::1] zoi, real[:] ref_pt, real[:] clipping_pt):
    cdef real_point ref_pt_c
    cdef real_point clipping_pt_c
    cdef Polygon *zoi_polygon = NULL

    with nogil:
        # Convert Python to C
        ref_pt_c[XX] = ref_pt[XX]
        ref_pt_c[YY] = ref_pt[YY]

        clipping_pt_c[XX] = clipping_pt[XX]
        clipping_pt_c[YY] = clipping_pt[YY]

        zoi_polygon = polygon_new_from_memview(zoi)

        fast_clip_zoi(zoi_polygon, ref_pt_c, clipping_pt_c, NULL)

    zoi_py = polygon_as_array(zoi_polygon)

    polygon_destroy(zoi_polygon)

    return zoi_py
