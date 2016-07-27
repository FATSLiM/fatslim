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
#cython: initializedcheck=False
from __future__ import print_function

# Cython processors DEFs
from email.mime.application import MIMEApplication

DEF DIM = 3
DEF XX = 0
DEF YY = 1
DEF ZZ = 2

DEF RET_OK = 1
DEF RET_YES = 1
DEF RET_ERROR = 0
DEF RET_NO = 0

DEF NOTSET=-12345
DEF EPSILON=1e-5
DEF EPSILON2=1e-10

DEF COS_45 = 0.70710678
DEF PI = 3.14159265359

DEF MAX_INTERLEAFLET_DISTANCE = 10 # nm
DEF MAX_LEAFLET_COM_DISTANCE = 5 # nm

DEF STACK_ALLOCATION_INCREMENT = 100

DEF MIN_LEAFLET_SIZE = 15
DEF MIN_NEIGHBORHOOD_SIZE = 4

DEF APL_DEFAULT_CUTOFF = 3.0
DEF APL_DEFAULT_AREA_LIMIT = 10
DEF APL_INTERACTION_ADJUSTEMENT = 0.7 # This constant was roughly estimated from lipid-only bilayer

DEF THICKNESS_DEFAULT_CUTOFF = 6.0
DEF THICKNESS_MAXIMUM_MINMAX_RATIO = 1.5
DEF THICKNESS_MIN_COS_DX = 0.98480775301 # Cos 10 deg
DEF THICKNESS_MIN_COS_NORMAL = 0.98480775301 # Cos 10 deg
DEF THICKNESS_DEBUG_BEADID = 160

DEF OUTPUT_RESOLUTION_DEFAULT = 600
DEF OUTPUT_DPI_DEFAULT = 96


DEF AGGREGATE_PROCESSING_RESET = 0
DEF AGGREGATE_PROCESSING_DONE = 1
DEF AGGREGATE_PROCESSING_REJECTED = 2

# Cython C imports (no Python here!)
from cpython.exc cimport PyErr_CheckSignals
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.math cimport fabs, exp, pow, sqrt
from libc.stdlib cimport malloc, realloc, free, abort
from libc.stdio cimport fprintf, stderr
from libc.string cimport strcmp
from cython.parallel cimport prange

from .typedefs cimport real, real_abs, fsl_int
from .typedefs cimport rvec, rvec_copy, rvec_clear, rvec_norm2, rvec_norm, rvec_dprod, rvec_smul, \
    rvec_normalize, rvec_inc, rvec_cprod, rvec_sub
from .typedefs cimport matrix, mat_from_rvec, invert_mat

cdef extern from "openmp_wrapper.h" nogil:
    cdef fsl_int omp_get_thread_num()


from .core_base cimport Frame, PBCBox, topol_atom_t, topol_residue_t
from .core_base cimport OPENMP_NUM_THREADS
from .core_base cimport fft_unified, fft_allatom, fft_coarse

from .core_ns cimport ns_neighborhood, ns_neighborhood_holder, \
    fast_neighbor_search, free_neighborhood_holder

from .core_geometry cimport real_point, Polygon, polygon_destroy, polygon_new_from_polygon, \
    polygon_get_area, fast_get_clipped_polygon, polygon_as_array, polygon_is_inside, polygon_new, \
    polygon_copy, fast_get_zoi, fast_clip_zoi

# Python imports
import numpy as np

###################################################################################################
#
# Useful functions
#
###################################################################################################
def test_parallelism():
    cdef fsl_int sum=0, i

    print("Testing parallelism with %i threads... " % OPENMP_NUM_THREADS, end="")
    for i in prange(OPENMP_NUM_THREADS, nogil=True, num_threads=OPENMP_NUM_THREADS):
        sum += 1

    if sum == OPENMP_NUM_THREADS:
        print("OK")
        return True
    else:
        print("FAILED!")
        return False

def test_put_atoms_on_plane(Aggregate leaflet, real cutoff=2.0):
    cdef ns_neighborhood_holder *holder = NULL
    cdef fsl_int size = len(leaflet)
    cdef real [:, ::1] coords = leaflet.coords
    cdef real [:, ::1] normals = leaflet.normals
    cdef real_point *coords_2d
    cdef fsl_int neighborhood_size
    holder = fast_neighbor_search(coords,
                                  coords,
                                  leaflet.frame.box,
                                  cutoff)

    coords_2d = <real_point *> malloc(sizeof(real_point) * size)
    py_2d_coords = []
    for i in range(size):
        neighborhood_size = holder.neighborhoods[i].size

        put_atoms_on_plane(&coords[i, XX],
                           &normals[i, XX],
                           holder.neighborhoods[i],
                           coords,
                           leaflet.frame.box,
                           coords_2d)

        tmp_2d_coords = np.empty((neighborhood_size, 2))

        for j in range(neighborhood_size):
            tmp_2d_coords[j, XX] = coords_2d[j][XX]
            tmp_2d_coords[j, YY] = coords_2d[j][YY]
        py_2d_coords.append(tmp_2d_coords)

    free(coords_2d)

    return py_2d_coords

cdef real fast_average(real[:] vals) nogil:
    cdef fsl_int i, size = vals.shape[0]
    cdef real avg = 0

    if size == 0:
        return avg

    for i in prange(size, schedule="static", num_threads=OPENMP_NUM_THREADS):
        avg += vals[i]

    return avg / size

cdef inline void rvec_to_basis(rvec v, matrix conv_mat, rvec result) nogil:
    cdef real v1, v2, v3

    v1 = conv_mat[XX][XX] * v[XX] + conv_mat[XX][YY] * v[YY] + conv_mat[XX][ZZ] * v[ZZ]
    v2 = conv_mat[YY][XX] * v[XX] + conv_mat[YY][YY] * v[YY] + conv_mat[YY][ZZ] * v[ZZ]
    v3 = conv_mat[ZZ][XX] * v[XX] + conv_mat[ZZ][YY] * v[YY] + conv_mat[ZZ][ZZ] * v[ZZ]

    result[XX] = v1
    result[YY] = v2
    result[ZZ] = v3


cdef real rvec_cprod_2d(rvec a, rvec b) nogil:
    return a[XX] * b[YY] - a[YY] * b[XX]

cdef inline fsl_int left_of(rvec v1, rvec v2, rvec v3) nogil:
    cdef rvec tmp1, tmp2
    cdef real x

    tmp1[XX] = v2[XX] - v1[XX]
    tmp1[YY] = v2[YY] - v1[YY]

    tmp2[XX] = v3[XX] - v2[XX]
    tmp2[YY] = v3[YY] - v2[YY]

    x = rvec_cprod_2d(tmp1, tmp2)

    if x != 0:
        return -1 if x < 0 else 1
    else:
        return 0

cdef inline bint same_side(rvec l0, rvec l1, rvec pt1, rvec pt2) nogil:
    cdef real a = l0[YY] - l1[YY]
    cdef real b = l1[XX] - l0[XX]
    return (a * (pt1[XX] - l0[XX]) + b * (pt1[YY] - l0[YY])) * (a * (pt2[XX] - l0[XX]) + b * (pt2[YY] - l0[YY])) > 0

cdef inline bint line_segment_intersect_2d(rvec l0, rvec l1, rvec s1, rvec s2, rvec inter) nogil except *:
    # See http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect for details
    cdef rvec l_dir, s_dir
    cdef real cprod, s

    # Get line and segment direction
    rvec_sub(l1, l0, l_dir)
    rvec_sub(s2, s1, s_dir)

    cprod = l_dir[XX] * s_dir[YY] - l_dir[YY] * s_dir[XX]

    if cprod * cprod < EPSILON2: # line and segment are parallel
        return False

    s = (- l_dir[YY] * (l0[XX] - s1[XX]) + l_dir[XX] * (l0[YY] - s1[YY])) / cprod

    if 0 <= s <= 1: # Intersection point is inside the segment
        inter[XX] = s1[XX] + s * s_dir[XX]
        inter[YY] = s1[YY] + s * s_dir[YY]
        return True
    return False

cdef inline void complete_basis(rvec z_axis, rvec x_axis, rvec y_axis) nogil:
    cdef real a, b
    rvec_normalize(z_axis)

    if z_axis[ZZ] < -0.9999999:
        x_axis[XX] = 0.0
        x_axis[YY] = -1.0
        x_axis[ZZ] = 0.0

        y_axis[XX] = -1.0
        y_axis[YY] = 0.0
        y_axis[ZZ] = 0.0

    else:
        a = 1.0  / (1.0 + z_axis[ZZ])
        b = -z_axis[XX] * z_axis[YY] * a

        x_axis[XX] = 1.0 - z_axis[XX] * z_axis[XX] * a
        x_axis[YY] = b
        x_axis[ZZ] = -z_axis[XX]

        y_axis[XX] = b
        y_axis[YY] = 1.0 - z_axis[YY] * z_axis[YY] * a
        y_axis[ZZ] = -z_axis[YY]

#    if (real_abs(z_axis[XX]) >= real_abs(z_axis[YY])) and real_abs(z_axis[XX]) > EPSILON:
#        x_axis[XX] = -z_axis[ZZ]
#        x_axis[YY] = 0
#        x_axis[ZZ] = z_axis[XX]
#    else:
#        x_axis[XX] = 0
#        x_axis[YY] = z_axis[ZZ]
#        x_axis[ZZ] = -z_axis[YY]

#    rvec_normalize(x_axis)
#    rvec_cprod(z_axis, x_axis, y_axis)

cdef bint check_leaflet_compatibility(Frame frame, Aggregate leaflet1, Aggregate leaflet2) nogil:
    cdef real[:,::1] coords = frame.fast_get_bead_coords_bbox()

    # We are seeking for a planar membrane
    if leaflet1.is_planar:
        # So leaflet 2 must also be planar
        if not leaflet2.is_planar:
            return False

        # And distance between the leaflets centers of mass must be small enough
        if frame.box.fast_leaflet_distance(leaflet1.xcm, leaflet2.xcm, leaflet1.avg_normal)\
                > MAX_INTERLEAFLET_DISTANCE:
            return False

        # And average normal should be roughly antiparallel
        if rvec_dprod(leaflet1.avg_normal, leaflet2.avg_normal) > -COS_45:
            return False
    else: # We are seeking for a vesicle

        # So leaflet 2 must NOT be planar
        if leaflet2.is_planar:
            return False

        # And distance between the leaflets centers of mass must be small
        if frame.box.fast_distance(leaflet1.xcm, leaflet2.xcm) > MAX_LEAFLET_COM_DISTANCE:
           return False

    # If we are here, it may be safe to say that the leaflets may belong to the same membrane
    return True


cdef void qsort(real *a, fsl_int start, fsl_int end) nogil:
    if (end - start) < 17:
        insertion_sort(a, start, end)
        return
    cdef fsl_int boundary = partition(a, start, end)
    qsort(a, start, boundary)
    qsort(a, boundary+1, end)

cdef fsl_int partition(real *a, fsl_int start, fsl_int end) nogil:
    cdef fsl_int i = start, j = end-1
    cdef real pivot = a[j]
    while True:
        # assert all(x < pivot for x in a[start:i])
        # assert all(x >= pivot for x in a[j:end])

        while a[i] < pivot:
            i += 1
        while i < j and pivot <= a[j]:
            j -= 1
        if i >= j:
            break
        #assert a[j] < pivot <= a[i]
        swap(a, i, j)
        #assert a[i] < pivot <= a[j]
    #assert i >= j and i < end
    swap(a, i, end-1)
    #assert a[i] == pivot
    # assert all(x < pivot for x in a[start:i])
    # assert all(x >= pivot for x in a[i:end])
    return i

cdef inline void swap(real *a, fsl_int i, fsl_int j) nogil:
    a[i], a[j] = a[j], a[i]

cdef void insertion_sort(real *a, fsl_int start, fsl_int end) nogil:
    cdef fsl_int i, j
    cdef real v
    for i in range(start, end):
        #invariant: [start:i) is sorted
        v = a[i]; j = i-1
        while j >= start:
            if a[j] <= v: break
            a[j+1] = a[j]
            j -= 1
        a[j+1] = v

cdef void xcm_from_neighborhood(rvec ref,
                                rvec ref_normal,
                                real[:, ::1]coords,

                                ns_neighborhood *neighborhood,
                                PBCBox box,
                                rvec xcm) nogil:
    cdef fsl_int i, d
    cdef rvec dx
    cdef real total_weight=0.0
    cdef real dx_norm, weight

    # Initialize xcm
    rvec_clear(xcm)
    #rvec_smul(neighborhood.size, xcm, xcm)

    # fprintf(stderr,
    #         "DEBUG XCM neighborhood: ref(%.3f, %.3f, %.3f) -> %i neighbors\n",
    #         ref[XX], ref[YY], ref[ZZ],
    #         <int> neighborhood.size)

    # Step 1: Get XCM
    for i in range(neighborhood.size):
        box.fast_pbc_dx(ref, &coords[neighborhood.beadids[i], XX], dx)

        dx_norm = rvec_norm(dx)

        weight = 1 - dx_norm / neighborhood.cutoff
        weight = 1
        total_weight += weight

        # fprintf(stderr,
        #     "    -> #%i coords(%.3f, %.3f, %.3f) -> dx: (%.3f, %.3f, %.3f) -> norm: %.3f\n",
        #     <int> i +1,
        #     coords[neighborhood.beadids[i], XX], coords[neighborhood.beadids[i], YY], coords[neighborhood.beadids[i], ZZ],
        #     dx[XX], dx[YY], dx[ZZ],
        #         dx_norm)

        rvec_smul(weight, dx, dx)
        rvec_inc(xcm, dx)

    rvec_smul(1.0/total_weight, xcm, xcm)

    # fprintf(stderr,
    #         "DEBUG XCM neighborhood: ref(%.3f, %.3f, %.3f) -> avg_dx: (%.3f, %.3f, %.3f)\n",
    #         ref[XX], ref[YY], ref[ZZ],
    #         xcm[XX], xcm[YY], xcm[ZZ])
    rvec_inc(xcm, ref)

    # Step 2: Make sure it is inside the brick-shaped box
    for i in range(DIM - 1, -1, -1):
        while xcm[i] < 0:
            for d in range(i+1):
                xcm[d] += box.c_pbcbox.box[i][d]
        while xcm[i] >= box.c_pbcbox.box[i][i]:
            for d in range(i+1):
                xcm[d] -= box.c_pbcbox.box[i][d]

cdef real thickness_from_neighborhood(rvec ref, rvec ref_normal,
                                      real[:, ::1]other_coords,
                                      real[:, ::1]same_coords,
                                      real[:, ::1]other_normals,
                                      real[:, ::1]same_normals,
                                      ns_neighborhood *neighborhood_other_leaflet,
                                      ns_neighborhood *neighborhood_same_leaflet,
                                      PBCBox box) nogil:
    cdef real ref_thickness = NOTSET, other_thickness = NOTSET
    cdef rvec dx
    cdef real ref_max_cos = NOTSET, cos_trial
    cdef real other_max_cos = NOTSET
    cdef real dprod_normal, dprod_dx, dx_norm
    cdef real total_weight, weight
    cdef fsl_int i

    cdef fsl_int n=0
    cdef real dprod_normals
    cdef rvec other_normal, other_coord

    cdef real avg_thickness = 0.0
    cdef real thickness
    cdef real min_thickness = 1e3
    cdef real max_thickness = -1

    cdef rvec ref_xcm
    cdef rvec avg_dx
    cdef int n_used = 0


    # First: get neighborhood xcm
    rvec_clear(ref_xcm)

    # Step 1: Get XCM
    for i in range(neighborhood_same_leaflet.size):
        dprod_normal = rvec_dprod(&same_normals[neighborhood_same_leaflet.beadids[i], XX],
                                  ref_normal)

        if dprod_normal < THICKNESS_MIN_COS_NORMAL:
            continue

        box.fast_pbc_dx(ref, &same_coords[neighborhood_same_leaflet.beadids[i], XX], dx)
        dx_norm = rvec_norm(dx)

        if dx_norm < EPSILON:
            continue

        weight = weight = (real_abs(dprod_normal) - THICKNESS_MIN_COS_NORMAL) / (1.0 - THICKNESS_MIN_COS_NORMAL)
        total_weight += weight

        rvec_smul(weight, dx, dx)
        rvec_inc(ref_xcm, dx)

    if total_weight > EPSILON:
        rvec_smul(1.0/total_weight, ref_xcm, ref_xcm)

    rvec_inc(ref_xcm, ref)

    # Step 2: Make sure it is inside the brick-shaped box
    for i in range(DIM - 1, -1, -1):
        while ref_xcm[i] < 0:
            for d in range(i+1):
                ref_xcm[d] += box.c_pbcbox.box[i][d]
        while ref_xcm[i] >= box.c_pbcbox.box[i][i]:
            for d in range(i+1):
                ref_xcm[d] -= box.c_pbcbox.box[i][d]


    box.fast_pbc_dx(ref_xcm, ref,dx)
    rvec_clear(avg_dx)

    total_weight = 0.0
    for i in range(neighborhood_other_leaflet.size):
        rvec_copy(&other_coords[neighborhood_other_leaflet.beadids[i], XX], other_coord)
        rvec_copy(&other_normals[neighborhood_other_leaflet.beadids[i], XX], other_normal)

        dprod_normal = rvec_dprod(other_normal, ref_normal)

        # Make sure that both beads belong to different leaflet
        if dprod_normal > -COS_45:
            continue

        # Get the distance (through the bilayer) between the ref bead and its twin
        box.fast_pbc_dx_leaflet(ref_xcm, other_coord,
                                dx,
                                ref_normal)
        dx_norm = rvec_norm(dx)

        # we check dx because the image which is on the right side of the bilayer may not be a good twin (too far)
        if dx_norm > neighborhood_other_leaflet.cutoff:
            continue

        dprod_dx = real_abs(rvec_dprod(dx, ref_normal))
        cos_trial = dprod_dx / dx_norm

        if cos_trial < THICKNESS_MIN_COS_DX:
            continue

        weight = 1 - dx_norm / neighborhood_other_leaflet.cutoff
        weight = (real_abs(dprod_normal) - THICKNESS_MIN_COS_NORMAL) / (1.0 - THICKNESS_MIN_COS_NORMAL)
        #weight = 1
        if weight > 0:
            rvec_smul(weight, dx, dx)
            total_weight += weight

            rvec_inc(avg_dx, dx)
            n_used += 1

    if total_weight < EPSILON:
        avg_thickness = NOTSET
    else:
        rvec_smul(1.0/total_weight, avg_dx, avg_dx)

        dprod_dx = real_abs(rvec_dprod(avg_dx, ref_normal))
        dx_norm = rvec_norm(avg_dx)

        cos_trial = dprod_dx / dx_norm

        avg_thickness = dprod_dx

    return avg_thickness


cdef void put_atoms_on_plane(rvec ref, rvec ref_normal,
                             ns_neighborhood *neighborhood,
                             real[:,::1] coords,
                             PBCBox box,
                             real_point *coords_2d) nogil:
    cdef fsl_int size = neighborhood.size
    cdef fsl_int i
    cdef rvec plane_z, plane_y, plane_x
    cdef matrix conv_mat_revert, conv_mat
    cdef rvec dx, proj

    # Build plane basis
    rvec_copy(ref_normal, plane_z)
    complete_basis(plane_z, plane_x, plane_y)

    # Get conversion matrix
    mat_from_rvec(plane_x, plane_y, plane_z, conv_mat_revert)
    invert_mat(conv_mat_revert, conv_mat)

    for i in range(size):
        # Get the closest image of the neighbor
        box.fast_pbc_dx(ref, &coords[neighborhood.beadids[i], XX], dx)

        # Get the coordinates into the new basis (plane_x, plane_y, plane_z)
        rvec_to_basis(dx, conv_mat, proj)

        # Store the projection onto the (plane_x, plane_y) plane
        coords_2d[i][XX] = proj[XX]
        coords_2d[i][YY] = proj[YY]



cdef real apl_from_neighborhoods(real_point *coords_2d,
                                 fsl_int size,
                                 real_point *interacting_coords_2d,
                                 real *interacting_weights,
                                 fsl_int interacting_size,
                                 real area_limit) nogil:
    cdef fsl_int i, n, j, n_interactions
    cdef real lipid_area, cell_area
    cdef Polygon *polygon_lipids
    cdef Polygon *interacting_polygon = NULL
    cdef Polygon *intersecting_polygon = NULL
    cdef real area_ratio, intersecting_area, cur_intersecting_area
    cdef real delta_area
    cdef Polygon *cell
    cdef Polygon *clipped_cell=NULL
    cdef real interaction_factor
    cdef real area, small_cutoff, big_cutoff, clipped_area
    cdef bint keep, need_interaction
    cdef real total_weight
    cdef real weight, avg_dist2, weighting_factor

    cdef real_point ref_2d, xcm_interacting
    cdef rvec dx

    if size < 4:
        return -1.0

    ref_2d[XX] = 0
    ref_2d[YY] = 0

    cell = fast_get_zoi(ref_2d, coords_2d, size)
    area = polygon_get_area(cell)

    # Handle interacting atoms
    interaction_factor = 1.0
    if interacting_size > 0:
        # First: put interacting coords on the 2D plane
        #interacting_coords_2d = put_atoms_on_plane(ref, ref_normal, interacting_neighbors, interacting_coords, box)

        total_weight = 0
        xcm_interacting[XX] = 0
        xcm_interacting[YY] = 0
        for i in range(interacting_size):
            if polygon_is_inside(cell, interacting_coords_2d[i]):
                total_weight += interacting_weights[i]

                xcm_interacting[XX] += interacting_weights[i] * interacting_coords_2d[i][XX]
                xcm_interacting[YY] += interacting_weights[i] * interacting_coords_2d[i][YY]

        if total_weight > 0:
            xcm_interacting[XX] /= total_weight
            xcm_interacting[YY] /= total_weight

            clipped_cell = polygon_new_from_polygon(cell)

            fast_clip_zoi(clipped_cell, ref_2d, xcm_interacting,NULL)

            clipped_area = polygon_get_area(clipped_cell)

            weighting_factor = 1
            interaction_factor = weighting_factor * clipped_area/area + (1 - weighting_factor)

            polygon_destroy(clipped_cell)

    polygon_destroy(cell)

    if area < 0: # No valid cell, no area!
        #fprintf(stderr, "WARNING: No valid area from %i neighbors\n", size)
        return -1.0
    elif area > area_limit: # Probably no a valid area
        #fprintf(stderr, "WARNING: Area is probably wrong: %.3f nm^2 (>%.3f nm^2) -> skipped\n", area, area_limit)
        return -1.0

    area *= interaction_factor

    return area

########################################################################################################################
#
# Classes
#
########################################################################################################################
cdef class Aggregate(object):
    def __init__(self, Frame frame not None, fsl_int[:] beadids not None):
        cdef fsl_int size = beadids.shape[0]
        cdef fsl_int i,j, current_id

        self.frame = frame

        # Make sure that bead IDs are sorted in ascending order
        tmp_beadids = np.asarray(beadids.copy())
        tmp_beadids.sort()
        self.beadids = tmp_beadids

        # Get headgroup atomids
        hg_atomids = []
        for i in range(size):
            current_id = self.beadids[i]
            for j in range(self.frame.trajectory.hg_bead_atomids_offsets[current_id], self.frame.trajectory.hg_bead_atomids_offsets[current_id+1]):
                hg_atomids.append(self.frame.trajectory.hg_bead_atomids[j])
        self.hg_atomids = np.array(hg_atomids, dtype=np.int64)

        # Get lipid atomids
        lipid_atomids = []
        for i in range(size):
            current_id = self.beadids[i]
            for j in range(self.frame.trajectory.lipid_atomids_offsets[current_id], self.frame.trajectory.lipid_atomids_offsets[current_id+1]):
                lipid_atomids.append(self.frame.trajectory.lipid_atomids[j])
        self.lipid_atomids = np.array(lipid_atomids, dtype=np.int64)

        # Group lipid by type
        lipids_by_type = {}
        nlipids_by_type = {}
        self.resids = np.empty(size, dtype=np.int64)
        for i in range(size):
            resname = self.get_resname(i)
            self.resids[i] = self.get_resid(i)
            try:
                lipids_by_type[resname].append(i)
            except KeyError:
                lipids_by_type[resname] = [i, ]
        for key, val in lipids_by_type.items():
            lipids_by_type[key] = np.array(val, dtype=np.int64)
            nlipids_by_type[key] = len(val)
        self.lipid_types = list(lipids_by_type.keys())
        self.beadids_by_type = list(lipids_by_type.values())
        self.nlipids_by_type = nlipids_by_type

        # Allocate memory for coords, normals and directions
        self.coords = np.empty([size,3])
        self.directions = np.empty([size,3])
        self.normals = np.empty([size,3])
        self.neighborhood_normals = np.empty([size,3])
        for i in range(DIM):
            self.avg_normal[i] = NOTSET

        self.neighborhood_cutoff = self.frame.proximity_cutoff
        self.neighborhoods = NULL

        # APL-related
        self.apl_values = np.empty(size)
        self.apl_by_types = []
        for val in self.beadids_by_type:
            self.apl_by_types.append(np.empty(len(val)))
        self.area_by_types = np.empty(len(self.lipid_types))
        self.area = NOTSET
        self.apl_min = NOTSET
        self.apl_max = NOTSET
        self.apl_avg = NOTSET
        self.apl_cutoff = NOTSET

        # thickness-related
        self.thickness_values = np.empty(size)
        self.thickness_avg = NOTSET
        self.thickness_cutoff = NOTSET
        self.thickness_min = NOTSET
        self.thickness_max = NOTSET


        with nogil:
            self.refresh_cache()

            if size > MIN_LEAFLET_SIZE and rvec_norm2(self.avg_normal) > 0.1:
                self.is_planar = True
            else:
                self.is_planar = False


    def __dealloc__(self):
        free_neighborhood_holder(self.neighborhoods)
        self.neighborhoods = NULL

    cdef void refresh_cache(self) nogil:
        cdef fsl_int i, size=self.beadids.shape[0]
        cdef real[:, ::1] coords
        cdef real[:, ::1] normals
        cdef real[:, ::1] directions
        cdef real[:, ::1] self_coords = self.coords
        cdef real[:, ::1] self_normals = self.normals
        cdef real[:, ::1] self_directions = self.directions
        cdef fsl_int[:] self_beadids = self.beadids
        cdef rvec dx
        cdef fsl_int beadid

        coords = self.frame.fast_get_bead_coords_bbox()
        normals = self.frame.fast_get_normals()
        directions = self.frame.fast_get_directions()
        if self.neighborhood_cutoff < 0:
            self.neighborhood_cutoff = self.frame.proximity_cutoff

        # Compute center of mass & average_normal
        rvec_clear(self.xcm)
        for i in range(size):
            beadid = self_beadids[i]
            self.frame.box.fast_pbc_dx(self.frame.box.center, &coords[beadid, XX], dx)
            rvec_inc(self.xcm, dx)

            rvec_copy(&coords[beadid, XX], &self_coords[i, XX])
            rvec_copy(&normals[beadid, XX], &self_normals[i, XX])
            rvec_copy(&directions[beadid, XX], &self_directions[i, XX])
        rvec_smul(1.0/size, self.xcm, self.xcm)
        rvec_inc(self.xcm, self.frame.box.center)

        # Update the smoothed normals
        self.set_neighborhood_cutoff(NOTSET, force=True)

    cdef fsl_int fast_get_atomid(self, fsl_int internalid) nogil:
        cdef fsl_int current_id = self.beadids[internalid]
        cdef fsl_int offset = self.frame.trajectory.hg_bead_atomids_offsets[current_id]
        return self.frame.trajectory.hg_bead_atomids[offset]

    cpdef fsl_int get_atomid(self, fsl_int internalid):
        return self.fast_get_atomid(internalid)

    def get_residue(self, fsl_int internalid):
        cdef fsl_int atomid = self.fast_get_atomid(internalid)
        return self.frame.trajectory.topology.get_residue_from_atomid(atomid)

    cpdef str get_resname(self, fsl_int internalid):
        cdef fsl_int atomid = self.fast_get_atomid(internalid)
        cdef topol_residue_t *residue = self.frame.trajectory.topology.fast_get_residue_from_atomid(atomid)

        if residue == NULL:
            raise RuntimeError

        return str(self.frame.trajectory.topology.resnames[residue.name_id].decode())

    cpdef fsl_int get_resid(self, fsl_int internalid):
        cdef fsl_int atomid = self.fast_get_atomid(internalid)
        return self.frame.trajectory.topology.fast_get_resid_from_atomid(atomid)

    cdef bint fast_same_restype(self, fsl_int internalid1, fsl_int internalid2) nogil:
        cdef fsl_int atomid1 = self.fast_get_atomid(internalid1)
        cdef fsl_int atomid2 = self.fast_get_atomid(internalid2)
        cdef topol_residue_t *residue1 = self.frame.trajectory.topology.fast_get_residue_from_atomid(atomid1)
        cdef topol_residue_t *residue2 = self.frame.trajectory.topology.fast_get_residue_from_atomid(atomid2)

        if residue1 == NULL or residue2 == NULL:
            return False

        return residue1.name_id == residue2.name_id

    cpdef bint same_restype(self, fsl_int internalid1, fsl_int internalid2):
        return self.fast_same_restype(internalid1, internalid2)

    cdef void compute_neighborhood_normal(self, fsl_int index) nogil:
        cdef fsl_int i, j, curid
        cdef rvec neighborhood_normal, ref_coords
        cdef ns_neighborhood *neighborhood
        cdef real[:, ::1] self_normals = self.normals
        cdef real[:, ::1] self_coords = self.coords

        rvec_copy(&self_coords[index, XX], ref_coords)

        # Initialize smoothed normal with the raw value
        rvec_copy(&self_normals[index, XX], neighborhood_normal)

        # Get neighborhood
        neighborhood = self.neighborhoods.neighborhoods[index]

        for i in range(neighborhood.size):
            curid = neighborhood.beadids[i]

            for j in range(DIM):
                neighborhood_normal[j] += self_normals[curid, j]

        rvec_normalize(neighborhood_normal)

        # Store the smoothed normal
        rvec_copy(neighborhood_normal, &self.neighborhood_normals[index, XX])


    cdef void set_neighborhood_cutoff(self, real cutoff, bint force=False) nogil:
        cdef fsl_int i
        cdef fsl_int size = self.neighborhood_normals.shape[0]
        cdef real[:, ::1] self_normals = self.neighborhood_normals

        if cutoff < 0:
            cutoff = self.neighborhood_cutoff

        # Delete old ref if needed
        if self.neighborhoods != NULL:
            free_neighborhood_holder(self.neighborhoods)
            self.neighborhoods = NULL

        # Update Neighborhoods
        self.neighborhoods = fast_neighbor_search(self.coords,
                                                  self.coords,
                                                  self.frame.box,
                                                  self.neighborhood_cutoff)

        for i in prange(size, schedule="dynamic", num_threads=OPENMP_NUM_THREADS):
            self.compute_neighborhood_normal(i)

        rvec_clear(self.avg_normal)
        for i in range(size):
            rvec_inc(self.avg_normal, &self_normals[i,XX])
        rvec_smul(1.0/size, self.avg_normal, self.avg_normal)

    cdef real fix_apl(self,
                      fsl_int refid,
                      real *tmp_apl,
                      bint onlyfix=False) nogil:
        cdef real fixed_apl, cur_apl
        cdef real total_weight, weight
        cdef real distance2
        cdef fsl_int i, curid
        cdef real[:, ::1] coords = self.coords
        cdef ns_neighborhood *neighborhood = self.neighborhoods.neighborhoods[refid]
        cdef real cutoff2 = neighborhood.cutoff
        cutoff2 *= cutoff2

        # Initialize fixed apl to current value, if it is set
        cur_apl = tmp_apl[refid]
        if cur_apl > 0: # APL is set
            fixed_apl = cur_apl
            total_weight = 1.0
            if onlyfix:
                return fixed_apl
        else:
            total_weight = 0.0
            fixed_apl = 0.0

        # Loop over all the neighbors
        for i in range(neighborhood.size):
            curid = neighborhood.beadids[i]
            cur_apl = tmp_apl[curid]

            if cur_apl < 0: # Don't use an APL which is not set
                continue

            if not self.fast_same_restype(refid, curid): # Use only beads that have similar type
                continue

            distance2 = self.frame.box.fast_distance2(&coords[refid, XX], &coords[curid, XX])
            #weight = 1 - pow(distance2/cutoff2, 0.25)
            #weight = 1 - distance2/cutoff2
            #weight = exp(-5 * distance2/cutoff2)
            weight = 1 - sqrt(distance2/cutoff2)
            #weight = 1

            fixed_apl += weight * cur_apl
            total_weight += weight

        if total_weight > 0:
            return fixed_apl / total_weight
        else:
            # We have no choice...
            # Ugly solution to get an APL: loop over all lipids of the same type and use the average APL as default value
            total_weight = 0.0
            fixed_apl = 0.0
            for curid in range(self.fast_size()):
                cur_apl = tmp_apl[curid]

                if cur_apl < 0: # Don't use an APL which is not set
                    continue

                if not self.fast_same_restype(refid, curid): # Use only beads that have similar type
                    continue

                fixed_apl += cur_apl
                total_weight += 1

            if total_weight > 0:
                fixed_apl /= total_weight
            else:
                fixed_apl = -1
            return fixed_apl

    cdef real fix_thickness(self,
                      fsl_int refid,
                      real *tmp_thickness,
                      bint onlyfix=True) nogil:
        cdef real fixed_thickness, cur_thickness
        cdef real total_weight, weight
        cdef real distance2
        cdef fsl_int i, curid
        cdef real[:, ::1] coords = self.coords
        cdef ns_neighborhood *neighborhood = self.neighborhoods.neighborhoods[refid]
        cdef real cutoff2 = neighborhood.cutoff
        cutoff2 *= cutoff2

        # Initialize fixed thickness to current value, if it is set
        cur_thickness = tmp_thickness[refid]
        if cur_thickness > 0: # Thickness is set
            fixed_thickness = cur_thickness
            total_weight = 1.0
            if onlyfix:
                return fixed_thickness
        else:
            total_weight = 0.0
            fixed_thickness = 0.0

        # Loop over all the neighbors
        for i in range(neighborhood.size):
            curid = neighborhood.beadids[i]
            cur_thickness = tmp_thickness[curid]

            if cur_thickness < 0: # Don't use an thickness which is not set
                continue

            distance2 = self.frame.box.fast_distance2(&coords[refid, XX], &coords[curid, XX])
            weight = 1 - pow(distance2/cutoff2, 0.25)

            fixed_thickness += weight * cur_thickness
            total_weight += weight


        if total_weight > EPSILON:
            return fixed_thickness / total_weight
        else:
            fprintf(stderr,
                    "WARNING: Could not fix thickness for lipid #%i! "
                    "(%i neighbors - total weight: %.3f)\n",
                    <int> refid, <int> neighborhood.size, total_weight)
            return -1.0


    cdef void compute_apl(self,
                          real cutoff,
                          real area_limit=APL_DEFAULT_AREA_LIMIT,
                          bint force=False) nogil except *:
        cdef fsl_int i, j, tid, beadid
        cdef fsl_int size = self.fast_size()
        cdef ns_neighborhood_holder *holder_interacting = NULL
        cdef ns_neighborhood_holder *self_neighborhoods = NULL
        cdef real[:, ::1] interacting_coords = self.frame.fast_get_interacting_group_coords_bbox()
        cdef ns_neighborhood *interacting_neighbors = NULL
        cdef real *tmp_apl
        cdef real[:, ::1] self_coords = self.coords
        cdef real[:, ::1] self_neighborhoods_normals
        cdef real[:] self_apl = self.apl_values
        cdef real apl_avg, apl_min, apl_max, area
        cdef real_point *lipid_coords_2d = NULL
        cdef fsl_int lipid_buffer_size = 0, lipid_offset
        cdef real_point *interacting_coords_2d = NULL
        cdef real *interacting_weights = NULL
        cdef fsl_int interacting_buffer_size = 0, interacting_size, interacting_offset
        cdef rvec *dx = NULL


        # Do nothing if the thickness is already computed and the force flag is not set.
        if self.apl_avg != NOTSET and not force and (fabs(self.apl_cutoff - cutoff) < EPSILON or cutoff < 0):
            return

        # Avoid breaking anything because the cutoff is crap
        if cutoff < 0:
            cutoff = APL_DEFAULT_CUTOFF

        # normals will only be updated if necessary so it is safe to call the method without further test
        self.set_neighborhood_cutoff(cutoff)
        self_neighborhoods_normals = self.neighborhood_normals
        self_neighborhoods = self.neighborhoods

        # Allocate memory for lipid buffer
        lipid_buffer_size = 0
        for i in range(size):
            if self_neighborhoods.neighborhoods[i].size > lipid_buffer_size:
                lipid_buffer_size = self_neighborhoods.neighborhoods[i].size
        lipid_coords_2d = <real_point *> malloc(lipid_buffer_size * sizeof(real_point) * OPENMP_NUM_THREADS)
        if lipid_coords_2d == NULL:
            abort()

        # Retrieve interacting neighbors if they exist
        if interacting_coords.shape[0] > 0:
            holder_interacting = fast_neighbor_search(self_coords,
                                      interacting_coords,
                                      self.frame.box,
                                      cutoff)
            interacting_buffer_size = 0
            for i in range(size):
                if holder_interacting.neighborhoods[i].size > interacting_buffer_size:
                    interacting_buffer_size = holder_interacting.neighborhoods[i].size

            interacting_coords_2d = <real_point *> malloc(interacting_buffer_size * sizeof(real_point) * OPENMP_NUM_THREADS)
            if interacting_coords_2d == NULL:
                abort()

            interacting_weights = <real *> malloc(interacting_buffer_size * sizeof(real) * OPENMP_NUM_THREADS)
            if interacting_weights == NULL:
                abort()

            dx = <rvec *> malloc(OPENMP_NUM_THREADS * sizeof(rvec))

        # Allocate memory for temp arrays
        tmp_apl = <real *> malloc(size * sizeof(real))
        if tmp_apl == NULL:
            abort()

        # Calculate raw area per lipid
        for i in prange(size, schedule="dynamic", num_threads=OPENMP_NUM_THREADS):
            lipid_offset = omp_get_thread_num() * lipid_buffer_size
            # Put lipid atoms on 2D plane
            put_atoms_on_plane(&self_coords[i, XX],
                               &self_neighborhoods_normals[i, XX],
                               self_neighborhoods.neighborhoods[i],
                               self_coords,
                               self.frame.box,
                               &lipid_coords_2d[lipid_offset])

            if holder_interacting != NULL:
                interacting_offset = omp_get_thread_num() * interacting_buffer_size
                interacting_size = holder_interacting.neighborhoods[i].size

                # Put interacting atoms on 2D plane
                put_atoms_on_plane(&self_coords[i, XX],
                                   &self_neighborhoods_normals[i, XX],
                                   holder_interacting.neighborhoods[i],
                                   interacting_coords,
                                   self.frame.box,
                                   &interacting_coords_2d[interacting_offset])

                # Compute weights of neighbors
                for j in range(interacting_size):
                    self.frame.box.fast_pbc_dx(&self_coords[i, XX],
                                               &interacting_coords[holder_interacting.neighborhoods[i].beadids[j], XX],
                                               &dx[omp_get_thread_num()][XX])

                    interacting_weights[interacting_offset + j] = 1 -\
                                                        fabs(dx[omp_get_thread_num()][XX]
                                                             * self_neighborhoods_normals[i, XX] +
                                                             dx[omp_get_thread_num()][YY]
                                                             * self_neighborhoods_normals[i, YY] +
                                                             dx[omp_get_thread_num()][ZZ]
                                                             * self_neighborhoods_normals[i, ZZ]) / cutoff

                tmp_apl[i] = apl_from_neighborhoods(&lipid_coords_2d[lipid_offset],
                                                    self_neighborhoods.neighborhoods[i].size,
                                                    &interacting_coords_2d[interacting_offset],
                                                    &interacting_weights[interacting_offset],
                                                    interacting_size,
                                                    area_limit)
            else:
                tmp_apl[i] = apl_from_neighborhoods(&lipid_coords_2d[omp_get_thread_num() * lipid_buffer_size],
                                                    self_neighborhoods.neighborhoods[i].size,
                                                    NULL,
                                                    NULL,
                                                    0,
                                                    area_limit)

        # Fix not valid area per lipid using neighbors
        area = 0
        apl_avg = 0
        for i in prange(size, schedule="dynamic", num_threads=OPENMP_NUM_THREADS):
        #for i in range(size):
            self_apl[i] = self.fix_apl(i, tmp_apl)

            area += self_apl[i]
            apl_avg += self_apl[i]
        apl_avg /= size

        # Find extrema
        apl_max = NOTSET
        apl_min = -NOTSET
        for i in range(size):
            if self.apl_values[i] > apl_max:
               apl_max = self.apl_values[i]
            if self.apl_values[i] < apl_min:
               apl_min = self.apl_values[i]

        # Free memory
        if holder_interacting != NULL:
            free_neighborhood_holder(holder_interacting)
            free(dx)
            free(interacting_weights)
            free(interacting_coords_2d)
        free(lipid_coords_2d)
        free(tmp_apl)

        # Store results
        self.apl_cutoff = cutoff
        self.apl_avg = apl_avg
        self.area = area
        self.apl_max = apl_max
        self.apl_min = apl_min
        with gil:
            for tid, val in enumerate(self.beadids_by_type):
                self.area_by_types[tid] = 0
                for i in range(len(val)):
                    beadid = val[i]
                    self.apl_by_types[tid][i] = self.apl_values[beadid]
                    self.area_by_types[tid] += self.apl_values[beadid]




    cdef void compute_thickness(self,
                                Aggregate other,
                                real interleaflet_cutoff=THICKNESS_DEFAULT_CUTOFF,
                                bint force=False) nogil except*:
        cdef ns_neighborhood_holder *holder_other_leaflet
        cdef fsl_int size = self.fast_size()
        cdef real *tmp_thickness
        cdef fsl_int i
        cdef real[:, ::1] self_coords = self.coords
        cdef real[:, ::1] self_neighborhood_normals = self.neighborhood_normals
        cdef real[:, ::1] other_coords = other.coords
        cdef real[:, ::1] other_neighborhood_normals = other.neighborhood_normals
        cdef real avg_thickness, thickness_min, thickness_max

        # Do nothing if the thickness is already computed and the force flag is not set.
        if self.thickness_avg != NOTSET and\
                not force and\
                (fabs(self.thickness_avg - interleaflet_cutoff) < EPSILON
                 or interleaflet_cutoff < 0):
            return

        # Avoid breaking anything because the cutoff is crap
        if interleaflet_cutoff < 0:
            return

        # Retrieve potential twins
        holder_other_leaflet = fast_neighbor_search(self_coords,
                                      other_coords,
                                      self.frame.box,
                                      interleaflet_cutoff)

        # Allocate memory for temp arrays
        tmp_thickness = <real *> malloc(size * sizeof(real))
        if tmp_thickness == NULL:
            abort()

        for i in prange(size, schedule="dynamic", num_threads=OPENMP_NUM_THREADS):
            tmp_thickness[i] = thickness_from_neighborhood(&self_coords[i, XX],
                                                           &self_neighborhood_normals[i, XX],
                                                           other_coords,
                                                           self_coords,
                                                           other_neighborhood_normals,
                                                           self_neighborhood_normals,
                                                           holder_other_leaflet.neighborhoods[i],
                                                           self.neighborhoods.neighborhoods[i],
                                                           self.frame.box)


        # Fix not valid thickness using neighbors
        for i in prange(size, schedule="dynamic", num_threads=OPENMP_NUM_THREADS):
        #for i in range(size):
            self.thickness_values[i] = self.fix_thickness(i, tmp_thickness)

        # Free memory
        free_neighborhood_holder(holder_other_leaflet)
        free(tmp_thickness)

        # Compute average thickness
        avg_thickness = 0
        thickness_max = NOTSET
        thickness_min = -NOTSET
        for i in range(size):
            avg_thickness += self.thickness_values[i]
            if self.thickness_values[i] < thickness_min:
                thickness_min = self.thickness_values[i]
            if self.thickness_values[i] > thickness_max:
                thickness_max = self.thickness_values[i]
        avg_thickness /= size

        # Store results
        self.thickness_avg = avg_thickness
        self.thickness_cutoff = interleaflet_cutoff
        self.thickness_min = thickness_min
        self.thickness_max = thickness_max

    # Python API
    property avg_normal:
        def __get__(self):
            return np.array([self.avg_normal[XX], self.avg_normal[YY], self.avg_normal[ZZ]],
                            dtype=np.float64)

    property xcm:
        def __get__(self):
            return np.array([self.xcm[XX], self.xcm[YY], self.xcm[ZZ]],
                            dtype=np.float64)

    property is_planar:
        def __get__(self):
            return self.is_planar

    cdef fsl_int fast_size(self) nogil:
        return  self.beadids.shape[0]

    property beadids:
        def __get__(self):
            return np.asarray(self.beadids)

    property resids:
        def __get__(self):
            return np.asarray(self.resids)

    property hg_atomids:
        def __get__(self):
            return np.asarray(self.hg_atomids)

    property lipid_atomids:
        def __get__(self):
            return np.asarray(self.lipid_atomids)

    property coords:
        def __get__(self):
            return np.asarray(self.coords)

    property normals:
        def __get__(self):
            return np.asarray(self.neighborhood_normals)

    property directions:
        def __get__(self):
            return np.asarray(self.directions)

    property raw_normals:
        def __get__(self):
            return np.asarray(self.normals)

    def __len__(self):
        return self.fast_size()

    def __repr__(self):
        if self.is_planar:
            str_type = "Planar"
        else:
            str_type = "Non planar"
        return "%s aggregate made of %i lipids (XCM: %.3f, %.3f, %.3f)" % \
               (str_type, self.fast_size(),
                self.xcm[XX], self.xcm[YY], self.xcm[ZZ])



cdef enum MEMBRANE_TYPE:
    MT_UNKNOWN, MT_PLANAR, MT_VESICLE

cdef class Membrane(object):
    def __init__(self, Frame frame not None, Aggregate l1 not None, Aggregate l2):
        self.frame = frame

        if l1.is_planar:
            membrane_type = MT_PLANAR
        else:
            membrane_type = MT_VESICLE


        if membrane_type == MT_VESICLE:
            # Make sure that leaflet 1 always corresponds to the outer leaflet (i.e. bigger)
            if l1.beadids.shape[0] >= l2.beadids.shape[0]:
                self.leaflet1 = l1
                self.leaflet2 = l2
            else:
                self.leaflet1 = l2
                self.leaflet2 = l1
        else:
            # Make sure that leaflet 1 always corresponds to the lower leaflet

            ref_dir = XX
            ref_val = real_abs(l1.avg_normal[ref_dir])
            for direction in [YY, ZZ]:
                val = real_abs(l1.avg_normal[direction])
                if val > ref_val:
                    ref_dir = direction
                    ref_val = val

            if l1.avg_normal[ref_dir] < 0:
                self.leaflet1 = l1
                self.leaflet2 = l2
            else:
                self.leaflet1 = l2
                self.leaflet2 = l1

        self.type = membrane_type

        # Thickness related
        self.thickness = NOTSET

        # APL related
        self.apl = NOTSET


    cdef void fast_compute_thickness(self, real interleaflet_cutoff=THICKNESS_DEFAULT_CUTOFF, bint force=False) nogil except*:
        cdef fsl_int l1_size = self.leaflet1.fast_size()
        cdef fsl_int l2_size = self.leaflet2.fast_size()

        # Compute the thickness for the leaflet 1
        self.leaflet1.compute_thickness(self.leaflet2, interleaflet_cutoff, force)

        # Compute the thickness for the leaflet 2
        self.leaflet2.compute_thickness(self.leaflet1, interleaflet_cutoff, force)

        # Store results
        self.thickness = (self.leaflet1.thickness_avg * l1_size + self.leaflet2.thickness_avg * l2_size) / (l1_size + l2_size)


    cdef void fast_compute_apl(self, real cutoff=APL_DEFAULT_CUTOFF, real area_limit=APL_DEFAULT_AREA_LIMIT, bint force=False) nogil except*:
        cdef fsl_int l1_size = self.leaflet1.fast_size()
        cdef fsl_int l2_size = self.leaflet2.fast_size()

        # Compute the APL for the leaflet 1
        self.leaflet1.compute_apl(cutoff, area_limit, force)

        # Compute the APL for the leaflet 2
        self.leaflet2.compute_apl(cutoff, area_limit, force)

        # Compute average APL and store results
        self.apl = (self.leaflet1.apl_avg * l1_size + self.leaflet2.apl_avg * l2_size) / (l1_size + l2_size)

    # Python API
    property beadids:
        def __get__(self):
            return np.asarray(self.leaflet1.beadids), np.asarray(self.leaflet2.beadids)

    property is_planar:
        def __get__(self):
            return self.type == MT_PLANAR

    property is_vesicle:
        def __get__(self):
            return self.type == MT_VESICLE

    def get_thickness(self,
                      real interleaflet_cutoff=THICKNESS_DEFAULT_CUTOFF,
                      bint only_average=True,
                      bint force=False):

        # First compute thickness
        with nogil:
            self.fast_compute_thickness(interleaflet_cutoff, force)

        if only_average:
            l1_val = (self.leaflet1.thickness_avg,
                      self.leaflet1.thickness_min,
                      self.leaflet1.thickness_max)
            l2_val = (self.leaflet2.thickness_avg,
                      self.leaflet2.thickness_min,
                      self.leaflet2.thickness_max)
        else:
            l1_val = np.asarray(self.leaflet1.thickness_values)
            l2_val = np.asarray(self.leaflet2.thickness_values)

        return self.thickness, l1_val, l2_val

    def get_apl(self,
                real cutoff=APL_DEFAULT_CUTOFF,
                real area_limit=APL_DEFAULT_AREA_LIMIT,
                bint by_type=True,
                bint only_average=True,
                bint force=False):
        cdef fsl_int i, tid
        cdef real apl_val
        cdef bytes resname

        # First compute APL
        with nogil:
            self.fast_compute_apl(cutoff, area_limit, force)

        # Then sort APL by restype
        if by_type:
            l1_apl = {}
            for tid, key in enumerate(self.leaflet1.lipid_types):
                if only_average:
                    l1_apl[key] = (self.leaflet1.nlipids_by_type[key],
                                   fast_average(self.leaflet1.apl_by_types[tid]),
                                   self.leaflet1.apl_by_types[tid].min(),
                                   self.leaflet1.apl_by_types[tid].max(),
                                   self.leaflet1.area_by_types[tid])
                else:
                    l1_apl[key] = self.leaflet1.apl_by_types[tid].copy()

            l2_apl = {}
            for tid, key in enumerate(self.leaflet2.lipid_types):
                if only_average:
                    l2_apl[key] = (self.leaflet2.nlipids_by_type[key],
                                   fast_average(self.leaflet2.apl_by_types[tid]),
                                   self.leaflet2.apl_by_types[tid].min(),
                                   self.leaflet2.apl_by_types[tid].max(),
                                   self.leaflet2.area_by_types[tid])
                else:
                    l2_apl[key] = self.leaflet2.apl_by_types[tid].copy()


        else:
            if only_average:
                l1_apl = (self.leaflet1.apl_avg,
                          self.leaflet1.apl_min,
                          self.leaflet1.apl_max,
                          self.leaflet1.area)
            else:
                l1_apl = np.asarray(self.leaflet1.apl_values)

            if only_average:
                l2_apl = (self.leaflet2.apl_avg,
                          self.leaflet2.apl_min,
                          self.leaflet2.apl_max,
                          self.leaflet2.area)
            else:
                l2_apl = np.asarray(self.leaflet2.apl_values)

        return self.apl, l1_apl, l2_apl

    def __getitem__(self, item):
        self_as_list = [self.leaflet1, self.leaflet2]
        return self_as_list[item]

    def __len__(self):
        cdef fsl_int l1_size = self.leaflet1.fast_size()
        cdef fsl_int l2_size = self.leaflet2.fast_size()

        return l1_size + l2_size


########################################################################################################################
#
# Analysis Stuff
#
########################################################################################################################
cdef list retrieve_aggregates(Frame frame, real cutoff):
    if cutoff < 0:
        raise ValueError("Cutoff MUST be > 0!")

    cdef fsl_int *stack
    cdef fsl_int stack_size
    cdef real[:, ::1] coords = frame.fast_get_bead_coords_bbox()
    cdef real[:, ::1] lipid_coords = frame.fast_get_lipid_coords_bbox()
    cdef real[:, ::1] directions = frame.fast_get_directions()
    cdef fsl_int[:] lipid_offsets = frame.lipid_atomids_offsets
    cdef fsl_int[:] processed_state
    cdef real[:, ::1] normals
    cdef fsl_int ncoords = coords.shape[0]
    cdef fsl_int[:] aggregates_id
    cdef fsl_int[:] cur_aggregate
    cdef fsl_int[:] rejected
    cdef fsl_int n_rejected
    cdef PBCBox box = frame.box
    cdef ns_neighborhood *neighborhood
    cdef fsl_int i, j, nid, ref_nid
    cdef fsl_int cur_aggregate_id, cur_leaflet_id, aggregate_size, aggregate_size_sum = 0

    cdef Aggregate ref_aggregate, aggregate, merged_aggregate

    # Aggregates buffer
    cdef fsl_int aggid, best_aggid, ref_aggid
    cdef fsl_int[:, ::1] aggregates_buffer
    cdef fsl_int[:] aggregate_sizes

    # Further tests
    cdef fsl_int[:] needs_check
    cdef fsl_int absorber, absorbed
    cdef fsl_int n_needs_check
    cdef fsl_int tmp_aggregate_size
    cdef fsl_int best_nid
    cdef real test_dprod, max_dprod, test_orientation, test_insider
    cdef real closest_dist, test_dist
    cdef fsl_int closest_nid, closest_aggid
    cdef real[:] xcm
    cdef real[:] local_normal
    cdef rvec xcm_rvec, dx
    cdef fsl_int bid


    # Allocate memory
    ncoords = coords.shape[0]
    stack = <fsl_int *> PyMem_Malloc(ncoords * sizeof(fsl_int))
    stack_size = 0
    aggregates_id = np.empty(ncoords, dtype=np.int64)
    rejected = np.empty(ncoords, dtype=np.int64)
    n_rejected = 0
    processed_state = np.zeros(ncoords, dtype=np.int64)
    cur_aggregate = np.empty(ncoords, dtype=np.int64)


    # Initialization
    for i in range(ncoords):
        aggregates_id[i] = NOTSET

    aggregates_checked = {}
    aggregates_tocheck = {}

    #print("Cutoff used for aggregation: %.2f" % cutoff)
    #print("Box: %s" % box)
    with nogil:
        #fprintf(stderr, "Cutoff used for aggregation: %.2f\n", cutoff)
        normals = frame.fast_get_normals(cutoff)
        cur_aggregate_id = -1
        cur_leaflet_id = -1

        # Step 1: Loop over the coordinates to find the "safe" aggregates
        for i in range(ncoords):

            # Only do something if the current bead has not been already accepted
            if processed_state[i] == AGGREGATE_PROCESSING_DONE:
                continue

            # Build a new aggregate
            cur_aggregate_id += 1
            aggregate_size = 0

            # Reset the stack to current bead id
            stack_size = 1
            stack[stack_size - 1] = i
            processed_state[i] = AGGREGATE_PROCESSING_DONE

            # Add the first ref bead to current leaflet
            aggregates_id[i] = cur_aggregate_id
            cur_aggregate[aggregate_size] = i
            aggregate_size += 1

            # Parse the neighbors' neighborhood to build the aggregate
            while stack_size > 0:

                # Pop stack
                ref_nid = stack[stack_size - 1]
                stack_size -= 1

                # Retrieve neighborhood
                neighborhood = frame.neighbors.neighborhoods[ref_nid]

                for j in range(neighborhood.size):
                    nid = neighborhood.beadids[j]

                    # Check if the bead already processed
                    if processed_state[nid] != 0:
                        continue

                    # Check if the normals are oriented towards the same direction
                    if rvec_dprod(&normals[ref_nid, XX], &normals[nid, XX]) < COS_45:
                        processed_state[nid] = AGGREGATE_PROCESSING_REJECTED
                        rejected[n_rejected] = nid
                        n_rejected += 1
                        continue

                    # If still here, add bead to current aggregate
                    aggregates_id[nid] = cur_aggregate_id
                    cur_aggregate[aggregate_size] = nid
                    aggregate_size += 1
                    processed_state[nid] = AGGREGATE_PROCESSING_DONE

                    # Append bead to stack
                    stack_size += 1
                    stack[stack_size - 1] = nid

            # Store aggregate
            with gil:
                if aggregate_size > MIN_LEAFLET_SIZE:
                    aggregates_checked[cur_aggregate_id] = cur_aggregate[:aggregate_size].copy()
                else:
                    #print("Small aggregate found: %s - ID: %i\n" % (np.asarray(cur_aggregate[:aggregate_size]).tolist(), cur_aggregate_id))
                    aggregates_tocheck[cur_aggregate_id] = cur_aggregate[:aggregate_size].copy()

            # Reset rejected
            for j in range(n_rejected):
                processed_state[rejected[j]] = AGGREGATE_PROCESSING_RESET
            n_rejected = 0

    # Step 2: Try to merge small aggregate to bigger ones

    changed = False
    for ref_aggid, beadids in aggregates_tocheck.items():
        #print("Inspecting agg#%i: %s" % (ref_aggid, np.asarray(beadids)))
        if len(beadids) == 0:
            continue

        max_dprod = -1.1
        best_aggid = NOTSET
        best_nid = NOTSET
        same_aggid = False
        #print("Stuff reset -> beadids: %s\n" % np.asarray(beadids).tolist())

        for ref_nid in beadids:
            # Retrieve neighborhood
            neighborhood = frame.neighbors.neighborhoods[ref_nid]

            # Compute the avg distances to big aggregates (if possible!)
            d2s = {}
            local_normals = {}
            local_xcms = {}
            ns = {}
            for j in range(neighborhood.size):
                nid = neighborhood.beadids[j]
                aggid = aggregates_id[nid]

                if aggid in aggregates_checked:
                    if len(aggregates_checked[aggid]) > MIN_LEAFLET_SIZE:
                        try:
                            d2s[aggid] += box.fast_distance2(&coords[ref_nid, XX], &coords[nid, XX])
                            local_normals[aggid] += np.asarray(normals[nid])
                            local_xcms[aggid].append(np.asarray(coords[nid]))
                            ns[aggid] += 1
                        except KeyError:
                            d2s[aggid] = box.fast_distance2(&coords[ref_nid, XX], &coords[nid, XX])
                            local_normals[aggid] = np.asarray(normals[nid])
                            local_xcms[aggid] = [np.asarray(coords[nid])]
                            ns[aggid] = 1

            # Handles the unlikely case where we could find a neighbor inside a checked aggregate
            if len(ns) == 0: # No choice, we create a new checked aggregate
                #fprintf(stderr, "NO active neighbor for bead #%i\n", ref_nid)
                try:
                    aggregates_checked[ref_aggid] = np.append(aggregates_checked[ref_aggid], np.array([ref_nid], dtype=np.int64), axis=0)
                except KeyError:
                    aggregates_checked[ref_aggid] = np.array([ref_nid], dtype=np.int64)
                continue

            # Get avg
            for aggid, n in ns.items():
                d2s[aggid] /= n
                local_normals[aggid] /= n
                box.fast_pbc_xcm(np.array(local_xcms[aggid], dtype=np.float64), xcm_rvec)
                local_xcms[aggid] = np.array([xcm_rvec[XX], xcm_rvec[YY], xcm_rvec[ZZ]],
                                             dtype=np.float64)

            # Get best local normal
            closest_aggid = NOTSET
            closest_dist = 1e10
            for aggid, d2 in d2s.items():
                if d2 < closest_dist:
                    closest_dist = d2
                    closest_aggid = aggid


            local_normal = local_normals[closest_aggid]
            xcm = local_xcms[closest_aggid]


            box.fast_pbc_dx(&coords[ref_nid, XX], &xcm[XX], dx)

            test_orientation = rvec_dprod(&local_normal[XX], &directions[ref_nid, XX])
            test_insider = rvec_dprod(&local_normal[XX], dx)

            # Check if the bead is located on the internal side of the potential leaflet
            if test_orientation > 0 or test_insider > 0: # That's the case
                aggregates_checked[closest_aggid] = np.append(aggregates_checked[closest_aggid], np.array([ref_nid], dtype=np.int64), axis=0)
                aggregates_id[ref_nid] = closest_aggid
            else:
                try:
                    aggregates_checked[ref_aggid] = np.append(aggregates_checked[ref_aggid], np.array([ref_nid], dtype=np.int64), axis=0)
                except KeyError:
                    aggregates_checked[ref_aggid] = np.array([ref_nid], dtype=np.int64)


    #fprintf(stderr, "%i small agg left\n", len(aggregates_tocheck))
    # Finally, build aggregates
    aggregates_list = list()

    for beadids in aggregates_checked.values():
        aggregates_list.append(Aggregate(frame, beadids))

    # Sort aggregates according to their populations
    aggregates_list.sort(key = len, reverse=True)

    # Free Memory
    PyMem_Free(<void *> stack)

    return aggregates_list


cdef list retrieve_membranes(Frame frame, real cutoff):
    cdef Aggregate ref_leaflet, leaflet
    cdef fsl_int membrane_type = MT_UNKNOWN

    # Sort aggregates according to their populations
    aggregates_py = retrieve_aggregates(frame, cutoff)

    # First, clean the aggregates
    clean_aggregates_py = list()
    for ref_leaflet in aggregates_py:
        if ref_leaflet.fast_size() > MIN_LEAFLET_SIZE:
            clean_aggregates_py.append(ref_leaflet)

    # Second, build the membranes
    membranes_py = list()

    while len(clean_aggregates_py) > 0:
        ref_leaflet = clean_aggregates_py.pop(0)

        compatibles = []
        for leaflet in clean_aggregates_py:
            if check_leaflet_compatibility(frame, ref_leaflet, leaflet):
                compatibles.append(leaflet)

        companion = None
        ref_dist = -1
        for leaflet in compatibles:
            if companion is None:
                companion = leaflet
                ref_dist = frame.box.fast_leaflet_distance(ref_leaflet.xcm,
                                                           leaflet.xcm,
                                                           ref_leaflet.avg_normal)
                continue

            cur_dist = frame.box.fast_leaflet_distance(ref_leaflet.xcm,
                                                       leaflet.xcm,
                                                       ref_leaflet.avg_normal)

            if cur_dist < ref_dist:
                companion = leaflet
                ref_dist = cur_dist

        if companion is not None:
            membranes_py.append(Membrane(frame, ref_leaflet, companion))
            clean_aggregates_py.remove(companion)


    return membranes_py
