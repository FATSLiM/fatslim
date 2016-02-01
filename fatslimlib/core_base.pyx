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
#cython: wraparound=False
#cython: profile=False
from __future__ import print_function

# Preprocessor constants

DEF NOTSET = -12345
DEF DIM = 3
DEF XX = 0
DEF YY = 1
DEF ZZ = 2
DEF RET_OK = 1
DEF RET_YES = 1
DEF RET_ERROR = 0
DEF RET_NO = 0
DEF NAME_SIZE = 6 # 5 chracters + NULL

DEF NAME_ALLOCATION_INCREMENT = 10
DEF RESIDUE_ALLOCATION_INCREMENT = 20
DEF RESIDUES_ALLOCATION_INCREMENT = 50
DEF ATOMS_ALLOCATION_INCREMENT = 1000

DEF EPSILON = 1e-6

DEF MAX_NTRICVEC=12
DEF BOX_MARGIN=1.0010

DEF DEFAULT_PROXIMITY_CUTOFF = 2.0

DEF FFT_ALLATOM_THRESHOLD = 0.015
DEF FFT_UNIFIED_THRESHOLD = 0.025
DEF FFT_COARSE_THRESHOLD = 0.25

#from libc.stdio cimport fprintf, stderr
from libc.math cimport sqrt, log10, floor, fabs
from libc.string cimport strncmp, strncpy
from libc.stdlib cimport malloc, free, abort
from cpython.mem cimport PyMem_Realloc, PyMem_Free
cimport cython
from cython.parallel cimport prange

from typedefs cimport real, fsl_int, real_min, real_abs
from typedefs cimport rvec, rvec_copy, rvec_inc, rvec_clear, rvec_smul, rvec_dprod, \
    rvec_normalize, rvec_dec, rvec_norm2
from typedefs cimport matrix, mat_clear

cdef extern from "eig_mat.h" nogil:
    void eigen_33_sym(matrix a, matrix eig_vec, rvec eig_val)

cdef extern from "openmp_wrapper.h" nogil:
    cdef fsl_int omp_get_max_threads()
    cdef fsl_int omp_get_num_procs()
    cdef fsl_int omp_get_thread_num()

from .core_ns cimport ns_neighborhood, fast_neighbor_search, free_neighborhood_holder
from .core_analysis cimport Aggregate, Membrane
cimport core_analysis

import numpy as np
from time import time
import os
import sys


# Constants
__cython_version__ = cython.__version__

# OpenMP-related stuff
cdef fsl_int OPENMP_MAX_THREADS = omp_get_max_threads()
cdef fsl_int OPENMP_NUM_THREADS = OPENMP_MAX_THREADS

cpdef fsl_int get_num_threads():
    return OPENMP_NUM_THREADS

cpdef fsl_int get_max_threads():
    return OPENMP_MAX_THREADS

cpdef fsl_int get_num_procs():
    return omp_get_num_procs()

cpdef set_num_threads(fsl_int num):
    global OPENMP_NUM_THREADS
    if not 0 < num <= OPENMP_MAX_THREADS:
        num = OPENMP_MAX_THREADS

    OPENMP_NUM_THREADS = num

    if OPENMP_NUM_THREADS != num:
        raise RuntimeError("Could not change number of OpenMP threads "
                           "(requested value: %i - actual value: %i)" % (num, OPENMP_NUM_THREADS))



cdef real norm2(real[:] a) nogil:
    return a[XX] * a[XX] + a[YY] * a[YY] + a[ZZ] * a[ZZ]

cdef void mat_to_memview(matrix src, real[:, ::1] dest) nogil:
    dest[XX, XX] = src[XX][XX]
    dest[XX, YY] = src[XX][YY]
    dest[XX, ZZ] = src[XX][ZZ]
    dest[YY, XX] = src[YY][XX]
    dest[YY, YY] = src[YY][YY]
    dest[YY, ZZ] = src[YY][ZZ]
    dest[ZZ, XX] = src[ZZ][XX]
    dest[ZZ, YY] = src[ZZ][YY]
    dest[ZZ, ZZ] = src[ZZ][ZZ]

cdef void rvec_to_memview(rvec src, real[:] dest) nogil:
    dest[XX] = src[XX]
    dest[YY] = src[YY]
    dest[ZZ] = src[ZZ]

cpdef pretty_delta(begin, end):
    delta = end - begin

    return pretty_duration(delta)

cpdef pretty_duration(duration):
    if duration != 0:
        expo = int(floor(log10(duration)))
        unit = {-1 : (3, "ms"), -2:(6, "us")}
        ex, unit_str = unit.setdefault(expo // 3, (0, "s"))
    else:
        ex = 0.0
        unit_str = "us"
    return u"%.0f %s" % (duration * 10 ** ex, unit_str)


cpdef real get_memory_usage():
    memory = -1
    with open('/proc/self/status', "rb") as fp:
        for line in fp:
            if line.startswith("VmRSS:"):
                memory = int(line.split()[1])
                break
    return memory

def verbose_print(msg, bint verbose, end="\n"):
    if verbose:
        print(msg, end=end)
        sys.stdout.flush()

cdef class Atom(object):

    def __init__(self,
                  fsl_int atomid=-1,
                  str name="unk",
                  real x=0,
                  real y=0,
                  real z=0):
        self.c_atom.atomid = atomid
        self.name = name.encode() # Used to keep at least one ref of the original Python string
        self.c_atom.name = self.name
        self.c_atom.coords[XX] = x
        self.c_atom.coords[YY] = y
        self.c_atom.coords[ZZ] = z

    # Python API
    property x:
        def __get__(self):
            return self.c_atom.coords[XX]

    property y:
        def __get__(self):
            return self.c_atom.coords[YY]

    property z:
        def __get__(self):
            return self.c_atom.coords[ZZ]

    property atomid:
        def __get__(self):
            return self.c_atom.atomid

    property name:
        def __get__(self):
            return str(self.name.decode())


cdef void calculate_lipid_direction(rvec ref_bead,
                                    cPBCBox_t box,
                                    real[:,::1] lipid_coords,
                                    fsl_int offset_start,
                                    fsl_int offset_end,
                                    rvec direction) nogil:
    cdef rvec xcm

    # Get lipid coords
    #lipid_coords = all_lipid_coords[offset_start:offset_end]

    # Get lipid center of mass
    fast_pbc_xcm_from_ref_c(box, lipid_coords, ref_bead, xcm,
                            offset_start, offset_end)

    # Compute direction
    fast_pbc_dx_c(box, xcm, ref_bead, direction)

    # Normalize it:
    rvec_normalize(direction)


cdef void calculate_neighborhood_normal(fsl_int refid, real[:, ::1] coords,
                                        real[:,::1] directions,
                                        ns_neighborhood *neighborhood,
                                        PBCBox box,
                                        rvec normal,
                                        fsl_int *useful_neighbors) nogil:
    cdef fsl_int i, nid
    cdef fsl_int neighborhood_size = neighborhood.size, useful_size = 0
    cdef rvec dx, xcm
    cdef rvec neighborhood_direction
    cdef rvec ref_direction
    cdef matrix cov_mat
    cdef rvec tmp_vec
    cdef rvec eig_vals
    cdef matrix eig_vecs

    # Don't compute anything if less than 3 neighbors
    # Assume the normal as the reverse of the lipid directions
    if neighborhood_size < 3:
        normal[XX] = directions[refid, XX]
        normal[YY] = directions[refid, YY]
        normal[ZZ] = directions[refid, ZZ]
        return

    # Compute center of masses
    rvec_clear(xcm)
    rvec_clear(neighborhood_direction)

    rvec_copy(&directions[refid, XX], ref_direction)

    for i in range(neighborhood_size):
        nid = neighborhood.beadids[i]

        # Only take into account the neighbors that are oriented toward the same direction (< 90°)
        if rvec_dprod(ref_direction, &directions[nid, XX]) > 0:
            useful_neighbors[useful_size] = nid
            useful_size += 1

            # HG xcm
            box.fast_pbc_dx(&coords[refid, XX], &coords[nid, XX], dx)
            rvec_inc(xcm, dx)

            rvec_inc(neighborhood_direction, &directions[nid, XX])

    # Don't compute anything if less than 3 useful neighbors
    # Instead, set it to notset as the flag for further computations
    if useful_size < 3:
        normal[XX] = ref_direction[XX]
        normal[YY] = ref_direction[YY]
        normal[ZZ] = ref_direction[ZZ]
        return

    rvec_smul(1.0/useful_size, xcm, xcm)
    rvec_inc(xcm, &coords[refid, XX])
    rvec_smul(1.0/useful_size, neighborhood_direction, neighborhood_direction)

    # Build covariance matrix
    mat_clear(cov_mat)
    for i in range(useful_size):
        nid = useful_neighbors[i]

        # Retrieve neighbor and its image
        box.fast_pbc_dx(xcm, &coords[nid, XX], tmp_vec)

        cov_mat[YY][YY] += tmp_vec[YY] * tmp_vec[YY]
        cov_mat[YY][ZZ] += tmp_vec[YY] * tmp_vec[ZZ]
        cov_mat[ZZ][ZZ] += tmp_vec[ZZ] * tmp_vec[ZZ]

        rvec_smul(tmp_vec[XX], tmp_vec, tmp_vec)

        cov_mat[XX][XX] += tmp_vec[XX]
        cov_mat[XX][YY] += tmp_vec[YY]
        cov_mat[XX][ZZ] += tmp_vec[ZZ]
    cov_mat[YY][XX] = cov_mat[XX][YY]
    cov_mat[ZZ][XX] = cov_mat[XX][ZZ]
    cov_mat[ZZ][YY] = cov_mat[YY][ZZ]

    cov_mat[XX][XX] /= useful_size
    cov_mat[XX][YY] /= useful_size
    cov_mat[XX][ZZ] /= useful_size
    cov_mat[YY][XX] /= useful_size
    cov_mat[YY][YY] /= useful_size
    cov_mat[YY][ZZ] /= useful_size
    cov_mat[ZZ][XX] /= useful_size
    cov_mat[ZZ][YY] /= useful_size
    cov_mat[ZZ][ZZ] /= useful_size

    # Get eigenvalues
    rvec_clear(eig_vals)
    mat_clear(eig_vecs)
    eigen_33_sym(cov_mat, eig_vecs, eig_vals)


    # Retrieve the normal from the chosen eigen vector
    rvec_clear(normal)
    if rvec_dprod(eig_vecs[0], &neighborhood_direction[XX]) < 0:
        rvec_dec(normal, eig_vecs[0])
    else:
        rvec_inc(normal, eig_vecs[0])
    rvec_normalize(normal)

cdef class AtomGroup(object):

    # initialization
    def __cinit__(self,
                  str name="group",
                  fsl_int groupid=0,
                  object atoms=()):
        # Store name
        self.name = name.encode()
        self.c_atomgroup.name = self.name

        # Store group id
        self.c_atomgroup.groupid = groupid

        # Check and store atoms
        self.c_atomgroup.n_atoms = 0
        self.c_atomgroup.atoms = NULL
        for atom in atoms:
            self.append(<Atom> atom)

    def append(self, Atom atom):
        if not isinstance(atom, Atom):
            raise TypeError("%s is not a %s instance" % (atom, Atom))
        self.fast_append(atom)

    cdef fast_append(self, Atom atom):
        self.c_atomgroup.n_atoms += 1
        self.c_atomgroup.atoms = <cAtom_t *> PyMem_Realloc(<void *> self.c_atomgroup.atoms,
                                             self.c_atomgroup.n_atoms * sizeof(cAtom_t))
        if self.c_atomgroup.atoms == NULL:
            raise MemoryError
        self.c_atomgroup.atoms[self.c_atomgroup.n_atoms - 1] = atom.c_atom


    def __len__(self):
        return self.c_atomgroup.n_atoms

    # C destructor
    def __dealloc__(self):
        PyMem_Free(self.c_atomgroup.atoms)


cdef void fast_pbc_dx_c(cPBCBox_t pbcbox, rvec ref, rvec other, rvec dx) nogil:
        cdef fsl_int i, j
        cdef rvec dx_start, trial

        for i in range(DIM):
            dx[i] = other[i] - ref[i]


        for i in range (DIM-1, -1, -1):
            while dx[i] > pbcbox.hbox_diag[i]:
                for j in range (i, -1, -1):
                    dx[j] -= pbcbox.box[i][j]

            while dx[i] <= pbcbox.mhbox_diag[i]:
                for j in range (i, -1, -1):
                    dx[j] += pbcbox.box[i][j]


cdef void fast_pbc_xcm_from_ref_c(cPBCBox_t pbcbox, real[:, ::1] coords, rvec ref, rvec xcm,
                                  fsl_int first=0, fsl_int last=NOTSET) nogil:
        cdef fsl_int i
        cdef fsl_int size
        cdef rvec dx
        cdef real[:,::1] bbox_coords

        if first < 0:
            first = 0

        if last == NOTSET:
            last = coords.shape[0]

        size = last - first

        if size < 1:
            return

        rvec_copy(ref, xcm)
        rvec_smul(size, xcm, xcm)

        # Step 1: Get XCM
        for i in range(first, last):
            fast_pbc_dx_c(pbcbox, ref, &coords[i, XX], dx)
            rvec_inc(xcm, dx)
        rvec_smul(1.0/size, xcm, xcm)

        # Step 2: Make sure it is inside the brick-shaped box
        for i in range(DIM - 1, -1, -1):
            while xcm[i] < 0:
                for d in range(i+1):
                    xcm[d] += pbcbox.box[i][d]
            while xcm[i] >= pbcbox.box[i][i]:
                for d in range(i+1):
                    xcm[d] -= pbcbox.box[i][d]

# noinspection PyNoneFunctionAssignment
cdef class PBCBox(object):

    def __init__(self, real[:,::1] box):
        self.update(box)

    cdef void fast_update(self, real[:,::1] box) nogil:
        cdef fsl_int i, j, k, d, jc, kc, shift
        cdef real d2old, d2new, d2new_c
        cdef rvec trial, pos
        cdef fsl_int ii, jj ,kk
        cdef fsl_int *order = [0, -1, 1, -2, 2]
        cdef bint use
        cdef real min_hv2, min_ss

        rvec_clear(self.center)
        # Update matrix
        for i in range(DIM):
            for j in range(DIM):
                self.c_pbcbox.box[i][j] = box[i, j]
                self.center[j] += 0.5 * box[i, j]
            self.bbox_center[i] = 0.5 *  box[i, i]

        # Update diagonals
        for i in range(DIM):
            self.c_pbcbox.fbox_diag[i]  =  box[i, i]
            self.c_pbcbox.hbox_diag[i]  =  self.c_pbcbox.fbox_diag[i] * 0.5
            self.c_pbcbox.mhbox_diag[i] = - self.c_pbcbox.hbox_diag[i]

        # Update maximum cutoff

        # Physical limitation of the cut-off
        # by half the length of the shortest box vector.
        min_hv2 = real_min(0.25 * rvec_norm2(&box[XX, XX]), 0.25 * rvec_norm2(&box[YY, XX]))
        min_hv2 = real_min(min_hv2, 0.25 * rvec_norm2(&box[ZZ, XX]))

        # Limitation to the smallest diagonal element due to optimizations:
        # checking only linear combinations of single box-vectors (2 in x)
        # in the grid search and pbc_dx is a lot faster
        # than checking all possible combinations.
        min_ss = real_min(box[XX, XX], real_min(box[YY, YY]- real_abs(box[ZZ, YY]), box[ZZ, ZZ]))

        self.c_pbcbox.max_cutoff2 = real_min(min_hv2, min_ss * min_ss)

        # Update shift vectors
        self.c_pbcbox.ntric_vec = 0
        # We will only use single shifts, but we will check a few
        # more shifts to see if there is a limiting distance
        # above which we can not be sure of the correct distance.
        for kk in range(5):
            k = order[kk]

            for jj in range(5):
                j = order[jj]

                for ii in range(5):
                    i = order[ii]

                    # A shift is only useful when it is trilinic
                    if j != 0 or k != 0:
                        d2old = 0
                        d2new = 0

                        for d in range(DIM):
                            trial[d] = i*box[XX, d] + j*box[YY, d] + k*box[ZZ, d]

                            # Choose the vector within the brick around 0,0,0 that
                            # will become the shortest due to shift try.

                            if d == DIM:
                                trial[d] = 0
                                pos[d] = 0
                            else:
                                if trial[d] < 0:
                                    pos[d] = real_min(self.c_pbcbox.hbox_diag[d], -trial[d])
                                else:
                                    pos[d] = max(-self.c_pbcbox.hbox_diag[d], -trial[d])

                            d2old += sqrt(pos[d])
                            d2new += sqrt(pos[d] + trial[d])

                        if BOX_MARGIN*d2new < d2old:
                            if  not (j < -1 or j > 1 or k < -1 or k > 1):
                                use = True

                                for dd in range(DIM):
                                    if dd == 0:
                                        shift = i
                                    elif dd == 1:
                                        shift = j
                                    else:
                                        shift = k

                                    if shift:
                                        d2new_c = 0

                                        for d in range(DIM):
                                            d2new_c += sqrt(pos[d] + trial[d] - shift*box[dd, d])

                                        if d2new_c <= BOX_MARGIN*d2new:
                                            use = False

                                if use: # Accept this shift vector.
                                    if self.c_pbcbox.ntric_vec >= MAX_NTRICVEC:
                                        with gil:
                                            print("\nWARNING: Found more than %d triclinic "
                                                  "correction vectors, ignoring some."
                                                  % MAX_NTRICVEC)
                                            print("  There is probably something wrong with "
                                                  "your box.")
                                            print(box)
                                    else:
                                        for d in range(DIM):
                                            self.c_pbcbox.tric_vec[self.c_pbcbox.ntric_vec][d] = \
                                                trial[d]
                                        self.c_pbcbox.tric_shift[self.c_pbcbox.ntric_vec][XX] = i
                                        self.c_pbcbox.tric_shift[self.c_pbcbox.ntric_vec][YY] = j
                                        self.c_pbcbox.tric_shift[self.c_pbcbox.ntric_vec][ZZ] = k
                                        self.c_pbcbox.ntric_vec += 1


    def update(self, real[:,::1] box):
        if box.shape[0] != DIM or box.shape[1] != DIM:
            raise ValueError("Box must be a %i x %i matrix. (shape: %i x %i)" %
                             (DIM, DIM, box.shape[0], box.shape[1]))
        if (box[XX, XX] == 0) or (box[YY, YY] == 0) or (box[ZZ, ZZ] == 0):
            raise ValueError("Box does not correspond to PBC=xyz")
        self.fast_update(box)

    property box:
        def __get__(self):
            box_py = np.empty((DIM, DIM))
            cdef real[:,::1] box_memview = box_py

            mat_to_memview(self.c_pbcbox.box, box_memview)

            return box_py

    property center:
        def __get__(self):
            return np.array([self.center[XX], self.center[YY], self.center[ZZ]], dtype=np.float64)

    property bbox_center:
        def __get__(self):
            return np.array([self.bbox_center[XX], self.bbox_center[YY], self.bbox_center[ZZ]],
                            dtype=np.float64)

    def to_gro(self):
        ret_str = "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f" % (
            self.c_pbcbox.box[XX][XX], self.c_pbcbox.box[YY][YY], self.c_pbcbox.box[ZZ][ZZ],
            self.c_pbcbox.box[XX][YY], self.c_pbcbox.box[XX][ZZ],
            self.c_pbcbox.box[YY][XX], self.c_pbcbox.box[YY][ZZ],
            self.c_pbcbox.box[ZZ][XX], self.c_pbcbox.box[ZZ][YY],
        )

        return ret_str

    def __str__(self):
        return str(getattr(self, "box"))


    cdef void fast_pbc_dx(self, rvec ref, rvec other, rvec dx) nogil:
        fast_pbc_dx_c(self.c_pbcbox, ref, other, dx)

    cdef void fast_pbc_dx_leaflet(self, rvec ref, rvec other, rvec dx, rvec ref_normal) nogil:
        cdef fsl_int i, j, xoffset, yoffset, zoffset, notset=1
        cdef rvec dx_raw, dx_try, offset, offset_x, offset_y, offset_z
        cdef real d2min, d2try
        cdef real coord_val, tmp
        cdef fsl_int shift_axis, sign

        # First: get the actual pbc dx
        self.fast_pbc_dx(ref, other, dx)

        # If dprod(dx, normal) < 0, dx goes through the bilayer
        # because dx = vector from ref to other
        if rvec_dprod(dx, ref_normal) < 0:
            return # Nothing to do!

        # Second: If the pbc dx is not right, we need to find the correct direction to go
        shift_axis = ZZ
        coord_val = fabs(ref_normal[ZZ])
        for i in range(DIM-1):
            tmp = fabs(ref_normal[i])
            if tmp > coord_val:
                shift_axis = i

        # Third: once we found the correct direction to shift,
        # make one shift in that direction to get the right distance
        sign = 1
        if ref_normal[shift_axis] > 0:
            sign = -1
        for i in range (shift_axis, -1, -1):
            dx[i] += sign * self.c_pbcbox.box[shift_axis][i]


    def dx(self, real[:] a_memview, real[:] b_memview):
        cdef rvec dx
        cdef real[:] dx_memview = np.empty(DIM)

        self.fast_pbc_dx(<rvec> &a_memview[0], <rvec> &b_memview[0], dx)
        rvec_to_memview(dx, dx_memview)

        return np.asarray(dx_memview)

    def dx_leaflet(self, real[:] a_memview, real[:] b_memview, real[:] normal_memview):
        cdef rvec dx
        cdef real[:] dx_memview = np.empty(DIM)

        self.fast_pbc_dx_leaflet(<rvec> &a_memview[0], <rvec> &b_memview[0],
                                 dx, <rvec> &normal_memview[0])
        rvec_to_memview(dx, dx_memview)

        return np.asarray(dx_memview)

    cdef real fast_distance2(self, rvec a, rvec b) nogil:
        cdef rvec dx
        self.fast_pbc_dx(a, b, dx)
        return rvec_norm2(dx)

    cdef real fast_distance(self, rvec a, rvec b) nogil:
        return sqrt(self.fast_distance2(a,b))

    def distance2(self, real[:] a_memview, real[:] b_memview):
        return self.fast_distance2(<rvec> &a_memview[0], <rvec> &b_memview[0])

    def distance(self, real[:] a_memview, real[:] b_memview):
        return sqrt(self.distance2(a_memview, b_memview))

    cdef void fast_pbc_xcm(self, real[:, ::1] coords, rvec xcm, fsl_int limit=NOTSET) nogil:
        self.fast_pbc_xcm_from_ref(coords, &coords[0, XX], xcm, limit)

    cdef void fast_pbc_xcm_from_ref(self, real[:, ::1] coords, rvec ref,
                                    rvec xcm, fsl_int limit=NOTSET) nogil:
        fast_pbc_xcm_from_ref_c(self.c_pbcbox, coords, ref, xcm, limit)

    def get_xcm(self, real[:, ::1] coords):
        cdef rvec xcm
        self.fast_pbc_xcm(coords, xcm)

        xcm_py = np.empty(DIM)
        rvec_to_memview(xcm, xcm_py)

        return xcm_py

    def asarray(self):
        cdef fsl_int i, j
        box = np.empty((DIM, DIM))

        for i in range(DIM):
            for j in range(DIM):
                box[i][j] = self.c_pbcbox.box[i][j]

        return box

    cdef real[:, ::1]fast_put_atoms_in_bbox(self, real[:,::1] coords) nogil:
        cdef fsl_int i, m, d, natoms, wd = 0
        cdef real[:,::1] bbox_coords

        natoms = coords.shape[0]
        with gil:
            if natoms == 0:
                bbox_coords = np.empty((0, DIM))
            else:
                bbox_coords = coords.copy()

        for i in range(natoms):
            for m in range(DIM - 1, -1, -1):
                while bbox_coords[i, m] < 0:
                    for d in range(m+1):
                        bbox_coords[i, d] += self.c_pbcbox.box[m][d]
                while bbox_coords[i, m] >= self.c_pbcbox.box[m][m]:
                    for d in range(m+1):
                        bbox_coords[i, d] -= self.c_pbcbox.box[m][d]
        return bbox_coords

    def put_atoms_in_bbox(self, real[:,::1] coords):
        return np.asarray(self.fast_put_atoms_in_bbox(coords))


cdef class TopologyGroup(object):

    def __init__(self, Topology parent, str resname="group", fsl_int resid=1):
        self.parent = parent
        self.resname = resname.encode()
        self.resid = resid
        self.name = self.resname + b"%i" % resid
        self.atoms = NULL
        self.size = 0
        self.allocated_size = 0

    cdef bint fast_contains(self, topol_atom_t atom) nogil except *:
        cdef fsl_int i
        for i in range(self.size):
            if self.atoms[i].atomid == atom.atomid:
                return True
        return False

    cdef void add_atom(self, topol_atom_t atom) nogil except *:
        if self.fast_contains(atom):
            return

        # Allocate memory
        with gil:
            while self.size + 1 > self.allocated_size:
                self.allocated_size += RESIDUE_ALLOCATION_INCREMENT

                self.atoms = \
                    <topol_atom_t *> PyMem_Realloc(self.atoms,
                                                   sizeof(topol_atom_t ) * self.allocated_size)
                if self.atoms == NULL:
                    raise MemoryError

        self.atoms[self.size] = atom
        self.size += 1

    cdef void merge_other_group(self, TopologyGroup other) nogil:
        cdef fsl_int i
        for i in range(other.size):
            self.add_atom(other.atoms[i])

    cdef fsl_int[:] fast_get_atomids(self) nogil:
        cdef fsl_int[:] atomids
        cdef fsl_int i
        cdef topol_atom_t atom

        # Create numpy array
        with gil:
            atomids = np.empty(self.size, dtype=np.int64)

        for i in range(self.size):
            atomids[i] = self.atoms[i].atomid

        return atomids

    cdef object group_by_resid(self):
        cdef fsl_int i, internal_id
        cdef fsl_int atomid
        cdef fsl_int last_resid = -1
        cdef topol_atom_t atom
        cdef TopologyGroup new_group

        groups = []
        for i in range(self.size):
            atom = self.atoms[i]

            if atom.resid != last_resid:
                new_group = TopologyGroup(self.parent)
                groups.append(new_group)
                last_resid = atom.resid
            new_group.add_atom(atom)

        return groups

    def __len__(self):
        return self.size

    def __dealloc__(self):
        PyMem_Free(self.atoms)

    property atomids:
        def __get__(self):
            return np.asarray(self.fast_get_atomids())

    property residue_atomids:
        def __get__(self):
            cdef fsl_int i, internal_id
            cdef fsl_int atomid
            cdef fsl_int last_resid = -1
            cdef topol_atom_t atom

            atomids = []
            for i in range(self.size):
                atom = self.atoms[i]

                if atom.resid != last_resid:
                    atomids.extend(list(self.parent._residues[atom.resid].atomids))
                    last_resid = atom.resid
            return np.array(atomids, dtype=np.int64)

    def __str__(self):
        return "TopologyGroup '%s': %i atoms" % (self.name, self.size)


cdef list topol_residue_atomids_as_list(topol_residue_t *residue):
    cdef fsl_int i

    obj = []
    for i in range(residue.size):
        obj.append(residue.atomids[i])

    return obj

cdef fsl_int topol_residue_atomid_index(topol_residue_t *residue, fsl_int atomid):
    cdef fsl_int i
    for i in range(residue.size):
        if residue.atomids[i] == atomid:
            return i
    return -1

cdef dict topol_residue_to_dict(topol_residue_t *residue, bytes resname):
    atomids = []
    for i in range(residue.size):
        atomids.append(residue.atomids[i])
    atomids.sort()

    resdict = {
        "resname": str(resname.decode()),
        "resid": residue.resid,
        "name": str(resname.decode()) + "%i" % residue.resid,
        "atomids": np.array(atomids, dtype=np.int64),
        "size": residue.size}

    return resdict

cdef void topol_residue_append_atom(topol_residue_t *residue, topol_atom_t *atom):
    if residue == NULL:
        raise RuntimeError

    if residue.size == residue.allocated_size:
        residue.allocated_size += RESIDUE_ALLOCATION_INCREMENT
        residue.atom_internal_ids = \
            <fsl_int *> PyMem_Realloc(<void *> residue.atom_internal_ids,
                                  sizeof(fsl_int) * residue.allocated_size)
        if residue.atom_internal_ids == NULL:
            raise MemoryError
        residue.atomids = \
            <fsl_int *> PyMem_Realloc(<void *> residue.atomids, sizeof(fsl_int) * residue.allocated_size)
        if residue.atomids == NULL:
            raise MemoryError

    atom.resid = residue.resid
    atom.residue_internal_id = residue.internal_id
    residue.atom_internal_ids[residue.size] = atom.internal_id
    residue.atomids[residue.size] = atom.atomid
    residue.size += 1

cdef class Topology(object):

    def __init__(self):
        self.atomnames = \
            <strname_t *> PyMem_Realloc(NULL, sizeof(strname_t) * NAME_ALLOCATION_INCREMENT)
        if self.atomnames == NULL:
            raise MemoryError
        self.atomnames_size = 0
        self.atomnames_allocated_size = NAME_ALLOCATION_INCREMENT

        self.resnames = \
            <strname_t *> PyMem_Realloc(NULL, sizeof(strname_t) * NAME_ALLOCATION_INCREMENT)
        if self.resnames == NULL:
            raise MemoryError
        self.resnames_size = 0
        self.resnames_allocated_size = NAME_ALLOCATION_INCREMENT

        self.atoms = NULL
        self.atoms_size = 0
        self.atoms_allocated_size = 0
        self.set_atoms_allocation(ATOMS_ALLOCATION_INCREMENT)

        self.atomids_to_internalids = \
            <fsl_int *> PyMem_Realloc(NULL, sizeof(fsl_int) * ATOMS_ALLOCATION_INCREMENT)
        if self.atomids_to_internalids == NULL:
            raise MemoryError
        for i in range(ATOMS_ALLOCATION_INCREMENT):
            self.atomids_to_internalids[i] = NOTSET
        self.atomids_to_internalids_size = 0
        self.atomids_to_internalids_allocated_size = ATOMS_ALLOCATION_INCREMENT

        self.residues = NULL
        self.residues_size = 0
        self.residues_allocated_size = 0
        self.set_residues_allocation(RESIDUES_ALLOCATION_INCREMENT)

        self.resids_to_internalids = \
            <fsl_int *> PyMem_Realloc(NULL, sizeof(fsl_int) * RESIDUES_ALLOCATION_INCREMENT)
        if self.resids_to_internalids == NULL:
            raise MemoryError
        for i in range(RESIDUES_ALLOCATION_INCREMENT):
            self.resids_to_internalids[i] = NOTSET

        self.resids_to_internalids_size = 0
        self.resids_to_internalids_allocated_size = RESIDUES_ALLOCATION_INCREMENT

    def __dealloc__(self):
        cdef fsl_int i
        PyMem_Free(self.atomnames)
        PyMem_Free(self.resnames)

        for i in range(self.residues_size):
            PyMem_Free(self.residues[i].atom_internal_ids)
            PyMem_Free(self.residues[i].atomids)
        PyMem_Free(self.residues)
        PyMem_Free(self.resids_to_internalids)

        PyMem_Free(self.atoms)
        PyMem_Free(self.atomids_to_internalids)

    cdef void set_atoms_allocation(self, fsl_int size):
        cdef fsl_int i
        self.atoms = <topol_atom_t *> PyMem_Realloc(<void *> self.atoms,
                                                    sizeof(topol_atom_t) * size)
        if self.atoms == NULL:
            raise MemoryError

        for i in range(self.atoms_allocated_size, size):
            self.atoms[i].atomid = NOTSET
            self.atoms[i].internal_id = NOTSET
            self.atoms[i].resid = NOTSET
            self.atoms[i].residue_internal_id = NOTSET
            self.atoms[i].name_id = NOTSET

        self.atoms_allocated_size = size

    cdef void set_residues_allocation(self, fsl_int size):
        cdef fsl_int i
        self.residues = <topol_residue_t *> PyMem_Realloc(<void *>self.residues,
                                                          sizeof(topol_residue_t) * size)
        if self.residues == NULL:
            raise MemoryError

        # Initialize residues
        for i in range(self.residues_allocated_size, size):

            self.residues[i].name_id = NOTSET
            self.residues[i].internal_id = NOTSET
            self.residues[i].resid = NOTSET
            self.residues[i].atom_internal_ids = \
                <fsl_int *> PyMem_Realloc(NULL, sizeof(fsl_int) * RESIDUE_ALLOCATION_INCREMENT)
            self.residues[i].atomids = \
                <fsl_int *> PyMem_Realloc(NULL, sizeof(fsl_int) * RESIDUE_ALLOCATION_INCREMENT)
            self.residues[i].allocated_size = RESIDUE_ALLOCATION_INCREMENT
            self.residues[i].size = 0

        self.residues_allocated_size = size

    def append(self, fsl_int resid, str resname_str, str atomname_str, fsl_int atomid, bint strip_name=True):
        cdef bytes resname, atomname

        if strip_name:
            resname = resname_str.strip().encode()
            atomname = atomname_str.strip().encode()
        else:
            resname = resname_str.encode()
            atomname = atomname_str.encode()

        self.internal_append(resid, resname, atomname, atomid)

    cdef void internal_append(self, fsl_int resid, bytes resname, bytes atomname, fsl_int atomid):
        cdef fsl_int i
        cdef fsl_int atomname_id, resname_id
        cdef fsl_int residue_internal_id, atom_internal_id
        cdef fsl_int old_size
        cdef bint residue_created = False
        cdef topol_residue_t *residue
        cdef topol_atom_t *atom

        # Handles atomname
        atomname_id = -1
        for i in range(self.atomnames_size):
            if strncmp(self.atomnames[i], atomname, NAME_SIZE-1) == 0:
                atomname_id = i
                break
        if atomname_id == -1: # We need to add the new atomname
            # Allocate memory if necessary
            if self.atomnames_size == self.atomnames_allocated_size:
                self.atomnames_allocated_size += NAME_ALLOCATION_INCREMENT
                self.atomnames = <strname_t *> PyMem_Realloc(<void *> self.atomnames,
                                                             sizeof(strname_t) *
                                                             self.atomnames_allocated_size)
                if self.atomnames == NULL:
                    raise MemoryError

            atomname_id = self.atomnames_size
            self.atomnames_size += 1
            strncpy(self.atomnames[atomname_id], atomname, NAME_SIZE)
            self.atomnames[atomname_id][NAME_SIZE-1] = b'\0'

        # Handles resname
        resname_id = -1
        for i in range(self.resnames_size):
            if strncmp(self.resnames[i], resname, NAME_SIZE-1) == 0:
                resname_id = i
                break
        if resname_id == -1: # We need to add the new resname
            # Allocate memory if necessary
            if self.resnames_size == self.resnames_allocated_size:
                self.resnames_allocated_size += NAME_ALLOCATION_INCREMENT
                self.resnames = <strname_t *> PyMem_Realloc(<void *> self.resnames,
                                                            sizeof(strname_t) *
                                                            self.resnames_allocated_size)
                if self.resnames == NULL:
                    raise MemoryError
            resname_id = self.resnames_size
            self.resnames_size += 1
            strncpy(self.resnames[resname_id], resname, NAME_SIZE)
            self.resnames[resname_id][NAME_SIZE-1] = b'\0'

        # Handles resid

        # First we need to handle the translation resid->internalid
        if resid >= self.resids_to_internalids_allocated_size:
            # Allocate enough memory to store the new resid plus extra for future resids
            self.resids_to_internalids_allocated_size += \
                (resid - self.resids_to_internalids_allocated_size) + RESIDUES_ALLOCATION_INCREMENT
            self.resids_to_internalids = \
                <fsl_int *> PyMem_Realloc(<void *>self.resids_to_internalids,
                                      sizeof(fsl_int) * self.resids_to_internalids_allocated_size)
            if self.resids_to_internalids == NULL:
                raise MemoryError

            # We also need to initialize the newly allocated memory:
            for i in range(self.resids_to_internalids_size,
                           self.resids_to_internalids_allocated_size):
                self.resids_to_internalids[i] = NOTSET
            self.resids_to_internalids_size = resid

        # Watch out! resids start from 1 but are stored from 0
        residue_internal_id = self.resids_to_internalids[resid-1]

        # Second we need to handle the internal id
        if residue_internal_id == NOTSET: # The residue is not registered yet
            # Allocate memory if necessary
            if self.residues_size == self.residues_allocated_size:
                self.set_residues_allocation(self.residues_allocated_size
                                             + RESIDUES_ALLOCATION_INCREMENT)

            residue_internal_id = self.residues_size
            self.residues_size += 1
            residue_created = True
            self.resids_to_internalids[resid - 1] = residue_internal_id
            self.resids_to_internalids_size = resid


        # Handles atomid

        # First we need to handle the translation atomid->internalid
        if atomid >= self.atomids_to_internalids_allocated_size:
            # Allocate enough memory to store the new resid plus extra for future resids
            self.atomids_to_internalids_allocated_size += \
                (atomid - self.atomids_to_internalids_allocated_size) + ATOMS_ALLOCATION_INCREMENT
            self.atomids_to_internalids = \
                <fsl_int *> PyMem_Realloc(<void *>self.atomids_to_internalids,
                                      sizeof(fsl_int) * self.atomids_to_internalids_allocated_size)
            if self.atomids_to_internalids == NULL:
                raise MemoryError

            # We also need to initialize the newly allocated memory:
            for i in range(self.atomids_to_internalids_size,
                           self.atomids_to_internalids_allocated_size):
                self.atomids_to_internalids[i] = NOTSET

        # Watch out! atomids start from 1 but are stored from 0
        atom_internal_id = self.atomids_to_internalids[atomid-1]

        # Second we create the atom
        if atom_internal_id != NOTSET: # The atom is not registered yet
            raise IndexError("Atom with atomid #%i is already registered" % atomid)

        # Allocate memory if necessary
        if self.atoms_size == self.atoms_allocated_size:
            self.set_atoms_allocation(self.atoms_allocated_size + ATOMS_ALLOCATION_INCREMENT)

        atom_internal_id = self.atoms_size
        self.atoms_size += 1

        # Create the atom
        atom = &self.atoms[atom_internal_id]
        atom.name_id = atomname_id
        atom.internal_id = atom_internal_id
        atom.atomid = atomid
        self.atomids_to_internalids[atomid-1] = atom_internal_id
        self.atomids_to_internalids_size = atomid

        # Populate the residue
        residue = &self.residues[residue_internal_id]

        if residue_created:
            residue.name_id = resname_id
            residue.internal_id = residue_internal_id
            residue.resid = resid

        # Add the atom to the residue
        topol_residue_append_atom(residue, atom)

    cdef fsl_int fast_get_atom_internal_id(self, fsl_int atomid) nogil:
        if not 0 < atomid <= self.atomids_to_internalids_size:
            return -1
        return self.atomids_to_internalids[atomid-1]

    cdef topol_atom_t *fast_get_atom(self, fsl_int atomid) nogil:
        cdef fsl_int internal_id = self.fast_get_atom_internal_id(atomid)

        if internal_id < 0:
            return NULL

        return &self.atoms[internal_id]

    def get_atom(self, fsl_int atomid):
        cdef topol_atom_t atom
        cdef topol_residue_t residue
        cdef fsl_int internal_id = self.fast_get_atom_internal_id(atomid)
        cdef bytes resname, atomname

        if internal_id < 0:
            for i in range(self.atomids_to_internalids_size):
                print("Atomid #%i -> internal id #%i" % (i+1, self.atomids_to_internalids[i]))
            raise KeyError("atomid not valid: %i (max atomid:%i)" %
                           (atomid, self.atomids_to_internalids_size))

        atom = self.atoms[internal_id]
        residue = self.residues[atom.residue_internal_id]

        return atom.resid, \
               str(self.resnames[residue.name_id].decode()), \
               str(self.atomnames[atom.name_id].decode()), \
               atom.atomid

    cdef fsl_int fast_get_resid_from_atomid(self, fsl_int atomid) nogil:
        cdef topol_atom_t atom
        cdef fsl_int internal_id = self.fast_get_atom_internal_id(atomid)

        if internal_id < 0:
            return -1

        atom = self.atoms[internal_id]
        return atom.resid

    cdef topol_residue_t *fast_get_residue_from_atomid(self, fsl_int atomid) nogil:
        cdef topol_atom_t atom
        cdef fsl_int internal_id = self.fast_get_atom_internal_id(atomid)

        if internal_id < 0:
            return NULL
        atom = self.atoms[internal_id]

        return &self.residues[atom.residue_internal_id]

    def get_residue_from_atomid(self, fsl_int atomid):
        cdef topol_residue_t *residue = self.fast_get_residue_from_atomid(atomid)

        if residue == NULL:
            raise KeyError("atomid not valid: %i" % atomid)

        return topol_residue_to_dict(residue, self.resnames[residue.name_id])

    cdef fsl_int fast_get_residue_internal_id(self, fsl_int resid) nogil:
        if not 0 < resid <= self.resids_to_internalids_size:
            return -1
        return self.resids_to_internalids[resid-1]

    cdef topol_residue_t *fast_get_residue(self, fsl_int resid) nogil:
        cdef fsl_int internal_id = self.fast_get_residue_internal_id(resid)

        if internal_id < 0:
            return NULL

        return &self.residues[internal_id]

    def get_residue(self, resid):
        cdef topol_residue_t *residue = self.fast_get_residue(resid)

        if residue == NULL:
            raise KeyError("resid not valid: %i" % resid)

        return topol_residue_to_dict(residue, self.resnames[residue.name_id])

    property nresidues:
        def __get__(self):
            return self.residues_size

    property natoms:
        def __get__(self):
            return self.atoms_size

    def __len__(self):
        return self.atoms_size

cdef class TopologyReader(object):

    def __init__(self, str filename, bint verbose=True):
        self.filename = filename.encode()
        self.topology = Topology()

        # load topology file
        verbose_print("Loading topology from '%s'... " %
                      os.path.basename(self.filename.decode()),
                      verbose,
                      end="")
        begin = time()
        self.load()
        verbose_print("%i topology atoms (from %i residues) loaded in %s" %
                      (self.topology.natoms,
                       self.topology.nresidues,
                       pretty_delta(begin, time())),
                      verbose)

    cdef load(self):
        raise NotImplementedError



cdef class CoordinateReader(object):

    def __init__(self, str filename, bint verbose=True):
        self.filename = filename.encode()

        # Preload file
        verbose_print("Preloading coordinates from '%s'... " %
                      os.path.basename(self.filename.decode()),
                      verbose,
                      end="")
        begin = time()
        self.preload()
        verbose_print("%i frames preloaded in %s" % (self.nframes, pretty_delta(begin, time())),
                      verbose)

    cdef preload(self):
        raise NotImplementedError

    cdef PBCBox load_box(self, fsl_int frame_id):
        raise NotImplementedError

    cdef assert_frame_id(self, fsl_int frame_id):
        if not 0 <= frame_id < self.nframes:
            raise IndexError("Frame index error: %i (%i frames)" % (frame_id, self.nframes))

    cdef real[:,::1] load_coords(self, fsl_int frame_id, fsl_int[:] atomids) nogil except *:
        with gil:
            raise NotImplementedError


cdef class IndexReader(object):

    def __init__(self, str filename, bint verbose=True):
        self.filename = filename.encode()
        self.loaded = False

        self.groups = {}


        verbose_print("Loading groups from '%s'... " %
                      os.path.basename(self.filename),
                      verbose,
                      end="")
        begin = time()
        self.load()
        self.loaded = True
        verbose_print("%i groups loaded in %s" % (len(self), pretty_delta(begin, time())),
                      verbose)

    def load(self):
        if self.loaded:
            return
        self.fast_load()

    cdef fast_load(self):
        raise NotImplementedError

    def __len__(self):
        return len(self.groups)

    def __getitem__(self, str item):
        return self.groups[item.encode()]


cdef class Frame(object):
    def __init__(self, Trajectory trajectory, fsl_int index, real timestep):
        self.trajectory = trajectory
        self.topology = trajectory.topology
        self.coords_reader = trajectory.coords_reader
        self.index = index
        self.timestep = timestep
        self.box = self.coords_reader.load_box(index)
        # -1 line below => there is one more offset compared to actual number of lipids!
        self.size = self.trajectory.lipid_atomids_offsets.shape[0] - 1

        # Headgroups
        self.hg_group_coords = None
        self.hg_coords_bbox = None
        # Lipids
        self.lipid_coords = None
        self.lipid_coords_aslist = None
        self.lipid_coords_bbox = None
        self.lipids_coords_bbox_aslist = None
        self.lipid_atomids_offsets = self.trajectory.lipid_atomids_offsets
        # Simplified lipids
        self.bead_coords = None
        self.bead_coords_bbox = None
        # Interacting groups
        self.interacting_group_coords = None
        self.interacting_group_coords_bbox = None

        # NS-related
        self.proximity_cutoff = NOTSET
        self.neighbors = NULL

        # Directions & Normals
        self.directions = None
        self.normals = None

        # Aggregates
        self.aggregates = []
        self.aggregates_retrieved = False
        self.aggregates_cutoff = NOTSET

        # Leaflets & Membranes
        self.membranes = []
        self.membranes_retrieved = False
        self.membranes_cutoff = NOTSET

    def __dealloc__(self):
        free_neighborhood_holder(self.neighbors)
        self.neighbors = NULL

    def __len__(self):
        return self.size

    cpdef fsl_int get_forcefield_type(self):
        return self.trajectory.get_forcefield_type()

    cdef real[:,::1] fast_get_hg_group_coords(self) nogil:
        if self.hg_group_coords is None:
            with gil:
                self.hg_group_coords = self.coords_reader.\
                    load_coords(self.index,self.trajectory.hg_group_atomids)
        return self.hg_group_coords

    # Associated python property
    property hg_group_coords:
        def __get__(self):
            cdef real[:, ::1] memview
            with nogil:
                memview = self.fast_get_hg_group_coords()
            return np.asarray(memview)

    cdef real[:,::1] fast_get_hg_group_coords_bbox(self) nogil:
        if self.hg_coords_bbox is None:
            self.hg_coords_bbox = self.box.fast_put_atoms_in_bbox(self.fast_get_hg_group_coords())
        return self.hg_coords_bbox

    # Associated python property
    property hg_group_coords_bbox:
        def __get__(self):
            cdef real[:, ::1] memview
            with nogil:
                memview = self.fast_get_hg_group_coords_bbox()
            return np.asarray(memview)

    cdef real[:, ::1] fast_get_lipid_coords(self) nogil:
        if self.lipid_coords is None:
            self.lipid_coords = self.coords_reader.\
                load_coords(self.index, self.trajectory.lipid_atomids)
        return self.lipid_coords

    cdef list get_lipid_coords_aslist(self):
        cdef fsl_int[:] offsets = self.lipid_atomids_offsets
        cdef real[:, ::1] coords = self.fast_get_lipid_coords()
        cdef fsl_int i
        if self.lipid_coords_aslist is None:
            coords_list = []

            for i in range(offsets.shape[0] - 1):
                coords_list.append(np.asarray(coords[offsets[i]:offsets[i+1]]))

            self.lipid_coords_aslist = coords_list

        return self.lipid_coords_aslist

    # Associated python property
    property lipid_coords:
        def __get__(self):
            return self.get_lipid_coords_aslist()

    cdef real[:, ::1] fast_get_lipid_coords_bbox(self) nogil:
        if self.lipid_coords_bbox is None:
            self.lipid_coords_bbox = self.box.fast_put_atoms_in_bbox(self.fast_get_lipid_coords())
        return self.lipid_coords_bbox

    cdef list get_lipid_coords_bbox_aslist(self):
        cdef fsl_int[:] offsets = self.lipid_atomids_offsets
        cdef real[:, ::1] coords = self.fast_get_lipid_coords_bbox()
        cdef fsl_int i
        if self.lipids_coords_bbox_aslist is None:
            coords_list = []

            for i in range(offsets.shape[0] - 1):
                coords_list.append(np.asarray(coords[offsets[i]:offsets[i+1]]))

            self.lipids_coords_bbox_aslist = coords_list

        return self.lipids_coords_bbox_aslist

    # Associated python property
    property lipid_coords_bbox:
        def __get__(self):
            return self.get_lipid_coords_bbox_aslist()

    cdef real[:,::1] fast_get_interacting_group_coords(self) nogil:
        if self.interacting_group_coords is None:
            with gil:
                if self.trajectory.interacting_atomids.shape[0] > 0:
                    self.interacting_group_coords = self.coords_reader.\
                        load_coords(self.index, self.trajectory.interacting_atomids)
                else:
                    self.interacting_group_coords = np.empty((0, DIM))
        return self.interacting_group_coords

    # Associated python property
    property interacting_group_coords:
        def __get__(self):
            cdef real[:, ::1] memview
            with nogil:
                memview = self.fast_get_interacting_group_coords()
            return np.asarray(memview)


    cdef real[:,::1] fast_get_interacting_group_coords_bbox(self) nogil:
        cdef real[:, ::1] coords_bbox = self.fast_get_interacting_group_coords()
        if self.interacting_group_coords_bbox is None:
            self.interacting_group_coords_bbox = self.box.fast_put_atoms_in_bbox(coords_bbox)
        return self.interacting_group_coords_bbox

    # Associated python property
    property interacting_group_coords_bbox:
        def __get__(self):
            cdef real[:, ::1] memview
            with nogil:
                memview = self.fast_get_interacting_group_coords_bbox()
            return np.asarray(memview)


    cdef real[:,::1] fast_get_bead_coords(self) nogil:
        cdef real[:, ::1] bead_coords, coords
        cdef real[:,::1] tmp_coords
        cdef fsl_int[:] hg_offsets = self.trajectory.hg_bead_atomids_offsets
        cdef fsl_int[:] hg_indices = self.trajectory.lipid_hg_indices
        cdef fsl_int i, j, residue_hg_size
        cdef rvec xcm

        if self.bead_coords is not None:
            return self.bead_coords

        # Load coords & allocate memory
        with gil:
            bead_coords = np.empty((self.trajectory.hg_bead_atomids_offsets.shape[0] - 1, DIM))
            tmp_coords = np.empty((self.trajectory.hg_bead_atomids_offsets.shape[0] - 1, DIM))
            coords = self.fast_get_lipid_coords()

        for i in range(hg_offsets.shape[0] - 1):
            residue_hg_size = hg_offsets[i+1] - hg_offsets[i]
            for j in range(residue_hg_size):
                tmp_coords[j, XX] = coords[hg_indices[hg_offsets[i] + j], XX]
                tmp_coords[j, YY] = coords[hg_indices[hg_offsets[i] + j], YY]
                tmp_coords[j, ZZ] = coords[hg_indices[hg_offsets[i] + j], ZZ]
            self.box.fast_pbc_xcm(tmp_coords[:residue_hg_size], xcm)
            rvec_to_memview(xcm, bead_coords[i])

        self.bead_coords = bead_coords
        return bead_coords

    # Associated python property
    property bead_coords:
        def __get__(self):
            cdef real[:, ::1] memview
            with nogil:
                memview = self.fast_get_bead_coords()
            return np.asarray(memview)


    cdef real[:,::1] fast_get_bead_coords_bbox(self) nogil:
        if self.bead_coords_bbox is None:
            self.bead_coords_bbox = self.box.fast_put_atoms_in_bbox(self.fast_get_bead_coords())
        return self.bead_coords_bbox

    # Associated python property
    property bead_coords_bbox:
        def __get__(self):
            cdef real[:, ::1] memview
            with nogil:
                memview = self.fast_get_bead_coords_bbox()
            return np.asarray(memview)

    cdef real[:, ::1] fast_get_directions(self) nogil except *:
        if self.directions is not None:
            return self.directions

        cdef fsl_int i
        cdef real[:,::1] bead_coords
        cdef fsl_int beads
        cdef real[:,::1] lipid_coords, all_lipid_coords
        cdef fsl_int[:] offsets = self.lipid_atomids_offsets
        cdef cPBCBox_t box = self.box.c_pbcbox
        cdef real[:, ::1] directions

        # Loads coordinates
        all_lipid_coords = self.fast_get_lipid_coords_bbox()
        bead_coords = self.fast_get_bead_coords_bbox()
        nbeads = bead_coords.shape[0]

        with gil:
            directions = np.empty((nbeads, DIM))

        #for i in range(nbeads):
        for i in prange(nbeads, schedule='dynamic', num_threads=OPENMP_NUM_THREADS):
            calculate_lipid_direction(&bead_coords[i, XX],
                                      box,
                                      all_lipid_coords,
                                      offsets[i],
                                      offsets[i+1],
                                      &directions[i, XX])

        with gil:
            self.directions = directions
        return directions

    # Associated python property
    property directions:
        def __get__(self):
            cdef real[:, ::1] memview
            with nogil:
                memview = self.fast_get_directions()
            return np.asarray(memview)


    cdef real[:, ::1] fast_get_normals(self, real proximity_cutoff=-1) nogil:
        cdef fsl_int ncoords = self.size
        cdef real[:, ::1] coords = self.fast_get_bead_coords_bbox()
        cdef real[:, ::1] directions = self.fast_get_directions()
        cdef real[:, ::1] normals
        cdef fsl_int *neighborhood_buffer = NULL
        cdef fsl_int buffer_size = 0
        cdef fsl_int nid, i

        if self.proximity_cutoff == proximity_cutoff and self.normals is not None:
            return self.normals

        if proximity_cutoff < 0:
            if self.normals is not None:
                return self.normals
            else:
                self.proximity_cutoff = self.trajectory.normal_cutoff
        else:
            self.proximity_cutoff = proximity_cutoff

        self.trajectory.normal_cutoff = self.proximity_cutoff

        # Delete old ref if needed
        if self.neighbors != NULL:
            free_neighborhood_holder(self.neighbors)
            self.neighbors = NULL

        # Retrieve neighborhoods
        self.neighbors = fast_neighbor_search(coords, coords,
                                  self.box,
                                  self.proximity_cutoff)

        with gil:
            normals = np.empty((ncoords, DIM))

        for i in range(ncoords):
            if self.neighbors.neighborhoods[i].size > buffer_size:
                buffer_size = self.neighbors.neighborhoods[i].size

        neighborhood_buffer = <fsl_int *> malloc(buffer_size * sizeof(fsl_int) * OPENMP_NUM_THREADS)
        if neighborhood_buffer == NULL:
            abort()

        for nid in prange(ncoords, schedule='dynamic', num_threads=OPENMP_NUM_THREADS):
        #for nid in range(ncoords):
            calculate_neighborhood_normal(nid,
                                          coords,
                                          directions,
                                          self.neighbors.neighborhoods[nid],
                                          self.box,
                                          &normals[nid, XX],
                                          &neighborhood_buffer[omp_get_thread_num()
                                                               * buffer_size])

        free(neighborhood_buffer)

        self.normals = normals
        return self.normals

    # Associated python property
    property normals:
        def __get__(self):
            cdef real[:, ::1] memview
            with nogil:
                memview = self.fast_get_normals()
            #memview = np.empty((self.size, DIM))
            return np.asarray(memview)


    cpdef list get_aggregates(self, real cutoff=DEFAULT_PROXIMITY_CUTOFF, bint update=False):
        if self.aggregates_cutoff == cutoff and self.aggregates_retrieved and not update:
            return self.aggregates

        self.aggregates_cutoff = cutoff
        if not update \
                and self.trajectory.aggregates is not None \
                and self.trajectory.normal_cutoff == cutoff:
            traj_aggregates = self.trajectory.aggregates
            self.aggregates = []
            for aggregate in traj_aggregates:
                self.aggregates.append(Aggregate(self, aggregate))
        else:
            self.aggregates = core_analysis.retrieve_aggregates(self, cutoff)
            self.trajectory.aggregates = []
            for aggregate in self.aggregates:
                self.trajectory.aggregates.append(np.asarray(aggregate.beadids))
        self.aggregates_retrieved = True

        return self.aggregates

    cpdef list get_membranes(self, real cutoff=DEFAULT_PROXIMITY_CUTOFF, bint update=False):
        if self.membranes_cutoff == cutoff and self.membranes_retrieved and not update:
            return self.membranes

        self.proximity_cutoff = cutoff
        if not update and self.trajectory.membranes is not None:
            traj_membranes = self.trajectory.membranes
            self.membranes = []
            for membrane in traj_membranes:
                l1 = Aggregate(self, membrane[0])
                l2 = Aggregate(self, membrane[1])
                self.membranes.append(Membrane(self, l1, l2))
        else:
            self.membranes = core_analysis.retrieve_membranes(self, cutoff)
            traj_membranes = []
            for membrane in self.membranes:
                memb_beadids = [np.asarray(membrane.leaflet1.beadids),
                                np.asarray(membrane.leaflet2.beadids)]

                traj_membranes.append(memb_beadids)
            self.trajectory.membranes = traj_membranes
        self.membranes_retrieved = True

        return self.membranes

    # Additional python properties
    property hg_atomids:
        def __get__(self):
            return np.asarray(self.trajectory.hg_group_atomids)

    property lipid_atomids:
        def __get__(self):
            cdef fsl_int i
            cdef fsl_int[:] offsets = self.trajectory.lipid_atomids_offsets

            lipid_atomids = []
            for i in range(offsets.shape[0] - 1):
                lipid_atomids.append(np.asarray(self.trajectory.
                                                lipid_atomids[offsets[i]:offsets[i+1]]))
            return lipid_atomids

    property index:
        def __get__(self):
            return self.index


FORCEFIELDTYPE_NOTSET = fft_notset
FORCEFIELDTYPE_ALLATOM = fft_allatom
FORCEFIELDTYPE_UNIFIED = fft_unified
FORCEFIELDTYPE_COARSE = fft_coarse
FORCEFIELDTYPE_UNKNOWN = fft_unknown

cdef class Trajectory(object):
    def __init__(self,
                 TopologyReader top_reader,
                 IndexReader index_reader,
                 CoordinateReader coords_reader,
                 bint be_verbose=True):
        self.be_verbose = be_verbose
        self.bead_group_name = None
        self.interacting_group_name = None

        assert isinstance(top_reader, TopologyReader)
        self.topol_reader = top_reader

        assert isinstance(index_reader, IndexReader)
        self.index_reader = index_reader

        assert isinstance(coords_reader, CoordinateReader)
        self.coords_reader = coords_reader

        # Load atomids & group
        self.hg_group_atomids = None
        self.lipid_atomids = None
        self.lipid_atomids_offsets = None
        self.interacting_atomids = None
        self.hg_bead_atomids = None
        self.hg_bead_atomids_offsets = None
        self.forcefield_type = fft_notset
        cur_frame = 0

        # Cache related
        self.normal_cutoff = DEFAULT_PROXIMITY_CUTOFF
        self.membranes = None
        self.aggregates = None

        self.topology = None
        self.initialized = False


    def initialize(self, str hg_group="headgroups", str interacting_group="protein"):
        cdef fsl_int size, i, resid_size
        cdef fsl_int last_resid, resid
        cdef fsl_int internal_id
        cdef fsl_int[:] hg_group_atomids, hg_bead_atomids_offsets, lipid_atomids_offsets, \
            lipid_atomids, resids, lipid_hg_indices
        cdef Topology topology
        cdef topol_residue_t *residue
        cdef topol_atom_t *atom
        cdef dict residues_dict
        cdef dict resdict
        cdef fsl_int atom_index

        # Preload file
        verbose_print("Initializing trajectory using groups: '%s' and '%s'... "
                      % (hg_group, interacting_group),
                      self.be_verbose,
                      end="")
        sys.stdout.flush()
        begin = time()

        topology = self.topol_reader.topology
        self.topology = topology

        # Handle headgroups
        self.bead_group_name = hg_group.encode()
        hg_group_atomids = self.index_reader[hg_group]
        self.hg_group_atomids = hg_group_atomids

        # Group by resid and retrieve offsets
        size = hg_group_atomids.shape[0]
        resids = np.empty(size, dtype=np.int64)
        hg_bead_atomids_offsets = np.empty(size+1, dtype=np.int64)
        lipid_atomids_offsets = np.empty(size+1, dtype=np.int64)
        lipid_hg_indices = np.empty(size, dtype=np.int64)
        lipid_atomids_list = []
        last_resid = -1
        lipid_atomids_offsets[0] = 0
        resid_size = 0

        residues_dict = {}
        for i in range(size):
            atom = topology.fast_get_atom(hg_group_atomids[i])
            if atom == NULL:
                raise KeyError("Atomid not registered in topology: %i. Does group '%s' from '%s'"
                                   " file contains H?" % (hg_group_atomids[i],
                                                          hg_group,
                                                          self.index_reader.filename))
            residue = &topology.residues[atom.residue_internal_id]

            resid = residue.resid

            if resid > last_resid: # WARNING: It is assumed that atomids (hence resids) are sorted!
                resids[resid_size] = resid
                hg_bead_atomids_offsets[resid_size] = i
                lipid_atomids_offsets[resid_size + 1] = lipid_atomids_offsets[resid_size] \
                                                        + residue.size
                lipid_atomids_list.extend(topol_residue_atomids_as_list(residue))
                last_resid = resid
                resid_size += 1


            try:
                atom_index = residues_dict[residue.name_id][atom.name_id]
            except KeyError:
                atom_index = topol_residue_atomid_index(residue, atom.atomid)
                try:
                    residues_dict[residue.name_id][atom.name_id] = atom_index
                except KeyError:
                    residues_dict[residue.name_id] = {atom.name_id: atom_index}

            lipid_hg_indices[i] = atom_index + lipid_atomids_offsets[resid_size-1]

        hg_bead_atomids_offsets[resid_size] = size
        lipid_atomids = np.array(lipid_atomids_list, dtype=np.int64)

        self.hg_bead_atomids = np.array(hg_group_atomids, dtype=np.int64)
        self.hg_bead_atomids_offsets = hg_bead_atomids_offsets[:resid_size + 1]

        self.lipid_atomids = lipid_atomids
        self.lipid_atomids_offsets = lipid_atomids_offsets[:resid_size + 1]
        self.lipid_hg_indices = lipid_hg_indices


        # Handle interacting group (if any)
        self.interacting_group_name = interacting_group.encode()
        try:
            self.interacting_atomids = self.index_reader[interacting_group]
        except KeyError:
            self.interacting_atomids = np.array([], dtype=np.int64)

        # Sanity check the interacting atomids
        for i in range(self.interacting_atomids.shape[0]):
            if topology.fast_get_atom_internal_id(self.interacting_atomids[i]) < 0:
                raise KeyError("Atomid not registered in topology: %i. Does group '%s' from '%s'"
                               " file contains H?" % (self.interacting_atomids[i],
                                                      interacting_group,
                                                      self.index_reader.filename))

        # Identify the forcefield type
        self.get_forcefield_type()

        self.initialized = True

        verbose_print("Done in %s" % (pretty_delta(begin, time())),
                      self.be_verbose)
        sys.stdout.flush()


    cpdef fsl_int get_forcefield_type(self):
        cdef real rvdw
        if self.forcefield_type != fft_notset: # Already done
            return self.forcefield_type

        rvdw = self.get_rvdw()

        if rvdw < FFT_ALLATOM_THRESHOLD:
            self.forcefield_type = fft_allatom
        elif rvdw < FFT_UNIFIED_THRESHOLD:
            self.forcefield_type = fft_unified
        elif rvdw < FFT_COARSE_THRESHOLD:
            self.forcefield_type = fft_coarse
        else:
            self.forcefield_type = fft_unknown

    cpdef real get_rvdw(self):
        cdef real[:, ::1] coords
        cdef fsl_int i
        cdef real dist, min_dist = 1.0e5
        cdef fsl_int lastid

        if self.lipid_atomids is None:
            raise ValueError("Trajectory is not initialized correctly!")

        # Load the first lipid from the first frame
        coords = self.coords_reader.load_coords(0, np.arange(self.lipid_atomids_offsets[0]+1,
                                                             self.lipid_atomids_offsets[1]+1,
                                                             dtype=np.int64))


        # We take the last atom as a ref, it is fair to assume that it corresponds to a carbon
        # from the lipid tail or a non-polar hydrogen (if we are dealing with an all-atom FF).
        # Then we get the minimal distance between this atom and the other. If this min dist
        # corresponds to a C-H bond it means the FF is all-atom (because the H is non polar).
        # If the distance is C-C, wa have unified atom. If greater, it must be a coarse grained FF.
        lastid = coords.shape[0] - 1
        for i in range(lastid):
            dist = (coords[i, XX] - coords[lastid, XX]) * (coords[i, XX] - coords[lastid, XX]) +\
                (coords[i, YY] - coords[lastid, YY]) * (coords[i, YY] - coords[lastid, YY]) +\
                (coords[i, ZZ] - coords[lastid, ZZ]) * (coords[i, ZZ] - coords[lastid, ZZ])

            if dist < min_dist:
                min_dist = dist


        return min_dist


    cdef assert_initialization(self):
        if not self.initialized:
            raise ValueError("Trajectory not initialized")

    def __getitem__(self, fsl_int item):
        self.assert_initialization()
        cdef real timestep
        item = int(item)
        if not 0 <= item < self.coords_reader.nframes:
            raise IndexError("Frame index out of range (%i frames in trajectory)" %
                             self.coords_reader.nframes)
        timestep = self.coords_reader.timesteps[item]
        return Frame(self, item, timestep)

    def __len__(self):
        return self.coords_reader.nframes


    def get_topology(self):
        return self.topology






