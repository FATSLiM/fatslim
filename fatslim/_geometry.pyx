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
DEF DIM = 3
DEF XX = 0
DEF YY = 1
DEF ZZ = 2
DEF EPSILON = 1e-6

cimport cython
from libc.math cimport fabs, sqrt, atan2, sin, cos

cimport numpy as np
import numpy as np

from ._typedefs cimport real, fsl_int
from ._typedefs cimport rvec_norm2, rvec_smul, rvec_copy, rvec_inc, rvec_clear, rvec_dprod, rvec_dec, rvec_normalize,\
    rvec_cprod_norm, rvec_cprod
from ._typedefs cimport mat_clear, mat_copy

###############################
# Utility class to handle PBC #
###############################


# Class to handle PBC calculations
@cython.initializedcheck(False)
@cython.boundscheck(False)
cdef class PBCBox(object):
    """
    Cython implementation of
    `PBC-related <https://en.wikipedia.org/wiki/Periodic_boundary_conditions>`_
    operations. This class is used by classes :class:`FastNS`
    and :class:`_NSGrid` to put all particles inside a brick-shaped box
    and to compute PBC-aware distance. The class can also handle
    non-PBC aware distance evaluations through ``periodic`` argument.

    .. warning::
        This class is not meant to be used by end users.

    .. warning::
        Even if MD triclinic boxes can be handled by this class,
        internal optimization is made based on the assumption that
        particles are inside a brick-shaped box. When this is not
        the case, calculated distances are not
        warranted to be exact.
    """

    def __init__(self, real[:, ::1] box):
        """
        Parameters
        ----------
        box : numpy.ndarray
            box vectors of shape ``(3, 3)`` or
            as returned by ``MDAnalysis.lib.mdamath.triclinic_vectors``
            ``dtype`` must be ``numpy.float32``
        periodic : boolean
            ``True`` for PBC-aware calculations
            ``False`` for non PBC aware calculations
        """

        self.is_triclinic = False
        self.update(box)

    cdef void fast_update(self, real[:, ::1] box) nogil:
        """
        Updates the internal box parameters for
        PBC-aware distance calculations. The internal
        box parameters are used to define the brick-shaped
        box which is eventually used for distance calculations.

        """
        cdef fsl_int i, j
        cdef dreal min_hv2, min_ss, tmp

        # Update matrix
        self.is_triclinic = False
        for i in range(DIM):
            for j in range(DIM):
                self.c_pbcbox.box[i][j] = box[i, j]

                if i != j:
                    # mdamath.triclinic_vectors explicitly sets the off-diagonal
                    # elements to zero if the box is orthogonal, so we can
                    # safely check floating point values for equality here
                    if box[i, j] != 0.0:
                        self.is_triclinic = True

        # Update diagonals
        for i in range(DIM):
            self.c_pbcbox.fbox_diag[i] = box[i, i]
            self.c_pbcbox.hbox_diag[i] = self.c_pbcbox.fbox_diag[i] * 0.5
            self.c_pbcbox.mhbox_diag[i] = - self.c_pbcbox.hbox_diag[i]

        # Update maximum cutoff

        # Physical limitation of the cut-off
        # by half the length of the shortest box vector.
        min_hv2 = min(0.25 * rvec_norm2(&box[XX, XX]), 0.25 * rvec_norm2(&box[YY, XX]))
        min_hv2 = min(min_hv2, 0.25 * rvec_norm2(&box[ZZ, XX]))

        # Limitation to the smallest diagonal element due to optimizations:
        # checking only linear combinations of single box-vectors (2 in x)
        # in the grid search and pbc_dx is a lot faster
        # than checking all possible combinations.
        tmp = box[YY, YY]
        if box[ZZ, YY] < 0:
            tmp -= box[ZZ, YY]
        else:
            tmp += box[ZZ, YY]

        min_ss = min(box[XX, XX], min(tmp, box[ZZ, ZZ]))
        self.c_pbcbox.max_cutoff2 = min(min_hv2, min_ss * min_ss)

    def update(self, real[:, ::1] box):
        """
        Updates internal MD box representation and parameters used for calculations.

        Parameters
        ----------
        box : numpy.ndarray
            Describes the MD box vectors as returned by
            :func:`MDAnalysis.lib.mdamath.triclinic_vectors`.
            `dtype` must be :class:`numpy.float32`

        Note
        ----
        Call to this method is only needed when the MD box is changed
        as it always called when class is instantiated.

        """

        if box.shape[0] != DIM or box.shape[1] != DIM:
            raise ValueError("Box must be a {} x {} matrix. Got: {} x {})".format(
                DIM, DIM, box.shape[0], box.shape[1]))
        if (box[XX, XX] < EPSILON) or (box[YY, YY] < EPSILON) or (box[ZZ, ZZ] < EPSILON):
            raise ValueError("Box does not correspond to PBC=xyz")
        self.fast_update(box)

    cdef void fast_pbc_dx(self, rvec ref, rvec other, rvec dx) nogil:
        """Dislacement between two points for both
        PBC and non-PBC conditions

        Modifies the displacement vector between two points based
        on the minimum image convention for PBC aware calculations.

        For non-PBC aware distance evaluations, calculates the
        displacement vector without any modifications
        """

        cdef fsl_int i, j

        for i in range(DIM):
            dx[i] = other[i] - ref[i]

        for i in range(DIM-1, -1, -1):
            while dx[i] > self.c_pbcbox.hbox_diag[i]:
                for j in range(i, -1, -1):
                    dx[j] -= self.c_pbcbox.box[i][j]

            while dx[i] <= self.c_pbcbox.mhbox_diag[i]:
                for j in range(i, -1, -1):
                    dx[j] += self.c_pbcbox.box[i][j]

    cdef dreal fast_distance2(self, rvec a, rvec b) nogil:
        """Distance calculation between two points
        for both PBC and non-PBC aware calculations

        Returns the distance obeying minimum
        image convention if periodic is set to ``True`` while
        instantiating the :class:`_PBCBox` object.
        """

        cdef rvec dx
        self.fast_pbc_dx(a, b, dx)
        return rvec_norm2(dx)

    cdef void fast_put_atoms_in_bbox(self, real[:, ::1] coords, real[:, ::1] bbox_coords) nogil:
        """Shifts all ``coords`` to an orthogonal brick shaped box

        All the coordinates are brought into an orthogonal
        box. The box vectors for the brick-shaped box
        are defined in ``fast_update`` method.

        """

        cdef fsl_int i, m, d, natoms

        natoms = coords.shape[0]

        if self.is_triclinic:
            for i in range(natoms):
                for m in range(DIM - 1, -1, -1):
                    while bbox_coords[i, m] < 0:
                        for d in range(m+1):
                            bbox_coords[i, d] += self.c_pbcbox.box[m][d]
                    while bbox_coords[i, m] >= self.c_pbcbox.box[m][m]:
                        for d in range(m+1):
                            bbox_coords[i, d] -= self.c_pbcbox.box[m][d]
        else:
            for i in range(natoms):
                for m in range(DIM):
                    while bbox_coords[i, m] < 0:
                        bbox_coords[i, m] += self.c_pbcbox.box[m][m]
                    while bbox_coords[i, m] >= self.c_pbcbox.box[m][m]:
                        bbox_coords[i, m] -= self.c_pbcbox.box[m][m]


    cdef void fast_pbc_xcm(self, real[:, ::1] coords, rvec xcm, fsl_int[:] indices) nogil:
        if indices.shape[0] < 1:
            return
        self.fast_pbc_xcm_from_ref(coords, &coords[indices[0], XX], xcm, indices)

    @cython.cdivision(True)
    cdef void fast_pbc_xcm_from_ref(self, real[:, ::1] coords, rvec ref,
                                    rvec xcm, fsl_int[:] indices) nogil:
        cdef fsl_int i
        cdef fsl_int size
        cdef rvec dx
        cdef real[:,::1] bbox_coords
        cdef fsl_int actual_index

        size = indices.shape[0]

        if size < 1:
            return

        rvec_copy(ref, xcm)
        rvec_smul(size, xcm, xcm)

        # Step 1: Get XCM
        for i in range(size):
            actual_index = indices[i]
            self.fast_pbc_dx(ref, &coords[actual_index, XX], dx)
            rvec_inc(xcm, dx)
        rvec_smul(1.0/size, xcm, xcm)

        # Step 2: Make sure it is inside the brick-shaped box
        for i in range(DIM - 1, -1, -1):
            while xcm[i] < 0:
                for d in range(i+1):
                    xcm[d] += self.c_pbcbox.box[i][d]
            while xcm[i] >= self.c_pbcbox.box[i][i]:
                for d in range(i+1):
                    xcm[d] -= self.c_pbcbox.box[i][d]

@cython.initializedcheck(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef bint normal_from_neighbours(real[:, ::1] positions,
                                 real[:, ::1] directions,
                                 fsl_int refid,
                                 fsl_int[:] neighbours_ids,
                                 PBCBox box,
                                 rvec normal) nogil:
    cdef fsl_int i, nid
    cdef fsl_int neighborhood_size = neighbours_ids.shape[0], useful_size = 0
    cdef rvec dx, xcm
    cdef rvec neighborhood_direction
    cdef rvec ref_direction
    cdef matrix cov_mat
    cdef rvec tmp_vec
    cdef rvec eig_vals
    cdef matrix eig_vecs
    #cdef fsl_int[:] useful_neighbors = np.empty_like(neighbours_ids)
    cdef fsl_int uselful_size

    # Don't compute anything if less than 3 neighbors
    # Assume the normal as the reverse of the lipid directions
    if neighborhood_size < 3:
        normal[XX] = directions[refid, XX]
        normal[YY] = directions[refid, YY]
        normal[ZZ] = directions[refid, ZZ]
        return False

    # Compute center of masses
    rvec_clear(xcm)
    rvec_clear(neighborhood_direction)

    rvec_copy(&directions[refid, XX], ref_direction)

    for i in range(neighborhood_size):
        nid = neighbours_ids[i]

        # Only take into account the neighbors that are oriented toward the same direction (< 90°)
        if rvec_dprod(ref_direction, &directions[nid, XX]) > 0:
            #useful_neighbors[useful_size] = nid
            useful_size += 1

            # HG xcm
            box.fast_pbc_dx(&positions[refid, XX], &positions[nid, XX], dx)
            rvec_inc(xcm, dx)

            rvec_inc(neighborhood_direction, &directions[nid, XX])

    # Don't compute anything if less than 3 useful neighbors
    # Instead, set it to notset as the flag for further computations
    if useful_size < 3:
        normal[XX] = ref_direction[XX]
        normal[YY] = ref_direction[YY]
        normal[ZZ] = ref_direction[ZZ]
        return False

    rvec_smul(1.0/useful_size, xcm, xcm)
    rvec_inc(xcm, &positions[refid, XX])
    rvec_smul(1.0/useful_size, neighborhood_direction, neighborhood_direction)

    # Build covariance matrix
    mat_clear(cov_mat)
    for i in range(neighborhood_size):
        nid = neighbours_ids[i]
        if rvec_dprod(ref_direction, &directions[nid, XX]) > 0:

            # Retrieve neighbor and its image
            box.fast_pbc_dx(xcm, &positions[nid, XX], tmp_vec)

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

    return True

@cython.cdivision(True)
cdef void eigen_33_sym(matrix a, matrix eig_vec, rvec eig_val) nogil:
    cdef matrix scaled_mat, m0, m1
    cdef real max_val
    cdef rvec roots
    cdef int i, j, rank0, rank1

    # scale the matrix
    max_val = a[XX][XX]
    for i in range(DIM):
        for j in range(DIM):
            if fabs(a[i][j]) > max_val:
                max_val = fabs(a[i][j])

    for i in range(DIM):
        for j in range(DIM):
          scaled_mat[i][j] = a[i][j] / max_val


    # Get the eigenvalues
    compute_roots(scaled_mat, roots)

    eig_val[XX] = roots[XX]
    eig_val[YY] = roots[YY]
    eig_val[ZZ] = roots[ZZ]

    # Compute A - eig_val[i] * I
    mat_copy(scaled_mat, m0)
    m0[XX][XX] -= eig_val[0]
    m0[YY][ZZ] -= eig_val[0]
    m0[ZZ][ZZ] -= eig_val[0]

    rank0 = compute_rank(m0)

    if rank0 == 0:  # eigenvalues are identical
        rvec_clear(eig_vec[XX])
        eig_vec[XX][XX] = 1.0
        rvec_clear(eig_vec[YY])
        eig_vec[YY][YY] = 1.0
        rvec_clear(eig_vec[ZZ])
        eig_vec[ZZ][ZZ] = 1.0
        return
    elif rank0 == 1:  #eig_val[0] = eig_val[1] < eig_val[2]
        complete_basis(m0[0], eig_vec[XX],eig_vec[YY])
        rvec_cprod(eig_vec[XX], eig_vec[YY], eig_vec[ZZ])
        return
    else:  # rank0 == 2
        rvec_cprod_norm(m0[0], m0[1], eig_vec[XX])

    # Compute A - eig_val[1] * I
    mat_copy(scaled_mat, m1)
    m1[XX][XX] -= eig_val[1]
    m1[YY][ZZ] -= eig_val[1]
    m1[ZZ][ZZ] -= eig_val[1]

    rank1 = compute_rank(m1)  # zero rank detected earlier, rank1 must be positive
    if rank1 == 1:  #eig_val[0] < eig_val[1] = eig_val[2]
        complete_basis(eig_vec[XX], eig_vec[YY], eig_vec[ZZ])
        return
    else:  # rank1 == 2
        rvec_cprod_norm(m1[0], m1[1], eig_vec[YY])

    # eigenvalues must be distinct at this point, rank2 must be 2
    rvec_cprod(eig_vec[XX], eig_vec[YY], eig_vec[ZZ])


cdef void complete_basis(rvec z_axis, rvec x_axis, rvec y_axis) nogil:

    rvec_normalize(z_axis)

    if fabs(z_axis[XX]) >= fabs(z_axis[YY]) :
        x_axis[XX] = -z_axis[ZZ]
        x_axis[YY] = 0
        x_axis[ZZ] = z_axis[XX]
    else:
        x_axis[XX] = 0
        x_axis[YY] = z_axis[ZZ]
        x_axis[ZZ] = -z_axis[YY]

    rvec_normalize(x_axis)
    rvec_cprod(z_axis, x_axis, y_axis)


cdef void compute_roots(matrix mat, rvec roots) nogil:
    cdef real c0, c1, c2
    cdef real one_over_3 = <real> 1 / <real> 3
    cdef real sqrt_3 = sqrt(3)
    cdef real c2_over_3, a_over_3
    cdef real half_mb, q, magnitude, angle, cs, sn, root0, root1, root2
    cdef real mat_XX = mat[XX][XX]
    cdef real mat_XY = mat[XX][YY]
    cdef real mat_XZ = mat[XX][ZZ]
    cdef real mat_YY = mat[YY][YY]
    cdef real mat_YZ = mat[YY][ZZ]
    cdef real mat_ZZ = mat[ZZ][ZZ]

    #
    # The characteristic equation is x^3 - c2*x^2 + c1*x - c0 = 0.  The
    # eigenvalues are the roots to this equation, all guaranteed to be
    # real-valued, because the matrix is symmetric.
    #
    c0 = mat_XX * mat_YY * mat_ZZ + 2.0 * mat_XY * mat_XZ * mat_YZ - mat_XX * mat_YZ * mat_YZ - \
            mat_YY * mat_XZ * mat_XZ - mat_ZZ * mat_XY * mat_XY

    c1 = mat_XX * mat_YY - mat_XY * mat_XY + mat_XX * mat_ZZ - mat_XZ * mat_XZ +\
            mat_YY * mat_ZZ - mat_YZ * mat_YZ

    c2 = mat_XX + mat_YY + mat_ZZ

    #
    # Construct the parameters used in classifying the roots of the equation
    # and in solving the equation for the roots in closed form.
    #
    c2_over_3 = c2 * one_over_3
    a_over_3 = (c1 - c2 * c2_over_3) * one_over_3

    if a_over_3 > 0.0:
        a_over_3 = 0.0

    half_mb = 0.5 * (c0 + c2_over_3 * (2.0 * c2_over_3 * c2_over_3 - c1))

    q = half_mb * half_mb + a_over_3 * a_over_3 * a_over_3

    if q > 0.0:
        q = 0.0

    #
    # Compute the eigenvalues by solving for the roots of the polynomial.
    #

    magnitude = sqrt( - a_over_3)
    angle = atan2(sqrt(-q), half_mb) * one_over_3
    cs = cos(angle)
    sn = sin(angle)

    root0 = c2_over_3 + 2.0 * magnitude * cs
    root1 = c2_over_3 - magnitude * (cs + sqrt_3 * sn)
    root2 = c2_over_3 - magnitude * (cs - sqrt_3 * sn)


    #
    # Sort roots so root0 < root1 < root2
    #
    if root1 >= root0:
        roots[0] = root0
        roots[1] = root1
    else:
        roots[0] = root1
        roots[1] = root0

    if root2 >= roots[1]:
        roots[2] = root2
    else:
        roots[2] = roots[1]
        if root2 >= roots[0]:
            roots[1] = root2
        else:
            roots[1] = roots[0]
            roots[0] = root2

@cython.cdivision(True)
cdef int compute_rank(matrix m) nogil:
    # Compute the maximum magnitude matrix entry.
    cdef dreal abs, save, max = -1, inv_max
    cdef int row, col, maxRow = -1, maxCol = -1

    for row in range(DIM):
        for col in range(row, DIM):
            abs = fabs(m[row][col])
            if abs > max:
                max = abs
                maxRow = row
                maxCol = col

    if max < EPSILON:
        #
        # The rank is 0. The eigenvalue has multiplicity 3.
        #
        return 0

    #
    # The rank is at least 1. Swap the row containing the maximum-magnitude
    # entry with row 0.
    #

    if maxRow != 0:
        for col in range(DIM):
            save = m[0][col]
            m[XX][col] = m[maxRow][col]
            m[maxRow][col] = save

    #
    # Row-reduce the matrix...
    # Scale row 0 to generate a 1-valued pivot.
    #
    inv_max = 1/m[XX][maxCol]
    m[XX][XX] *= inv_max
    m[XX][YY] *= inv_max
    m[XX][ZZ] *= inv_max

    #
    # Eliminate the maxCol column entries in rows 1 and 2.
    #
    if maxCol == XX:
        m[YY][YY] -= m[YY][XX] * m[XX][YY]
        m[YY][ZZ] -= m[YY][XX] * m[XX][ZZ]
        m[ZZ][YY] -= m[ZZ][XX] * m[XX][YY]
        m[ZZ][ZZ] -= m[ZZ][XX] * m[XX][ZZ]
        m[YY][XX] = 0
        m[ZZ][XX] = 0
    elif maxCol == YY:
        m[YY][XX] -= m[YY][YY] * m[XX][XX]
        m[YY][ZZ] -= m[YY][YY] * m[XX][ZZ]
        m[ZZ][XX] -= m[ZZ][YY] * m[XX][XX]
        m[ZZ][ZZ] -= m[ZZ][YY] * m[XX][ZZ]
        m[YY][YY] = 0
        m[ZZ][YY] = 0
    else:
        m[YY][XX] -= m[YY][ZZ] * m[XX][XX]
        m[YY][YY] -= m[YY][ZZ] * m[XX][YY]
        m[ZZ][XX] -= m[ZZ][ZZ] * m[XX][XX]
        m[ZZ][YY] -= m[ZZ][ZZ] * m[XX][YY]
        m[YY][ZZ] = 0
        m[ZZ][ZZ] = 0


    #
    # Compute the maximum-magnitude entry of the last two rows of the
    # row-reduced matrix.
    #
    max = -1
    maxRow = -1
    maxCol = -1
    for row in range(1, DIM):
        for col in range(DIM):
            abs = fabs(m[row][col])
            if abs > max:
                max = abs
                maxRow = row
                maxCol = col

    if max == XX:
        #
        # The rank is 1. The eigenvalue has multiplicity 2.
        #
        return 1

    #
    # If row 2 has the maximum-magnitude entry, swap it with row 1.
    #
    if maxRow == ZZ:
        for col in range(DIM):
            save = m[YY][col]
            m[YY][col] = m[ZZ][col]
            m[ZZ][col] = save

    #
    # Scale row 1 to generate a 1-valued pivot.
    #
    inv_max = 1 / m[YY][maxCol]
    m[YY][XX] *= inv_max
    m[YY][YY] *= inv_max
    m[YY][ZZ] *= inv_max

    #
    # The rank is 2. The eigenvalue has multiplicity 1.
    #
    return 2