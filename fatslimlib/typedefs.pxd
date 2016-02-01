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

cdef extern from "typedefs.h" nogil:
    # Types
    ctypedef int fsl_int
    ctypedef unsigned int fsl_uint
    ctypedef double real
    ctypedef real rvec[DIM]
    ctypedef real matrix[DIM][DIM]

    # Functions
    cdef real real_min(real a, real b)
    cdef real real_max(real a, real b)
    cdef real real_abs(real a)
    cdef real rvec_norm2(const rvec a)
    cdef real rvec_norm(const rvec a)
    cdef void rvec_normalize(rvec a)
    cdef void rvec_add(const rvec a,const rvec b,rvec c)
    cdef void rvec_inc(rvec a,const rvec b)
    cdef void rvec_dec(rvec a,const rvec b)
    cdef void rvec_sub(const rvec a,const rvec b,rvec c)
    cdef void rvec_cprod(const rvec a, const rvec b, rvec c)
    cdef void rvec_cprod_norm(const rvec a, const rvec b, rvec c)
    cdef real rvec_dprod(const rvec a, const rvec b)
    cdef void rvec_copy(const rvec src,rvec dest)
    cdef void rvec_clear(rvec a)
    cdef void mat_clear(matrix a)
    cdef void mat_copy(matrix src,matrix dest)
    cdef void rvec_smul(real a,const rvec v1,rvec v2)
    cdef void invert_mat(matrix mat, matrix inverted_mat)
    cdef void mat_from_rvec(rvec v1, rvec v2, rvec v3, matrix mat)
