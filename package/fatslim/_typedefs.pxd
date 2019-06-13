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

# Cython directives
# cython: language_level=3


cimport numpy as np

DEF DIM=3

ctypedef np.int_t fsl_int
ctypedef np.float32_t real
ctypedef np.float64_t dreal
ctypedef real rvec[DIM]
ctypedef fsl_int ivec[DIM]
ctypedef real matrix[DIM][DIM]

cdef void rvec_clear(rvec a) nogil
cdef real rvec_norm2(const rvec a) nogil
cdef void rvec_copy(rvec src, rvec dest) nogil
cdef void rvec_smul(real a, const rvec v1, rvec v2) nogil
cdef void rvec_inc(rvec a,const rvec b) nogil
cdef void rvec_dec(rvec a,const rvec b) nogil
cdef real rvec_norm(const rvec a) nogil
cdef void rvec_normalize(rvec a) nogil
cdef real rvec_dprod(const rvec a, const rvec b) nogil
cdef void mat_clear(matrix a) nogil
cdef void mat_copy(matrix src,matrix dest) nogil
cdef void rvec_cprod(const rvec a, const rvec b, rvec c) nogil
cdef void rvec_cprod_norm(const rvec a, const rvec b, rvec c) nogil
cdef void invert_mat(matrix mat, matrix inverted_mat) nogil
cdef void mat_from_rvec(rvec v1, rvec v2, rvec v3, matrix mat) nogil