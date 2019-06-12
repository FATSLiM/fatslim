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

from libc.math cimport sqrt

# Cython directives
# cython: language_level=3

DEF XX=0
DEF YY=1
DEF ZZ=2

cdef void rvec_clear(rvec a) nogil:
    a[XX] = 0
    a[YY] = 0
    a[ZZ] = 0

cdef real rvec_norm2(const rvec a) nogil:
    return a[XX]*a[XX] + a[YY]*a[YY] + a[ZZ]*a[ZZ]

cdef void rvec_copy(rvec src, rvec dest) nogil:
    dest[XX]=src[XX]
    dest[YY]=src[YY]
    dest[ZZ]=src[ZZ]

cdef void rvec_smul(real a, const rvec v1, rvec v2) nogil:
    v2[XX] = a * v1[XX]
    v2[YY] = a * v1[YY]
    v2[ZZ] = a * v1[ZZ]

cdef void rvec_inc(rvec a,const rvec b) nogil:
    cdef real x,y,z

    x=a[XX]+b[XX]
    y=a[YY]+b[YY]
    z=a[ZZ]+b[ZZ]

    a[XX]=x
    a[YY]=y
    a[ZZ]=z

cdef void rvec_dec(rvec a,const rvec b) nogil:
    cdef real x, y, z

    x=a[XX]-b[XX]
    y=a[YY]-b[YY]
    z=a[ZZ]-b[ZZ]

    a[XX]=x
    a[YY]=y
    a[ZZ]=z

cdef real rvec_norm(const rvec a) nogil:
    return sqrt(rvec_norm2(a))

cdef void rvec_normalize(rvec a) nogil:
    cdef real vec_norm = rvec_norm(a)
    a[XX] /= vec_norm
    a[YY] /= vec_norm
    a[ZZ] /= vec_norm

cdef real rvec_dprod(const rvec a, const rvec b) nogil:
    return a[XX] * b[XX] + a[YY] * b[YY] + a[ZZ] * b[ZZ]

cdef void mat_clear(matrix a) nogil:
    rvec_clear(a[XX])
    rvec_clear(a[YY])
    rvec_clear(a[ZZ])

cdef void mat_copy(matrix src,matrix dest) nogil:
    rvec_copy(src[XX],dest[XX])
    rvec_copy(src[YY],dest[YY])
    rvec_copy(src[ZZ],dest[ZZ])


cdef void rvec_cprod(const rvec a, const rvec b, rvec c) nogil:
    c[XX] = a[YY] * b[ZZ] - a[ZZ] * b[YY]
    c[YY] = a[ZZ] * b[XX] - a[XX] * b[ZZ]
    c[ZZ] = a[XX] * b[YY] - a[YY] * b[XX]


cdef void rvec_cprod_norm(const rvec a, const rvec b, rvec c) nogil:
    rvec_cprod(a, b, c)
    rvec_normalize(c)


cdef void invert_mat(matrix mat, matrix inverted_mat) nogil:
    cdef real det
    cdef real A, B, C, D, E, F, G, H, K
    cdef real a, b, c, d, e, f, g, h, k

    a = mat[XX][XX]
    b = mat[XX][YY]
    c = mat[XX][ZZ]
    d = mat[YY][XX]
    e = mat[YY][YY]
    f = mat[YY][ZZ]
    g = mat[ZZ][XX]
    h = mat[ZZ][YY]
    k = mat[ZZ][ZZ]

    A = e * k - f * h  # ek-fh
    B = f * g - d * k  # fg-dk
    C = d * h - e * g  # dh-eg
    D = c * h - b * k  # ch-bk
    E = a * k - c * g  # ak-cg
    F = g * b - a * h  # gb-ah
    G = b * f - c * e  # bf-ce
    H = c * d - a * f  # cd-af
    K = a * e - b * d  # ae-bd

    det = a * A + b * B + c * C
    det = 1 / det  # Beware: determinant MUST not be 0!

    inverted_mat[XX][XX] = det * A
    inverted_mat[XX][YY] = det * D
    inverted_mat[XX][ZZ] = det * G
    inverted_mat[YY][XX] = det * B
    inverted_mat[YY][YY] = det * E
    inverted_mat[YY][ZZ] = det * H
    inverted_mat[ZZ][XX] = det * C
    inverted_mat[ZZ][YY] = det * F
    inverted_mat[ZZ][ZZ] = det * K


cdef void mat_from_rvec(rvec v1, rvec v2, rvec v3, matrix mat) nogil:
    mat[XX][XX] = v1[XX]
    mat[XX][YY] = v2[XX]
    mat[XX][ZZ] = v3[XX]
    mat[YY][XX] = v1[YY]
    mat[YY][YY] = v2[YY]
    mat[YY][ZZ] = v3[YY]
    mat[ZZ][XX] = v1[ZZ]
    mat[ZZ][YY] = v2[ZZ]
    mat[ZZ][ZZ] = v3[ZZ]
