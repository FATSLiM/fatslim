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


cdef real rvec_norm(const rvec a) nogil:
    return sqrt(rvec_norm2(a))

cdef void rvec_normalize(rvec a) nogil:
    cdef real vec_norm = rvec_norm(a)
    a[XX] /= vec_norm
    a[YY] /= vec_norm
    a[ZZ] /= vec_norm