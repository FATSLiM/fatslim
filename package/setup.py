#!/usr/bin/env python
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
import sys
import os
from setuptools import setup
from setuptools import Extension
from setuptools.command.test import test as TestCommand

# Handle cython modules
from Cython.Build import cythonize

import numpy

cython_coverage = 0
force_cythonize = False
cython_linetrace = False
if os.environ.get("CYTHON_COVERAGE") in ("1", 1, True):
    print("INFO: Coverage is enabled for Cython code (Expect low performances)")
    cython_linetrace = True
    force_cythonize = True
    cython_coverage = 1

if __name__ == "__main__":
    # MUST match fatslim.__version__
    VERSION = "0.99.0-dev"

    # Extensions
    typedefs_ext = Extension("fatslim._typedefs",
                             ["fatslim/_typedefs.pyx"],
                             include_dirs=[numpy.get_include()],
                             define_macros=[('CYTHON_TRACE_NOGIL', cython_coverage)]
                             )

    geometry_ext = Extension("fatslim._geometry",
                             ["fatslim/_geometry.pyx"],
                             include_dirs=[numpy.get_include()],
                             define_macros=[('CYTHON_TRACE_NOGIL', cython_coverage)]
                             )

    core_ext = Extension("fatslim._core",
                         ["fatslim/_core.pyx"],
                         include_dirs=[numpy.get_include()],
                         define_macros=[('CYTHON_TRACE_NOGIL', cython_coverage)]
                         )

    aggregate_ext = Extension("fatslim._aggregate",
                              ["fatslim/_aggregate.pyx"],
                              include_dirs=[numpy.get_include()],
                              define_macros=[('CYTHON_TRACE_NOGIL', cython_coverage)]
                              )

    membrane_ext = Extension("fatslim._membrane",
                             ["fatslim/_membrane.pyx"],
                             include_dirs=[numpy.get_include()],
                             define_macros=[('CYTHON_TRACE_NOGIL', cython_coverage)]
                             )

    setup(
        name='fatslim',
        version=VERSION,
        ext_modules=cythonize(
            [
                typedefs_ext,
                geometry_ext,
                core_ext,
                aggregate_ext,
                membrane_ext
            ],
            force=force_cythonize,
            compiler_directives={'linetrace': cython_linetrace, 'binding': True},
        ),
    )
