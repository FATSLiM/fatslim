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
from setuptools.command.build_ext import build_ext as _build_ext

class build_ext(_build_ext): # subclass the build_ext class to defer numpy and cython import to prevent setup.py from failing if run on a clean system

    def finalize_options(self):
        super().finalize_options()
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        for extension in self.distribution.ext_modules:
            extension.include_dirs.append(numpy.get_include())
            extension.define_macros.append(
                ('CYTHON_TRACE_NOGIL', cython_coverage))
        from Cython.Build import cythonize
        self.distribution.ext_modules = cythonize(self.distribution.ext_modules,
                                                  language_level=3)

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
                             )

    geometry_ext = Extension("fatslim._geometry",
                             ["fatslim/_geometry.pyx"],
                             )

    core_ext = Extension("fatslim._core",
                         ["fatslim/_core.pyx"],
                         )

    aggregate_ext = Extension("fatslim._aggregate",
                              ["fatslim/_aggregate.pyx"],
                              )

    membrane_ext = Extension("fatslim._membrane",
                             ["fatslim/_membrane.pyx"],
                             )

    setup(
        name='fatslim',
        version=VERSION,
        packages=['fatslim'],
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
        setup_requires=["cython", "numpy", "mdanalysis"],
        install_requires=["numpy", "mdanalysis"],
        cmdclass={"build_ext": build_ext},
    )
