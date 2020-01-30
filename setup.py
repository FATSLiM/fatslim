#!/usr/bin/env python
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
"""
FATSLiM setup script
"""
from __future__ import print_function
from os import chdir, environ
from os.path import realpath, dirname
import sys
if sys.version_info.major < 3:
    print("ERROR: FATSLiM needs python 3. Python 2 is not supported anymore.")
    sys.exit(1)

from setuptools import setup, find_packages
from setuptools import Extension
from setuptools.command.test import test as TestCommand
from io import StringIO  # python 3
from contextlib import contextmanager

@contextmanager
def std_redirector(stream):
    old_stderr = sys.stderr
    old_stdout = sys.stdout
    sys.stderr = stream
    sys.stdout = stream
    try:
        yield
    finally:
        sys.stderr = old_stderr
        sys.stdout = old_stdout

chdir(realpath(dirname(__file__)))

import numpy

include_dirs = [numpy.get_include()]

# Handle cython modules
try:
    import Cython
    from distutils.version import LooseVersion
    req_version = "0.29"

    if not LooseVersion(Cython.__version__) >= LooseVersion(req_version):
        print("INFO: Cython is present but installed version (%s) "
              "does not support parallelisation." % Cython.__version__)
        del Cython
        del LooseVersion
        raise ImportError
    from Cython.Build import cythonize

    use_cython = True
    print("INFO: Cython will be used.")
except ImportError:
    use_cython = False

    print("INFO: Cython will NOT be used.")

    def cythonize(obj, **kwargs):
        return obj


# Check for OpenMP support and get the compilation flags
# Derived from: http://stackoverflow.com/questions/16549893/
def get_compilation_flags():
    import os
    import tempfile
    from distutils.ccompiler import new_compiler, CompileError, LinkError
    import shutil
    print("INFO: Checking for OpenMP support...")
    cc = new_compiler()
    if cc.compiler_type == "msvc":
        compile_flags = ["/openmp"]
        link_flags = []
    else:
        compile_flags = ["-fopenmp"]
        link_flags = ["-fopenmp"]

    # see http://openmp.org/wp/openmp-compilers/
    omp_test = \
    b"""
    #include <omp.h>\n
    #include <stdio.h>\n
    int main() {\n
    #pragma omp parallel\n
    printf("Hello from thread %d, nthreads %d\\n", omp_get_thread_num(), omp_get_num_threads());\n
    }\n
    """
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    # Write test.c
    filename = 'test.c'
    with open(filename, 'wb') as fp:
        fp.write(omp_test)

    try:
        buf = StringIO()
        with std_redirector(buf):
            objects = cc.compile([filename], output_dir=tmpdir,
                                 extra_postargs=compile_flags)
            cc.link_executable(objects, os.path.join(tmpdir, "a.out"),
                               extra_postargs=link_flags)
    except (CompileError, LinkError):
        print("WARNING: Could not find OpenMP, parallelism will "
              "not be available!")
        compile_flags = []
        link_flags = []
    else:
        print("INFO: OpenMP is present and will be used!")

    os.chdir(curdir)

    # Clean up
    shutil.rmtree(tmpdir)

    # Add generic flags
    if cc.compiler_type == "msvc":
        compile_flags += ["/wd4244", "/wd4267", "/wd4018", "/wd4996"]
    else:
        compile_flags += ["-Wno-maybe-uninitialized", "-Wno-unused-function", "-Wno-cpp",
                          "-Wno-shorten-64-to-32", "-Wno-unneeded-internal-declaration"]

    return compile_flags, link_flags

cmdclass = {}

# Test command
class PyTest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        # Default py.test options (overriden if --pytest-args/-a is set)
        self.pytest_args = ["-v", "-x", "--durations=3", "--tb=short"]

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        # import here, cause outside the eggs aren't loaded
        import pytest
        print("py.test args: '%s'" % self.pytest_args)
        errno = pytest.main(self.pytest_args)
        if errno != 0:
            print("Ouch... Congratulations! you probably found a bug!")
            print("Please contact the FATSLiM devs with this output so they can work it out!")
        sys.exit(errno)

cmdclass['test'] = PyTest


# ===============================================================================
# Common extension setup
# ===============================================================================

# List of common include dirs
include_dirs += ["./include", "./src", ".", "./fatslimlib"]

# List of default extra compile args
extra_compile_args,  extra_link_args = get_compilation_flags()




debug_mode = (environ.get('FATSLIM_DEBUG') is not None) or ("-g" in sys.argv) or \
             ("--debug" in sys.argv)
if debug_mode:
    print("NOTE: Debug flag is set for compilation")
    extra_compile_args += ["-O0", "-DDEBUG", "-Wcpp", "-g"]
    try:
        sys.argv.remove("--debug")
    except ValueError:
        pass

    try:
        sys.argv.remove("-g")
    except ValueError:
        pass

# List of source files
core_base = Extension("fatslimlib.core_base",
                      ['fatslimlib/core_base.%s' % ("pyx" if use_cython else "c"),
                       'src/eig_mat.c',
                       'src/typeutil.c'],
                      include_dirs=include_dirs,
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args)

core_datareading = Extension("fatslimlib.core_datareading", ['fatslimlib/core_datareading.%s' %
                                                             ("pyx" if use_cython else "c"),
                                                             'src/xdrfile/xdrfile.c',
                                                             'src/typeutil.c'],
                             include_dirs=include_dirs,
                             extra_compile_args=extra_compile_args,
                             extra_link_args=extra_link_args)

core_geometry = Extension("fatslimlib.core_geometry",
                          ['fatslimlib/core_geometry.%s' % ("pyx" if use_cython else "c"),
                           'src/typeutil.c'],
                          include_dirs=include_dirs,
                          extra_compile_args=extra_compile_args,
                          extra_link_args=extra_link_args)

core_ns = Extension("fatslimlib.core_ns",
                    ['fatslimlib/core_ns.%s' % ("pyx" if use_cython else "c"),
                     'src/typeutil.c'],
                    include_dirs=include_dirs,
                    extra_compile_args=extra_compile_args,
                    extra_link_args=extra_link_args)

core_analysis = Extension("fatslimlib.core_analysis", ['fatslimlib/core_analysis.%s' %
                                                       ("pyx" if use_cython else "c"),
                                                       'src/typeutil.c'],
                          include_dirs=include_dirs,
                          extra_compile_args=extra_compile_args,
                          extra_link_args=extra_link_args)

# ===============================================================================
# Run setup
# ===============================================================================
if __name__ == "__main__":
    from fatslimlib import __version__, __full_desc__, __url__

    setup(name='fatslim',
          version=__version__,
          description='Python/MD Toolbox',
          author='Sébastien Buchoux',
          author_email='sebastien.buchoux@gmail.com',
          url=__url__,
          long_description=__full_desc__,
          classifiers=[
              "Intended Audience :: Science/Research",
              "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
              "Operating System :: MacOS :: MacOS X",
              "Operating System :: Microsoft :: Windows",
              "Operating System :: POSIX",
              "Programming Language :: C",
              "Programming Language :: Cython",
              "Programming Language :: Python",
              "Topic :: Scientific/Engineering :: Bio-Informatics",
              "Topic :: Scientific/Engineering :: Chemistry",
                       ],
          packages=find_packages(),
          cmdclass=cmdclass,
          include_package_data=True,
          package_data={"fatslimlib": ["test/data/*.gro", "test/data/*.ndx", "test/data/*.trr",
                                       "test/data/*.xtc", "test/data/*.csv", "test/data/*.xvg"]},
          data_files=[("share/fatslim", ["extra/fatslim-completion.bash"])],
          ext_modules=cythonize([core_base,
                                 core_datareading,
                                 core_geometry,
                                 core_ns,
                                 core_analysis],
                                compiler_directives={"language_level": "3"}),
          install_requires=["numpy>=1.5"],
          tests_require=['pytest'],
          scripts=['fatslim'])
