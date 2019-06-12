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
import os
from setuptools import setup, find_packages


if __name__ == "__main__":
    setup(
        name='fatslimtest',
        version="0.99.0-dev",
        packages=find_packages(),
        package_data={"fatslimtest":
                      [
                          "data/*.gro", "data/*.xtc", "data/*.ndx",
                          "data/*.pickle"
                      ],
                      },
    )
