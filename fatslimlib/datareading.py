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
from __future__ import print_function

# Global imports
import os

# Local imports
from .core_base import Trajectory
from .core_datareading import _Topology_readers, _Coordinates_readers, _Index_loaders


class UnknownFileType(Exception):
    pass


def get_readers(top_file, index_file, coords_file=None, verbose=True):

    # Handle topology file
    top_ext = os.path.splitext(top_file)[1].lower()
    try:
        top_reader = _Topology_readers[top_ext]
    except KeyError:
        raise UnknownFileType("'%s' is not a topology file. Known types: %s" % (top_file,
                                                                                ", ".join(_Topology_readers.keys())))

    # Handle index file
    index_ext = os.path.splitext(index_file)[1].lower()
    try:
        index_loader = _Index_loaders[index_ext]
    except KeyError:
        raise UnknownFileType("'%s' is not an index file. Known types: %s" % (index_file,
                                                                              ", ".join(_Index_loaders.keys())))

    # Handle and load coordinates file
    if coords_file is None:
        coords_file = top_file
    elif not os.path.isfile(coords_file):
        coords_file = top_file
    coords_ext = os.path.splitext(coords_file)[1].lower()
    try:
        coords_reader = _Coordinates_readers[coords_ext]
    except KeyError:
        raise UnknownFileType("'%s' is not a coordinates file. Known types: %s" % (coords_file,
                                                                                   ", ".join(_Coordinates_readers.keys())))

    return \
        top_reader(top_file, verbose=verbose), \
        index_loader(index_file, verbose=verbose), \
        coords_reader(coords_file, verbose=verbose)


def load_trajectory(top_file, index_file, coords_file=None, verbose=True):
    readers = get_readers(top_file, index_file, coords_file, verbose=verbose)

    return Trajectory(*readers, be_verbose=verbose)
