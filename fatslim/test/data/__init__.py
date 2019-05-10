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
import pickle


def _fullpath_from_basename(fname):
    return os.path.join(DATADIR, fname)


DATADIR = "%s" % os.path.dirname(__file__)

# Model systems
MODEL_FLAT_GRO = _fullpath_from_basename("model_flat.gro")  # Flat membrane
MODEL_FLAT_NDX = _fullpath_from_basename("model_flat.ndx")  # Flat membrane index file
MODEL_FLAT_BBOX_GRO = _fullpath_from_basename("model_flat_bbox.gro")  # Flat membrane in brick-shaped box
MODEL_VESICLE_GRO = _fullpath_from_basename("model_vesicle.gro")  # Vesicle
MODEL_BICELLE_GRO = _fullpath_from_basename("model_bicelle.gro")  # Bicelle
MODEL_CURVED_GRO = _fullpath_from_basename("model_curved.gro")  # Curved
MODEL_BULGED_GRO = _fullpath_from_basename("model_bulged.gro")  # Bulged

with open(_fullpath_from_basename("models_metadata.pickle"), "rb") as fp:
    MODELS_METADATA = pickle.load(fp)

# Real systems
BIG_DEFORMED_GRO = _fullpath_from_basename("deformed.gro")



