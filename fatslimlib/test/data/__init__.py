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

import os


def _fullpath_from_basename(fname):
    return os.path.join(DATADIR, fname)

DATADIR = "%s" % os.path.dirname(__file__)

# Model bilayer
MODEL_BILAYER_GRO = _fullpath_from_basename("model_bilayer.gro")
MODEL_BILAYER_NDX = _fullpath_from_basename("model_bilayer.ndx")
MODEL_BILAYER_THICKNESS_XVG = _fullpath_from_basename("model_bilayer_thickness.xvg")
MODEL_BILAYER_THICKNESS_CSV = _fullpath_from_basename("model_bilayer_thickness_raw.csv")
MODEL_BILAYER_APL_CSV = _fullpath_from_basename("model_bilayer_apl_raw.csv")


# Simple (yet experimental) bilayer
BILAYER_GRO = _fullpath_from_basename("bilayer.gro")
BILAYER_NDX = _fullpath_from_basename("bilayer.ndx")

# All-atom bilayer
BILAYER_ALLATOM_GRO = _fullpath_from_basename("bilayer_allatom.gro")
BILAYER_ALLATOM_NDX = _fullpath_from_basename("bilayer_allatom.ndx")

# Simple bilayer with small peptide
BILAYER_PEPTIDE_GRO = _fullpath_from_basename("bilayer_peptide.gro")
BILAYER_PEPTIDE_NDX = _fullpath_from_basename("bilayer_peptide.ndx")
BILAYER_PEPTIDE_H_NDX = _fullpath_from_basename("bilayer_peptide_hydrogen.ndx")

# Bilayer with protein
MODEL_BILAYER_PROT_GRO = _fullpath_from_basename("model_bilayer_prot.gro")
MODEL_BILAYER_PROT_NDX = _fullpath_from_basename("model_bilayer_prot.ndx")
BILAYER_PROT_GRO = _fullpath_from_basename("bilayer_prot.gro")
BILAYER_PROT_NDX = _fullpath_from_basename("bilayer_prot.ndx")

# 2 Model bilayers
MODEL_MULTIBILAYER_GRO = _fullpath_from_basename("model_multibilayer.gro")
MODEL_MULTIBILAYER_NDX = _fullpath_from_basename("model_multibilayer.ndx")

# 2 Bilayers
MULTIBILAYER_GRO = _fullpath_from_basename("bilayer_multi.gro")
MULTIBILAYER_NDX = _fullpath_from_basename("bilayer_multi.ndx")

# Model vesicle
MODEL_BILAYER_VESICLE_GRO = _fullpath_from_basename("model_bilayer_vesicle.gro")
MODEL_BILAYER_VESICLE_NDX = _fullpath_from_basename("model_bilayer_vesicle.ndx")

# Bilayer with cholesterol
BILAYER_CHOL_GRO = _fullpath_from_basename("bilayer_chol.gro")
BILAYER_CHOL_NDX = _fullpath_from_basename("bilayer_chol.ndx")
BILAYER_CHOL_APL_GROUPED = _fullpath_from_basename("bilayer_chol_apl_group.xvg")
BILAYER_CHOL_APL_NOGROUP = _fullpath_from_basename("bilayer_chol_apl_nogroup.xvg")


# Big Bilayer-related variables
MODEL_BIG_GRO = _fullpath_from_basename("model_big.gro")
MODEL_BIG_NDX = _fullpath_from_basename("model_big.ndx")

# Big and tricky bilayer
BIG_PROT_GRO = _fullpath_from_basename("deformed.gro")
BIG_PROT_NDX = _fullpath_from_basename("deformed.ndx")

# Vesicle-related variables
MODEL_VESICLE_GRO = _fullpath_from_basename("model_vesicle.gro")
MODEL_VESICLE_NDX = _fullpath_from_basename("model_vesicle.ndx")
VESICLE_GRO = _fullpath_from_basename("dppc_vesicle.gro")
VESICLE_NDX = _fullpath_from_basename("dppc_vesicle.ndx")
VESICLE_XTC = _fullpath_from_basename("dppc_vesicle.xtc")
VESICLE_TRR = _fullpath_from_basename("dppc_vesicle.trr")
VESICLE_APL_XVG = _fullpath_from_basename("dppc_vesicle_apl.xvg")
VESICLE_AREA_XVG = _fullpath_from_basename("dppc_vesicle_area.xvg")
VESICLE_THICKNESS_XVG = _fullpath_from_basename("dppc_vesicle_thickness.xvg")



