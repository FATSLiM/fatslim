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

DEF DIM=3
DEF MAX_NTRICVEC=12
DEF NAME_SIZE = 6 # 5 chracters + NULL

from .typedefs cimport real, rvec, matrix, fsl_int
from .core_ns cimport ns_neighborhood_holder

# Types
ctypedef char[NAME_SIZE] strname_t

# Constants
cdef enum ForcefieldType:
    fft_notset,
    fft_allatom, fft_unified, fft_coarse,
    fft_unknown

# OpenMP-related stuff
cdef fsl_int OPENMP_MAX_THREADS
cdef fsl_int OPENMP_NUM_THREADS
cpdef fsl_int get_num_threads()
cpdef fsl_int get_max_threads()
cpdef fsl_int get_num_procs()
cpdef set_num_threads(fsl_int num)

# Utility functions
cdef real norm2(real[:] a) nogil
cdef void mat_to_memview(matrix src, real[:, ::1] dest) nogil
cdef void rvec_to_memview(rvec src, real[:] dest) nogil
cpdef pretty_delta(begin, end)
cpdef pretty_duration(duration)

# Other objects
cdef struct cAtom_t:
    rvec coords
    fsl_int atomid
    char *name

cdef class Atom:
    # Attributes
    cdef cAtom_t c_atom
    cdef bytes name
    cdef object __weakref__

cdef struct cAtomgroup_t:
    cAtom_t *atoms
    fsl_int n_atoms
    char *name
    fsl_int groupid
    rvec xcm

cdef class AtomGroup:
    # Attributes
    cdef cAtomgroup_t c_atomgroup
    cdef object __weakref__
    cdef bytes name

    # Methods
    cdef fast_append(self, Atom atom)

cdef struct cPBCBox_t:
    matrix     box
    rvec       fbox_diag
    rvec       hbox_diag
    rvec       mhbox_diag
    real       max_cutoff2
    fsl_int        ntric_vec
    fsl_int[DIM]   tric_shift[MAX_NTRICVEC]
    real[DIM]  tric_vec[MAX_NTRICVEC]

cdef class PBCBox:
    # Attributes
    cdef cPBCBox_t c_pbcbox
    cdef rvec center
    cdef rvec bbox_center

    # Methods
    cdef void fast_update(self, real[:,::1] box) nogil
    cdef void fast_pbc_dx(self, rvec ref, rvec other, rvec dx) nogil
    cdef void fast_pbc_dx_leaflet(self, rvec ref, rvec other, rvec dx, rvec normal) nogil
    cdef real fast_distance2(self, rvec ref, rvec other) nogil
    cdef real fast_distance(self, rvec aref, rvec other) nogil
    cdef real fast_leaflet_distance2(self, rvec a, rvec b, rvec normal) nogil
    cdef real fast_leaflet_distance(self, rvec a, rvec b, rvec normal) nogil
    cdef void fast_pbc_xcm(self, real[:, ::1] coords, rvec xcm, fsl_int limit=*) nogil
    cdef void fast_pbc_xcm_from_ref(self, real[:, ::1] coords, rvec ref, rvec xcm, fsl_int limit=*) nogil
    cdef real[:, ::1]fast_put_atoms_in_bbox(self, real[:,::1] coords) nogil

cdef struct topol_atom_t:
    fsl_int name_id
    fsl_int internal_id
    fsl_int atomid
    fsl_int resid
    fsl_int residue_internal_id

cdef struct topol_residue_t:
    fsl_int name_id
    fsl_int internal_id
    fsl_int resid
    fsl_int size
    fsl_int allocated_size
    fsl_int *atom_internal_ids
    fsl_int *atomids

cdef class TopologyGroup(object):
    # Attrbutes
    cdef readonly Topology parent
    cdef bytes name
    cdef bytes resname
    cdef readonly fsl_int resid
    cdef topol_atom_t *atoms
    cdef readonly fsl_int size
    cdef fsl_int allocated_size

    # Methods
    cdef bint fast_contains(self, topol_atom_t atom) nogil except *
    cdef void add_atom(self, topol_atom_t atom) nogil except *
    cdef void merge_other_group(self, TopologyGroup other) nogil
    cdef fsl_int[:] fast_get_atomids(self) nogil
    cdef object group_by_resid(self)


cdef class Topology(object):
    # Attributes
    cdef strname_t *atomnames
    cdef fsl_int atomnames_size
    cdef fsl_int atomnames_allocated_size

    cdef strname_t *resnames
    cdef fsl_int resnames_size
    cdef fsl_int resnames_allocated_size

    cdef topol_atom_t *atoms
    cdef fsl_int atoms_size
    cdef fsl_int atoms_allocated_size

    cdef fsl_int *atomids_to_internalids
    cdef fsl_int atomids_to_internalids_size
    cdef fsl_int atomids_to_internalids_allocated_size

    cdef topol_residue_t *residues
    cdef fsl_int residues_size
    cdef fsl_int residues_allocated_size

    cdef fsl_int *resids_to_internalids
    cdef fsl_int resids_to_internalids_size
    cdef fsl_int resids_to_internalids_allocated_size


    # Methods
    cdef void set_residues_allocation(self, fsl_int size)
    cdef void set_atoms_allocation(self, fsl_int size)
    cdef void internal_append(self, fsl_int resid, bytes resname, bytes atomname, fsl_int atomid)
    cdef fsl_int fast_get_atom_internal_id(self, fsl_int atomid) nogil
    cdef topol_atom_t *fast_get_atom(self, fsl_int atomid) nogil
    cdef fsl_int fast_get_resid_from_atomid(self, fsl_int atomid) nogil
    cdef fsl_int fast_get_residue_internal_id(self, fsl_int resid) nogil
    cdef topol_residue_t *fast_get_residue(self, fsl_int resid) nogil
    cdef topol_residue_t *fast_get_residue_from_atomid(self, fsl_int atomid) nogil

cdef class TopologyReader(object):
    # Attributes
    cdef bytes filename
    cdef readonly Topology topology

    # Methods
    cdef load(self)


cdef class CoordinateReader(object):
    # Attributes
    cdef bytes filename
    cdef readonly fsl_int nframes
    cdef fsl_int[:] coordinate_offsets
    cdef fsl_int[:] box_offsets
    cdef readonly real[:] timesteps

    # Methods
    cdef preload(self)
    cdef PBCBox load_box(self, fsl_int frame_id)
    cdef assert_frame_id(self, fsl_int frame_id)
    cdef real[:, ::1] load_coords(self, fsl_int frame_id, fsl_int[:] atomids) nogil except *


cdef class IndexReader(object):
    # Attributes
    cdef bint loaded
    cdef bytes filename
    cdef readonly object groups

    # Methods
    cdef fast_load(self)

cdef class Frame:
    # Attributes
    cdef fsl_int index
    cdef readonly real timestep
    cdef readonly Trajectory trajectory
    cdef readonly Topology topology
    cdef CoordinateReader coords_reader
    cdef readonly PBCBox box
    cdef readonly fsl_int size
    # Headgroups
    cdef real[:, ::1] hg_group_coords, hg_coords_bbox
    # Lipids
    cdef real[:, ::1] lipid_coords, lipid_coords_bbox
    cdef list lipid_coords_aslist, lipids_coords_bbox_aslist
    cdef fsl_int[:] lipid_atomids_offsets
    # Simplified lipids
    cdef real[:, ::1] bead_coords, bead_coords_bbox
    # Interacting groups
    cdef real[:, ::1] interacting_group_coords, interacting_group_coords_bbox

    # Related to NS
    cdef real proximity_cutoff
    cdef ns_neighborhood_holder *neighbors

    cdef real[:, ::1] directions
    cdef real[:, ::1] normals

    cdef list aggregates
    cdef bint aggregates_retrieved
    cdef real aggregates_cutoff

    cdef list membranes
    cdef bint membranes_retrieved
    cdef real membranes_cutoff

    # Methods
    cpdef fsl_int get_forcefield_type(self)
    cdef real[:, ::1] fast_get_hg_group_coords(self) nogil
    cdef real[:, ::1] fast_get_hg_group_coords_bbox(self) nogil
    cdef real[:, ::1] fast_get_interacting_group_coords(self) nogil
    cdef real[:, ::1] fast_get_interacting_group_coords_bbox(self) nogil
    cdef real[:, ::1] fast_get_bead_coords(self) nogil
    cdef real[:, ::1] fast_get_bead_coords_bbox(self) nogil
    cdef real[:, ::1] fast_get_lipid_coords(self) nogil
    cdef list get_lipid_coords_aslist(self)
    cdef real[:, ::1] fast_get_lipid_coords_bbox(self) nogil
    cdef list get_lipid_coords_bbox_aslist(self)
    cdef real[:, ::1] fast_get_directions(self) nogil except *
    cdef real[:, ::1] fast_get_normals(self, real proximity_cutoff=*) nogil
    cpdef list get_aggregates(self, real cutoff=*, bint update=*)
    cpdef list get_membranes(self, real cutoff=*, bint update=*)


cdef class Trajectory(object):
    # Attributes
    cdef bytes bead_group_name
    cdef bytes interacting_group_name
    cdef Topology topology
    cdef TopologyReader topol_reader
    cdef IndexReader index_reader
    cdef CoordinateReader coords_reader
    cdef bint be_verbose

    cdef fsl_int cur_frame
    cdef bint initialized

    cdef fsl_int[:] hg_group_atomids, interacting_atomids, lipid_atomids, hg_bead_atomids
    cdef fsl_int[:] lipid_atomids_offsets, hg_bead_atomids_offsets
    cdef fsl_int[:] lipid_hg_indices
    cdef fsl_int forcefield_type

    # Cache-related
    cdef real normal_cutoff
    cdef list membranes
    cdef list aggregates

    # Methods
    cdef assert_initialization(self)
    cpdef fsl_int get_forcefield_type(self)
    cpdef real get_rvdw(self)



