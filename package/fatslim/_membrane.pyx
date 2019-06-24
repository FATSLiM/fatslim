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

# Cython directives
# cython: language_level=3

cimport cython

DEF NOTSET = -12345
DEF DIM = 3
DEF XX = 0
DEF YY = 1
DEF ZZ = 2
DEF EPSILON = 1e-6
DEF PI = 3.141592653589793

import numpy as np
cimport numpy as np

import warnings

from ._typedefs cimport real, rvec, fsl_int, matrix
from ._typedefs cimport rvec_dprod, rvec_copy

from ._typedefs cimport rvec_norm, rvec_normalize, rvec_clear, rvec_inc
from ._typedefs cimport mat_from_rvec, invert_mat

from ._core cimport SimplifiedLipid, LipidRegistry, FixedQueue

from ._aggregate cimport LipidAggregate

from ._geometry cimport PBCBox
from ._geometry cimport complete_basis, rvec_to_basis
from ._geometry cimport real_point, Polygon, polygon_new, polygon_get_area, fast_clip_zoi, polygon_destroy, polygon_append, polygon_empty

from libc.math cimport acos, fabs


@cython.boundscheck(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
def identify_membranes(LipidRegistry system, bint force_update=False):
    # Constants
    DEF MIN_LEAFLET_SIZE = 36 # 6x6 lipids

    DEF UNKNOWN = 0
    DEF LEAFLET1 = 1
    DEF LEAFLET2 = 2
    DEF FRONTIER = 3
    DEF EDGE = 4
    DEF REJECTED = 5
    DEF ADDED_TO_LEAFLET1 = 6
    DEF ADDED_TO_LEAFLET2 = 7
    DEF ADDED_TO_BORDER1 = 8
    DEF ADDED_TO_BORDER2 = 9

    DEF UNPROCESSED = 0
    DEF ADDED_TO_STACK = 1
    DEF PROCESSED = 2

    # Cdef variables
    cdef LipidAggregate aggregate
    cdef fsl_int i, count
    cdef fsl_int parent_index, nid
    cdef fsl_int[:] lipid_status = np.empty(system._nlipids, dtype=int)
    cdef fsl_int[:] l1_members = np.empty(system._nlipids, dtype=int)
    cdef fsl_int[:] l2_members = np.empty(system._nlipids, dtype=int)
    cdef FixedQueue l1_neighbours = FixedQueue(system._nlipids)
    cdef FixedQueue l2_neighbours = FixedQueue(system._nlipids)
    cdef FixedQueue candidates = FixedQueue(system._nlipids)
    cdef fsl_int[:] processing_status = np.empty(system._nlipids, dtype=int)
    cdef real[:] sum_normals = np.empty(system._nlipids, dtype=np.float32)
    cdef rvec sum_normal
    cdef real edger_limit
    cdef fsl_int n_leaflets, nlipids_l1, nlipids_l2
    cdef FixedQueue stack = FixedQueue(system._nlipids)
    cdef rvec dx, l1_neighbours_cog, l2_neighbours_cog
    cdef rvec dx1, dx2
    cdef rvec interleaflet_axis, l1_neighbours_normal, l2_neighbours_normal
    cdef real d, d_l1, d_l2, dot_value

    cdef FixedQueue l1_core = FixedQueue(system._nlipids)
    cdef FixedQueue l2_core = FixedQueue(system._nlipids)
    cdef FixedQueue l1_border = FixedQueue(system._nlipids)
    cdef FixedQueue l2_border = FixedQueue(system._nlipids)


    membranes = []

    # Make sure everything is updated
    system.update()

    # Step 1: Get "naive" aggregates (i.e. aggregates based only on distance)
    system.find_aggregates(force_update)

    # Step 2: Iterate over the known aggregates to find if they can be splitted into potential leaflets
    potential_leaflets = []
    for aggregate in system.aggregates:
        # Make sure the aggregate is updated
        aggregate.update()

        with nogil:
            # If aggregate is small, do not bother to check if it could be a leaflet
            if aggregate._size < MIN_LEAFLET_SIZE:
                continue

            # (Re)initialize limit
            edger_limit = 0

            # Try to identify edges using normals
            for parent_index in range(aggregate._size):
                rvec_clear(sum_normal)

                rvec_inc(sum_normal, &aggregate._lipid_normals[parent_index, XX])

                count = 1
                for i in range(aggregate._lipid_neighbours.nneighbours[parent_index]):
                    nid = aggregate._lipid_neighbours.neighbours[parent_index][i]
                    if rvec_dprod(&aggregate._lipid_normals[parent_index, XX], &aggregate._lipid_normals[nid, XX]) > 0:
                        count += 1
                        rvec_inc(sum_normal, &aggregate._lipid_normals[nid, XX])

                sum_normals[parent_index] = rvec_norm(sum_normal) / count
                edger_limit += sum_normals[parent_index]


            edger_limit /= aggregate._size
            edger_limit -= 0.1

            if edger_limit < 0.6:
                edger_limit = 0.6

            #print("DEBUG: edger limit={:.3f}".format(edger_limit))

            # Loop over lipids to flag the ones that
            # 1. should be rejected from any potential leaflets -> REJECTED
            # 2. might be part of an edge between two leaflets (eg in a bicelle or a pore) -> EDGE
            # 3. lipids which are close to edges but not part of it -> FRONTIER
            lipid_status[:] = UNKNOWN
            for parent_index in range(aggregate._size):

                if sum_normals[parent_index] < 0.5:
                    lipid_status[parent_index] = REJECTED

                elif sum_normals[parent_index] < edger_limit:
                    lipid_status[parent_index] = EDGE

                    for i in range(aggregate._lipid_neighbours.nneighbours[parent_index]):
                        nid = aggregate._lipid_neighbours.neighbours[parent_index][i]

                        if lipid_status[nid] != EDGE and lipid_status[nid] != REJECTED:
                            lipid_status[nid] = FRONTIER

            #print("DEBUG: Flagging lipids -> {} rejected, {} edgers, {} frontierers".format(
            #    len(np.argwhere(np.asarray(lipid_status) == REJECTED).flatten()),
            #    len(np.argwhere(np.asarray(lipid_status) == EDGE).flatten()),
            #    len(np.argwhere(np.asarray(lipid_status) == FRONTIER).flatten()),
            #))

            # (Re)initialize process status
            processing_status[:] = UNPROCESSED

            current_leaflet = -1

            # Build coherent leaflets from aggregate. Three possibility here:
            # 1. no leaflet can be identified (eg the aggregate is just a lipid cluster rather than a proper leaflet
            # 2. one single leaflet is identified (most likely case...)
            # 3. two leaflets can be identified (eg bicelle or classical bilayer with "bridges" between both leaflets such
            #    as a pore or flip-floping lipids)
            n_leaflets = 0
            nlipids_l1 = 0
            nlipids_l2 = 0
            for parent_index in range(aggregate._size):

                # if lipid status is already clear, no need to do anything
                if lipid_status[parent_index] != UNKNOWN:
                    continue

                if n_leaflets == 0:
                    current_leaflet = LEAFLET1
                    n_leaflets += 1
                    nlipids_leaflet1 = 0
                elif n_leaflets == 1: # Already one leaflet defined

                    # Check the leaflet size, if too small, it is erased
                    if nlipids_l1 < MIN_LEAFLET_SIZE:
                        for i in range(nlipids_l1):
                            lipid_status[l1_members[i]] = REJECTED
                            processing_status[l1_members[i]] = UNPROCESSED
                    else:
                        current_leaflet = LEAFLET2
                        n_leaflets += 1
                else:
                    # Check the leaflet size, if too small, it is erased
                    if nlipids_l2 < MIN_LEAFLET_SIZE:
                        for i in range(nlipids_l2):
                            lipid_status[l2_members[i]] = REJECTED
                            processing_status[l2_members[i]] = UNPROCESSED
                    else:
                        with gil:
                            warnings.warn("Fatslim only support splitting and aggregate into 2 leaflets")
                        break

                # Reset stack to hold only the leaflet seed
                stack.fast_empty()
                stack.fast_add(parent_index)

                while stack.size > 0:

                    ref_index = stack.fast_pop()

                    processing_status[ref_index] = PROCESSED
                    lipid_status[ref_index] = current_leaflet

                    if current_leaflet == LEAFLET1:
                        l1_members[nlipids_l1] = ref_index
                        nlipids_l1 += 1
                    else:
                        l2_members[nlipids_l2] = ref_index
                        nlipids_l2 += 1

                    for i in range(aggregate._lipid_neighbours.nneighbours[ref_index]):
                        nid = aggregate._lipid_neighbours.neighbours[ref_index][i]

                        if processing_status[nid] != UNPROCESSED:
                            continue

                        if lipid_status[nid] != UNKNOWN:
                            continue

                        if rvec_dprod(&aggregate._lipid_normals[ref_index, XX],
                                      &aggregate._lipid_normals[nid, XX]) < 0.0:
                            lipid_status[nid] = REJECTED
                            continue

                        stack.fast_add(nid)
                        processing_status[nid] = ADDED_TO_STACK

            if n_leaflets == 2 and nlipids_l2 < MIN_LEAFLET_SIZE:
                for i in range(nlipids_l2):
                    lipid_status[l2_members[i]] = REJECTED

                n_leaflets = 1

            # print("DEBUG: About to handle frontier lipids -> l1:{}, l2:{} ->  rejected, {} edgers, {} frontierers".format(
            #     len(np.argwhere(np.asarray(lipid_status) == LEAFLET1).flatten()),
            #     len(np.argwhere(np.asarray(lipid_status) == LEAFLET2).flatten()),
            #     len(np.argwhere(np.asarray(lipid_status) == REJECTED).flatten()),
            #     len(np.argwhere(np.asarray(lipid_status) == EDGE).flatten()),
            #     len(np.argwhere(np.asarray(lipid_status) == FRONTIER).flatten()),
            # ))

            # Add frontiers to leaflets
            for parent_index in range(aggregate._size):
                if lipid_status[parent_index] not in [FRONTIER, REJECTED]:
                    continue

                d_l1 = 1e6
                d_l2 = 1e6
                flag = False
                for i in range(aggregate._lipid_neighbours.nneighbours[parent_index]):
                    nid = aggregate._lipid_neighbours.neighbours[parent_index][i]

                    if lipid_status[nid] == LEAFLET1:
                        system.box.fast_pbc_dx(
                            &aggregate._lipid_positions[parent_index, XX],
                            &aggregate._lipid_positions[nid, XX],
                            dx
                        )
                        d = rvec_norm(dx)

                        dot_value = rvec_dprod(
                            &aggregate._lipid_normals[parent_index, XX],
                            &aggregate._lipid_normals[nid, XX]
                        )

                        if d < d_l1 and dot_value > 0.866:
                            d_l1 = d
                            flag = True

                    elif lipid_status[nid] == LEAFLET2:
                        system.box.fast_pbc_dx(
                            &aggregate._lipid_positions[parent_index, XX],
                            &aggregate._lipid_positions[nid, XX],
                            dx
                        )
                        d = rvec_norm(dx)

                        dot_value = rvec_dprod(
                            &aggregate._lipid_normals[parent_index, XX],
                            &aggregate._lipid_normals[nid, XX]
                        )

                        if d < d_l2 and dot_value > 0.866:
                            d_l2 = d
                            flag = True

                if flag:
                    if d_l1 < d_l2:
                        lipid_status[parent_index] = ADDED_TO_LEAFLET1
                    else:
                        lipid_status[parent_index] = ADDED_TO_LEAFLET2


            # print("DEBUG: after dealing with frontier lipids -> {} rejected, {} edgers, {} frontierers,"
            #       "{} added to l1, {} adder to l2".format(
            #     len(np.argwhere(np.asarray(lipid_status) == REJECTED).flatten()),
            #     len(np.argwhere(np.asarray(lipid_status) == EDGE).flatten()),
            #     len(np.argwhere(np.asarray(lipid_status) == FRONTIER).flatten()),
            #     len(np.argwhere(np.asarray(lipid_status) == ADDED_TO_LEAFLET1).flatten()),
            #     len(np.argwhere(np.asarray(lipid_status) == ADDED_TO_LEAFLET2).flatten()),
            # ))

            if n_leaflets == 1:
                # Check if make sense to add the rejected lipid the unique leaflet
                for seed_index in range(aggregate._size):
                    if lipid_status[seed_index] not in [EDGE, FRONTIER]:
                        continue

                    l1_neighbours.fast_empty()

                    processing_status[:] = UNPROCESSED

                    stack.fast_empty()
                    stack.fast_add(seed_index)

                    exhaust_stack = False


                    while stack.size > 0:
                        ref_index = stack.fast_pop()

                        processing_status[ref_index] = PROCESSED

                        for i in range(aggregate._lipid_neighbours.nneighbours[ref_index]):
                            nid = aggregate._lipid_neighbours.neighbours[ref_index][i]

                            if processing_status[nid] != UNPROCESSED:
                                continue

                            if lipid_status[nid] == LEAFLET1:
                                l1_neighbours.fast_add(nid)
                                exhaust_stack = True

                            if not exhaust_stack:
                                stack.fast_add(nid)
                                processing_status[nid] = ADDED_TO_STACK

                    if l1_neighbours.size == 0:
                        lipid_status[seed_index] = REJECTED
                    else:
                        lipid_status[seed_index] = ADDED_TO_BORDER1

                        for i in range(l1_neighbours.size):
                            nid = l1_neighbours.elements[i]

                            dot_value = rvec_dprod(
                                &aggregate._lipid_normals[seed_index, XX],
                                &aggregate._lipid_normals[nid, XX]
                            )

                            system.box.fast_pbc_dx(
                                &aggregate._lipid_positions[seed_index, XX],
                                &aggregate._lipid_positions[nid, XX],
                                dx
                            )

                            d = fabs(
                                rvec_dprod(
                                    &aggregate._lipid_normals[seed_index, XX],
                                    dx
                                ))

                            if dot_value < 0.707 or d > 10:
                                lipid_status[seed_index] = REJECTED
                                break
            else:
                for seed_index in range(aggregate._size):
                    if lipid_status[seed_index] not in [EDGE, FRONTIER]:
                        continue

                    l1_neighbours.fast_empty()
                    l2_neighbours.fast_empty()

                    candidates.fast_empty()
                    candidates.fast_add(seed_index)

                    processing_status[:] = UNPROCESSED


                    stack.fast_empty()
                    stack.fast_add(seed_index)

                    exhaust_stack = False

                    while stack.size > 0:

                        ref_index = stack.fast_pop()

                        processing_status[ref_index] = PROCESSED

                        for i in range(aggregate._lipid_neighbours.nneighbours[ref_index]):
                            nid = aggregate._lipid_neighbours.neighbours[ref_index][i]

                            if processing_status[nid] != UNPROCESSED:
                                continue

                            if lipid_status[nid] == LEAFLET1:
                                l1_neighbours.fast_add_unique(nid)

                            elif lipid_status[nid] == LEAFLET2:
                                l2_neighbours.fast_add_unique(nid)

                            elif lipid_status[nid] in (EDGE, FRONTIER):
                                candidates.fast_add_unique(nid)

                            if l1_neighbours.size > 0 and l2_neighbours.size > 0:
                                exhaust_stack = True

                            if not exhaust_stack:
                                stack.fast_add(nid)
                                processing_status[nid] = ADDED_TO_STACK

                    system.box.fast_pbc_centroid(
                        aggregate._lipid_positions,
                        l1_neighbours_cog,
                        l1_neighbours.elements[:l1_neighbours.size]
                    )

                    system.box.fast_pbc_centroid(
                        aggregate._lipid_positions,
                        l2_neighbours_cog,
                        l2_neighbours.elements[:l2_neighbours.size]
                    )

                    rvec_clear(l1_neighbours_normal)
                    for i in range(l1_neighbours.size):
                        rvec_inc(l1_neighbours_normal, &aggregate._lipid_normals[l1_neighbours.elements[i], XX])

                    rvec_clear(l2_neighbours_normal)
                    for i in range(l2_neighbours.size):
                        rvec_inc(l2_neighbours_normal, &aggregate._lipid_normals[l2_neighbours.elements[i], XX])

                    for i in range(DIM):
                        interleaflet_axis[i] = l1_neighbours_normal[i] - l2_neighbours_normal[i]
                    rvec_normalize(interleaflet_axis)

                    system.box.fast_pbc_dx(l1_neighbours_cog, l1_neighbours_cog, dx)
                    d = fabs(
                            rvec_dprod(
                                dx,
                                interleaflet_axis
                            ))

                    #print("\nDEBUG: L1 COG: [{:.3f}, {:.3f}, {:.3f}],  L2 COG: [{:.3f}, {:.3f}, {:.3f}],  "
                    #      "interleaflet axis: [{:.3f}, {:.3f}, {:.3f}]".format(
                    #    l1_neighbours_cog[XX], l1_neighbours_cog[YY], l1_neighbours_cog[ZZ],
                    #    l2_neighbours_cog[XX], l2_neighbours_cog[YY], l2_neighbours_cog[ZZ],
                    #    interleaflet_axis[XX], interleaflet_axis[YY], interleaflet_axis[ZZ]
                    #))

                    while candidates.size > 0:
                        seed_index = candidates.fast_pop()

                        system.box.fast_pbc_dx(l1_neighbours_cog, &aggregate._lipid_positions[seed_index, XX], dx1)
                        d_l1 = fabs(rvec_dprod(dx1, interleaflet_axis))

                        system.box.fast_pbc_dx(&aggregate._lipid_positions[seed_index, XX], l2_neighbours_cog, dx2)
                        d_l2 = fabs(rvec_dprod(dx2, interleaflet_axis))


                        if d_l2 > d_l1:
                            lipid_status[seed_index] = ADDED_TO_BORDER1
                            #print("Candidate #{} added to border 1 (dprod1: {:.3} vs dpprod2:{:.3f})".format(
                            #    seed_index, d_l1, d_l2))
                        else:
                            lipid_status[seed_index] = ADDED_TO_BORDER2
                            #print("Candidate #{} added to border 2 (dprod1: {:.3} vs dpprod2:{:.3f})".format(
                            #    seed_index, d_l1, d_l2))

            # print("DEBUG: aftermath -> {} rejected, {} edgers, {} frontierers,"
            #       "{} added to l1, {} adder to l2, {} added to l1 border, {} added to l2 border".format(
            #     len(np.argwhere(np.asarray(lipid_status) == REJECTED).flatten()),
            #     len(np.argwhere(np.asarray(lipid_status) == EDGE).flatten()),
            #     len(np.argwhere(np.asarray(lipid_status) == FRONTIER).flatten()),
            #     len(np.argwhere(np.asarray(lipid_status) == ADDED_TO_LEAFLET1).flatten()),
            #     len(np.argwhere(np.asarray(lipid_status) == ADDED_TO_LEAFLET2).flatten()),
            #     len(np.argwhere(np.asarray(lipid_status) == ADDED_TO_BORDER1).flatten()),
            #     len(np.argwhere(np.asarray(lipid_status) == ADDED_TO_BORDER2).flatten()),
            # ))

            l1_core.fast_empty()
            l1_border.fast_empty()
            l2_core.fast_empty()
            l2_border.fast_empty()

            for i in range(aggregate._size):
                if lipid_status[i] in [LEAFLET1, ADDED_TO_LEAFLET1]:
                    l1_core.fast_add(aggregate._lipid_ids[i])
                elif lipid_status[i] == ADDED_TO_BORDER1:
                    l1_border.fast_add(aggregate._lipid_ids[i])
                elif lipid_status[i] in [LEAFLET2, ADDED_TO_LEAFLET2]:
                    l2_core.fast_add(aggregate._lipid_ids[i])
                elif lipid_status[i] == ADDED_TO_BORDER2:
                    l2_border.fast_add(aggregate._lipid_ids[i])

        if l1_core.size > 0:
            core_indices = np.asarray(l1_core.elements[:l1_core.size]).copy()
            border_indices = np.asarray(l1_border.elements[:l1_border.size]).copy()

            potential_leaflets.append(Monolayer(core_indices, system, border_indices))

        if l2_core.size > 0:
            core_indices = np.asarray(l2_core.elements[:l2_core.size]).copy()
            border_indices = np.asarray(l2_border.elements[:l2_border.size]).copy()

            potential_leaflets.append(Monolayer(core_indices, system, border_indices))

    membranes = []
    if len(potential_leaflets) < 2:
        warnings.warn("Only {} potential leaflets found: No membrane".format(len(potential_leaflets)))
    elif len(potential_leaflets) == 2:
        try:
            membranes.append(Bilayer(potential_leaflets[0], potential_leaflets[1]))
        except ValueError:
            warnings.warn("The only two leaflets found are incompatible: No membrane")
    else:
        raise NotImplementedError

    return membranes


cdef class Monolayer(LipidAggregate):
    def __init__(self, fsl_int[:] core_ids, LipidRegistry system, fsl_int[:] border_ids=None):
        cdef fsl_int[:] lipids_ids
        cdef fsl_int i

        if border_ids is None:
            border_ids = np.empty(0, dtype=int)

        lipids_ids = np.empty(core_ids.shape[0] + border_ids.shape[0], dtype=int)
        for i in range(core_ids.shape[0]):
            lipids_ids[i] = core_ids[i]
        for i in range(border_ids.shape[0]):
            lipids_ids[i + core_ids.shape[0]] = border_ids[i]

        super().__init__(lipids_ids, system)

        self._lipid_core_ids = np.sort(core_ids)
        self._core_size = self._lipid_core_ids.shape[0]

        self._lipid_border_ids = np.sort(border_ids)
        self._border_size = self._lipid_border_ids.shape[0]

        self._parent = None

        # APL-related
        self._lastupdate_apl = -1
        self._area = 0
        self._apl = 0
        self._lipid_apls = np.empty(self._size, dtype=np.float32)

    def __repr__(self):
        return "<Monolayer: {} lipids (core: {} - border: {})>".format(self._size, self._core_size, self._border_size)

    def is_compatible(self, other):
        if not isinstance(other, Monolayer):
            print("DEBUG: Bad object")
            return False

        if self.is_planar:
            if not other.is_planar:
                print("DEBUG: incoherent planarity (other not planar)")
                return False

            if np.dot(self.normal, other.normal) > -0.707:
                print("DEBUG: normals too close (dprod: {:.3f})".format(np.dot(self.normal, other.normal)))
                return False

        else:
            if other.is_planar:
                print("DEBUG: incoherent planarity (other is planar)")
                return False

            if self.system.box.pbc_distance(self.position, other.position) > 1.0:
                print("DEBUG: COG too far (d: {:.3f}".format(
                    self.system.box.pbc_distance(self.position, other.position)
                ))
                return False

        return True

    @property
    def lipid_thicknesses(self):
        if self._parent is None:
            raise AttributeError("Not a membrane leaflet: thickness is not available")

        if self._parent._leaflet1 == self:
            start = 0
        else:
            start = self._parent._leaflet1._size
        end = start + self._size

        return self._parent.lipid_thicknesses[start:end].copy()

    @property
    def thickness(self):
        return np.nanmean(self.lipid_thicknesses)


    cdef compute_apl(self, bint force_update=False):
        cdef fsl_int i, j
        cdef fsl_int ref_beadid, beadid
        cdef rvec plane_x, plane_y, plane_z
        cdef real[:] ref_normal, ref_position
        cdef real[:] current_position
        cdef matrix conv_mat_revert, conv_mat
        cdef rvec dx, proj

        cdef Polygon *cell
        cdef real_point tmp_pt
        cdef Polygon *buffer = polygon_new()
        cdef real_point ref_point_2d
        cdef real_point current_position_2d


        # No need to do anything if the membranes are already identified for the current frame
        if self._lastupdate_apl == self.system._lastupdate and not force_update:
            return self._apl, self._lipid_apls, self._area

        warnings.warn("Current implementation of APL calculation supports only the core lipids and ignores all the "
                      " interacting atoms (eg protein)")


        self.update(force_update)
        box = self.system.box

        cell = polygon_new()

        ref_point_2d[XX] = 0
        ref_point_2d[YY] = 0

        self._lipid_apls[:] = np.nan

        for i in range(self._core_size):
            ref_beadid = self._lipid_core_ids[i]

            ref_normal = self.system._lipid_normals[ref_beadid]
            ref_position = self.system._lipid_positions[ref_beadid]


            # Compute projection matrix
            # Build plane basis
            rvec_copy(&ref_normal[XX], plane_z)
            complete_basis(plane_z, plane_x, plane_y)

            # Get conversion matrix
            mat_from_rvec(plane_x, plane_y, plane_z, conv_mat_revert)
            invert_mat(conv_mat_revert, conv_mat)

            # print("Resid {} -> Convertion matrix:\n{:.3f}, {:.3f}, {:.3f}\n{:.3f}, {:.3f}, {:.3f}\n{:.3f}, {:.3f}, {:.3f}\n".format(
            #     self.system.lipids[ref_beadid].resid,
            #     conv_mat_revert[XX][XX], conv_mat_revert[XX][YY], conv_mat_revert[XX][ZZ],
            #     conv_mat_revert[YY][XX], conv_mat_revert[YY][YY], conv_mat_revert[YY][ZZ],
            #     conv_mat_revert[ZZ][XX], conv_mat_revert[ZZ][YY], conv_mat_revert[ZZ][ZZ],
            # ))

            # (Re)initialize cell
            polygon_empty(cell)

            # Start with a big square
            tmp_pt[XX] = -2 * self.system.ns_cutoff
            tmp_pt[YY] = -2 * self.system.ns_cutoff
            polygon_append(cell, tmp_pt)
            tmp_pt[XX] = 2 * self.system.ns_cutoff
            polygon_append(cell, tmp_pt)
            tmp_pt[YY] = 2 * self.system.ns_cutoff
            polygon_append(cell, tmp_pt)
            tmp_pt[XX] = -2 * self.system.ns_cutoff
            polygon_append(cell, tmp_pt)


            for j in range(self.system._lipid_neighbours.nneighbours[ref_beadid]):
                beadid = self.system._lipid_neighbours.neighbours[ref_beadid][j]

                current_position = self.system._lipid_positions[beadid]

                # Get the closest image of the neighbor
                box.fast_pbc_dx(&ref_position[XX], &current_position[XX], dx)

                # Get the coordinates into the new basis (plane_x, plane_y, plane_z)
                rvec_to_basis(dx, conv_mat, proj)

                # Store the projection onto the (plane_x, plane_y) plane
                current_position_2d[XX] = proj[XX]
                current_position_2d[YY] = proj[YY]

                # print("Resid {} neighbour {}: [{:.3f}, {:.3f}]".format(
                #     self.system.lipids[ref_beadid].resid,
                #     self.system.lipids[beadid].resid,
                #     current_position_2d[XX], current_position_2d[YY]
                # ))

                fast_clip_zoi(cell, ref_point_2d, current_position_2d, buffer)

            self._lipid_apls[i] = polygon_get_area(cell)


        polygon_destroy(cell)
        polygon_destroy(buffer)

        self._lastupdate_apl = self.system._lastupdate
        self._apl = np.nanmean(self._lipid_apls)
        self._area = np.nansum(self._lipid_apls)

        return self._apl, self._lipid_apls, self._area

    @property
    def apl(self):
        self.compute_apl()
        return self._apl

    @property
    def lipid_apls(self):
        self.compute_apl()
        return np.asarray(self._lipid_apls).copy()

    @property
    def area(self):
        self.compute_apl()
        return self._area



cdef class Bilayer:
    def __init__(self, Monolayer leaflet1, Monolayer leaflet2):
        if not leaflet1.is_compatible(leaflet2):
            raise ValueError("Can not build bilayer from incompatible monolayers")

        self.system = leaflet1.system
        self._is_planar = leaflet1._is_planar

        if leaflet1._is_planar:
            if leaflet1._position[ZZ] > leaflet2._position[ZZ]:
                self._leaflet1 = leaflet1
                self._leaflet2 = leaflet2
            else:
                self._leaflet1 = leaflet2
                self._leaflet2 = leaflet1
        elif leaflet1._size > leaflet2._size:
            self._leaflet1 = leaflet1
            self._leaflet2 = leaflet2
        else:
            self._leaflet1 = leaflet2
            self._leaflet2 = leaflet1

        self._leaflet1._parent = self
        self._leaflet2._parent = self

        self._size = self._leaflet1._size + self._leaflet2._size

        # Thickness-related
        self._lastupdate_thickness = -1
        self._thickness = 0
        self._lipid_thicknesses = np.empty(self._size, dtype=np.float32)
        self._lipid_interleaflet_gaps = np.empty(self._size, dtype=np.float32)

    def __repr__(self):
        if self._leaflet1._is_planar:
            membrane_type = "Planar membrane"
            l1_type = "upper"
            l2_type = "lower"
        else:
            membrane_type = "Vesicle"
            l1_type = "outer"
            l2_type = "inner"
        return "<{}: {} lipids on {} leaflet - {} lipids on {} leaflet>".format(
            membrane_type,
            self._leaflet1._size,
            l1_type,
            self._leaflet2._size,
            l2_type
        )

    def __getitem__(self, item):
        if item in (-2, 0):
            return self._leaflet1
        elif item in (-1, 1):
            return self._leaflet2
        else:
            raise IndexError("A *bi*layer contains two leaflets")

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.cdivision(True)
    cdef compute_thickness(self, bint force_update=False):
        # Cdef variables
        cdef PBCBox box
        cdef Monolayer ref_leaflet, other_leaflet
        cdef fsl_int internal_id, ref_internal_id, ref_beadid, beadid
        cdef real dprod_best, dprod_val, norm, dprod_val2, d_min
        cdef rvec ref_position, ref_normal
        cdef rvec dx
        cdef fsl_int locnorm_beadid
        cdef rvec locnorm_position, locnorm_normal, locnorm_dx

        cdef fsl_int closest_beadid
        cdef rvec closest_position, closest_normal, closest_dx

        cdef fsl_int leafnorm_beadid
        cdef rvec leafnorm_position, leafnorm_normal, leafnorm_dx

        cdef fsl_int other_beadid
        cdef rvec other_position, other_normal

        # No need to do anything if thicknessed are already calculated for the current frame
        if self._lastupdate_thickness == self.system._lastupdate and not force_update:
            return self._thickness, self._lipid_thicknesses

        # Update leaflets (and system)
        self._leaflet1.update(force_update=force_update)
        self._leaflet2.update(force_update=force_update)

        # Load useful stuff
        box = self.system.box
        universe_positions = self.system.universe_coords_bbox

        for internal_id in range(self._size):
            if internal_id < self._leaflet1._size:
                ref_leaflet = self._leaflet1
                other_leaflet = self._leaflet2
                ref_internal_id = internal_id
            else:
                ref_leaflet = self._leaflet2
                other_leaflet = self._leaflet1
                ref_internal_id = internal_id - self._leaflet1._size


            if ref_internal_id > ref_leaflet._core_size:
                self._lipid_thicknesses[ref_internal_id] = NOTSET
                continue

            # Get actual beadid
            ref_beadid = ref_leaflet._lipid_ids[ref_internal_id]

            # Get average normal and position for the reference bead
            self.system.compute_weighted_average(ref_beadid, ref_position, ref_normal)


            # Get the bead from the other leaflet which can be used to compute thickness
            # The candidate is selected by checking the dot product of dx and reference normal
            beadid_best = -1
            dprod_best = 0

            for j in range(other_leaflet._size):
                beadid = other_leaflet._lipid_ids[j]


                box.fast_pbc_dx_leaflet(
                    &self.system._lipid_positions[ref_beadid, XX],
                    &self.system._lipid_positions[beadid, XX],
                    dx,
                    &ref_normal[XX])

                dprod_val = rvec_dprod(&ref_normal[XX],
                                       dx)

                norm = rvec_norm(dx)

                dprod_val2 = rvec_dprod(&ref_normal[XX],
                                       &self.system._lipid_normals[beadid, XX])

                if dprod_val/norm < dprod_best and dprod_val2 < -0.707:
                    dprod_best = dprod_val/norm
                    beadid_best = beadid


            box.fast_pbc_dx_leaflet(
                    &self.system._lipid_positions[ref_beadid, XX],
                    &self.system._lipid_positions[beadid_best, XX],
                    dx,
                    &ref_normal[XX])

            # print("\nRef bead resid {}, Best bead resid {}: dx:[{}, {}, {}], ref normal:{}, dprod:{}".format(
            #     self.system.lipids[ref_beadid].resid,
            #     self.system.lipids[beadid_best].resid,
            #     dx[XX], dx[YY], dx[ZZ],
            #     np.asarray(ref_normal),
            #     dprod_best
            # ))

            box.fast_pbc_dx_leaflet(
                    &self.system._lipid_positions[ref_beadid, XX],
                    &self.system._lipid_positions[beadid_best, XX],
                    locnorm_dx,
                    &ref_normal[XX])

            self.system.compute_weighted_average(beadid_best, locnorm_position, locnorm_normal)
            locnorm_beadid = beadid_best


            dprod_val = rvec_dprod(&ref_leaflet._normal[XX],
                                   &ref_normal[XX])

            rvec_copy(locnorm_normal, other_normal)
            rvec_copy(locnorm_position, other_position)
            other_beadid = beadid_best

            if self._is_planar and acos(dprod_val) > 10 / 180 * PI:
                    # print("lipid resid {} normal is off (by {:.1f}°): {} vs avg normal: {}".format(
                    #     self.system.lipids[ref_beadid].resid,
                    #     acos(dprod_val) / np.pi * 180,
                    #     np.asarray(ref_normal),
                    #     np.asarray(self._normal)
                    # ))


                    # Get closest bead
                    d_min = 100000
                    beadid_best = -1

                    for j in range(other_leaflet._size):
                        beadid = other_leaflet._lipid_ids[j]


                        box.fast_pbc_dx_leaflet(
                            &self.system._lipid_positions[ref_beadid, XX],
                            &self.system._lipid_positions[beadid, XX],
                            dx,
                            &ref_normal[XX])

                        norm = rvec_norm(dx)

                        if norm < d_min:
                            beadid_best = beadid
                            d_min = norm

                    box.fast_pbc_dx_leaflet(
                        &self.system._lipid_positions[ref_beadid, XX],
                        &self.system._lipid_positions[beadid_best, XX],
                        closest_dx,
                        &ref_normal[XX])

                    self.system.compute_weighted_average(beadid_best, closest_position, closest_normal)
                    closest_beadid = beadid_best


                    # Get best using leaflet normal
                    beadid_best = -1
                    dprod_best = 0

                    for j in range(other_leaflet._size):
                        beadid = other_leaflet._lipid_ids[j]


                        box.fast_pbc_dx_leaflet(
                            &self.system._lipid_positions[ref_beadid, XX],
                            &self.system._lipid_positions[beadid, XX],
                            dx,
                            &ref_leaflet._normal[XX])

                        dprod_val = rvec_dprod(&ref_leaflet._normal[XX],
                                               dx)

                        norm = rvec_norm(dx)

                        if dprod_val/norm < dprod_best:
                            dprod_best = dprod_val/norm
                            beadid_best = beadid

                    box.fast_pbc_dx_leaflet(
                    &self.system._lipid_positions[ref_beadid, XX],
                    &self.system._lipid_positions[beadid_best, XX],
                    leafnorm_dx,
                    &ref_normal[XX])


                    self.system.compute_weighted_average(beadid_best, leafnorm_position, leafnorm_normal)
                    leafnorm_beadid = beadid_best

                    delta = fabs(rvec_norm(closest_dx) - rvec_norm(locnorm_dx))/rvec_norm(closest_dx)


                    # print("-> Candidates for resid {} (p:[{:.1f},{:.1f},{:.1f}], n: [{:.1f},{:.1f},{:.1f}]):\n  ->locnorm ({}): p:[{:.1f},{:.1f},{:.1f}], n:[{:.1f},{:.1f},{:.1f}], d:{:.1f}\n  ->"
                    #       "leafnorm ({}): p:[{:.1f},{:.1f},{:.1f}], n:[{:.1f},{:.1f},{:.1f}], d:{:.1f}\n  ->"
                    #       "closest ({}): p:[{:.1f},{:.1f},{:.1f}], n:[{:.1f},{:.1f},{:.1f}], d:{:.1f} => delta: {:.3f}\n".format(
                    #     self.system.lipids[ref_beadid].resid,
                    #     ref_position[XX], ref_position[YY], ref_position[ZZ],
                    #     ref_normal[XX], ref_normal[YY], ref_normal[ZZ],
                    #     self.system.lipids[locnorm_beadid].resid,
                    #     locnorm_position[XX], locnorm_position[YY], locnorm_position[ZZ],
                    #     locnorm_normal[XX], locnorm_normal[YY], locnorm_normal[ZZ],
                    #     rvec_norm(locnorm_dx),
                    #     self.system.lipids[leafnorm_beadid].resid,
                    #     leafnorm_position[XX], leafnorm_position[YY], leafnorm_position[ZZ],
                    #     leafnorm_normal[XX], leafnorm_normal[YY], leafnorm_normal[ZZ],
                    #     rvec_norm(leafnorm_dx),
                    #     self.system.lipids[closest_beadid].resid,
                    #     closest_position[XX], closest_position[YY], closest_position[ZZ],
                    #     closest_normal[XX], closest_normal[YY], closest_normal[ZZ],
                    #     rvec_norm(closest_dx),
                    #     delta
                    # ))

                    if delta > 0.2:
                        if delta > 5:
                            self._lipid_thicknesses[internal_id] = NOTSET
                            continue
                        else:
                            rvec_copy(leafnorm_normal, other_normal)
                            rvec_copy(leafnorm_position, other_position)
                            rvec_copy(&ref_leaflet._normal[XX], ref_normal)
                            other_beadid = leafnorm_beadid

            box.fast_pbc_dx_leaflet(
            ref_position,
            other_position,
            dx,
            ref_normal)

            norm = rvec_norm(dx)

            dprod_val = fabs(rvec_dprod(dx, ref_normal))
            other_dprod_val = fabs(rvec_dprod(dx, other_normal))


            self._lipid_thicknesses[internal_id] = 0.5 * (dprod_val + other_dprod_val)

        for internal_id in range(self._size):
            if self._lipid_thicknesses[internal_id] < 0:
                self._lipid_thicknesses[internal_id] = np.nan

        self._thickness = np.nanmean(self._lipid_thicknesses)
        self._lastupdate_thickness = self.system._lastupdate

        return self._thickness, self._lipid_thicknesses

    @property
    def thickness(self):
        self.compute_thickness()
        return self._thickness

    @property
    def lipid_thicknesses(self):
        self.compute_thickness()
        return np.asarray(self._lipid_thicknesses)
