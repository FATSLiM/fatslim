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

from ._typedefs cimport real, rvec, fsl_int, matrix
from ._typedefs cimport rvec_dprod, rvec_copy

from ._typedefs cimport rvec_norm
from ._typedefs cimport mat_from_rvec, invert_mat

from ._core cimport SimplifiedLipid

from ._aggregate cimport LipidAggregate

from ._geometry cimport PBCBox
from ._geometry cimport complete_basis, rvec_to_basis
from ._geometry cimport real_point, Polygon, polygon_new, polygon_get_area, fast_clip_zoi, polygon_destroy, polygon_append, polygon_empty

from libc.math cimport acos



cdef class Membrane(object):

    def __init__(self, leaflet1, leaflet2):
        self.system = leaflet1.system

        if leaflet1.system != leaflet2.system:
            raise ValueError("Leaflets do not belong to the same system")

        self._leaflets = list()

        if (leaflet1.is_planar and leaflet1.position[2] > leaflet2.position[2]) or not leaflet1.is_planar:
            self._leaflets.append(Leaflet(leaflet1, self, 0))
            self._leaflets.append(Leaflet(leaflet2, self, 1))
        else:
            self._leaflets.append(Leaflet(leaflet2, self, 0))
            self._leaflets.append(Leaflet(leaflet1, self, 1))

    def __getitem__(self, item):
        if item in (-2, 0):
            return self._leaflets[0]
        elif item in (-1, 1):
            return self._leaflets[1]
        else:
            raise IndexError("A membrane contains two leaflets")



cdef class Leaflet(LipidAggregate):

    def __init__(self, ids_or_aggregate, Membrane membrane, fsl_int leaflet_id):
        cdef SimplifiedLipid lipid

        if isinstance(ids_or_aggregate, LipidAggregate):
            ids_or_aggregate = ids_or_aggregate.indices

        system = membrane.system

        super().__init__(ids_or_aggregate, system)

        for i, val in enumerate(self.indices):
            lipid = self.system.lipids[val]

            lipid._regid_leaflet = i
            lipid._leaflet = self

        self._membrane = membrane
        self._leaflet_id = leaflet_id

        # Thickness
        self._lastupdate_thickness = -1
        self._thickness = 0
        self._lipid_thicknesses = np.empty(self._size, dtype=np.float32)
        self._lipid_interleaflet_gaps = np.empty(self._size, dtype=np.float32)

        # APL
        self._lastupdate_apl = -1
        self._apl = 0
        self._lipid_apls = np.empty(self._size, dtype=np.float32)
        self._area = 0


    cdef compute_thickness(self, bint force_update=False):
        cdef fsl_int i, j
        cdef fsl_int ref_beadid, beadid
        cdef Leaflet other_leaflet
        cdef PBCBox box
        cdef rvec dx
        cdef rvec ref_normal, ref_position

        cdef fsl_int locnorm_beadid
        cdef rvec locnorm_position, locnorm_normal, locnorm_dx

        cdef fsl_int closest_beadid
        cdef rvec closest_position, closest_normal, closest_dx

        cdef fsl_int leafnorm_beadid
        cdef rvec leafnorm_position, leafnorm_normal, leafnorm_dx

        cdef fsl_int other_beadid
        cdef rvec other_position, other_normal

        cdef real[:,::] universe_positions


        if self._membrane is None:
            raise ValueError("Can not calculate thickness if leaflet does not belong to a membrane")

        # No need to do anything if the membranes are already identified for the current frame
        if self._lastupdate_thickness == self.system._lastupdate and not force_update:
            return self._thickness, self._lipid_thicknesses


        self.update(force_update)
        box = self.system.box
        universe_positions = self.system.universe_coords_bbox

        if self._leaflet_id == 0:
            other_leaflet = self._membrane._leaflets[1]
        else:
            other_leaflet = self._membrane._leaflets[0]

        for i in range(self._size):
            ref_beadid = self._lipid_ids[i]

            # Get average normal and position for the reference bead
            self.system.compute_weighted_average(ref_beadid, ref_position, ref_normal)

            # print("Weighted position for resid {}: [{:.3f}, {:.3f}, {:.3f}], weighted normal: [{:.3f}, {:.3f}, {:.3f}] - actual position: {}".format(
            #     self.system.lipids[ref_beadid].resid,
            #     ref_position[XX], ref_position[YY], ref_position[ZZ],
            #     ref_normal[XX], ref_normal[YY], ref_normal[ZZ],
            #     np.asarray(self.system._lipid_positions[ref_beadid])
            # ))

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


            dprod_val = rvec_dprod(&self._normal[XX],
                                   &ref_normal[XX])

            rvec_copy(locnorm_normal, other_normal)
            rvec_copy(locnorm_position, other_position)
            other_beadid = beadid_best

            if self._isplanar and acos(dprod_val) > 10 / 180 * np.pi:
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
                            &self._normal[XX])

                        dprod_val = rvec_dprod(&self._normal[XX],
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

                    delta = abs(rvec_norm(closest_dx) - rvec_norm(locnorm_dx))/rvec_norm(closest_dx)


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
                            self._lipid_thicknesses[i] = np.nan
                            continue
                        else:
                            rvec_copy(leafnorm_normal, other_normal)
                            rvec_copy(leafnorm_position, other_position)
                            rvec_copy(&self._normal[XX], ref_normal)
                            other_beadid = leafnorm_beadid

            box.fast_pbc_dx_leaflet(
            ref_position,
            other_position,
            dx,
            ref_normal)

            norm = rvec_norm(dx)

            dprod_val = abs(rvec_dprod(dx, ref_normal))
            other_dprod_val = abs(rvec_dprod(dx, other_normal))


            self._lipid_thicknesses[i] = 0.5 * (dprod_val + other_dprod_val)

            # print("Thickness for resid {}: ref {:.3f}, dx norm: {:.3f}, other: {:.3f}, avg:{:.3f}".format(
            #     self.system.lipids[ref_beadid].resid,
            #     dprod_val,
            #     norm,
            #     abs(rvec_dprod(dx, other_normal)),
            #     0.5 * (dprod_val + abs(rvec_dprod(dx, other_normal)))
            # ))

            # Compute interleaflet gap
            ref_d = 0

            for ix in self.system.lipids[ref_beadid]._ix:
                box.fast_pbc_dx(
                    &self.system._lipid_positions[ref_beadid, XX],
                    &universe_positions[ix, XX],
                    dx)

                dprod_val = rvec_dprod(dx, ref_normal)

                # print("d ref atom-atom {}: {:.3f} (dz={:.3f})".format(
                #     self.system.universe.atoms[ix].name,
                #     dprod_val,
                #     self.system._lipid_positions[ref_beadid, ZZ] - self.system.universe.atoms[ix].position[ZZ]
                # ))

                if dprod_val < ref_d:
                    ref_d = dprod_val
            # print("")

            other_d = 0
            for ix in self.system.lipids[other_beadid]._ix:
                box.fast_pbc_dx_leaflet(
                    &self.system._lipid_positions[ref_beadid, XX],
                    &universe_positions[ix, XX],
                    dx,
                ref_normal)

                dprod_val = rvec_dprod(dx, ref_normal)

                if dprod_val < other_d:
                    other_d = dprod_val

            self._lipid_interleaflet_gaps[i] = 0.5 * self._lipid_thicknesses[i] - np.abs(other_d)


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
        return np.asarray(self._lipid_thicknesses).copy()

    @property
    def lipid_interleaflet_gaps(self):
        self.compute_thickness()
        return np.asarray(self._lipid_interleaflet_gaps).copy()


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

        if self._membrane is None:
            raise ValueError("Can not calculate APL if leaflet does not belong to a membrane")

        # No need to do anything if the membranes are already identified for the current frame
        if self._lastupdate_apl == self.system._lastupdate and not force_update:
            return self._apl, self._lipid_apls


        self.update(force_update)
        box = self.system.box

        cell = polygon_new()

        ref_point_2d[XX] = 0
        ref_point_2d[YY] = 0

        for i in range(self._size):
            ref_beadid = self._lipid_ids[i]

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
        self._apl = np.mean(self._lipid_apls)
        self._area = np.sum(self._lipid_apls)

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