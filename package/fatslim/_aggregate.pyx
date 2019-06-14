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


from ._typedefs cimport real, rvec, fsl_int
from ._typedefs cimport rvec_clear

from ._typedefs cimport rvec_norm


cdef class LipidAggregate:

    def __init__(self, fsl_int[:] lipid_ids, LipidRegistry system):
        self.system = system

        self._lipid_ids = np.sort(lipid_ids)
        self._is_lipid_id_used = np.zeros(self.system._nlipids, dtype=int)
        for beadid in self._lipid_ids:
            self._is_lipid_id_used[beadid] = 1

        self._size = self._lipid_ids.shape[0]

        self._lastupdate = -1

        # planarity and average normal and position
        self._isplanar = False
        self._normal = np.empty(DIM, dtype=np.float32)
        self._position = np.empty(DIM, dtype=np.float32)

        # clusterization
        self._clustered = np.empty(self.system._nlipids, dtype=int)
        self._positions_clustered_buffer = np.empty((self.system._nlipids, DIM), dtype=np.float32)
        self._positions_clustered = np.empty((self._size, DIM), dtype=np.float32)
        self._cluster_stack = np.empty((self.system._nlipids, 2), dtype=int)

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.cdivision(True)
    cdef update(self, force_update=False):
        cdef fsl_int i, j, bead_id
        cdef real norm

        # No need to do anything if the membranes are already identified for the current frame
        if self._lastupdate == self.system._lastupdate and not force_update:
            return

        self.system.update(force_update)

        with nogil:
            self._normal[:] = 0

            for i in range(self._size):
                bead_id = self._lipid_ids[i]

                # print("Bead resid {}: position: {}, normal: {}".format(
                #     self.system.lipids[bead_id].resid,
                #     np.asarray(self.system._lipid_positions[bead_id]),
                #     np.asarray(self.system._lipid_normals[bead_id])
                # ))

                for j in range(DIM):
                    self._normal[j] += self.system._lipid_normals[bead_id, j]
                    self._position[j] += self.system._lipid_positions[bead_id, j]


            self.fast_clusterize()


            norm = rvec_norm(&self._normal[XX])

            # print("\nNormal: {}, norm: {}, size: {} - position: {}".format(
            #     np.asarray(self._normal),
            #     norm,
            #     self._size,
            #     np.asarray(self._position)
            # ))

            if norm > 0.75 * self._size:
                self._isplanar = True
            else:
                self._isplanar = False

            for j in range(DIM):
                self._normal[j] /= norm

            self._lastupdate = self.system._lastupdate

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.cdivision(True)
    cdef void fast_clusterize(self, bint force_update=False) nogil:
        cdef fsl_int i, j
        cdef fsl_int beadid, ref_beadid, nid
        cdef fsl_int stack_size = 0, stack_index = 0
        cdef rvec ref_position
        cdef fsl_int counter
        cdef rvec dx, cluster_centroid, raw_dx

        # No need to do anything if self.update was already called
        if self._lastupdate == self.system._lastupdate and not force_update:
            return

        # Reinitialize clutered-related variables
        for i in range(self._size):
            beadid = self._lipid_ids[i]
            self._clustered[beadid] = 0

            for j in range(DIM):
                self._positions_clustered_buffer[beadid, j] = self.system._lipid_positions[beadid, j]


        self._cluster_stack[stack_size][0] = -1
        self._cluster_stack[stack_size][1] = self._lipid_ids[0]
        stack_size = 1
        stack_index = 0

        counter = 0
        while stack_size > stack_index:

            counter += 1

            ref_beadid = self._cluster_stack[stack_index][0]
            beadid  = self._cluster_stack[stack_index][1]
            stack_index += 1

            if self._clustered[beadid] == 2:
                continue

            if ref_beadid < 0:
                for i in range(DIM):
                    ref_position[i] = self._positions_clustered_buffer[beadid, i]
            else:
                for i in range(DIM):
                    ref_position[i] = self._positions_clustered_buffer[ref_beadid, i]

            self.system.box.fast_pbc_dx(ref_position,
                                        &self._positions_clustered_buffer[beadid, XX],
                                        dx)

            for j in range(DIM):
                self._positions_clustered_buffer[beadid, j] = ref_position[j] + dx[j]

            self._clustered[beadid] = 2

            for j in range(self.system._lipid_neighbours.nneighbours[beadid]):
                nid = self.system._lipid_neighbours.neighbours[beadid][j]

                if self._is_lipid_id_used[nid] == 1:
                    if self._clustered[nid] == 0:
                        self._cluster_stack[stack_size][0] = beadid
                        self._cluster_stack[stack_size][1] = nid
                        stack_size += 1

                        self._clustered[nid] = 1

        # Get cluster centroid
        rvec_clear(cluster_centroid)
        for i in range(self._size):
            beadid = self._lipid_ids[i]

            for j in range(DIM):
                cluster_centroid[j] += self._positions_clustered_buffer[beadid, j]
                self._positions_clustered[i][j] = self._positions_clustered_buffer[beadid][j]
        for j in range(DIM):
            cluster_centroid[j] /= self._size

        for i in range(DIM):
            self._position[i] = cluster_centroid[i]


    def __len__(self):
        return self._size

    @property
    def indices(self):
        return np.asarray(self._lipid_ids, dtype=int)

    @property
    def lipids(self):
        lipids = []
        for i in self.indices:
            lipids.append(self.system.lipids[i])
        return lipids

    @property
    def size(self):
        return self._size

    def __getitem__(self, item):
        return self.indices[item]

    @property
    def is_planar(self):
        self.update()
        return self._isplanar

    @property
    def normal(self):
        self.update()
        return np.asarray(self._normal)

    @property
    def centroid(self):
        self.update()
        return np.asarray(self._position)

    @property
    def position(self):
        return self.centroid

    @property
    def positions(self):
        self.update()
        return np.array(self._positions_clustered)

    @property
    def positions_raw(self):
        self.update()
        return np.asarray(self.system.lipid_positions[self._lipid_ids])

    def __str__(self):
        return "Lipid aggregate made of {} lipids".format(self._size)

    def __repr__(self):
        return "<LipidAggregate with {} lipids>".format(self._size)