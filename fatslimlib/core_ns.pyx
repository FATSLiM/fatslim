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
#cython: cdivision=True
#cython: boundscheck=False

# Preprocessor DEFs
DEF DIM = 3
DEF XX = 0
DEF YY = 1
DEF ZZ = 2
DEF RET_OK = 1
DEF RET_ERROR = 0
DEF EPSILON = 1e-5
DEF NEIGHBORHOOD_ALLOCATION_INCREMENT = 10
DEF GRID_ALLOCATION_INCREMENT = 10

from libc.stdio cimport fprintf, stderr

# Cython C imports (no Python here!)
from cython.parallel cimport prange
from libc.stdlib cimport malloc, realloc, free, abort

from .typedefs cimport real, rvec, matrix, rvec_norm2, fsl_int
from .core_ns cimport PBCBox
from .core_base cimport OPENMP_NUM_THREADS

# Python imports
import numpy as np

########################################################################################################################
#
# Neighbor Search Stuff
#
########################################################################################################################
cdef struct ns_grid:
    fsl_int size
    fsl_int[DIM] ncells
    real[DIM] cellsize
    fsl_int *nbeads
    fsl_int **beadids

cdef ns_grid initialize_nsgrid(matrix box,
                               float cutoff) nogil:
    cdef ns_grid grid
    cdef fsl_int i

    for i in range(DIM):
        grid.ncells[i] = <fsl_int> (box[i][i] / cutoff)
        if grid.ncells[i] == 0:
            grid.ncells[i] = 1
        grid.cellsize[i] = box[i][i] / grid.ncells[i]

    grid.size = grid.ncells[XX] * grid.ncells[YY] * grid.ncells[ZZ]
    return grid

cdef fsl_int populate_grid(ns_grid *grid,
                       real[:,::1] coords) nogil:
    cdef fsl_int ncoords = coords.shape[0]
    cdef bint ret_val

    ret_val = populate_grid_array(grid,
                                  <rvec *> &coords[0, 0],
                                  ncoords)

    return ret_val

cdef fsl_int populate_grid_array(ns_grid *grid,
                             rvec *coords,
                             fsl_int ncoords) nogil:
    cdef fsl_int i, cellindex = -1
    cdef fsl_int grid_size = grid.size
    cdef fsl_int *allocated_size = NULL

    if grid_size != grid.ncells[XX] * grid.ncells[YY] * grid.ncells[ZZ]: # Grid not initialized
        return RET_ERROR

    # Allocate memory
    grid.nbeads = <fsl_int *> malloc(sizeof(fsl_int) * grid_size)
    if grid.nbeads == NULL:
        fprintf(stderr,"FATAL: Could not allocate memory for NS grid.nbeads (requested: %i bytes)\n",
                sizeof(fsl_int) * grid_size)
        abort()

    allocated_size = <fsl_int *> malloc(sizeof(fsl_int) * grid_size)
    if allocated_size == NULL:
        fprintf(stderr,"FATAL: Could not allocate memory for NS allocated_size (requested: %i bytes)\n",
                sizeof(fsl_int) * grid_size)
        abort()

    for i in range(grid_size):
        grid.nbeads[i] = 0
        allocated_size[i] = GRID_ALLOCATION_INCREMENT

    grid.beadids = <fsl_int **> malloc(sizeof(fsl_int *) * grid_size)
    if grid.beadids == NULL:
        fprintf(stderr,"FATAL: Could not allocate memory for NS grid.beadids (requested: %i bytes)\n",
                sizeof(fsl_int *) * grid_size)
        abort()

    for i in range(grid_size):
        grid.beadids[i] = <fsl_int *> malloc(sizeof(fsl_int) * allocated_size[i])
        if grid.beadids[i] == NULL:
            fprintf(stderr,"FATAL: Could not allocate memory for NS grid.beadids[i] (requested: %i bytes)\n",
                sizeof(fsl_int) * allocated_size[i])
            abort()

    # Get cell indices for coords
    for i in range(ncoords):
        cellindex = <fsl_int> (coords[i][ZZ] / grid.cellsize[ZZ]) * (grid.ncells[XX] * grid.ncells[YY]) +\
                    <fsl_int> (coords[i][YY] / grid.cellsize[YY]) * grid.ncells[XX] + \
                    <fsl_int> (coords[i][XX] / grid.cellsize[XX])

        grid.beadids[cellindex][grid.nbeads[cellindex]] = i
        grid.nbeads[cellindex] += 1

        if grid.nbeads[cellindex] >= allocated_size[cellindex]:
            allocated_size[cellindex] += GRID_ALLOCATION_INCREMENT
            grid.beadids[cellindex] = <fsl_int *> realloc(<void *> grid.beadids[cellindex], sizeof(fsl_int) * allocated_size[cellindex])

    free(allocated_size)
    return RET_OK

cdef void destroy_nsgrid(ns_grid *grid) nogil:
    cdef fsl_int i
    if grid.nbeads != NULL:
        free(grid.nbeads)

    for i in range(grid.size):
        if grid.beadids[i] != NULL:
            free(grid.beadids[i])
    free(grid.beadids)


cdef ns_neighborhood_holder *create_neighborhood_holder() nogil:
    cdef ns_neighborhood_holder *holder

    holder = <ns_neighborhood_holder *> malloc(sizeof(ns_neighborhood_holder))

    return holder

cdef void free_neighborhood_holder(ns_neighborhood_holder *holder) nogil:
    cdef fsl_int i

    if holder == NULL:
        return

    for i in range(holder.size):
        if holder.neighborhoods[i].beadids != NULL:
            free(holder.neighborhoods[i].beadids)
        free(holder.neighborhoods[i])
    free(holder.neighborhoods)
    free(holder)

cdef ns_neighborhood *retrieve_neighborhood(rvec current_coords, real[:, ::1]neighborcoords, ns_grid *grid, PBCBox box, real cutoff2) nogil:
    cdef fsl_int d, m
    cdef fsl_int xi, yi, zi, bid
    cdef real d2
    cdef rvec shifted_coord, dx, neighbor_coord, corrected_coords

    cdef bint already_checked[27]
    cdef bint skip
    cdef fsl_int nchecked = 0, icheck
    cdef fsl_int cell_index

    cdef ns_neighborhood *neighborhood = <ns_neighborhood *> malloc(sizeof(ns_neighborhood))
    if neighborhood == NULL:
        abort()

    neighborhood.size = 0
    neighborhood.allocated_size = NEIGHBORHOOD_ALLOCATION_INCREMENT
    neighborhood.beadids = <fsl_int *> malloc(NEIGHBORHOOD_ALLOCATION_INCREMENT * sizeof(fsl_int))
    if neighborhood.beadids == NULL:
        abort()

    for zi in range(3):
        for yi in range(3):
            for xi in range(3):
                # Calculate and/or reinitialize shifted coordinates
                shifted_coord[XX] = current_coords[XX] + (xi - 1) * grid.cellsize[XX]
                shifted_coord[YY] = current_coords[YY] + (yi - 1) * grid.cellsize[YY]
                shifted_coord[ZZ] = current_coords[ZZ] + (zi - 1) * grid.cellsize[ZZ]

                # Make sure the shifted coordinates is inside the brick-shaped box
                for m in range(DIM - 1, -1, -1):

                    while shifted_coord[m] < 0:
                        for d in range(m+1):
                            shifted_coord[d] += box.c_pbcbox.box[m][d]


                    while shifted_coord[m] >= box.c_pbcbox.box[m][m]:
                        for d in range(m+1):
                            shifted_coord[d] -= box.c_pbcbox.box[m][d]

                # Get the cell index corresponding to the coord
                cell_index = <fsl_int> (shifted_coord[ZZ] / grid.cellsize[ZZ]) * (grid.ncells[XX] * grid.ncells[YY]) +\
                             <fsl_int> (shifted_coord[YY] / grid.cellsize[YY]) * grid.ncells[XX] + \
                             <fsl_int> (shifted_coord[XX] / grid.cellsize[XX])

                # Just a safeguard
                if cell_index >= grid.size:
                    continue

                # Check the cell index was not already selected
                skip = False
                for icheck in range(nchecked):
                    if already_checked[icheck] == cell_index:
                        skip = True
                        break
                if skip:
                    continue

                # Search for neighbors inside this cell
                for i_bead in range(grid.nbeads[cell_index]):
                    bid = grid.beadids[cell_index][i_bead]

                    box.fast_pbc_dx(current_coords, &neighborcoords[bid, XX], dx)

                    d2 = rvec_norm2(dx)

                    if d2 < cutoff2:
                        if d2 < EPSILON: # Don't add the current bead as its own neighbor!
                            continue

                        # Update neighbor lists
                        neighborhood.beadids[neighborhood.size] = bid
                        neighborhood.size += 1
                        if neighborhood.size >= neighborhood.allocated_size:
                            neighborhood.allocated_size += NEIGHBORHOOD_ALLOCATION_INCREMENT
                            neighborhood.beadids = <fsl_int *> realloc(<void*> neighborhood.beadids, neighborhood.allocated_size * sizeof(fsl_int))
                            if neighborhood.beadids == NULL:
                                abort()

                # Register the cell as checked
                already_checked[nchecked] = cell_index
                nchecked += 1

    return neighborhood


cdef ns_neighborhood_holder *ns_core(real[:, ::1] refcoords,
                                     real[:, ::1] neighborcoords,
                                     ns_grid *grid,
                                     PBCBox box,
                                     real cutoff) nogil:
    cdef fsl_int coordid, i, j
    cdef fsl_int ncoords = refcoords.shape[0]
    cdef fsl_int ncoords_neighbors = neighborcoords.shape[0]
    cdef real cutoff2 = cutoff * cutoff
    cdef ns_neighborhood_holder *holder

    cdef fsl_int *neighbor_buf
    cdef fsl_int buf_size, ibuf

    holder = create_neighborhood_holder()
    if holder == NULL:
        fprintf(stderr,"FATAL: Could not allocate memory for NS holder\n",
                sizeof(fsl_int) * ncoords)
        abort()

    holder.size = ncoords
    holder.neighborhoods = <ns_neighborhood **> malloc(sizeof(ns_neighborhood *) * ncoords)
    if holder.neighborhoods == NULL:
        fprintf(stderr,"FATAL: Could not allocate memory for NS holder.neighborhoods (requested: %i bytes)\n",
                sizeof(ns_neighborhood) * ncoords)
        abort()

    # Here starts the real core and the iteration over coordinates
    for coordid in prange(ncoords, schedule='dynamic', num_threads=OPENMP_NUM_THREADS):
        holder.neighborhoods[coordid] = retrieve_neighborhood(&refcoords[coordid, XX],
                                                              neighborcoords,
                                                              grid,
                                                              box,
                                                              cutoff2)
        holder.neighborhoods[coordid].cutoff = cutoff

    return holder

cdef ns_neighborhood_holder *fast_neighbor_search(real[:, ::1] ref_coords,
                                                  real[:, ::1] neighbor_coords,
                                                  PBCBox box,
                                                  real cutoff) nogil:
    cdef ns_grid grid_coords
    cdef ns_neighborhood_holder *holder


    # Check number of coordinates
    if neighbor_coords.shape[0] == 0 or ref_coords.shape[0] == 0:
        return NULL

    # Check cutoff
    if cutoff <= 0:
        return NULL

    grid_coords = initialize_nsgrid(box.c_pbcbox.box, cutoff)

    if populate_grid(&grid_coords, neighbor_coords) != RET_OK:
        destroy_nsgrid(&grid_coords)
        return NULL

    # Retrieve neighbors from grid
    holder = ns_core(ref_coords, neighbor_coords, &grid_coords, box, cutoff)

    # Free memory
    destroy_nsgrid(&grid_coords)

    return holder

def neighbor_search(PBCBox box, real[:, ::1] ref_coords, real[:, ::1] neighbor_coords=None, real cutoff=1.0):
    cdef real[:, ::1] ref_coords_bbox, neighbor_coords_box
    cdef fsl_int nid, i
    cdef ns_neighborhood_holder *holder
    cdef ns_neighborhood *neighborhood

    print "Performing a neighbor search on %i beads using a cutoff of %.2f (%i threads used)" % (ref_coords.shape[0],
    cutoff, OPENMP_NUM_THREADS)

    # Make sure atoms are inside the brick-shaped box
    ref_coords_bbox = box.fast_put_atoms_in_bbox(ref_coords)

    if neighbor_coords is None:
        neighbor_coords_box = ref_coords_bbox
    else:
        neighbor_coords_box = box.fast_put_atoms_in_bbox(neighbor_coords)

    with nogil:
        holder = fast_neighbor_search(ref_coords_bbox,
                                      neighbor_coords_box,
                                      box,
                                      cutoff)

    neighbors = []
    for nid in range(holder.size):
        neighborhood = holder.neighborhoods[nid]
        neighborhood_py = np.empty(neighborhood.size, dtype=np.int64)
        for i in range(neighborhood.size):
            neighborhood_py[i] = neighborhood.beadids[i]
        neighbors.append(neighborhood_py)


    # Free Memory
    free_neighborhood_holder(holder)
    return neighbors
