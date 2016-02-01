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
#cython: profile=False

from __future__ import print_function

# C libraries
from cython.parallel cimport prange, parallel
from libc.stdio cimport FILE, fopen, fclose, fseek, ftell, SEEK_SET, fgets, printf, fread
from libc.string cimport strncpy, strlen, strcmp
from libc.stdlib cimport atof, atoi
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from libc.stdlib cimport free, malloc

from . import core_base
cimport core_base
from .core_base cimport OPENMP_NUM_THREADS
from typedefs cimport real, rvec, fsl_int, fsl_uint
cimport cython
import numpy as np

cdef extern from "xdrfile/xdrfile.h":
    ctypedef struct XDRFILE:
        pass

    # XDRFile I/O
    XDRFILE *xdrfile_open(const char *path, const char *mode) nogil
    int xdrfile_close(XDRFILE *xfp) nogil
    int xdrfile_getpos(XDRFILE* xfp) nogil
    int xdrfile_setpos(XDRFILE* xfp, unsigned int pos) nogil

    # XDRFile Reading
    int xdrfile_decompress_coord_double_partial(double *ptr,
                                             int *size,
                                             int *ncoords,
                                             double *precision,
                                             XDRFILE* xfp) nogil
    int xdrfile_read_float(float *           ptr,
					       int               ndata,
                           XDRFILE *         xfp) nogil
    int xdrfile_read_double(double *           ptr,
					        int               ndata,
                            XDRFILE *         xfp) nogil
    int xdrfile_read_int(int *         ptr,
					     int           ndata,
					     XDRFILE *     xfp) nogil



# Python libraries
import numpy
import os

# Constants
DEF XTC_MAGIC = 1995
DEF TRR_MAGIC = 1993
DEF TRR_HEADER_SIZE = 76
DEF XTC_HEADER_SIZE = 56
DEF DIM = 3
DEF XX = 0
DEF YY = 1
DEF ZZ = 2
DEF RET_OK = 1
DEF RET_YES = 1
DEF RET_ERROR = 0
DEF RET_NO = 0

ctypedef float fvec[DIM]

cdef bint is_solvent(bytes resname):

    if resname == b"W":
        return True

    if resname == b"PW":
        return True

    if resname == b"SOL":
        return True

    if resname == b"WAT":
        return True

    if resname == b"HOH":
        return True

    if resname == b"OHH":
        return True

    if resname == b"TIP":
        return True

    if resname == b"T3P":
        return True

    if resname == b"T4P":
        return True

    if resname == b"T5P":
        return True

    if resname == b"T3H":
        return True
    return False

cdef class GroReaderTopol(core_base.TopologyReader):
    cdef load(self):
        cdef int lino = -1, natoms = -1, last_resid = -1
        cdef int resid, atomid, resid_offset = 0
        cdef bytes resname
        cdef str atomname
        cdef bint skip = False

        topol = self.topology
        with open(self.filename, "r") as fp:
            for line in fp:
                lino += 1
                if lino == 0:
                    continue
                elif lino == 1:
                    natoms = int(line)
                else:
                    if lino < natoms + 2:
                        resid = atoi(line[:5].encode())

                        # Correct the resid if necessary because resids are modulo 10000 in .gro
                        if (resid+resid_offset * 10000) < last_resid:
                            resid_offset += 1
                        resid += resid_offset * 10000

                        if resid == last_resid and skip:
                            if skip: # Already selected to be skipped
                                continue
                        else:
                            if resid != last_resid:
                                skip = False

                            # Update the last_resid counter
                            last_resid = resid

                            # Read residue name
                            resname = line[5:10].strip().encode()

                            if is_solvent(resname): # Classified as solvent, we skip it
                                residue = None
                                skip = True
                                #print("Skipping %s - resid:%i" % (resname, resid))
                                continue

                            # Read atom information
                            atomname = line[10:15].strip()
                            #atomname = strip_name(line[10:15])

                            if atomname[0:1] == "H": # Labelled as hydrogen -> skipped
                                continue

                            atomid = atoi(line[15:20].encode()) + ((lino - 1)// 100000) * 100000

                            topol.internal_append(resid, resname, atomname.encode(), atomid)




_Topology_readers = {".gro": GroReaderTopol}

cdef class NdxReader(core_base.IndexReader):
    cdef fast_load(self):
        cdef str line
        cdef bytes last_group = b""

        if self.filename.split(b".")[-1] != b"ndx":
            raise ValueError("%s is not a supported index file. supported: .ndx" % self.filename)
        groups = {}

        with open(self.filename, "r") as fp:

            # Use Python object to parse line... Much more easier!
            for line in fp:
                if line[0] == '[':
                    if last_group != b"":
                        groups[last_group] = numpy.fromstring(groups[last_group], dtype=np.int64, sep=" ")
                    last_group = line.strip(" []\n").encode()
                    groups[last_group] = b""
                else:
                    groups[last_group] += line.strip().encode() + b" "

        groups[last_group] = numpy.fromstring(groups[last_group], dtype=np.int64, sep=" ")

        self.groups = groups

_Index_loaders = {".ndx": NdxReader}


cdef class GroReaderCoords(core_base.CoordinateReader):
    cdef fsl_int coordline_length

    cdef preload(self):
        cdef int lino = -1, natoms = -1
        cdef bint need_read = True

        self.nframes = 1
        self.timesteps = numpy.array([0], dtype=np.float64)

        if self.filename.decode().split(".")[-1] != "gro":
            raise ValueError("%s is not a supported gro file. supported: .gro" % self.filename)

        with open(self.filename, "r") as fp:
            while need_read:
                lino += 1
                line = fp.readline()
                if lino == 0:
                    continue
                elif lino == 1:
                    natoms = int(line)
                    self.coordinate_offsets = numpy.array((fp.tell(),), dtype=np.int64)
                elif lino == 2:
                    self.coordline_length = len(line)
                elif lino == (natoms + 1):
                    self.box_offsets = numpy.array((fp.tell(),), dtype=np.int64)
                    need_read = False



    cdef core_base.PBCBox load_box(self, int frame_id):
        cdef int n_coords, i, j
        self.assert_frame_id(frame_id)

        box = numpy.zeros((3,3), np.float64)

        with open(self.filename, "rb") as fp:
            with cython.boundscheck(True):
                fp.seek(self.box_offsets[frame_id])
            box_str = fp.read().strip(b'\n\r').split()

        n_coords = len(box_str)

        box[0, 0] = float(box_str[0])
        box[1, 1] = float(box_str[1])
        box[2, 2] = float(box_str[2])

        if n_coords == 9:
            box[0, 1] = float(box_str[3])
            box[0, 2] = float(box_str[4])
            box[1, 0] = float(box_str[5])
            box[1, 2] = float(box_str[6])
            box[2, 0] = float(box_str[7])
            box[2, 1] = float(box_str[8])

        return core_base.PBCBox(box)

    cdef real[:,::1] load_coords(self, int frame_id, fsl_int[:] atomids) nogil except *:
        cdef int natoms = atomids.shape[0]
        cdef int i
        cdef FILE* cfile
        cdef real[:,::1] coords_memview
        cdef char* cfilename
        cdef fsl_int coords_offset
        cdef fsl_int coord_offset
        cdef int min_offset
        cdef int max_offset
        cdef int offset_range
        cdef char *coords_buffer
        cdef char *coord_buffer
        cdef double val

        with gil:
            coords_memview = numpy.empty((natoms, 3))
            cfilename = self.filename
            coords_offset = self.coordinate_offsets[frame_id]

        # Do nothing if there is no atom to read!
        if natoms == 0:
             return coords_memview


        min_offset = atomids[0]
        max_offset = atomids[0]
        for i in range(natoms):
            if atomids[i] < min_offset:
                min_offset = atomids[i]
            elif atomids[i] > max_offset:
                max_offset = atomids[i]

        offset_range = (max_offset - min_offset + 1) * self.coordline_length # +1 is because we also the last line!


        with gil:
            coords_buffer = <char *> PyMem_Malloc(offset_range * sizeof(char))

        # Open coords file
        cfile = fopen(cfilename, "rb")
        if cfile == NULL:
            with gil:
                raise IOError("No such file or directory: '%s'" % self.filename)

        # Seek the proper position
        fseek(cfile, coords_offset + (min_offset - 1) * self.coordline_length, SEEK_SET) # NOTE: atomid starts @ 1

        # Read buffer & close the file
        if fread(coords_buffer, offset_range, 1, cfile) != 1:
            fclose(cfile)
            with gil:
                raise RuntimeError("EOF reached!")
        fclose(cfile)

        with parallel(num_threads=OPENMP_NUM_THREADS):
            coord_buffer = <char *> malloc(9 * sizeof(char))
            # Initialize buffer end!
            coord_buffer[8] = b'\0'

            for i in prange(natoms, schedule='static'):
                coord_offset = (atomids[i] - min_offset) * self.coordline_length

                # Extract coordinates from buffer
                strncpy(coord_buffer, coords_buffer+coord_offset+20, 8)
                coords_memview[i, XX] = atof(coord_buffer)

                strncpy(coord_buffer, coords_buffer+coord_offset+28, 8)
                coords_memview[i, YY] = atof(coord_buffer)

                strncpy(coord_buffer, coords_buffer+coord_offset+36, 8)
                coords_memview[i, ZZ] = atof(coord_buffer)

            free(coord_buffer)

        with gil:
            PyMem_Free(coords_buffer)

        return coords_memview


cdef class XtcReaderCoords(core_base.CoordinateReader):
    cdef preload(self):
        cdef fsl_uint offset=0
        cdef fsl_uint max_offset = os.path.getsize(self.filename)
        cdef int magic=-1, num_atoms_frame, frame_size, num_coords_bytes, nframes = 0
        cdef float timestep_float
        cdef XDRFILE *xfp

        frame_offsets = []
        timesteps = []

        # Open trajectory
        xfp = xdrfile_open(self.filename, "rb")
        if xfp == NULL:
            raise IOError("Could not read file: %s" % self.filename)

        # Loop over all available frames
        while offset < max_offset:
            # Get to the proper offset
            if xdrfile_setpos(xfp, offset) == 0:
                raise IOError("Could not set position in file: %s" % self.filename)

            # Initialize the frame length
            frame_size = XTC_HEADER_SIZE

            # Check the magic number
            xdrfile_read_int(&magic, 1, xfp)
            if magic != XTC_MAGIC:
                raise IOError("Corrupted XTC file: '%s'" % self.filename)

            # Read timestep
            if xdrfile_setpos(xfp, offset + 12) == 0:
                raise IOError("Could not set position in file: %s" % self.filename)
            xdrfile_read_float(&timestep_float, 1, xfp)

            # Read number of atoms
            if xdrfile_setpos(xfp, offset + 52) == 0:
                raise IOError("Could not set position in file: %s" % self.filename)
            xdrfile_read_int(&num_atoms_frame, 1, xfp)

            if num_atoms_frame <= 9: # No compression used
                frame_size += num_atoms_frame * 3 * 4
            else:
                # Read number of bytes used to store coordinates
                if xdrfile_setpos(xfp, offset + 88) == 0:
                    raise IOError("Could not set position in file: %s" % self.filename)
                xdrfile_read_int(&num_coords_bytes, 1, xfp)

                # As coords bytes might not a factor of 4, we need to make sure of that.
                # NOTE: this is due to xdr file spec, not related to gmx stuff.
                if num_coords_bytes % 4 != 0:
                    num_coords_bytes = 4 * (num_coords_bytes // 4 + 1)

                frame_size += 36 + num_coords_bytes

            # Store data
            frame_offsets.append(offset)
            timesteps.append(timestep_float)

            # Update frame counter and offset
            nframes += 1
            offset += frame_size

        # Close file and release memory
        xdrfile_close(xfp)

        frame_offsets = numpy.array(frame_offsets, dtype=np.int64)
        timesteps = numpy.array(timesteps, dtype=np.float64)
        self.timesteps = timesteps
        self.coordinate_offsets = frame_offsets + 52
        self.box_offsets = frame_offsets + 16
        self.nframes = nframes

    cdef core_base.PBCBox load_box(self, int frame_id):
        cdef XDRFILE *xfp
        cdef fsl_uint pos = self.box_offsets[frame_id]
        cdef real[:,::1] box = numpy.empty((DIM, DIM), dtype=np.float64)
        cdef float float_box[DIM][DIM]
        cdef int i,j

        xfp = xdrfile_open(self.filename, "rb")
        if xfp == NULL:
            raise IOError("Could not read file: %s" % self.filename)
        if xdrfile_setpos(xfp, pos) == 0:
            raise IOError("Could not set position in file: %s" % self.filename)
        if xdrfile_read_float(float_box[0], DIM*DIM, xfp) != DIM*DIM:
            raise IOError("Could not read matrix")
        # Convert float to real
        for i in range(DIM):
            for j in range(DIM):
                box[i, j] = <real> float_box[i][j]
        # Close file and release memory
        xdrfile_close(xfp)

        return core_base.PBCBox(box)



    cdef real[:, ::1] load_coords(self, int frame_id, fsl_int[:] atomids) nogil except *:
        cdef XDRFILE *xfp
        cdef fsl_uint pos
        cdef int size = atomids.shape[0], total_size
        cdef int i
        cdef fsl_int tmpid, maxid = -1
        cdef int natoms, ndecompressed
        cdef real[:, ::1] coords
        cdef rvec *real_coords
        cdef real prec = 0.0
        cdef char *filename

        with gil:
            coords = numpy.empty((size, DIM))
            pos = self.coordinate_offsets[frame_id]
            filename = self.filename

        # Retrieve the max atomid to get how many coordinates needs to be retrieved
        for i in range(size):
            tmpid = atomids[i]
            if tmpid > maxid:
                maxid = tmpid
        natoms = maxid

        # Open traj file
        xfp = xdrfile_open(filename, "rb")
        if xfp == NULL:
            with gil:
                raise IOError("Could not read file: %s" % self.filename)
        if xdrfile_setpos(xfp, pos) == 0:
            with gil:
                raise IOError("Could not set position in file: %s" % self.filename)

        # Get the number of coordinates in frame
        xdrfile_read_int(&total_size, 1, xfp)

        # Allocate memory
        with gil:
            real_coords = <rvec *> PyMem_Malloc(sizeof(rvec) * total_size)

        # Reset position to the beginning of the frame
        if xdrfile_setpos(xfp, pos) == 0:
            with gil:
                raise IOError("Could not set position in file: %s" % self.filename)

        # Load coordinates
        if xdrfile_decompress_coord_double_partial(real_coords[0], &total_size, &natoms, &prec, xfp) < -1:
            with gil:
                raise IOError("Could not decompress coordinates")

        # Close file
        xdrfile_close(xfp)

        # Build the output array
        for i in range(size):
            tmpid = atomids[i] - 1 # Reminder: atomids start at 1

            coords[i, XX] = real_coords[tmpid][XX]
            coords[i, YY] = real_coords[tmpid][YY]
            coords[i, ZZ] = real_coords[tmpid][ZZ]

        # Free memory
        with gil:
            PyMem_Free(real_coords)

        return coords

cdef class TrrReaderCoords(core_base.CoordinateReader):
    cdef bint use_double
    cdef int natoms
    cdef preload(self):
        cdef fsl_uint offset=0
        cdef fsl_uint max_offset = os.path.getsize(self.filename)
        cdef int magic=-1, num_atoms_frame, frame_size, num_coords_bytes, nframes = 0
        cdef float timestep_float
        cdef real timestep
        cdef XDRFILE *xfp
        cdef char *version = "GMX_trn_file"
        cdef int ir_size, e_size, box_size, vir_size, pres_size, top_size, sym_size, x_size, v_size, f_size
        cdef int float_size = 0
        cdef bint use_double

        frame_offsets = []
        timesteps = []

        self.natoms = 0

        # Open trajectory
        xfp = xdrfile_open(self.filename, "rb")
        if xfp == NULL:
            raise IOError("Could not read file: %s" % self.filename)

        # Loop over all available frames
        while offset < max_offset:
            # Get to the proper offset
            if xdrfile_setpos(xfp, offset) == 0:
                raise IOError("Could not set position in file: %s" % self.filename)

            # Check the magic number
            xdrfile_read_int(&magic, 1, xfp)
            if magic != TRR_MAGIC:
                raise IOError("Corrupted TRR file: '%s'" % self.filename)

            # Check version string
            xdrfile_read_int(&frame_size, 1, xfp)
            if frame_size != strlen(version)+1:
                raise IOError("Corrupted TRR file: '%s'" % self.filename)

            if xdrfile_setpos(xfp, offset + 8 + 16) == 0: # Skip the version string
                raise IOError("Could not set position in file: %s" % self.filename)
            xdrfile_read_int(&ir_size, 1, xfp)
            xdrfile_read_int(&e_size, 1, xfp)
            xdrfile_read_int(&box_size, 1, xfp)
            xdrfile_read_int(&vir_size, 1, xfp)
            xdrfile_read_int(&pres_size, 1, xfp)
            xdrfile_read_int(&top_size, 1, xfp)
            xdrfile_read_int(&sym_size, 1, xfp)
            xdrfile_read_int(&x_size, 1, xfp)
            xdrfile_read_int(&v_size, 1, xfp)
            xdrfile_read_int (&f_size, 1, xfp)
            xdrfile_read_int(&num_atoms_frame, 1, xfp)

            frame_size = TRR_HEADER_SIZE + ir_size + e_size + box_size + vir_size + pres_size + top_size + sym_size + \
                         x_size + v_size + f_size

            if x_size == 0:
                raise ValueError("Frames without coordinates are not supported. Please correct your trajectory")

            # Check trajectory homogeneity
            if offset == 0:
                self.natoms = num_atoms_frame
            else:
                if num_atoms_frame != self.natoms:
                    raise ValueError("Only trajectory with same number of atoms per frame is supported.\n\
                                     Please correct your trajectory file.")

            # Check float size
            if box_size != 0:
                float_size = box_size / 9
            else:
                float_size = x_size / (3 * num_atoms_frame)

            if float_size == 4:
                use_double = False
            elif float_size == 8:
                use_double = True
            else:
                raise ValueError("Float size is not correct (value: %i), trajectory must be corrupted!" % float_size)

            if offset == 0:
                self.use_double = use_double
            else:
                if self.use_double != use_double:
                    raise ValueError("Only trajectory with homogeneous float size is supported.")

            if xdrfile_setpos(xfp, offset + TRR_HEADER_SIZE) == 0: # Skip the step and nre
                raise IOError("Could not set position in file: %s" % self.filename)
            if use_double:
                xdrfile_read_double(&timestep, 1, xfp)
            else:
                xdrfile_read_float(&timestep_float, 1, xfp)
                timestep = timestep_float

            # Store data
            frame_offsets.append(offset)
            timesteps.append(timestep)

            # Update frame counter and offset
            nframes += 1
            offset += frame_size + 2 * float_size # the 2 floats are timestep and lamba_val


        xdrfile_close(xfp)
        frame_offsets = numpy.array(frame_offsets, dtype=np.int64)
        timesteps = numpy.array(timesteps, dtype=np.float64)
        self.timesteps = timesteps
        self.box_offsets = frame_offsets + TRR_HEADER_SIZE + 2 * float_size
        self.coordinate_offsets = frame_offsets + TRR_HEADER_SIZE + 2 * float_size + box_size + vir_size + pres_size
        self.nframes = nframes

    cdef core_base.PBCBox load_box(self, int frame_id):
        cdef XDRFILE *xfp
        cdef fsl_uint pos = self.box_offsets[frame_id]
        cdef real[:,::1] box = numpy.empty((DIM, DIM), dtype=np.float64)
        cdef float float_box[DIM][DIM]
        cdef real double_box[DIM][DIM]
        cdef int i,j

        xfp = xdrfile_open(self.filename, "rb")
        if xfp == NULL:
            raise IOError("Could not read file: %s" % self.filename)
        if xdrfile_setpos(xfp, pos) == 0:
            raise IOError("Could not set position in file: %s" % self.filename)
        if self.use_double:
            if xdrfile_read_double(double_box[0], DIM*DIM, xfp) != DIM*DIM:
                raise IOError("Could not read matrix")
            for i in range(DIM):
                for j in range(DIM):
                    box[i, j] = double_box[i][j]
        else:
            if xdrfile_read_float(float_box[0], DIM*DIM, xfp) != DIM*DIM:
                raise IOError("Could not read matrix")
            # Convert float to real
            for i in range(DIM):
                for j in range(DIM):
                    box[i, j] = <real> float_box[i][j]

        # Close file and release memory
        xdrfile_close(xfp)

        return core_base.PBCBox(box)



    cdef real[:, ::1] load_coords(self, int frame_id, fsl_int[:] atomids) nogil except *:
        cdef XDRFILE *xfp
        cdef fsl_uint pos
        cdef int size = atomids.shape[0]
        cdef int i
        cdef fsl_int atomid
        cdef int natoms, ndecompressed
        cdef real[:, ::1] coords
        cdef rvec real_coords
        cdef fvec float_coords
        cdef int float_size = sizeof(real)
        cdef char *filename

        with gil:
            pos = self.coordinate_offsets[frame_id]
            coords = numpy.empty((size, DIM))
            filename = self.filename

        # Open traj file
        xfp = xdrfile_open(filename, "rb")
        if xfp == NULL:
            with gil:
                raise IOError("Could not read file: %s" % self.filename)

        if self.use_double:
            float_size = sizeof(real)

            # Build the output array
            for i in range(size):
                atomid = atomids[i] - 1 # Reminder: atomids start at 1

                if xdrfile_setpos(xfp, pos + atomid * float_size) == 0:
                    with gil:
                        raise IOError("Could not set position in file: %s" % self.filename)

                # Load coordinates
                if xdrfile_read_double(real_coords, DIM, xfp) < -1:
                    with gil:
                        raise IOError("Could not load coordinates")


                coords[i, XX] = real_coords[XX]
                coords[i, YY] = real_coords[YY]
                coords[i, ZZ] = real_coords[ZZ]

        else:
            float_size = sizeof(float)

            # Build the output array
            for i in range(size):
                atomid = atomids[i] - 1 # Reminder: atomids start at 1

                if xdrfile_setpos(xfp, pos + atomid * float_size * DIM) == 0:
                    with gil:
                        raise IOError("Could not set position in file: %s" % self.filename)

                # Load coordinates
                if xdrfile_read_float(float_coords, DIM, xfp) < -1:
                    with gil:
                        raise IOError("Could not load coordinates")

                coords[i, XX] = float_coords[XX]
                coords[i, YY] = float_coords[YY]
                coords[i, ZZ] = float_coords[ZZ]

        # Close file
        xdrfile_close(xfp)

        return coords

_Coordinates_readers = {".gro": GroReaderCoords,
                        ".xtc": XtcReaderCoords,
                        ".trr": TrrReaderCoords}
