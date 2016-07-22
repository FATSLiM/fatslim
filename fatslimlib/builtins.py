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
from __future__ import print_function, absolute_import, division
import os
import sys
import numpy as np
import time

from . import __version__
from .command import Command, AnalyticalCommand, MembraneAnalysisCommand, MembraneProperty
from .util import backup_file


class CmdVersion(Command):
    _name = "version"
    _help = "Displays FATSLiM version"

    def run(self):
        from . import print_full_version_info
        print_full_version_info()


class CmdHelp(Command):
    _name = "help"
    _help = "Displays help about FATSLiM's commands"

    def fill_parser(self):
        super(CmdHelp, self).fill_parser()
        self.parser.add_argument("command", nargs="?",
                                 help="Command to display information about",
                                 default=None)

    def run(self):
        if self.namespace.command is None:
            print(self.registry.usage)

            print("Available commands:")

            cmd_names = list(self.registry.commands.keys())
            cmd_names.sort()

            fmt_str = "  %%-%is%%s" % (max([len(val) for val in cmd_names]) + 4)

            for name in cmd_names:
                print(fmt_str % (name, self.registry.commands[name].desc))

            print("\nSee '%s help <command>' for more information on a specific command" %
                  self.registry.name)

        else:
            try:
                self.registry.commands[self.namespace.command].help()
            except KeyError:
                print("Unkown command: %s. See '%s --help" % (self.namespace.command,
                                                              self.registry.name))


class CmdTest(Command):
    _name = "self-test"
    _help = "Performs FATSLiM's sanity check up"

    def run(self):  # pragma: no cover
        from . import print_full_version_info
        import os

        print_full_version_info(with_header=False)

        try:
            # noinspection PyUnresolvedReferences
            import pytest
        except ImportError:
            self.print_error("pytest is required to run self-test. Please install it")
            return

        try:
            # noinspection PyUnresolvedReferences
            from . import test
        except ImportError:
            self.print_error("Could not find fatslim's test module. Are you sure it is installed?")
            return

        ret_code = pytest.main(["-rxs", "--tb=short", "%s" % os.path.dirname(test.__file__)])
        if ret_code == 0:
            print("\nOK buddy, I am ready to rock!")
        else:
            print("Ouch... Congratulations! you probably found a bug!")
            print("Please contact the FATSLiM devs with this output so they can work it out!")


class CmdAggregates(AnalyticalCommand):
    _name = "aggregates"
    _help = "Identifies and reports aggregates"

    def __init__(self, registry):
        super(CmdAggregates, self).__init__(registry)
        self.num_aggregates = None

    def fill_parser(self):
        super(CmdAggregates, self).fill_parser()

        self.parser_output_group.add_argument("--output", "-o", default=None, type=str,
                                              help="Plot number of aggregates over the selected "
                                                   "trajectory frames")

        self.parser_output_group.add_argument("--output-index-hg", default=None, type=str,
                                              help="Index group to store aggregates headgroups")

        self.parser_output_group.add_argument("--output-index", default=None, type=str,
                                              help="Index group to store aggregates")

        self.parser_analysis_group.add_argument("--cutoff", default=2.0, type=float,
                                                help="Cutoff distance for aggregate identication")

    def prepare_results(self, num_frames):
        self.num_aggregates = np.zeros(num_frames, dtype=int)
        self.timesteps = np.zeros(num_frames, dtype=int)

    def process_frame(self, frame):
        aggregates = frame.get_aggregates(self.namespace.cutoff)

        self.timesteps[self.processing_index] = frame.timestep

        self.num_aggregates[self.processing_index] = len(aggregates)

        n_aggregates = 0
        if self.namespace.output_index is not None:
            fname, ext = os.path.splitext(self.namespace.output_index)
            fname = "%s_%04i.ndx" % (fname, frame.timestep * .001)

            with open(fname, "w") as fp:
                for agg in aggregates:
                    n_aggregates += 1
                    fp.write("[ aggregate_%i ]\n" % n_aggregates)

                    slice_i = 0
                    atomids = agg.lipid_atomids[slice_i:slice_i + 15]

                    while len(atomids) > 0:
                        fp.write(" ".join(["%i" % val for val in atomids]))
                        fp.write("\n")
                        slice_i += 15
                        atomids = agg.lipid_atomids[slice_i:slice_i + 15]

                    fp.write("\n")

        n_aggregates = 0
        if self.namespace.output_index_hg is not None:
            fname, ext = os.path.splitext(self.namespace.output_index_hg)
            fname = "%s_%04i.ndx" % (fname, frame.timestep * .001)

            with open(fname, "w") as fp:
                for agg in aggregates:
                    n_aggregates += 1
                    fp.write("[ aggregate_%i ]\n" % n_aggregates)

                    slice_i = 0
                    atomids = agg.hg_atomids[slice_i:slice_i + 15]

                    while len(atomids) > 0:
                        fp.write(" ".join(["%i" % val for val in atomids]))
                        fp.write("\n")
                        slice_i += 15
                        atomids = agg.hg_atomids[slice_i:slice_i + 15]

                    fp.write("\n")

    def process_results(self):
        filepath = self.namespace.output
        if filepath is not None:
            fname, ext = os.path.splitext(filepath)

            if ext != ".xvg":
                filepath = "%s.xvg" % fname

            backup_file(filepath)

            with open(filepath, "w") as fp:
                fp.write(self.get_xvg_header())
                fp.write("@ title \"Number of aggregates\"\n")
                fp.write("@ xaxis label \"Time (ps)\"\n")
                fp.write("@ yaxis label \"Number of aggregates\"\n")
                fp.write("@TYPE xy\n")

                for i, val in enumerate(self.timesteps):
                    fp.write("% 8.3f    % 4i\n" % (val, self.num_aggregates[i]))


class CmdMembranes(AnalyticalCommand):
    _name = "membranes"
    _help = "Identifies and report membranes"

    def __init__(self, registry):
        super(CmdMembranes, self).__init__(registry)
        self.num_membranes = None

    def fill_parser(self):
        super(CmdMembranes, self).fill_parser()

        self.parser_output_group.add_argument("-o", "--output", default=None, type=str,
                                              help="Plot number of membranes over the "
                                                   "selected trajectory frames")

        self.parser_output_group.add_argument("--output-index-hg", default=None, type=str,
                                              help="Index group to store leaflet headgroups")

        self.parser_output_group.add_argument("--output-index", default=None, type=str,
                                              help="Index group to store leaflets")

        self.parser_analysis_group.add_argument("--cutoff", default=2.0, type=float,
                                                help="Cutoff distance (in nm) for leaflet "
                                                     "identification")

    def prepare_results(self, num_frames):
        self.num_membranes = np.zeros(num_frames, dtype=int)
        self.timesteps = np.zeros(num_frames, dtype=int)

    def process_frame(self, frame):
        membranes = frame.get_membranes(self.namespace.cutoff)

        self.timesteps[self.processing_index] = frame.timestep

        self.num_membranes[self.processing_index] = len(membranes)

        if self.namespace.output_index is not None:
            fname, ext = os.path.splitext(self.namespace.output_index)
            fname = "%s_%04i.ndx" % (fname, frame.timestep * .001)

            with open(fname, "w") as fp:
                for mid, membrane in enumerate(membranes):
                    for i, leaf in enumerate(membrane):
                        fp.write("[ membrane_%i_leaflet_%i ]\n" % (mid + 1, i + 1))

                        slice_i = 0
                        atomids = leaf.lipid_atomids[slice_i:slice_i + 15]

                        while len(atomids) > 0:
                            fp.write(" ".join(["%i" % val for val in atomids]))
                            fp.write("\n")
                            slice_i += 15
                            atomids = leaf.lipid_atomids[slice_i:slice_i + 15]

                        fp.write("\n")

        if self.namespace.output_index_hg is not None:
            fname, ext = os.path.splitext(self.namespace.output_index_hg)
            fname = "%s_%04i.ndx" % (fname, frame.timestep * .001)

            with open(fname, "w") as fp:
                for mid, membrane in enumerate(membranes):
                    for i, leaf in enumerate(membrane):
                        fp.write("[ membrane_%i_leaflet_%i ]\n" % (mid + 1, i + 1))

                        slice_i = 0
                        atomids = leaf.hg_atomids[slice_i:slice_i + 15]

                        while len(atomids) > 0:
                            fp.write(" ".join(["%i" % val for val in atomids]))
                            fp.write("\n")
                            slice_i += 15
                            atomids = leaf.hg_atomids[slice_i:slice_i + 15]

                        fp.write("\n")

    def process_results(self):
        filepath = self.namespace.output
        if filepath is not None:
            fname, ext = os.path.splitext(filepath)

            filepath = "%s.xvg" % fname

            backup_file(filepath)

            with open(filepath, "w") as fp:
                fp.write(self.get_xvg_header())
                fp.write("@ title \"Number of membranes\"\n")
                fp.write("@ xaxis label \"Time (ps)\"\n")
                fp.write("@ yaxis label \"Number of membranes\"\n")
                fp.write("@TYPE xy\n")

                for i, val in enumerate(self.timesteps):
                    fp.write("% 8.3f    % 4i\n" % (val, self.num_membranes[i]))


class PropertyThickness(MembraneProperty):
    _fullname = "thickness"
    _shortname = "thickness"
    _units = "nm"

    def fill_parsers(self, parser_analysis, parser_input, parser_output):
        super(PropertyThickness, self).fill_parsers(parser_analysis, parser_input, parser_output)

        parser_analysis.add_argument("--thickness-cutoff",
                                     dest="thickness_cutoff", default=6.0, type=float,
                                     help="Cutoff distance (in nm) used to identify inter-leaflet "
                                          "neighbors")

    def get_values(self, membrane):
        mem_val, l1_val, l2_val = membrane.get_thickness(interleaflet_cutoff=
                                                         self.parent.namespace.thickness_cutoff,
                                                         only_average=True,
                                                         force=False)

        self.membrane_avg_values[self.parent.processing_index] = mem_val
        self.leaflet1_avg_values[self.parent.processing_index] = l1_val[0]
        self.leaflet2_avg_values[self.parent.processing_index] = l2_val[0]

        return mem_val, l1_val[0], l2_val[0]

    def get_raw_values(self, membrane):
        mem_val, l1_val, l2_val = membrane.get_thickness(interleaflet_cutoff=
                                                         self.parent.namespace.thickness_cutoff,
                                                         only_average=False,
                                                         force=False)
        return l1_val, l2_val


class CmdThickness(MembraneAnalysisCommand):
    _name = "thickness"
    _help = "Retrieves bilayer thickness"
    _property_classes = [PropertyThickness]


class PropertyAPL(MembraneProperty):
    _fullname = "area per lipid"
    _shortname = "apl"
    _units = "nm^2"
    _groupable_by_type = True

    def fill_parsers(self, parser_analysis, parser_input, parser_output):
        super(PropertyAPL, self).fill_parsers(parser_analysis, parser_input, parser_output)

        parser_analysis.add_argument("--apl-cutoff",
                                     dest="apl_cutoff", default=3.0, type=float,
                                     help="Cutoff distance (in nm) used to approximate "
                                          "planar region")

        parser_analysis.add_argument("--apl-limit",
                                     dest="apl_limit", default=10.0, type=float,
                                     help="Upper limit (in nm2) considered when calculating "
                                          "individual APL values")

    def get_values(self, membrane):
        mem_val, l1_val_dict, l2_val_dict = membrane.get_apl(
            cutoff=self.parent.namespace.apl_cutoff,
            area_limit=self.parent.namespace.apl_limit,
            by_type=self.needs_grouping,
            only_average=True,
            force=False)

        if self.needs_grouping:
            mem_val_by_type = {}
            mem_val = [0.0, 0]
            l1_val = [0.0, 0]
            l2_val = [0.0, 0]

            for key, val in l1_val_dict.items():
                mem_val[0] += val[0] * val[1]
                mem_val[1] += val[0]

                l1_val[0] += val[0] * val[1]
                l1_val[1] += val[0]

                try:
                    self.leaflet1_avg_values_by_type[key][self.parent.processing_index] = val[1]
                except KeyError:
                    self.leaflet1_avg_values_by_type[key] = np.zeros(len(self.membrane_avg_values))
                    self.leaflet1_avg_values_by_type[key][self.parent.processing_index] = val[1]

                mem_val_by_type[key] = [val[1] * val[0], val[0]]

            for key, val in l2_val_dict.items():
                mem_val[0] += val[0] * val[1]
                mem_val[1] += val[0]

                l2_val[0] += val[0] * val[1]
                l2_val[1] += val[0]

                try:
                    self.leaflet2_avg_values_by_type[key][self.parent.processing_index] = val[1]
                except KeyError:
                    self.leaflet2_avg_values_by_type[key] = np.zeros(len(self.membrane_avg_values))
                    self.leaflet2_avg_values_by_type[key][self.parent.processing_index] = val[1]

                try:
                    mem_val_by_type[key][0] += val[1] * val[0]
                    mem_val_by_type[key][1] += val[0]
                except KeyError:  # pragma: no cover
                    mem_val_by_type[key] = [val[1] * val[0], val[0]]

            for key, val in mem_val_by_type.items():
                val = float(val[0]) / val[1]

                try:
                    self.membrane_avg_values_by_type[key][self.parent.processing_index] = val
                except KeyError:
                    self.membrane_avg_values_by_type[key] = np.zeros(len(self.membrane_avg_values))
                    self.membrane_avg_values_by_type[key][self.parent.processing_index] = val

            mem_val = float(mem_val[0]) / mem_val[1]
            l1_val = float(l1_val[0]) / l1_val[1]
            l2_val = float(l2_val[0]) / l2_val[1]

            self.membrane_avg_values[self.parent.processing_index] = mem_val
            self.leaflet1_avg_values[self.parent.processing_index] = l1_val
            self.leaflet2_avg_values[self.parent.processing_index] = l2_val

        else:
            l1_val = l1_val_dict[0]
            l2_val = l2_val_dict[0]

            self.membrane_avg_values[self.parent.processing_index] = mem_val
            self.leaflet1_avg_values[self.parent.processing_index] = l1_val
            self.leaflet2_avg_values[self.parent.processing_index] = l2_val

        return mem_val, l1_val, l2_val

    def get_raw_values(self, membrane):
        mem_val, l1_val, l2_val = membrane.get_apl(cutoff=self.parent.namespace.apl_cutoff,
                                                   area_limit=self.parent.namespace.apl_limit,
                                                   by_type=False,
                                                   only_average=False,
                                                   force=False)

        return l1_val, l2_val


class PropertyArea(MembraneProperty):
    _fullname = "area"
    _shortname = "area"
    _units = "nm^2"
    _groupable_by_type = False
    _raw_values_supported = False

    def get_values(self, membrane):
        mem_val, l1_val, l2_val = membrane.get_apl(
            cutoff=self.parent.namespace.apl_cutoff,
            area_limit=self.parent.namespace.apl_limit,
            by_type=self.needs_grouping,
            only_average=True,
            force=False)

        l1_val = l1_val[3]
        l2_val = l2_val[3]

        self.membrane_avg_values[self.parent.processing_index] = 0.5 * (l1_val + l2_val)
        self.leaflet1_avg_values[self.parent.processing_index] = l1_val
        self.leaflet2_avg_values[self.parent.processing_index] = l2_val

        return mem_val, l1_val, l2_val

    def get_raw_values(self, membrane):  # pragma: no cover
        pass


class CmdApl(MembraneAnalysisCommand):
    _name = "apl"
    _help = "Retrieves area per lipid"
    _property_classes = [PropertyAPL, PropertyArea]


class CmdBenchmark(AnalyticalCommand):
    _name = "benchmark"
    _help = "Benchmark internal routine (not intented to be a end user command!)"

    def run(self):  # pragma: no cover
        import time
        from .datareading import load_trajectory
        from .core_base import pretty_delta, pretty_duration
        import gc

        def get_memory_usage():
            memory = -1
            with open('/proc/self/status', "rb") as fp:
                for line in fp:
                    if line.startswith("VmRSS:"):
                        memory = int(line.split()[1])
                        break
            return memory

        def h_memory(memory):
            units = ["kB", "MB", "GB", "TB"]
            unitid = 0
            while memory > 1024:
                memory /= 1024.
                unitid += 1
            return "%.1f %s" % (memory, units[unitid])

        mem_before = get_memory_usage()
        time_before = time.time()
        print("Memory usage BEFORE starting: %s" % h_memory(mem_before))

        print("\nInitializing Trajectory...")
        mem_begin = get_memory_usage()
        time_begin = time.time()
        traj = load_trajectory(self.namespace.conf,
                               self.namespace.index,
                               self.namespace.trajectory,
                               self.verbose)
        time_end = time.time()
        mem_end = get_memory_usage()
        print("=> Memory usage after initialization: %s (increment: %s) - Time taken: %s" %
              (h_memory(mem_end),
               h_memory(mem_end - mem_begin),
               pretty_delta(time_begin, time_end)))

        print("\nPreloading Trajectory...")
        mem_begin = get_memory_usage()
        time_begin = time.time()
        traj.initialize(hg_group=self.namespace.hg_group)
        time_end = time.time()
        mem_end = get_memory_usage()
        print("=> Memory usage after preloading: %s (increment: %s) - Time taken: %s" %
              (h_memory(mem_end),
               h_memory(mem_end - mem_begin),
               pretty_delta(time_begin, time_end)))

        size = len(traj)
        first_frame = max(0, self.namespace.begin)
        last_frame = min(size, self.namespace.end)
        if last_frame < 0:
            last_frame = size
        num_frames = last_frame - first_frame

        if num_frames <= 0:
            return

        getting_times = np.zeros(num_frames)
        getting_mems = np.zeros(num_frames)
        coords_times = np.zeros(num_frames)
        coords_mems = np.zeros(num_frames)
        coords_bbox_times = np.zeros(num_frames)
        coords_bbox_mems = np.zeros(num_frames)
        directions_times = np.zeros(num_frames)
        directions_mems = np.zeros(num_frames)
        normals_times = np.zeros(num_frames)
        normals_mems = np.zeros(num_frames)
        membrane_times = np.zeros(num_frames)
        membrane_mems = np.zeros(num_frames)
        apl_times = np.zeros(num_frames)
        apl_mems = np.zeros(num_frames)
        thickness_times = np.zeros(num_frames)
        thickness_mems = np.zeros(num_frames)
        gc_times = np.zeros(num_frames)
        gc_mems = np.zeros(num_frames)
        total_frame_times = np.zeros(num_frames)

        full_details = False
        for i in range(first_frame, last_frame):
            if full_details:
                print("\nProcessing frame #%i (current memory usage: %s):" %
                      (
                          i - first_frame + 1,
                          h_memory(get_memory_usage())
                      ))
            else:
                print("Processing frame #%i (current memory usage: %s)" %
                      (
                          i - first_frame + 1,
                          h_memory(get_memory_usage())
                      ), end="")
            mem_begin = get_memory_usage()
            time_begin = time.time()
            frame = traj[i]
            time_end = time.time()
            mem_end = get_memory_usage()
            frame_begin = time_begin
            if full_details:
                print("  => Memory usage after getting frame: %s (increment: %s) - "
                      "Time taken: %s" %
                      (h_memory(mem_end),
                       h_memory(mem_end - mem_begin),
                       pretty_delta(time_begin, time_end)))
            getting_times[i - first_frame] = time_end - time_begin
            getting_mems[i - first_frame] = mem_end - mem_begin
            total_frame_times[i - first_frame] += time_end - time_begin

            # Getting coords
            mem_begin = get_memory_usage()
            time_begin = time.time()
            frame.bead_coords
            time_end = time.time()
            mem_end = get_memory_usage()
            if full_details:
                print("  => Memory usage after getting coords: %s (increment: %s) - "
                      "Time taken: %s" %
                      (h_memory(mem_end),
                       h_memory(mem_end - mem_begin),
                       pretty_delta(time_begin, time_end),))
            coords_times[i - first_frame] = time_end - time_begin
            coords_mems[i - first_frame] = mem_end - mem_begin
            total_frame_times[i - first_frame] += time_end - time_begin

            # Getting coords
            mem_begin = get_memory_usage()
            time_begin = time.time()
            frame.bead_coords_bbox
            time_end = time.time()
            mem_end = get_memory_usage()
            if full_details:
                print("  => Memory usage after getting coords in bbox : %s (increment: %s) - "
                      "Time taken: %s" %
                      (h_memory(mem_end),
                       h_memory(mem_end - mem_begin),
                       pretty_delta(time_begin, time_end),))
            coords_bbox_times[i - first_frame] = time_end - time_begin
            coords_bbox_mems[i - first_frame] = mem_end - mem_begin
            total_frame_times[i - first_frame] += time_end - time_begin

            # Getting directions
            mem_begin = get_memory_usage()
            time_begin = time.time()
            frame.directions
            time_end = time.time()
            mem_end = get_memory_usage()
            if full_details:
                print("  => Memory usage after getting directions: %s (increment: %s) - "
                      "Time taken: %s" %
                      (h_memory(mem_end),
                       h_memory(mem_end - mem_begin),
                       pretty_delta(time_begin, time_end),
                       ))
            directions_times[i - first_frame] = time_end - time_begin
            directions_mems[i - first_frame] = mem_end - mem_begin
            total_frame_times[i - first_frame] += time_end - time_begin

            # Getting normals
            mem_begin = get_memory_usage()
            time_begin = time.time()
            frame.normals
            time_end = time.time()
            mem_end = get_memory_usage()
            if full_details:
                print("  => Memory usage after getting normals: %s (increment: %s) - "
                      "Time taken: %s" %
                      (h_memory(mem_end),
                       h_memory(mem_end - mem_begin),
                       pretty_delta(time_begin, time_end)
                       ))
            normals_times[i - first_frame] = time_end - time_begin
            normals_mems[i - first_frame] = mem_end - mem_begin
            total_frame_times[i - first_frame] += time_end - time_begin

            # Getting membranes
            mem_begin = get_memory_usage()
            time_begin = time.time()
            membrane = frame.get_membranes()[0]
            time_end = time.time()
            mem_end = get_memory_usage()
            if full_details:
                print("  => Memory usage after getting membranes: %s (increment: %s) - "
                      "Time taken: %s" %
                      (h_memory(mem_end),
                       h_memory(mem_end - mem_begin),
                       pretty_delta(time_begin, time_end)))
            membrane_times[i - first_frame] = time_end - time_begin
            membrane_mems[i - first_frame] = mem_end - mem_begin
            total_frame_times[i - first_frame] += time_end - time_begin

            # Getting APL
            mem_begin = get_memory_usage()
            time_begin = time.time()
            membrane.get_apl()
            time_end = time.time()
            mem_end = get_memory_usage()
            if full_details:
                print("  => Memory usage after getting APL: %s (increment: %s) - "
                      "Time taken: %s" %
                      (h_memory(mem_end),
                       h_memory(mem_end - mem_begin),
                       pretty_delta(time_begin, time_end)))
            apl_times[i - first_frame] = time_end - time_begin
            apl_mems[i - first_frame] = mem_end - mem_begin
            total_frame_times[i - first_frame] += time_end - time_begin

            # Getting Thickness
            mem_begin = get_memory_usage()
            time_begin = time.time()
            membrane.get_thickness()
            time_end = time.time()
            mem_end = get_memory_usage()
            if full_details:
                print("  => Memory usage after getting thickness: %s (increment: %s) - "
                      "Time taken: %s" %
                      (h_memory(mem_end),
                       h_memory(mem_end - mem_begin),
                       pretty_delta(time_begin, time_end)))
            thickness_times[i - first_frame] = time_end - time_begin
            thickness_mems[i - first_frame] = mem_end - mem_begin
            total_frame_times[i - first_frame] += time_end - time_begin

            # Garbage collection
            mem_begin = get_memory_usage()
            time_begin = time.time()
            try:
                del frame
                del membrane
            except UnboundLocalError:
                pass
            gc.collect()
            time_end = time.time()
            mem_end = get_memory_usage()
            frame_end = time_end
            if full_details:
                print("  => Memory usage after garbage collection: %s (increment: %s) - "
                      "Time taken: %s" %
                      (h_memory(mem_end),
                       h_memory(mem_end - mem_begin),
                       pretty_delta(time_begin, time_end)))
            else:
                print(" Done in %s" % pretty_delta(frame_begin, frame_end))
            gc_times[i - first_frame] = time_end - time_begin
            gc_mems[i - first_frame] = mem_end - mem_begin
            total_frame_times[i - first_frame] += time_end - time_begin

        mem_after = get_memory_usage()
        time_after = time.time()
        print("\nGrand total summary: Memory usage: %s (increment: %s) - Time taken: %s" %
              (h_memory(mem_after),
               h_memory(mem_after - mem_before),
               pretty_delta(time_before, time_after))
              )
        print("Average values over %i frames:" % num_frames)
        print("Average processing time per frame: %s-+%s" %
              (pretty_duration(total_frame_times.mean()),
               pretty_duration(total_frame_times.std())))
        print("Getting frames: Memory usage: %s - Time taken: %s (%.1f%% frame processing)" %
              (
                  h_memory(getting_mems.mean()),
                  pretty_duration(getting_times.mean()),
                  (getting_times / total_frame_times * 100.0).mean()
              ))
        print("Getting coords: Memory usage: %s - Time taken: %s (%.1f%% frame processing)" %
              (
                  h_memory(coords_mems.mean()),
                  pretty_duration(coords_times.mean()),
                  (coords_times / total_frame_times * 100.0).mean()
              ))
        print("Getting coords_bbox: Memory usage: %s - Time taken: %s (%.1f%% frame processing)" %
              (
                  h_memory(coords_bbox_mems.mean()),
                  pretty_duration(coords_bbox_times.mean()),
                  (coords_bbox_times / total_frame_times * 100.0).mean()
              ))
        print("Total I/O: %.1f%% of frame processing" %
              ((getting_times + coords_times + coords_bbox_times) /
               total_frame_times * 100.0).mean())
        print("Getting directions: Memory usage: %s - Time taken: %s (%.1f%% frame processing)" %
              (
                  h_memory(directions_mems.mean()),
                  pretty_duration(directions_times.mean()),
                  (directions_times / total_frame_times * 100.0).mean()
              ))
        print("Getting normals: Memory usage: %s - Time taken: %s (%.1f%% frame processing)" %
              (
                  h_memory(normals_mems.mean()),
                  pretty_duration(normals_times.mean()),
                  (normals_times / total_frame_times * 100.0).mean()
              ))
        print("Getting membrane: Memory usage: %s - Time taken: %s (%.1f%% frame processing)" %
              (
                  h_memory(membrane_mems.mean()),
                  pretty_duration(membrane_times.mean()),
                  (membrane_times / total_frame_times * 100.0).mean()
              ))
        print("Getting apl: Memory usage: %s - Time taken: %s (%.1f%% frame processing)" %
              (
                  h_memory(apl_mems.mean()),
                  pretty_duration(apl_times.mean()),
                  (apl_times / total_frame_times * 100.0).mean()
              ))
        print("Getting thickness: Memory usage: %s - Time taken: %s (%.1f%% frame processing)" %
              (
                  h_memory(thickness_mems.mean()),
                  pretty_duration(thickness_times.mean()),
                  (thickness_times / total_frame_times * 100.0).mean()
              ))
        print("Garbage collection: Memory usage: %s - Time taken: %s (%.1f%% frame processing)" %
              (
                  h_memory(gc_mems.mean()),
                  pretty_duration(gc_times.mean()),
                  (gc_times / total_frame_times * 100.0).mean()
              ))
