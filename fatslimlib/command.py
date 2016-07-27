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

# Global imports
from __future__ import print_function
import argparse
import csv
import numpy as np
import os
import sys
import time
import gc

# Local imports
from .datareading import load_trajectory
from .core_base import pretty_delta, pretty_duration
from .util import backup_file

COMMAND_STATE_NOTSTARTED = 0
COMMAND_STATE_RUNNING = 1
COMMAND_STATE_SUCCESS = 2
COMMAND_STATE_INTERRUPTED = 3
COMMAND_STATE_FAILURE = 4
COMMAND_STATE_ABORTED = 5


class AlreadyRegisteredCommandError(Exception):
    pass


class UnknownArgumentError(Exception):
    pass


class CommandRegistry(object):
    def __init__(self, name, debug_mode=False, verbose=True):
        self._name = name
        self._commands = {}

        self.usage = "usage: %s [--help|-h] [--version] command [<args>]" \
                     % self.name
        self.debug_mode = debug_mode
        self.verbose = verbose
        self.warnings = []
        self.errors = []

    @property
    def name(self):
        return self._name

    @property
    def commands(self):
        return self._commands

    def print_debug(self, msg, *args, **kwargs):
        if self.debug_mode:
            msg = "DEBUG: %s" % msg
            kwargs["file"] = sys.stderr
            print(msg, *args, **kwargs)

    def print_warning(self, msg, *args, **kwargs):
        self.warnings.append(msg)
        msg = "WARNING: %s" % msg
        kwargs["file"] = sys.stderr
        print(msg, *args, **kwargs)

    def print_error(self, msg, *args, **kwargs):
        self.errors.append(msg)
        msg = "ERROR: %s" % msg
        kwargs["file"] = sys.stderr
        print(msg, *args, **kwargs)

    def add_command(self, cmd):

        """Add command to registry

        :param cmd: A subclass of `Command`
        :return: None
        """
        if not issubclass(cmd, Command):  # pragma: no cover
            raise TypeError("Only subclasses of `Command` can be added to CommandRegistry")

        # Do nothing if cmd is Command itself
        if cmd == Command:  # pragma: no cover
            return

        # Instantiate cmd
        cmd = cmd(self)

        # Verify the command (or one of its aliases) is not already registered
        for identifier in cmd.identifiers:
            if identifier in self._commands:
                raise AlreadyRegisteredCommandError(identifier)

        # Register the command
        for identifier in cmd.identifiers:
            self._commands[identifier] = cmd

    def run(self, args=None):
        exit_code = 0

        # First strip arguments to get global options, command and command arguments
        if args is None:  # pragma: no cover
            args = sys.argv[1:]

        global_options = []
        command = ''
        command_args = []

        try:
            while args[0].startswith("-"):
                val = args.pop(0)
                global_options.append(val)
            command = args.pop(0)
            command_args = args
        except IndexError:
            pass

        # Handle global --version
        try:
            global_options.remove("--version")
        except ValueError:
            pass
        else:
            if ("--help" not in global_options) and ("-h" not in global_options):
                command = "version"

        # Handle global --help
        try:
            global_options.remove("--help")
        except ValueError:
            pass
        else:
            if command == "":
                command_args = []
                command = "help"
            else:  # pragma: no cover
                command_args += ["--help"]

        # Handle global -h
        try:
            global_options.remove("-h")
        except ValueError:
            pass
        else:  # pragma: no cover
            if command == "":
                command_args = []
                command = "help"
            else:
                command_args += ["-h"]

        # Handle --debug
        if "--debug" in command_args:  # pragma: no cover
            self.debug_mode = True

        # Stop here if there are unknown global options
        if len(global_options) != 0:  # pragma: no cover
            self.print_error("Unknown global option: %s" % global_options[0])
            print(self.usage)
            return 1

        # Try to run the command
        try:
            # If command is not set yet, help is displayed
            if command == "":
                command = "help"
            cmd_obj = self._commands[command]
        except KeyError:
            self.print_error("Unknown command: %s. See '%s --help'" % (command, self.name))
            return 1
        else:

            if isinstance(cmd_obj, AnalyticalCommand) and not self.debug_mode:
                print("Running command: '%s'... This may take some time, be patient!" % command)

            if self.debug_mode:  # pragma: no cover
                import cProfile
                import pstats
                import io

                self.print_debug("Running command: '%s' with arguments: %s" % (command,
                                                                               command_args))
                self.print_debug("Watch out! Profiling is on, this may slow things down.")

                pr = cProfile.Profile()
                pr.enable()
            else:
                pr = None
                pstats = None

            cmd_obj.prerun()
            begin_time = time.time()

            cmd_obj.run_with_args(command_args)

            finish_time = time.time()

            if self.debug_mode:  # pragma: no cover
                pr.disable()

            duration = finish_time - begin_time
            duration_units = "s"

            if duration < 0.001:
                duration *= 1e6
                duration_units = "µs"
            elif duration < 1:
                duration *= 1e3
                duration_units = "ms"

            print("'%s' command executed in %.3f %s (CPU)" % (command, duration, duration_units))

            if cmd_obj.state == COMMAND_STATE_FAILURE:  # pragma: no cover
                print("")
                if cmd_obj.exc is not None:
                    import traceback
                    exc_type, exc_value, exc_tb = cmd_obj.exc
                    exc_string = traceback.format_exception_only(exc_type, exc_value)[-1]
                    exc_file, exc_lino = traceback.extract_tb(exc_tb)[-1][:2]
                    message = "%s If you think this may be a bug, please send the following " \
                              "information to the developpers:\n" % exc_string
                    message += "-" * 80
                    message += "\n  %s was raised in %s at line %i\n" % (exc_type.__name__,
                                                                         exc_file, exc_lino)
                    message += "  Complete traceback:\n    ->%s" % "    ->"\
                        .join(traceback.format_tb(exc_tb))
                    message += "-" * 80
                    self.print_error(message)
                else:
                    self.print_error("Unknown failure!")
            if cmd_obj.state in [COMMAND_STATE_ABORTED, COMMAND_STATE_FAILURE]:
                exit_code = 1

            if self.debug_mode:  # pragma: no cover
                ps = pstats.Stats(pr, stream=sys.stderr)
                ps.sort_stats('cumulative', 'time', 'nfl')
                self.print_debug("Profile stats:")
                ps.print_stats(20)

        return exit_code


class Command(object):
    """Abstract class for FATSLiM's commands"""
    _name = "Abstract command"
    _help = "Please override me!"
    _aliases = []

    def __init__(self, registry):
        if not isinstance(registry, CommandRegistry):
            raise TypeError("registry must be a CommandRegistry instance")
        self._registry = registry

        self.verbose = registry.verbose
        self.debug_mode = registry.debug_mode

        self.parser = None
        self.raw_args = ""
        self.namespace = None

        self.state = COMMAND_STATE_NOTSTARTED
        self.exc = None

        self.fill_parser()

    @property
    def registry(self):
        return self._registry

    @property
    def name(self):
        """Command name"""
        return self._name

    @property
    def fullname(self):
        """Full name"""
        return "%s-%s" % (self.registry.name, self.name)

    @property
    def aliases(self):
        """Aliases for command"""
        return self._aliases

    @property
    def identifiers(self):
        identifiers = [self._name]
        identifiers.extend(self._aliases)
        return tuple(identifiers)

    @property
    def desc(self):
        return self._help

    def help(self):
        print("command: %s %s" % (self.registry.name, self.name))
        print("purpose: %s" % self.desc)
        self.parser.print_help()

    def print_debug(self, msg, *args, **kwargs):
        if self.debug_mode:
            msg = "DEBUG: %s" % msg
            kwargs["file"] = sys.stderr
            print(msg, *args, **kwargs)

    def print_warning(self, msg, *args, **kwargs):
        self.registry.print_warning(msg, *args, **kwargs)

    def print_error(self, msg, *args, **kwargs):
        self.registry.print_error(msg, *args, **kwargs)

    def print_verbose(self, *args, **kwargs):
        if self.verbose:
            print(*args, **kwargs)
            sys.stdout.flush()

    def fill_parser(self):
        # Create and populate parser
        parser = argparse.ArgumentParser(prog="%s %s" % (self.registry.name, self._name),
                                         add_help=False,
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        group = parser.add_argument_group("General options")

        group.add_argument("--help", "-h", action="store_true", default=False,
                           help="See help about %s command" % self.name)
        group.add_argument("--debug", action="store_true", default=False,
                           help="Enable debug output")
        group.add_argument("--verbose", "-v", action="store_true", default=False,
                           dest="verbose",
                           help="Be loud and noisy")
        self.parser = parser

    def process_results(self):
        pass

    def _process_results(self):
        self.process_results()

    def _run(self):
        output = None
        try:
            self.state = COMMAND_STATE_RUNNING
            output = self.run()
        except KeyboardInterrupt:  # pragma: no cover
            self.state = COMMAND_STATE_INTERRUPTED
            print("\nInterruption requested!")
        except (ValueError, TypeError, IOError):
            self.state = COMMAND_STATE_FAILURE
            self.exc = sys.exc_info()
            raise
        else:
            self.state = COMMAND_STATE_SUCCESS
        return output

    def prerun(self):
        pass

    def run(self):  # pragma: no cover
        pass

    def run_with_args(self, args=()):
        self.raw_args = " ".join(args)
        namespace, unknown_args = self.parser.parse_known_args(args)

        if len(unknown_args) != 0:
            if len(unknown_args) > 1:  # pragma: no cover
                plural = "s"
            else:
                plural = ""
            self.print_error("'%s' command does not recognize the following argument%s: %s\n" %
                             (self.name, plural, ",".join(unknown_args)))
            self.parser.print_help()
            self.state = COMMAND_STATE_ABORTED
            return

        self.verbose = namespace.verbose
        self.debug_mode = namespace.debug

        del namespace.verbose
        del namespace.debug

        if namespace.help:
            self.help()
        else:
            del namespace.help

            self.namespace = namespace
            self._run()
            self._process_results()


class AnalyticalCommand(Command):
    def __init__(self, registry):
        self.parser_input_group = None
        self.parser_output_group = None
        self.parser_analysis_group = None
        self.timesteps = None
        self.processing_index = -1
        super(AnalyticalCommand, self).__init__(registry)

    def fill_parser(self):
        super(AnalyticalCommand, self).fill_parser()

        self.parser_input_group = self.parser.add_argument_group("Options to specify input files")

        self.parser_input_group.add_argument("--conf", "-c", default="conf.gro",
                                             help="Configuration file")
        self.parser_input_group.add_argument("--trajectory", "-t", default="traj.trr",
                                             help="Trajectory file")
        self.parser_input_group.add_argument("--index", "-n", default="index.ndx",
                                             help="Index file")

        self.parser_output_group = self.parser.add_argument_group("Options to specify output"
                                                                  " files")

        self.parser_analysis_group = self.parser.add_argument_group("Options related to "
                                                                    "analysis parameters")

        self.parser_analysis_group.add_argument("--hg-group", default="headgroups",
                                                help="Index group name used to define lipid "
                                                     "headgroups")

        self.parser_analysis_group.add_argument("--interacting-group", default="protein",
                                                help="Index group name used to define interacting "
                                                     "atoms (e.g. protein).")

        self.parser_analysis_group.add_argument("--nthreads", default=-1, type=int,
                                                help="Number of threads to use")

        self.parser_analysis_group.add_argument("--begin", "-b", default=-1, type=int,
                                                help="First timestep (ps) to use for analysis")

        self.parser_analysis_group.add_argument("--begin-frame", default=-1, type=int,
                                                help="First frame (index) to use for analysis")

        self.parser_analysis_group.add_argument("--end", "-e", default=-1, type=float,
                                                help="Last timestep (ps) to use for analysis")

        self.parser_analysis_group.add_argument("--end-frame", default=-1, type=float,
                                                help="Last frame (index) to use for analysis")

    def get_xvg_header(self):
        import time
        from . import __shortname__, __version__, __cmdlinename__
        header = "# This file was created on %s\n" % time.asctime()

        header += "# Created by: %s v.%s\n" % (__shortname__, __version__)

        header += "# Command line: %s %s %s\n" % (__cmdlinename__, self.name, self.raw_args)

        return header

    def _run(self):
        from .core_base import set_num_threads, get_num_threads
        nthreads = self.namespace.nthreads
        set_num_threads(nthreads)

        print("Analysis will be performed using %i threads." % (get_num_threads()))
        begin = time.time()
        nframes = super(AnalyticalCommand, self)._run()
        end = time.time()

        output = "Analysis done in %s" % pretty_delta(begin, end)
        if nframes is not None:
            if nframes > 0:
                output += " (%s per frame)" % pretty_duration((end - begin) / nframes)
        if output is not None and len(output) > 0:
            self.print_verbose(output)

    def initialize_trajetory(self, traj):
        traj.initialize(hg_group=self.namespace.hg_group,
                        interacting_group=self.namespace.interacting_group)

    def prepare_results(self, num_frames):
        """Initializes the numpy.ndarray used to store the results.

        :param num_frames: Number of frames selected for analysis
        :type num_frames: int
        :returns: None
        """
        self.timesteps = np.zeros(num_frames)

    def process_frame(self, frame):
        """Performs the analysis a one frame.

        :param frame: The frame to analyze
        :type frame: fatslimlib.core_base.pyx.Frame
        :returns: the output from the analysis
        :rtype: str"""
        self.timesteps[self.processing_index] = frame.timestep

        return ""

    def _process_results(self):
        begin = time.time()
        self.print_verbose("Processing results...", end="")
        output = self.process_results()

        end = time.time()
        self.print_verbose(" done in %s" % pretty_delta(begin, end))
        if output is not None and len(output) > 0:
            print("Results:")
            print(output)

    def run(self):
        traj = load_trajectory(self.namespace.conf,
                               self.namespace.index,
                               self.namespace.trajectory,
                               self.verbose)
        self.initialize_trajetory(traj)
        size = len(traj)

        first_frame = max(0, self.namespace.begin_frame)
        last_frame = min(size, self.namespace.end_frame)
        if last_frame < 0:
            last_frame = size

        if self.namespace.begin >= 0:
            first_frame = size
            for i, timestep in enumerate(traj.timesteps):
                if timestep >= self.namespace.begin:
                    first_frame = i
                    break

        if self.namespace.end >= 0:
            for i, timestep in enumerate(traj.timesteps):
                if timestep > self.namespace.end:
                    last_frame = i
                    break

        num_frames = last_frame - first_frame

        if num_frames <= 0:
            self.print_warning("No frames are selected for analysis!")
            return
        self.prepare_results(num_frames)

        self.print_verbose("Beginning the analysis of %i frame%s:" % (num_frames,
                                                                      "s"[num_frames == 1:]))

        grand_begin = -1
        last_line_length = 0
        for i in range(first_frame, last_frame):
            frame = traj[i]
            self.processing_index = (i - first_frame)

            output = "Analysing frame % 5i/% 5i (time: %i ps)..." % (self.processing_index + 1,
                                                                     num_frames,
                                                                     frame.timestep)
            line_length = len(output)
            print("\r%s" % output, end="")
            sys.stdout.flush()

            # Get current time
            begin = time.time()
            if grand_begin < 0:
                grand_begin = begin

            # Actually process frame
            frame_output = self.process_frame(frame)

            # Force garbage collection
            del frame
            gc.collect()

            # Get current time to evaluate how much time was needed to process the frame
            end = time.time()

            # Print the duration
            duration = end - grand_begin
            remaining = duration * num_frames / (self.processing_index + 1) - duration

            output = " done in %s (Remaining: %.0f s)" % (pretty_delta(begin, end),
                                                          remaining)

            line_length += len(output)

            print("%s%s" % (output, " " * max(0, last_line_length - line_length)),
                  end="")

            self.print_verbose("\n%s" % frame_output)

            last_line_length = line_length

        if not self.verbose:
            print("")

        return num_frames


class MembraneProperty(object):
    _fullname = "Unknown property"
    _shortname = "prop"
    _units = "A.U."
    _groupable_by_type = False
    _raw_values_supported = True

    def __init__(self, parent):
        """Membrane property class
        :param parent: The parent command of the property
        :type parent: fatslimlib.command.MembraneAnalysisCommand
        """
        if not isinstance(parent, MembraneAnalysisCommand):
            raise TypeError("Parent be a MembraneAnalysisCommand instance")
        self.parent = parent

        self.membrane_avg_values = None
        self.membrane_avg_values_by_type = None
        self.leaflet1_avg_values = None
        self.leaflet1_avg_values_by_type = None
        self.leaflet2_avg_values = None
        self.leaflet2_avg_values_by_type = None

        self.needs_grouping = False

    @property
    def fullname(self):
        return self._fullname

    @property
    def shortname(self):
        return self._shortname

    @property
    def units(self):
        return self._units

    @property
    def raw_values_supported(self):
        return self._raw_values_supported

    def fill_parsers(self, parser_analysis, parser_input, parser_output):
        """Fill the parser with options specific to the property.
        :param parser_analysis: parser used to store options related to analysis parameters
        :type parser_analysis: argparse.ArgumentParser
        :param parser_input: parser used to store options related to input files
        :type parser_input: argparse.ArgumentParser
        :param parser_output: parser used to store options related to output files
        :type parser_output: argparse.ArgumentParser
        """

        parser_output.add_argument("--plot-%s" % self._shortname,
                                   dest="plot_%s" % self._shortname,
                                   default=None, type=str,
                                   help="Plot %s over trajectory (.xvg file)" % self._fullname.lower())

        if self._raw_values_supported:
            parser_output.add_argument("--export-%s-raw" % self._shortname,
                                       dest="export_%s_raw" % self._shortname,
                                       default=None, type=str,
                                       help="Save %s raw values (one .csv file per frame)"
                                            % self._fullname.lower())

        if self._groupable_by_type:
            parser_analysis.add_argument("--%s-by-type" % self._shortname,
                                         dest="%s_by_type" % self._shortname,
                                         action="store_true",
                                         help="Group %s values by lipid types" % self._fullname)

    def prepare(self, num_frames):
        """Allocate array to store values.
        :param num_frames: The number of frames considered
        :type num_frames: int
        :returns: None
        """
        self.membrane_avg_values = np.zeros(num_frames)
        self.leaflet1_avg_values = np.zeros(num_frames)
        self.leaflet2_avg_values = np.zeros(num_frames)

        if self._groupable_by_type:
            self.needs_grouping = getattr(self.parent.namespace, "%s_by_type" % self._shortname)
            self.membrane_avg_values_by_type = {}
            self.leaflet1_avg_values_by_type = {}
            self.leaflet2_avg_values_by_type = {}

    def get_values(self, membrane):  # pragma: no cover
        """Retrieve property values from a membrane.
        :param membrane: A membrane instance
        :type membrane: fatslimlib.core_analysis.pyx.Membrane
        :returns: membrane value, leaflet #1 value, leaflet #2 value
        :rtype: tuple
        """
        raise NotImplementedError

    def get_results(self, limit=-1):
        """Computes and return average values.
        :param int limit: Specify the number of values to take into account (negative values are
                          ignored)
        :returns: A 6-item tuple that contains:
                  mean value for membrane, standard deviation for membrane value,
                  mean value for leaflet #1, standard deviation for leaflet #1 value,
                  mean value for leaflet #2, standard deviation for leaflet #2 value
        :rtype: tuple
        """
        if self.membrane_avg_values is None:
            return 0, 0, 0, 0, 0, 0

        if limit > 0:
            membrane_values = self.membrane_avg_values[:limit]
            leaflet1_values = self.leaflet1_avg_values[:limit]
            leaflet2_values = self.leaflet2_avg_values[:limit]
        else:
            membrane_values = self.membrane_avg_values
            leaflet1_values = self.leaflet1_avg_values
            leaflet2_values = self.leaflet2_avg_values

        return membrane_values.mean(), membrane_values.std(), \
               leaflet1_values.mean(), leaflet1_values.std(), \
               leaflet2_values.mean(), leaflet2_values.std()

    def save_values_to_xvg(self, filename):
        """Save values to xvg file.
        :param str filename: Name of the .xvg file.
        :returns: An output to be displayed for information
        :rtype: str"""

        fname, ext = os.path.splitext(filename)
        xvg_file = "%s.xvg" % fname

        output = backup_file(xvg_file)

        with open(xvg_file, "w") as fp:
            output_plot = self.parent.get_xvg_header()
            output_plot += "@TYPE xy\n"
            output_plot += "@ title \"%s\"\n" % self.fullname.capitalize()
            output_plot += "@ xaxis  label \"Time (ps)\"\n"
            output_plot += "@ yaxis  label \"%s (%s)\"\n" % (self.fullname.capitalize(),
                                                             self.units)
            output_plot += "@ legend on\n"
            output_plot += "@ legend box on\n"
            output_plot += "@ legend loctype view\n"
            output_plot += "@ legend 0.78, 0.8\n"
            output_plot += "@ legend length 2\n"

            if self.needs_grouping:
                counter = 0
                membrane_types = sorted(self.membrane_avg_values_by_type.keys())
                for name in membrane_types:
                    output_plot += "@ s%i legend \"Membrane %s\"\n" % (counter, name)
                    counter += 1

                leaflet1_types = sorted(self.leaflet1_avg_values_by_type.keys())
                for name in leaflet1_types:
                    output_plot += "@ s%i legend \"%s %s\"\n" % (counter,
                                                                 self.parent.leaflet1_type.
                                                                 capitalize(),
                                                                 name)
                    counter += 1

                leaflet2_types = sorted(self.leaflet2_avg_values_by_type.keys())
                for name in leaflet2_types:
                    output_plot += "@ s%i legend \"%s %s\"\n" % (counter,
                                                                 self.parent.leaflet2_type.
                                                                 capitalize(),
                                                                 name)
                    counter += 1

                for i, timestep in enumerate(self.parent.timesteps):
                    output_plot += "%8.3f " % timestep

                    for key in membrane_types:
                        output_plot += "%8.3f " % self.membrane_avg_values_by_type[key][i]

                    for key in leaflet1_types:
                        output_plot += "%8.3f " % self.leaflet1_avg_values_by_type[key][i]

                    for key in leaflet2_types:
                        output_plot += "%8.3f " % self.leaflet2_avg_values_by_type[key][i]

                    output_plot += "\n"

            else:
                output_plot += "@ s0 legend \"Membrane\"\n"
                output_plot += "@ s1 legend \"%s\"\n" % self.parent.leaflet1_type.capitalize()
                output_plot += "@ s2 legend \"%s\"\n" % self.parent.leaflet2_type.capitalize()

                for i, timestep in enumerate(self.parent.timesteps):
                    output_plot += "%8.3f " % timestep
                    output_plot += "%8.3f " % self.membrane_avg_values[i]
                    output_plot += "%8.3f " % self.leaflet1_avg_values[i]
                    output_plot += "%8.3f " % self.leaflet2_avg_values[i]
                    output_plot += "\n"

            fp.write(output_plot)

        output += "%s values saved to '%s'\n" % (self.fullname.capitalize(),
                                                 xvg_file)
        return output

    def get_raw_values(self, membrane):  # pragma: no cover
        """Returns the raw property values for both leaflets.
        :param fatslimlib.core_analysis.pyx.Membrane membrane: the membrane used to compute property.
        :returns: (leaflet #1 raw values, leaflet #2 raw values)
        :rtype: tuple
        """
        raise NotImplementedError

    def save_raw_values_to_csv(self, membrane, filename):
        """Compute and save raw values to csv file.
        :param fatslimlib.core_analysis.pyx.Membrane membrane: the membrane used to compute property.
        :param str filename: Name of the .csv file."""

        fname, ext = os.path.splitext(filename)
        csv_file = "%s.csv" % fname

        output = backup_file(csv_file)

        with open(csv_file, "w") as fp:
            writer = csv.writer(fp)
            writer.writerow(['resid', 'leaflet',
                             'X coords', 'Y coords', 'Z coords',
                             '%s' % self.fullname.capitalize()])

            for i, raw_values in enumerate(self.get_raw_values(membrane)):
                resids = membrane[i].resids
                coords = membrane[i].coords
                if i == 0:
                    leaflet_type = self.parent.leaflet1_type
                else:
                    leaflet_type = self.parent.leaflet2_type

                for j in range(len(resids)):
                    writer.writerow([resids[j], leaflet_type,
                                     "%.3f" % coords[j][0],
                                     "%.3f" % coords[j][1],
                                     "%.3f" % coords[j][2],
                                     "%.3f" % raw_values[j]])

        output += "%s raw values saved to '%s'\n" % (self.fullname.capitalize(),
                                                     csv_file)
        return output


class MembraneAnalysisCommand(AnalyticalCommand):
    _property_classes = []

    def __init__(self, registry):
        self.membrane_nlipids = None
        self.membrane_type = "unknown"

        self.leaflet1_nlipids = None
        self.leaflet1_type = "unknown"

        self.leaflet2_nlipids = None
        self.leaflet2_type = "unknown"

        self.properties = []
        self.property_needs_raw = {}
        for klass in self._property_classes:
            inst = klass(self)
            self.properties.append(inst)
            self.property_needs_raw[inst.shortname] = False


        super(MembraneAnalysisCommand, self).__init__(registry)

    def fill_parser(self):
        super(MembraneAnalysisCommand, self).fill_parser()

        self.parser_analysis_group.add_argument("--idfreq", default=0, type=int,
                                                help="Frequency used to update membrane "
                                                     "identification")

        self.parser_analysis_group.add_argument("--cutoff", default=2, type=float,
                                                help="Cutoff distance for leaflets identification")

        for prop in self.properties:
            prop.fill_parsers(self.parser_analysis_group,
                              self.parser_input_group,
                              self.parser_output_group)

    def prepare_results(self, num_frames):
        super(MembraneAnalysisCommand, self).prepare_results(num_frames)

        self.membrane_nlipids = np.zeros(num_frames, dtype=np.int)
        self.leaflet1_nlipids = np.zeros(num_frames, dtype=np.int)
        self.leaflet2_nlipids = np.zeros(num_frames, dtype=np.int)

        for prop in self.properties:
            prop.prepare(num_frames)

            if prop.raw_values_supported:
                fname = getattr(self.namespace, "export_%s_raw" % prop.shortname)
                if fname is not None:
                    self.property_needs_raw[prop.shortname] = os.path.splitext(fname)[0]
                else:
                    self.property_needs_raw[prop.shortname] = None

    def process_membrane(self, membrane):
        output = "  Membrane (%s - %i lipids) - Leaflet #1 (%s - %i lipids) - Leaflet #2 " \
                 "(%s - %i lipids)\n" % (self.membrane_type,
                                         self.membrane_nlipids[self.processing_index],
                                         self.leaflet1_type,
                                         self.leaflet1_nlipids[self.processing_index],
                                         self.leaflet2_type,
                                         self.leaflet2_nlipids[self.processing_index])

        for prop in self.properties:
            membrane_avg, leaflet1_avg, leaflet2_avg = prop.get_values(membrane)

            output += "  -> Average %s values (%s): Membrane=%.3f - Leaflet #1=%.3f - " \
                      "Leaflet #2=%.3f\n" % (prop.fullname.lower(),
                                             prop.units,
                                             membrane_avg,
                                             leaflet1_avg,
                                             leaflet2_avg)

            if self.property_needs_raw[prop.shortname]:
                fname = "%s_frame_%05i.csv" % (self.property_needs_raw[prop.shortname],
                                               self.processing_index)
                output += prop.save_raw_values_to_csv(membrane, fname)
        return output

    def process_frame(self, frame):
        output = super(MembraneAnalysisCommand, self).process_frame(frame)
        try:
            membranes = frame.get_membranes(self.namespace.cutoff)
        except (ValueError, IndexError):  # pragma: no cover
            output += "->No membrane found!"
        else:
            n_membranes = len(membranes)
            output = "->%i membrane%s found!\n" % (n_membranes, "s"[n_membranes == 1:])

            if len(membranes) > 1:
                self.print_warning("%i membranes found but only the first one will "
                                   "be analyzed!" % len(membranes))

            membrane = membranes[0]

            self.leaflet1_nlipids[self.processing_index] = len(membrane[0])
            self.leaflet2_nlipids[self.processing_index] = len(membrane[1])
            self.membrane_nlipids[self.processing_index] = len(membrane)

            if membrane.is_planar:
                self.membrane_type = "planar"
                self.leaflet1_type = "lower leaflet"
                self.leaflet2_type = "upper leaflet"
            elif membrane.is_vesicle:
                self.membrane_type = "vesicle"
                self.leaflet1_type = "outer leaflet"
                self.leaflet2_type = "inner leaflet"
            else:  # pragma: no cover
                self.membrane_type = "Unknown"
                self.leaflet1_type = "leaflet #1"
                self.leaflet1_type = "leaflet #2"

            output += self.process_membrane(membrane)

        return output

    def process_results(self):
        output = u"Average values over %i processed frames:\n" % (self.processing_index + 1)

        for prop in self.properties:
            m_mean, m_std, l1_mean, l1_std, l2_mean, l2_std = \
                prop.get_results(self.processing_index+1)

            output += u"%s: Membrane=%.3f\u00b1%.3f%s - " \
                      u"%s=%.3f\u00b1%.3f%s - " \
                      u"%s=%.3f\u00b1%.3f%s\n" % (prop.fullname.capitalize(),
                                                m_mean, m_std, prop.units,
                                                self.leaflet1_type.capitalize(),
                                                l1_mean, l1_std, prop.units,
                                                self.leaflet2_type.capitalize(),
                                                l2_mean, l2_std, prop.units,
                                                )

            xvg_file = getattr(self.namespace, "plot_%s" % prop.shortname)
            if xvg_file is not None:
                output += u"%s\n" % prop.save_values_to_xvg(xvg_file)

        return output
