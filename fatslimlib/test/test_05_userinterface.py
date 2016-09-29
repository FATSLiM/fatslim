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
from __future__ import print_function
from numpy.testing import assert_allclose
ATOL = 5e-3

try:
    from StringIO import StringIO  # python 2, not cStringIO due to unicode strings
except ImportError:
    from io import StringIO  # python 3
from contextlib import contextmanager
import sys
from unittest import TestCase
import pytest
import os
from tempfile import mkdtemp
import filecmp
import numpy as np

from .. import core_base
from .. import print_full_version_info, print_greetings, print_goodbye, main,\
    _format_version_tuple, __version__,\
    __description__, __copyright__
from . import data


@contextmanager
def stdout_redirector(stream):
    old_stdout = sys.stdout
    old_stderr = sys.stderr
    sys.stdout = stream
    sys.stderr = stream
    try:
        yield
    finally:
        sys.stdout = old_stdout
        sys.stderr = old_stderr

@contextmanager
def stderr_redirector(stream):
    old_stderr = sys.stderr
    sys.stderr = stream
    try:
        yield
    finally:
        sys.stderr = old_stderr


def test_version_formatting():
    assert _format_version_tuple((1, 0)) == "1.0"
    assert _format_version_tuple((1, 0, 0, 'final', 0)) == "1.0"
    assert _format_version_tuple((1, 0, 0, 'final')) == "1.0"
    assert _format_version_tuple((1, 0, 1, 'final', 0)) == "1.0.1"
    assert _format_version_tuple((1, 4, 0)) == "1.4"
    assert _format_version_tuple((1, 2, 0, 'dev', 0)) == "1.2.dev0"
    assert _format_version_tuple((1, 1, 1, 'candidate', 2)) == "1.1.1rc2"
    assert _format_version_tuple((2, 1, 0, 'alpha')) == "2.1a0"
    assert _format_version_tuple((2, 1, 0, 'beta', 1)) == "2.1b1"
    assert _format_version_tuple((2, 1, 0, 'final', 42)) == "2.1"
    assert _format_version_tuple((1, 4, 0, 'post', 0)) == "1.4.post0"


def test_greetings():
    header = u"FATSLiM - %s" % __description__
    fsl_version = u"version %s" % __version__
    copyright = u"%s" % __copyright__

    lines = [header, fsl_version, copyright, '']

    buf = StringIO()
    with stdout_redirector(buf):
        print_greetings()
    output = buf.getvalue()
    buf.close()


    lino = 0
    for line in output.splitlines():
        if len(line) == 0:
            continue
        assert line == lines[lino]
        lino += 1


def test_goodbye():
    ok = "Goodbye!"
    notok = "Sorry, I made a boo boo..."

    results = [ok, notok]

    for i in range(2):
        buf = StringIO()
        with stdout_redirector(buf):
            print_goodbye(i)
        output = buf.getvalue().strip()
        buf.close()

        assert output == results[i]


def build_version_string():
    import numpy
    import platform

    header = "FATSLiM - %s" % __description__
    fsl_version = "FATSLiM version: %s" % __version__
    py_version = "Python version: %s (%s)" % (platform.python_version(), sys.executable)
    cython_version = "Cython version (file generation): %s" % core_base.__cython_version__
    compiler_version = "Python compiler: %s" % platform.python_compiler()
    arch = "CPU architecture: %s" % platform.architecture()[0]
    openmp_info = "OpenMP: %i CPUs (default number of threads: %i - max: %i)" % (core_base.get_num_procs(),
                                                                         core_base.get_num_threads(),
                                                                         core_base.get_max_threads())
    numpy_version = "NumPy version: %s" % numpy.__version__

    return [header, fsl_version, py_version, cython_version, compiler_version, arch, openmp_info,
             numpy_version]


def test_full_version_string():
    lines = build_version_string()

    buf = StringIO()
    with stdout_redirector(buf):
        print_full_version_info()
    output = buf.getvalue()
    buf.close()

    lino = 0
    for line in output.splitlines():
        if len(line) == 0:
            continue
        assert line == lines[lino]
        lino += 1


def test_main_funct():

    buf = StringIO()
    with stdout_redirector(buf):
        exit_code = main(argv=[])
    buf.close()

    assert exit_code == 0

    buf = StringIO()
    with stdout_redirector(buf):
        exit_code = main(argv=["undefinedcommand"])
    buf.close()

    assert exit_code == 1


def test_command_no_registry():
    from ..command import Command

    with pytest.raises(TypeError):
        Command(None)


def test_registry_command():
    from ..command import CommandRegistry

    test = CommandRegistry("test")

    with pytest.raises(TypeError):
        test.add_command(None)


def test_property_parent():
    from ..command import MembraneProperty

    with pytest.raises(TypeError):
        MembraneProperty(None)


class TestCommands(TestCase):
    @classmethod
    def setUpClass(cls):
        from .. import builtins
        from .. import command

        cls._command_registry = command.CommandRegistry("test")

        for name in dir(builtins):
            if not name.startswith("_"):
                attr = getattr(builtins, name)
                try:
                    if issubclass(attr, command.Command) and name.startswith("Cmd"):
                        cls._command_registry.add_command(attr)
                except TypeError:
                    pass

        cls.tempdir = mkdtemp("fsl_test")

    @classmethod
    def tearDownClass(cls):
        for fname in os.listdir(cls.tempdir):
            print("Deleting: %s" % os.path.join(cls.tempdir, fname))
            os.remove(os.path.join(cls.tempdir, fname))
        os.rmdir(cls.tempdir)

    def run_command(self, argv):
        buf = StringIO()
        with stdout_redirector(buf):
            exit_code = self._command_registry.run(list(argv))
        output = buf.getvalue()
        buf.close()
        assert exit_code == 0, "Command '%s' crashed with args: %s\nCaptured output:\n%s" % \
                               (argv[0], ",".join(argv[1:]), output)
        return output

    def run_command_nocheck(self, argv):
        buf = StringIO()
        with stdout_redirector(buf):
            self._command_registry.run(list(argv))
        output = buf.getvalue()
        buf.close()
        return output

    def assert_filecmp(self, test_output, ref_output):
        if filecmp.cmp(test_output, ref_output):  # Files are identical
            return

        with open(ref_output) as fp:
            ref_content = fp.readlines()

        with open(test_output) as fp:
            lino = 0
            for line in fp:
                if line.startswith("#"):
                    continue
                else:
                    if os.path.splitext(ref_output)[1] == ".xvg" and not line.startswith("@"):
                        ref_values = np.array([float(val) for val in ref_content[lino].split()])
                        actual_values = np.array([float(val) for val in line.split()])

                        for i, ref_value in enumerate(ref_values):
                            assert_allclose(actual_values[i], ref_value, atol=ATOL)

                    else:
                        assert line.strip() == ref_content[lino].strip(), "Reference value is %s " \
                                                                          "but actual value is " \
                                                                          "%s" % \
                                                                          (repr(ref_content[lino]),
                                                                          repr(line))
                    lino += 1

    def test_util_backup(self):
        import shutil
        from .. import util

        tmpfile = os.path.join(self.tempdir, "test.gro")
        backupbackupfile1 = "%s.%02i.old" % (tmpfile, 1)
        backupbackupfile2 = "%s.%02i.old" % (tmpfile, 2)

        with open(tmpfile, 'w') as fp:
            fp.write("Alive")
        shutil.copy(tmpfile, backupbackupfile1)
        output = util.backup_file(tmpfile)

        assert output == "'%s' backed up to '%s'\n" % (tmpfile, backupbackupfile2)

        os.unlink(backupbackupfile1)
        os.unlink(backupbackupfile2)

    def test_registry_add_bad_command(self):
        with pytest.raises(TypeError):
            self._command_registry.add_command(None)

    def test_registry_add_registered_command(self):
        from .. import builtins
        from .. import command

        with pytest.raises(command.AlreadyRegisteredCommandError):
            self._command_registry.add_command(builtins.CmdVersion)

    def test_registry_print_debug(self):
        buf = StringIO()
        with stdout_redirector(buf):
            self._command_registry.print_debug("Test")
        output = buf.getvalue()
        buf.close()

        assert output == ""

        self._command_registry.debug_mode = True

        buf = StringIO()
        with stdout_redirector(buf):
            self._command_registry.print_debug("Test")
        output = buf.getvalue()
        buf.close()

        assert output == "DEBUG: Test\n"

        self._command_registry.debug_mode = False

    def test_registry_print_warnings(self):
        buf = StringIO()
        with stderr_redirector(buf):
            self._command_registry.print_warning("Test")
        output = buf.getvalue()
        buf.close()

        assert output == "WARNING: Test\n"

    def test_registry_print_error(self):
        buf = StringIO()
        with stderr_redirector(buf):
            self._command_registry.print_error("Test")
        output = buf.getvalue()
        buf.close()

        assert output == "ERROR: Test\n"

    def test_unknown_command(self):
        output = self.run_command_nocheck(["undefinedcommand"]).splitlines()

        assert output[0] == "ERROR: Unknown command: undefinedcommand. See 'test --help'"

    def test_unknown_argument(self):
        output = self.run_command_nocheck(["version", "whatever"]).splitlines()

        assert output[0] == "ERROR: 'version' command does not recognize the following argument:" \
                            " whatever"

    def test_command_failure(self):
        with pytest.raises(IOError):
            output = self.run_command(["apl"])

    def test_help_command_empty(self):
        output = self.run_command(["--help"])

        assert "usage: " in output
        assert "Available commands:" in output
        assert "See 'test help <command>' for more information on a specific command" in output
        assert "'help' command executed" in output

    def test_command_empty(self):
        output = self.run_command([""])

        assert "usage: " in output
        assert "Available commands:" in output
        assert "See 'test help <command>' for more information on a specific command" in output
        assert "'help' command executed" in output

    def test_help_command_unknown(self):
        output = self.run_command(["help", "unknown"]).splitlines()

        assert output[0] == "Unkown command: unknown. See 'test --help"

    def test_help_command_normal(self):
        output = self.run_command(["help", "help"])

        inside_posarg = False
        inside_options = False
        for line in output.splitlines():
            line = line.strip()
            if len(line) == 0:
                continue

            if line.startswith("command: "):
                assert line == "command: test help"
                continue
            elif line.startswith("purpose: "):
                assert line == "purpose: Displays help about FATSLiM's commands"
                continue
            elif line.startswith("usage: "):
                assert line == "usage: test help [--help] [--debug] [--verbose] [command]"
                continue
            elif line.startswith("positional arguments"):
                inside_options = False
                inside_posarg = True
                continue
            elif line.startswith("General options:"):
                inside_posarg = False
                inside_options = True
                continue

            if inside_options:
                if line.startswith("--"):
                    options = ['--help, -h     See help about help command (default: False)',
                               '--debug        Enable debug output (default: False)',
                               '--verbose, -v  Be loud and noisy (default: False)']
                    assert line in options, "unknown option line: '%s'" % line
                    continue

            elif inside_posarg:
                assert line == "command        Command to display information about " \
                               "(default: None)"
                continue

            if line.startswith("'help' command executed in"):
                continue

            assert False, "Unkown line content: '%s'" % line

    def test_help_command_alternate(self):
        output = self.run_command(["version", "--help"])

        assert "Displays FATSLiM version" in output

        assert "'version' command executed" in output

    def test_help_command_alternate2(self):
        output = self.run_command(["version", "-h"])

        assert "Displays FATSLiM version" in output

        assert "'version' command executed" in output

    def test_command_fullname(self):
        assert self._command_registry.commands["version"].fullname == "test-version"

    def test_command_aliases(self):
        assert self._command_registry.commands["version"].aliases == []

    def test_command_print_debug(self):
        cmd = self._command_registry.commands["version"]

        buf = StringIO()
        with stdout_redirector(buf):
            cmd.print_debug("Test")
        output = buf.getvalue()
        buf.close()

        assert output == ""

        cmd.debug_mode = True

        buf = StringIO()
        with stdout_redirector(buf):
            cmd.print_debug("Test")
        output = buf.getvalue()
        buf.close()

        assert output == "DEBUG: Test\n"

        cmd.debug_mode = False

    def test_command_print_verbose(self):
        cmd = self._command_registry.commands["version"]

        cmd.verbose = False

        buf = StringIO()
        with stdout_redirector(buf):
            cmd.print_verbose("Test")
        output = buf.getvalue()
        buf.close()

        assert output == ""

        cmd.verbose = True

        buf = StringIO()
        with stdout_redirector(buf):
            cmd.print_verbose("Test")
        output = buf.getvalue()
        buf.close()

        assert output == "Test\n"



    def test_cmd_print_warnings(self):
        buf = StringIO()
        with stderr_redirector(buf):
            self._command_registry.commands["version"].print_warning("Test")
        output = buf.getvalue()
        buf.close()

        assert output == "WARNING: Test\n"

    def test_version_command(self):
        lines = build_version_string()
        output = self.run_command(["version"])

        output = output.splitlines()

        for line in lines:
            output_line = ""
            while output_line == "":
                output_line = output.pop(0).strip()
            assert line == output_line, "'%s' != '%s'" % (line, output_line)

        exc_time_line = output.pop(0).strip()

        assert exc_time_line.startswith("'version' command executed in")

        if len(output) != 0:
            assert False, "Unkown remaining content: '%s'" % str(output)

    def test_version_command_alternate(self):
        output = self.run_command(["--version"])

        assert "'version' command executed" in output

    def test_version_help(self):
        output = self.run_command(["--version", "--help"])

        assert "'help' command executed" in output

    def test_aggregate_command_ok(self):
        cutoff = "2.0"
        agg_hg_ndx = "agg_hg"
        output_agg_hg = os.path.join(self.tempdir, "%s.ndx" % agg_hg_ndx)
        agg_lipid_ndx = "agg_lipids"
        output_agg_lipids = os.path.join(self.tempdir, "%s.ndx" % agg_lipid_ndx)
        output_agg_number = os.path.join(self.tempdir, "number_agg.bad")

        output = self.run_command(["aggregates",
                                   "-c", data.MODEL_BILAYER_GRO,
                                   "-n", data.MODEL_BILAYER_NDX,
                                   "-o", output_agg_number,
                                   "--cutoff", cutoff,
                                   "--output-index", output_agg_lipids,
                                   "--output-index-hg", output_agg_hg,
                                   ],)

        assert "'aggregates' command executed" in output

        # Check that the index group for hg is fine
        with open(os.path.join(self.tempdir, "%s_0000.ndx" % agg_hg_ndx)) as fp:
            in_agg1 = False
            in_agg2 = False
            agg1 = []
            agg2 = []
            for line in fp:
                line = line.strip()
                if len(line) == 0:
                    continue

                elif line == "[ aggregate_1 ]":
                    in_agg1 = True
                    in_agg2 = False
                    continue
                elif line == "[ aggregate_2 ]":
                    in_agg1 = False
                    in_agg2 = True
                    continue

                if in_agg1:
                    agg1.extend([int(val) for val in line.split()])
                    continue
                elif in_agg2:
                    agg2.extend([int(val) for val in line.split()])
                    continue

                assert False, "Unknown line: '%s'" % line

        for i in range(2):
            ref_agg = []

            for j in range(36):
                ref_agg.append(4 + (i*36+j) * 50)
                ref_agg.append(8 + (i*36+j) * 50)

            if i == 0:
                np.testing.assert_array_equal(agg1, ref_agg)
            else:
                np.testing.assert_array_equal(agg2, ref_agg)

        os.unlink(os.path.join(self.tempdir, "%s_0000.ndx" % agg_hg_ndx))

        # Check that the index group for lipids is fine
        with open(os.path.join(self.tempdir, "%s_0000.ndx" % agg_lipid_ndx)) as fp:
            in_agg1 = False
            in_agg2 = False
            agg1 = []
            agg2 = []
            for line in fp:
                line = line.strip()
                if len(line) == 0:
                    continue

                elif line == "[ aggregate_1 ]":
                    in_agg1 = True
                    in_agg2 = False
                    continue
                elif line == "[ aggregate_2 ]":
                    in_agg1 = False
                    in_agg2 = True
                    continue

                if in_agg1:
                    agg1.extend([int(val) for val in line.split()])
                    continue
                elif in_agg2:
                    agg2.extend([int(val) for val in line.split()])
                    continue

                assert False, "Unknown line: '%s'" % line

        for i in range(2):
            ref_agg = np.arange(36*50) + 1 + i*(36*50)

            if i == 0:
                np.testing.assert_array_equal(agg1, ref_agg)
            else:
                np.testing.assert_array_equal(agg2, ref_agg)

        os.unlink(os.path.join(self.tempdir, "%s_0000.ndx" % agg_lipid_ndx))

        # Check output xvg
        with open("%s.xvg" % os.path.splitext(output_agg_number)[0]) as fp:
            for line in fp:
                line = line.strip()

                if line == "":
                    continue
                if line.startswith("#"):
                    continue
                if line.startswith("@"):
                    assert line in ['@ title "Number of aggregates"',
                                    '@ xaxis label "Time (ps)"',
                                    '@ yaxis label "Number of aggregates"',
                                    '@TYPE xy']
                else:
                    assert line == "0.000       2"

        os.unlink("%s.xvg" % os.path.splitext(output_agg_number)[0])

    def test_aggregate_command_not_ok(self):
        cutoff = "0.5"
        output_agg_number = os.path.join(self.tempdir, "number_agg.bad")

        output = self.run_command(["aggregates",
                                   "-c", data.MODEL_BILAYER_GRO,
                                   "-n", data.MODEL_BILAYER_NDX,
                                   "-o", output_agg_number,
                                   "--cutoff", cutoff,
                                   ],)

        assert "'aggregates' command executed" in output

        # Check output xvg
        with open("%s.xvg" % os.path.splitext(output_agg_number)[0]) as fp:
            for line in fp:
                line = line.strip()

                if line == "" or line.startswith("#") or line.startswith("@"):
                    continue
                else:
                    assert line == "0.000      72"  # with used cutoff there should be no real agg

        os.unlink("%s.xvg" % os.path.splitext(output_agg_number)[0])

    def test_analysis_no_frame(self):
        output = self.run_command(["aggregates",
                                   "-c", data.MODEL_BILAYER_GRO,
                                   "-n", data.MODEL_BILAYER_NDX,
                                   "-b", "1000000"
                                   ],)

        print (output)

        assert "WARNING: No frames are selected for analysis!" in output

    def test_membrane_command_ok(self):
        cutoff = "2.0"
        membrane_hg = "membrane_hg"
        output_membrane_hg = os.path.join(self.tempdir, "%s.ndx" % membrane_hg)
        membrane_lipid_ndx = "membrane_lipids"
        output_membrane_lipids = os.path.join(self.tempdir, "%s.ndx" % membrane_lipid_ndx)
        output_membrane_number = os.path.join(self.tempdir, "number_membrane.xvg")

        output = self.run_command(["membranes",
                                   "-c", data.MODEL_BILAYER_GRO,
                                   "-n", data.MODEL_BILAYER_NDX,
                                   "-o", output_membrane_number,
                                   "--cutoff", cutoff,
                                   "--output-index", output_membrane_lipids,
                                   "--output-index-hg", output_membrane_hg,
                                   ],)

        assert "'membranes' command executed" in output

        # Check that the index group for hg is fine
        with open(os.path.join(self.tempdir, "%s_0000.ndx" % membrane_hg)) as fp:
            in_leaflet1 = False
            in_leaflet2 = False
            leaflet1 = []
            leaflet2 = []
            for line in fp:
                line = line.strip()
                if len(line) == 0:
                    continue

                elif line == "[ membrane_1_leaflet_1 ]":
                    in_leaflet1 = True
                    in_leaflet2 = False
                    continue
                elif line == "[ membrane_1_leaflet_2 ]":
                    in_leaflet1 = False
                    in_leaflet2 = True
                    continue

                if in_leaflet1:
                    leaflet1.extend([int(val) for val in line.split()])
                    continue
                elif in_leaflet2:
                    leaflet2.extend([int(val) for val in line.split()])
                    continue

                assert False, "Unknown line: '%s'" % line

        for i in range(2):
            ref_leaflet = []

            for j in range(36):
                ref_leaflet.append(4 + (i*36+j) * 50)
                ref_leaflet.append(8 + (i*36+j) * 50)

            if i == 1:
                np.testing.assert_array_equal(leaflet1, ref_leaflet)
            else:
                np.testing.assert_array_equal(leaflet2, ref_leaflet)

        os.unlink(os.path.join(self.tempdir, "%s_0000.ndx" % membrane_hg))

        # Check that the index group for lipids is fine
        with open(os.path.join(self.tempdir, "%s_0000.ndx" % membrane_lipid_ndx)) as fp:
            in_leaflet1 = False
            in_leaflet2 = False
            leaflet1 = []
            leaflet2 = []
            for line in fp:
                line = line.strip()
                if len(line) == 0:
                    continue

                elif line == "[ membrane_1_leaflet_1 ]":
                    in_leaflet1 = True
                    in_leaflet2 = False
                    continue
                elif line == "[ membrane_1_leaflet_2 ]":
                    in_leaflet1 = False
                    in_leaflet2 = True
                    continue

                if in_leaflet1:
                    leaflet1.extend([int(val) for val in line.split()])
                    continue
                elif in_leaflet2:
                    leaflet2.extend([int(val) for val in line.split()])
                    continue

                assert False, "Unknown line: '%s'" % line

        for i in range(2):
            ref_leaflet = np.arange(36*50) + 1 + i*(36*50)

            if i == 1:
                np.testing.assert_array_equal(leaflet1, ref_leaflet)
            else:
                np.testing.assert_array_equal(leaflet2, ref_leaflet)

        os.unlink(os.path.join(self.tempdir, "%s_0000.ndx" % membrane_lipid_ndx))

        # Check output xvg
        with open("%s.xvg" % os.path.splitext(output_membrane_number)[0]) as fp:
            for line in fp:
                line = line.strip()

                if line == "":
                    continue
                if line.startswith("#"):
                    continue
                if line.startswith("@"):
                    assert line in ['@ title "Number of membranes"',
                                    '@ xaxis label "Time (ps)"',
                                    '@ yaxis label "Number of membranes"',
                                    '@TYPE xy']
                else:
                    assert line == "0.000       1"

        os.unlink("%s.xvg" % os.path.splitext(output_membrane_number)[0])

    def test_membrane_command_vesicle(self):

        output = self.run_command(["membranes",
                                   "-c", data.MODEL_VESICLE_GRO,
                                   "-n", data.MODEL_VESICLE_NDX])

        assert "'membranes' command executed" in output

        assert os.listdir(self.tempdir) == []

    def test_id_freq(self):
        output = self.run_command(["thickness",
                                   "-c", data.VESICLE_GRO,
                                   "-n", data.VESICLE_NDX,
                                   "-t", data.VESICLE_XTC,
                                   "--verbose",
                                   "--idfreq", "0"])

        assert output.count("1 membrane identified!") == 1
        assert output.count("identification skipped!") == 10

        output = self.run_command(["thickness",
                                   "-c", data.VESICLE_GRO,
                                   "-n", data.VESICLE_NDX,
                                   "-t", data.VESICLE_XTC,
                                   "--verbose",
                                   "--idfreq", "5"])
        assert output.count("1 membrane identified!") == 3
        assert output.count("identification skipped!") == 8

    def test_thickness_command_planar(self):
        output_file = os.path.join(self.tempdir, "test_thickness.xvg")

        output = self.run_command(["thickness",
                                   "-c", data.MODEL_BILAYER_GRO,
                                   "-n", data.MODEL_BILAYER_NDX,
                                   "--plot-thickness", output_file,
                                   ],)

        assert "'thickness' command executed" in output

        self.assert_filecmp(output_file, data.MODEL_BILAYER_THICKNESS_XVG)

        os.unlink(output_file)

    def test_thickness_command_vesicle_xtc(self):
        output_file = os.path.join(self.tempdir, "test_thickness.xvg")

        output = self.run_command(["thickness",
                                   "-c", data.VESICLE_GRO,
                                   "-n", data.VESICLE_NDX,
                                   "-t", data.VESICLE_XTC,
                                   "--plot-thickness", output_file,
                                   ],)

        assert "'thickness' command executed" in output

        self.assert_filecmp(output_file, data.VESICLE_THICKNESS_XVG)

        os.unlink(output_file)

    def test_thickness_raw(self):
        output_file = os.path.join(self.tempdir, "test_thickness.csv")

        output = self.run_command(["thickness",
                                   "-c", data.MODEL_BILAYER_GRO,
                                   "-n", data.MODEL_BILAYER_NDX,
                                   "--export-thickness-raw", output_file,
                                   ],)

        assert "'thickness' command executed" in output
        output_file = "%s_frame_00000.csv" % os.path.splitext(output_file)[0]

        self.assert_filecmp(output_file, data.MODEL_BILAYER_THICKNESS_CSV)

        os.unlink(output_file)

    def test_thickness_multiple(self):
        output = self.run_command(["thickness",
                                   "-c", data.MODEL_MULTIBILAYER_GRO,
                                   "-n", data.MODEL_MULTIBILAYER_NDX,
                                   ],)

        assert "2 membranes found but only the first one will be analyzed!" in output

    def test_property_limit(self):
        cmd = self._command_registry.commands["thickness"]

        buf = StringIO()
        with stdout_redirector(buf):
            cmd.run_with_args(["-c", data.VESICLE_GRO,
                               "-n", data.VESICLE_NDX,
                               "-t", data.VESICLE_XTC])
        buf.close()

        prop = cmd.properties[0]
        assert prop.get_results(-1)[0] == prop.membrane_avg_values.mean()
        assert prop.get_results(2)[0] == prop.membrane_avg_values[:2].mean()

    def test_apl_command_nogroup(self):
        output_file = os.path.join(self.tempdir, "test_apl.xvg")

        output = self.run_command(["apl",
                                   "-c", data.BILAYER_CHOL_GRO,
                                   "-n", data.BILAYER_CHOL_NDX,
                                   "--plot-apl", output_file,
                                   ],)

        assert "'apl' command executed" in output

        self.assert_filecmp(output_file, data.BILAYER_CHOL_APL_NOGROUP)

        os.unlink(output_file)

    def test_apl_command_by_type(self):
        output_file = os.path.join(self.tempdir, "test_apl.xvg")

        output = self.run_command(["apl",
                                   "-c", data.BILAYER_CHOL_GRO,
                                   "-n", data.BILAYER_CHOL_NDX,
                                   "--plot-apl", output_file,
                                   "--apl-by-type"
                                   ],)

        assert "'apl' command executed" in output

        self.assert_filecmp(output_file, data.BILAYER_CHOL_APL_GROUPED)

        os.unlink(output_file)

    def test_apl_raw(self):
        output_file = os.path.join(self.tempdir, "test_apl.csv")

        output = self.run_command(["apl",
                                   "-c", data.MODEL_BILAYER_GRO,
                                   "-n", data.MODEL_BILAYER_NDX,
                                   "--export-apl-raw", output_file,
                                   ],)

        assert "'apl' command executed" in output
        output_file = "%s_frame_00000.csv" % os.path.splitext(output_file)[0]

        self.assert_filecmp(output_file, data.MODEL_BILAYER_APL_CSV)

        os.unlink(output_file)

    def test_apl_command_vesicle_xtc(self):
        output_file = os.path.join(self.tempdir, "test_apl.xvg")

        output = self.run_command(["apl",
                                   "-c", data.VESICLE_GRO,
                                   "-n", data.VESICLE_NDX,
                                   "-t", data.VESICLE_XTC,
                                   "--plot-apl", output_file,
                                   ],)

        assert "'apl' command executed" in output

        self.assert_filecmp(output_file, data.VESICLE_APL_XVG)

        os.unlink(output_file)

    def test_area_command_nogroup(self):
        output_file = os.path.join(self.tempdir, "test_area.xvg")
        output = self.run_command(["apl",
                                   "-c", data.VESICLE_GRO,
                                   "-n", data.VESICLE_NDX,
                                   "--plot-area", output_file,
                                   ],)

        assert "'apl' command executed" in output

        self.assert_filecmp(output_file, data.VESICLE_AREA_XVG)

        os.unlink(output_file)

    def test_area_protein(self):
        output_file = os.path.join(self.tempdir, "test_area.xvg")
        output = self.run_command(["apl",
                                   "-c", data.BILAYER_PROT_GRO,
                                   "-n", data.BILAYER_PROT_NDX,
                                   "--plot-area", output_file,
                                   ], )

        assert "'apl' command executed" in output

        self.assert_filecmp(output_file, data.BILAYER_PROT_AREA_XVG)

        os.unlink(output_file)

    def test_begin_frame_wrong(self):
        output = self.run_command(["thickness",
                                   "-c", data.VESICLE_GRO,
                                   "-n", data.VESICLE_NDX,
                                   "-t", data.VESICLE_XTC,
                                   "--begin-frame", "2000",
                                   ])

        assert "'thickness' command executed" in output

        assert "WARNING: No frames are selected for analysis" in output

    def test_begin_good(self):
        output = self.run_command(["thickness",
                                   "-c", data.VESICLE_GRO,
                                   "-n", data.VESICLE_NDX,
                                   "-t", data.VESICLE_XTC,
                                   "--begin", "2000",
                                   ])

        assert "'thickness' command executed" in output
        assert "7 processed frames" in output

    def test_end(self):
        output = self.run_command(["thickness",
                                   "-c", data.VESICLE_GRO,
                                   "-n", data.VESICLE_NDX,
                                   "-t", data.VESICLE_XTC,
                                   "--end", "2000",
                                   ])

        assert "'thickness' command executed" in output
        assert "5 processed frames" in output
