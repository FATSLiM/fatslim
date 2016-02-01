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
import sys

__authors__ = "Sebastien Buchoux <sebastien.buchoux@gmail.com>"
__copyright__ = "Copyright (c) 2013-2016 %s" % __authors__
__shortname__ = "FATSLiM"
__cmdlinename__ = "fatslim"
__description__ = "Fast Analysis Toolbox for Simulations of Lipid Membranes"
__full_desc__ = """FATSLiM (Fast Analysis Toolbox for Simulations of Lipid Membranes and its goal
 is to provide an efficient, yet robust, tool to extract physical parameters from MD trajectories.

The main objective of FATSLiM is to decrease as possible the amount of time needed to
analyze MD trajectories: the ultimate goal is to process a several-gigabyte-big file is just a
few minutes or less. This is why a rather important part of FATSLiM's development is
focused on code optimization and simplification in order to maximize its efficiency."""
__licence__ = "GNU Public License 3"
__url__ = "https://github.com/seb-buch/FATSLiM"

version_tuple = (0, 1, 0, 'final', 0)


def _format_version_tuple(version_tup=version_tuple):
    if len(version_tup) == 2:
        main_version = "%d.%d" % version_tup
    else:
        if version_tup[2] == 0:
            main_version = "%d.%d" % version_tup[:2]
        else:
            main_version = "%d.%d.%d" % version_tup[:3]
    if len(version_tup) <= 3:
        return main_version

    release_type = version_tup[3]
    try:
        sub = version_tup[4]
    except IndexError:
        sub = 0

    if release_type == "final":
        return main_version
    elif release_type == "alpha":
        release_type = "a"
    elif release_type == "beta":
        release_type = "b"
    elif release_type == "candidate":
        release_type = "rc"
    else:
        release_type = ".%s" % release_type
    sub_string = "%s%i" % (release_type, sub)
    return main_version + sub_string

__version__ = _format_version_tuple(version_tuple)


def print_greetings():
    print(u"%s - %s" % (__shortname__, __description__))
    print(u"version %s" % __version__)
    print(u"%s\n" % __copyright__)


def print_goodbye(exit_code):
    if exit_code == 0:
        print("Goodbye!")
    else:
        print("Sorry, I made a boo boo...")


def print_full_version_info(with_header=True):
    import platform
    import numpy
    from . import core_base

    if with_header:
        print("%s - %s\n" % (__shortname__, __description__))
    print("%s version: %s" % (__shortname__, __version__))
    print("Python version: %s (%s)" % (platform.python_version(), sys.executable))
    print("Cython version (file generation): %s" % core_base.__cython_version__)
    print("Python compiler: %s" % platform.python_compiler())
    print("CPU architecture: %s" % platform.architecture()[0])
    print("OpenMP: %i CPUs (default number of threads: %i - max: %i)" % (core_base.get_num_procs(),
                                                                         core_base.get_num_threads(),
                                                                         core_base.get_max_threads()))
    print("NumPy version: %s" % numpy.__version__)


_command_registry = None


def main(name=__cmdlinename__, argv=None):
    """Main entry point to FATSLiM user interface

    :param name: Name passed to the CommandRegistry
    :param argv: Command line arguments
    """
    from . import command
    from . import builtins
    import inspect

    try:
        print_greetings()

        global _command_registry
        if _command_registry is None:
            _command_registry = command.CommandRegistry(name)

            for name in dir(builtins):
                if not name.startswith("_"):
                    attr = getattr(builtins, name)
                    if inspect.isclass(attr) and \
                            issubclass(attr, command.Command) and \
                            name.startswith("Cmd"):
                            _command_registry.add_command(attr)
        exit_code = _command_registry.run(argv)
    except KeyboardInterrupt:  # pragma: no cover
        print("Interruption requested!")
        exit_code = 0
    except Exception as exc:  # pragma: no cover
        import traceback
        print("\nFATAL ERROR: %s\n" % exc)
        print("Well, this was unexpected... Congrats, you probably found a bug in %s!" %
              __shortname__)
        print("Please report the following traceback to the devs:\n")
        traceback.print_exc()
        print()
        exit_code = 1

    print_goodbye(exit_code)

    return exit_code


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
