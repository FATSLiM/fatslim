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


def backup_file(filepath):
    import os

    ret_val = ""

    if os.path.exists(filepath):
        num = 1

        while os.path.exists("%s.%02i.old" % (filepath, num)):
            num += 1

        backup_filepath = "%s.%02i.old" % (filepath, num)

        ret_val = "'%s' backed up to '%s'\n" % (filepath, backup_filepath)

        os.rename(filepath, backup_filepath)

    return ret_val
