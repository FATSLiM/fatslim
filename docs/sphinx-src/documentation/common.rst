User interface
##############

Command line interface
**********************

Once FATSLiM is correctly :ref:`installed <installation>` and :ref:`configured <post-installation>`,
it is available from the command line interface by launching the main (and only) executable: ``fatslim``.

Usage
=====

The user interface for FATSLiM is pretty straightforward: simply type ``fatslim`` followed by the command you want to use.
For example, ``fatslim version`` will display detailed information about FATSLiM installation such as versions used for dependencies.

.. _general_opt:


General commands and options
****************************

General options
===============

These options will take precendence over any subsequent command.
For instance ``fatslim --help apl`` will display help about the ``apl`` command and exit, area per lipid calculation will not be done.

.. _help_opt:

--help / -h
"""""""""""

``fatslim --help`` or ``fatslim -h`` without any command will display general help.

.. note::

    ``fatslim`` alone (i.e. without any command specified) will also default to displaying help.


If a command a specified (e.g. ``fatslim -h thickness``), the help message for this command will be displayed.

.. _version_opt:

--version
"""""""""

``fatslim --version`` will display detailed information about FATSLiM installation (Python, Cython or NumPy versions, Number of threads available, etc...).

.. note::

    This is an alias of the :ref:`version_cmd` command.


General commands
================

These commands do not refer to any analysis.

.. _help_cmd:

help
""""

``fatslim help`` is just an alias for :ref:`help_opt`. So ``fatslim help`` will display general help and ``fatslim help command`` will display help for a specific command.

.. _version_cmd:

version
"""""""

``fatslim version`` is just an alias for :ref:`version_opt`.


Analytical commands and related options
***************************************

Analysis-related paremeters / options
=====================================

For consistency purpose, all analytical commands share a common set of parameters and options which will be described here.

.. note::

    Similarly to GROMACS philosophy, you may omit any parameter or option when its value equals the default one.


Input files
"""""""""""

To perform any analysis, FATSLiM needs at least two input files: a :ref:`configuration file <conf_file_input>` and an :ref:`index file <index_file_input>`.
When working with a MD trajectory, you also need to specify a :ref:`trajectory file <traj_file_input>`.

.. _conf_file_input:

Configuration file
~~~~~~~~~~~~~~~~~~

- **Associated parameter:** ``--conf`` or ``-c``

- **Purpose:** This file provides the "topology" of your system. This is used by FATSLiM to identify lipid based on residue indices (atoms with the same residue index will be considered as belonging to the same lipid).
  This file is mandatory.

- **Accepted file extensions:** `.gro`_

- **Default value:** ``conf.gro``

.. _.gro: http://manual.gromacs.org/current/online/gro.html


.. _traj_file_input:

Trajectory file
~~~~~~~~~~~~~~~

- **Associated parameter:** ``--trajectory`` or ``-t``

- **Purpose:** This file provides the atom coordinates. If this file is provided FATSLiM will use atoms from it and not from the :ref:`conf file <conf_file_input>`.
  If omitted, FATSLiM will use coordinates from the :ref:`conf file <conf_file_input>`.

- **Accepted file extensions:** `.gro`_, `.trr`_, `.xtc`_.

- **Default value:** ``traj.trr``

.. _.trr: http://manual.gromacs.org/current/online/trr.html

.. _.xtc: http://manual.gromacs.org/current/online/xtc.html

.. _index_file_input:

Index file
~~~~~~~~~~

- **Associated parameter:** ``--index`` or ``-n``

- **Purpose:** This file provides a `GROMACS index file <http://manual.gromacs.org/current/online/ndx.html>`_ used to identify lipid headgroups and, if present and needed, interacting atoms.
  This file is mandatory. See :ref:`tuto_generate_ndx` for further details.

- **Accepted file extensions:** `.ndx`_.

- **Default value:** ``index.ndx``

.. _.ndx: http://manual.gromacs.org/current/online/ndx.html


Analysis options
""""""""""""""""

.. _hg_group:

Headgroup group
~~~~~~~~~~~~~~~

- **Associated parameter:** ``--hg-group``

- **Purpose:** ``--hg-group`` specifies the name of the group (from the :ref:` index file <index_file_input>`) which defines the atoms used to describe the lipid head groups.

- **Default value:** ``headgroups``


Number of threads
~~~~~~~~~~~~~~~~~

- **Associated parameter:** ``--nthreads``

- **Purpose:** Sets the number of openMP threads to use. a negative value means "use all available CPUs".

- **Default value:** ``-1``


.. _begin_end_timestep_opt:

Specify analysis window using timesteps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **Associated parameters:** ``--begin``/``-b`` and ``--end``/``-e``

- **Purpose:** ``--begin`` and ``--end`` allow users to specify the first and the last timestep (in ps) to be used for analysis.
  If the timestep value does not correspond to one stored in the trajectory, the closest actual timestep will be used. If the specified value is negative, it will be ignored.

- **Default values:** ``-1`` for both parameters

.. versionchanged:: 0.1.3

    Prior to version 0.1.3, ``--begin`` and ``--end`` was used to specify frame indices and not timesteps.

.. seealso::

    :ref:`--begin-frame and --end-frame <begin_end_frame_opt>`


.. _begin_end_frame_opt:

Specify analysis window using timesteps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **Associated parameters:** ``--begin-frame`` and ``--end-frame``

- **Purpose:** ``--begin-frame`` and ``--end-frame`` allow users to specify the first and the last frame index to be used for analysis.
  If the index does not correspond to an actual frame in the trajectory, the closest actual index will be used. If the specified value is negative, it will be ignored.

- **Default values:** ``-1`` for both parameters

.. seealso::

    :ref:`--begin and --end <begin_end_frame_opt>`


.. _cutoff_leaflet_opt:

Cutoff distance for leaflet identification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **Associated parameter:** ``--cutoff``

- **Purpose:** This option allows user to specify the cutoff distance (in nm) to be used when leaflet identification is performed.
  See :ref:`leaflets` for details.

- **Default value:** ``2.0``