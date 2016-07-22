Membrane area and Area per lipid calculation
############################################

How FATSLiM estimates membrane area and area per lipid
******************************************************

.. _algo_apl:

Step 0. Membrane identification
===============================

As no knowledge of the other leaflet is required to perform area calculation on one leaflet, membrane identification is not required *per se*.
Yet, FATSLiM performs such calculation on a fully-identified membrane.
Please refer to the :ref:`corresponding section <algo_membrane_id>` in the :ref:`previous chapter <chapter_leaflet_membrane>` for details.

As for :ref:`thickness <algo_thickness>`, membrane area and area per lipid area estimated for every single lipid:
each lipid is successively taken as a reference and its accessible area (area per lipid) is calculated as described in the following.


Step 1. Estimating reference lipid's accessible area
====================================================

Depending on if a protein (or any molecule) is embedded inside the membrane, the algorithm is slightly modulated but it can be resumed by the following steps:

.. figure:: images/apl_calculation.png
    :align: center


Step 2. Area per lipid and membrane area
========================================

Because the algorithm gives an accessible area for each lipid, area per lipid can obviously be approximated as the average value of the accessible areas.
Additionally, FATSLiM can also group results by lipid type (actually lipid names are used) and give as...

Calculation examples
====================

Here are a few examples of membrane area and area per lipid done with FATSLiM and other software for comparison.

.. seealso::

    Detailed description of these example systems is available :ref:`here <tuto_example_systems>`.

    Check :ref:`tutorial <tutorials>` section to learn how to make these calculation with FATSLiM.

Area per lipid
""""""""""""""

+-----------------+---------------+----------------------------------+--------------------------------------+--------------------------------------+-----------------------------------+---------------------------------+
|                                 | Flat membrane                                                                                                  | Vesicle                                                             |
+                                 +----------------------------------+--------------------------------------+--------------------------------------+-----------------------------------+---------------------------------+
|                                 | :ref:`lipid <tuto_lipid_system>` | :ref:`protein <tuto_protein_system>` | :ref:`peptide <tuto_peptide_system>` | :ref:`model <tuto_model_vesicle>` | :ref:`real <tuto_real_vesicle>` |
+=================+===============+==================================+======================================+======================================+===================================+=================================+
|                 | FATSLiM       | 49.0                             | 68.9                                 | 61.8                                 | 63.9 / 39.8                       | 79.5 / 49.9                     |
+                 +---------------+----------------------------------+--------------------------------------+--------------------------------------+-----------------------------------+---------------------------------+
|                 | `APL@Voro`_   | 48.6                             | 64.6                                 | 61.9                                 | |---|                             | |---|                           |
+                 +---------------+----------------------------------+--------------------------------------+--------------------------------------+-----------------------------------+---------------------------------+
| Area per lipid  | `GridMAT-MD`_ | 48.6                             | 65.1                                 | 61.8                                 | |---|                             | |---|                           |
+ (|ang|:sup:`2`) +---------------+----------------------------------+--------------------------------------+--------------------------------------+-----------------------------------+---------------------------------+
|                 | `MEMBPLUGIN`_ | 46.4                             | |---|                                | |---|                                | |---|                             | |---|                           |
+                 +---------------+----------------------------------+--------------------------------------+--------------------------------------+-----------------------------------+---------------------------------+
|                 | Manually [1]_ | 48.6                             | 67.8                                 | 61.2                                 | 64.0 / 40.0                       | 80.7 / 50.2                     |
+-----------------+---------------+----------------------------------+--------------------------------------+--------------------------------------+-----------------------------------+---------------------------------+

.. _APL@Voro: http://www.aplvoro.org/
.. _GridMAT-MD: http://www.bevanlab.biochem.vt.edu/GridMAT-MD/
.. _MEMBPLUGIN: https://sourceforge.net/projects/membplugin/
.. |---| unicode:: U+2014   .. em dash
.. |ang| unicode:: U+212B .. angstrom symbol

.. note::

    No area per lipid value means that the software is not able to work with such system.

Membrane area
"""""""""""""

+-----------------+---------------+----------------------------------+--------------------------------------+--------------------------------------+-----------------------------------+---------------------------------+
|                                 | Flat membrane                                                                                                  | Vesicle                                                             |
+                                 +----------------------------------+--------------------------------------+--------------------------------------+-----------------------------------+---------------------------------+
|                                 | :ref:`lipid <tuto_lipid_system>` | :ref:`protein <tuto_protein_system>` | :ref:`peptide <tuto_peptide_system>` | :ref:`model <tuto_model_vesicle>` | :ref:`real <tuto_real_vesicle>` |
+=================+===============+==================================+======================================+======================================+===================================+=================================+
|                 | FATSLiM       | 476.5                            | 38.6                                 | 38.9                                 | 1254 / 313                        | 1471 / 588                      |
+                 +---------------+----------------------------------+--------------------------------------+--------------------------------------+-----------------------------------+---------------------------------+
|                 | `APL@Voro`_   | 472.8                            | 36.1                                 | 39.0                                 | |---|                             | |---|                           |
+                 +---------------+----------------------------------+--------------------------------------+--------------------------------------+-----------------------------------+---------------------------------+
| Area            | `GridMAT-MD`_ | 472.7                            | 36.4                                 | 38.9                                 | |---|                             | |---|                           |
+ (nm\ :sup:`2`\ )+---------------+----------------------------------+--------------------------------------+--------------------------------------+-----------------------------------+---------------------------------+
|                 | `MEMBPLUGIN`_ | 453.5                            | |---|                                | |---|                                | |---|                             | |---|                           |
+                 +---------------+----------------------------------+--------------------------------------+--------------------------------------+-----------------------------------+---------------------------------+
|                 | Manually [1]_ | 472.7                            | 38.2                                 | 38.6                                 | 1256 / 314                        | 1493 / 592                      |
+-----------------+---------------+----------------------------------+--------------------------------------+--------------------------------------+-----------------------------------+---------------------------------+

.. note::

    No area value means that the software is not able to work with such system.



Associated command and parameters
*********************************

Command
=======

If you want FATSLiM to estimate membrane area and area per lipid, use the following command:

.. code-block:: bash

    fatslim apl

Parameters
==========

In addition to the common :ref:`analytical parameters <analytical_parameters>`,
Some parameters are specific to the ``apl`` command.

Analytical parameters
"""""""""""""""""""""

.. _parameter_apl_cutoff:

Cutoff distance for area per lipid calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **Associated parameter:** ``--apl-cutoff``

- **Purpose:** This option allows user to specify the cutoff distance (in nm) to be used when
  performing the neighbor search needed by the APL calculation algorithm.

- **Default value:** ``3.0``

.. _parameter_apl_limit:

Upper limit for area per lipid
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **Associated parameter:** ``--apl-limit``

- **Purpose:** This option allows user to specify the upper limit (in nm\ :sup:`2`\ ) for a valid
  area per lipid value.

- **Default value:** ``10.0``

Grouping per lipid type
~~~~~~~~~~~~~~~~~~~~~~~

- **Associated parameter:** ``--apl-by-type``

- **Purpose:** This option allows user to specify that area per lipid values should be grouped by lipid type (i.e. same lipid name).

- **Default value:** ``False``


Output files
""""""""""""

Plotting area per lipid
~~~~~~~~~~~~~~~~~~~~~~~

- **Associated parameter:** ``--plot-apl``

- **Purpose:** This option specifies the filename where FATSLiM should save the area per lipid average values (for membrane and both leaflets) over time (as a XY plot).

- **Accepted file extensions:** `.xvg`_

- **Default value:** None (no output file)

.. _.xvg: http://manual.gromacs.org/current/online/xvg.html

Plotting area per lipid
~~~~~~~~~~~~~~~~~~~~~~~

- **Associated parameter:** ``--plot-area``

- **Purpose:** This option specifies the filename where FATSLiM should save the area average values (for membrane and both leaflets) over time (as a XY plot).

- **Accepted file extensions:** `.xvg`_

- **Default value:** None (no output file)


Raw area per lipid values
~~~~~~~~~~~~~~~~~~~~~~~~~

- **Associated parameter:** ``--export-apl-raw``

- **Purpose:** This option specifies the filename where FATSLiM should save the raw area per lipid (as calculated by the algorithm |--| one value per lipid).
  These values are saved in a `comma separated values <.csv>`_ file.
  To ease further processing the file contains the following columns:

    * residue number (resid)
    * leaflet identifier (e.g. "lower leaflet")
    * lipid coordinates (three columns for x, y and z)
    * area per lipid (in nm\ :sup:`2`\ ).

- **Accepted file extensions:** `.csv`_

- **Default value:** None (no output file)

.. _.csv: https://en.wikipedia.org/wiki/Comma-separated_values

.. |--| unicode:: U+2013   .. en dash

.. [1] See FATSLiM's original paper for details