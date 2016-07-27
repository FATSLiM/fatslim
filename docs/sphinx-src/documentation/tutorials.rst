.. _tutorials:

********************
Tutorials & Examples
********************

This chapter contains a few information and examples to get familiarized with FATSLiM.

Introduction
============

.. _tuto_test_systems:

Test systems
------------

The following MD systems are used both in the tutorials and as tests to check FATSLiM accuracy (see
:ref:`thickness <thickness_accuracy>` and :ref:`area per lipid <apl_accuracy>` chapters)

.. _tuto_lipid_system:

Lipid-only system
"""""""""""""""""

.. figure:: images/bilayer_chol_nobox.png
    :align: center

This system is originally from `Lukat et al. <http://dx.doi.org/10.1021/ci400172g>`_, in which it is called system M2.
It consists of a coarse-grained (Martini force field) planar bilayer made of
828 DPPC, 540 DLPC and 576 cholesterol molecules. Details about this system can be found in the
original paper and the corresponding files are freely available from
the APL@Voro `website <aplvoro_downloads>`_.

- Input files used in tests and/or tutorials:
    * :download:`bilayer_chol.gro <tutorials/bilayer_chol.gro>`
    * :download:`bilayer_chol.ndx <tutorials/bilayer_chol.ndx>`

.. _aplvoro_downloads: http://www.aplvoro.org/index.php?section=downloads

.. _tuto_protein_system:

Protein system
""""""""""""""

.. figure:: images/bilayer_prot_nowater.png
    :align: center

This is originally from `Krüger and Fischer <http://dx.doi.org/10.1007/s00249-009-0487-0>`_ but is also mentioned in `Lukat et al. <http://dx.doi.org/10.1021/ci400172g>`_ as system M3.
It is made of Vpu pore embedded in a membrane made of 112 DOPC. This pore results from the aggregation of
five VPU1−32 WT peptides and is then not a protein *per se*. Yet, the size and morphology of
this pentameric assembly is comparable to a real trans-membrane protein and the analysis
is quite similar. Details about this system can be found in the original paper and the corresponding
files are freely available from the APL@Voro `website <aplvoro_downloads>`_.

- Input files used in tests and/or tutorials:
    * :download:`bilayer_prot.gro <tutorials/bilayer_prot.gro>`
    * :download:`bilayer_prot.ndx <tutorials/bilayer_prot.ndx>`

.. _tuto_peptide_system:

Peptide system
""""""""""""""

.. figure:: images/bilayer_kalp.png
    :align: center

This system is based on work by `Kandasamy and Larson <http://dx.doi.org/10.1529/biophysj.105.073395>`_ and consists of a single
trans-membrane KALP15 peptide (sequence: Ac-GKK(LA)\ :sub:`4` - LKKA-NH\ :sub:`2`\ ) embedded in
a bilayer made of 126 DPPC.

- Input files used in tests and/or tutorials:
    * :download:`bilayer_peptide.gro <tutorials/bilayer_peptide.gro>`
    * :download:`bilayer_peptide.ndx <tutorials/bilayer_peptide.ndx>`
    * :download:`bilayer_peptide.xtc <tutorials/bilayer_peptide.xtc>`

.. _tuto_model_vesicle:

Model vesicle
"""""""""""""

.. figure:: images/model_vesicle.png
    :align: center

The model vesicle made of 2748 DPPC. It was built *ex nihilo* by positioning and orienting
lipid molecules on two concentric spheres. Because
this system is artificially generated, its properties are known: the membrane thickness is
set to 5 nm and the area per lipid for both leaflets are set to 40 |ang|:sup:`2`
(inner) and 64 |ang|:sup:`2` (outer), respectively.

- Input files used in tests and/or tutorials:
    * :download:`model_vesicle.gro <tutorials/model_vesicle.gro>`
    * :download:`model_vesicle.ndx <tutorials/model_vesicle.ndx>`


.. _tuto_real_vesicle:

Real vesicle
""""""""""""

.. figure:: images/dppc_vesicle.png
    :align: center

This vesicle (3030 DPPC) was obtained from the self-aggregation of MARTINI lipids.

- Input files used in tests and/or tutorials:
    * :download:`dppc_vesicle.gro <tutorials/dppc_vesicle.gro>`
    * :download:`dppc_vesicle.ndx <tutorials/dppc_vesicle.ndx>`
    * :download:`dppc_vesicle.xtc <tutorials/dppc_vesicle.xtc>`

.. _tuto_generate_ndx:

Generating the index file
-------------------------

FATSLiM uses a `GROMACS index file <.ndx>`_ to identify the atoms corresponding to lipid headgroup so
every GROMACS users should be at ease and already know how to use the ``gmx make_ndx`` `utility <make_ndx>`_
provided by GROMACS.

.. _.ndx: http://manual.gromacs.org/current/online/ndx.html
.. _make_ndx: http://manual.gromacs.org/current/programs/gmx-make_ndx.html

In the following examples, the atom selection along with the index file will be provided to make things
as clear as possible. See :ref:`below <tuto_make_ndx>` for a complete and comprehensive example.


.. _tuto_membrane_identification:

Tutorial #1: Simple membrane identification
===========================================

- **Goal:** Identify leaflets from the :ref:`lipid-only system <tuto_lipid_system>` and save them to
  an index file
- **Configuration file:** :download:`bilayer_chol.gro <tutorials/bilayer_chol.gro>`

.. _tuto_make_ndx:

Index file
----------

- **Atoms selected as headgroups**: *PO4* for DUPC and DPPC residues  and *ROH* for CHOL residue.

.. important::

    FATSLiM uses these atoms to identify lipids assuming the selected atoms are lipid headgroups.
    You can select several atoms from the same residue (e.g. *PO4* and *NC3* beads for Martini phospholipids)
    without any problem (FATSLiM will just use the center of geometry as the lipid head group position)
    as long as all the atoms belong to the actual lipid headgroup.

When working with Martini lipids, the most common choice is to choose the phosphate moiety (*PO4* bead)
to describe phospholipid headgroups. For the cholesterol molecule, the alcool moiety (*ROH* bead) is usually used.

Before performing analysis, you must create an index file containing these atoms so FATSLiM can
identify lipids:

1. Launch ``gmx make_ndx`` using :download:`bilayer_chol.gro <tutorials/bilayer_chol.gro>` as input file:

.. code-block:: bash

    gmx make_ndx -f bilayer_chol.gro -o bilayer_chol.ndx

You should see the default groups GROMACS creates:

.. code-block:: bash

    0 System              : 33624 atoms
    1 Other               : 33624 atoms
    2 DPPC                :  9936 atoms
    3 DUPC                :  6480 atoms
    4 W                   : 12600 atoms
    5 CHOL                :  4608 atoms

2. Create the group to store headgroups. As the *PO4* bead and the *ROH* only belong to the phospholipids and the cholesterol,
   respectively, you can simply use the following selection string:

.. code-block:: bash

    a PO4 ROH

This will create a new group:

.. code-block:: bash

    6 PO4_ROH             :  1944 atoms

3. Optionally, you can rename this group to be more explicit:

.. code-block:: bash

    name 6 headgroups

4. Finally, you can delete the other groups which are not needed by FATSLiM (this is also completely optional):

.. code-block:: bash

    del 0-5

5. Quit ``gmx make_ndx`` and there should be an index file named ``bilayer_chol.ndx`` that will be
   used by FATSLiM

Voilà, you are now ready to use FATSLiM!

Analysis
--------

As, previously described, FATSLiM needs at least :ref:`two files <common_input_files>` to perform analysis. Here will be used:

- :download:`bilayer_chol.gro <tutorials/bilayer_chol.gro>` which will provide the atom coordinates and system topology

- bilayer_chol.ndx you created to provide lipid headgroups (a safe version of this file
  is available :download:`here <tutorials/bilayer_chol.ndx>`)

To identify the leaflets and store them into an index file use the following command:

.. code-block:: bash

    fatslim membranes -c bilayer_chol.gro -n bilayer_chol.ndx --output-index bilayer_chol_leaflet.ndx

Alternatively, if you did not rename the group in the previous section, you have to specify the name of the group
to use in the index file:

.. code-block:: bash

    fatslim membranes -c bilayer_chol.gro -n bilayer_chol.ndx --output-index bilayer_chol_leaflet.ndx --hg-group PO4_ROH


In both cases, this will create a file named :download:`bilayer_chol_leaflet_0000.ndx <tutorials/bilayer_chol_leaflet_0000.ndx>`
which can be used for further processing (e.g. extracting leaflets from your trajectory with ``gmx trjconv``)

.. note::

    The frame index (starting from 0) is appended to the filename as one index file is created per frame.
    Hence *bilayer_chol_leaflet_0000.ndx* instead of *bilayer_chol_leaflet.ndx*


Tutorial #2: Membrane thickness over trajectory
===============================================

- **Goal:** Plotting the membrane thickness over a trajectory
- **Configuration file:** :download:`bilayer_peptide.gro <tutorials/bilayer_peptide.gro>`
- **Trajectory file:** :download:`bilayer_peptide.xtc <tutorials/bilayer_peptide.xtc>`
- **Atoms selected as headgroups**: *P8* (Phosphorus atom)
- **Index file:** :download:`bilayer_peptide.ndx <tutorials/bilayer_peptide.ndx>`

.. note::

    Generating the (provided) index file using the atom selection is left as an exercise.

    Check the :ref:`above tutorial <tuto_make_ndx>` if needed.

FATSLiM is able to store results into a `.xvg <http://manual.gromacs.org/current/online/xvg.html>`_ file so you can plot them.
Here it is used to plot the thickness of a DPPC vesicle over a (pretty small |--| this is an example!) MD trajectory.
Run the following command:

.. code-block:: bash

    fatslim thickness -c bilayer_peptide.gro -n bilayer_peptide.ndx -t bilayer_peptide.xtc --plot-thickness bilayer_peptide_thickness.xvg

It should give an output similar to this one:

.. code-block:: bash

    FATSLiM - Fast Analysis Toolbox for Simulations of Lipid Membranes
    version 0.2.0
    Copyright (c) 2013-2016 Sébastien Buchoux

    Running command: 'thickness'... This may take some time, be patient!
    Analysis will be performed using 8 threads.
    Analysing frame    91/   91 (time: 10000 ps)... done in 5 ms (Remaining: 0 s)
    Results:
    Average values over 91 processed frames:
    Thickness: Membrane=3.916±0.061nm - Lower leaflet=3.907±0.060nm - Upper leaflet=3.926±0.066nm
    'bilayer_peptide_thickness.xvg' backed up to 'bilayer_peptide_thickness.xvg'
    Thickness values saved to 'bilayer_peptide_thickness.xvg'

    'thickness' command executed in 574.271 ms (CPU)
    Goodbye!

:download:`bilayer_peptide_thickness.xvg <tutorials/bilayer_peptide_thickness.xvg>` was created as a result.
This is a plain text XY plot that can be plot with a vast variety of tool such as `Grace <http://plasma-gate.weizmann.ac.il/Grace/>`_  or `Matplotlib <http://matplotlib.org/>`_:

.. figure:: images/bilayer_peptide_thickness.png
    :align: center


.. _tuto_apl_lipids:

Tutorial #3: Vesicle membrane area and area per lipid
=====================================================

- **Goal:** Plotting the membrane area and area per lipid of a vesicle over a trajectory
- **Configuration file:** :download:`dppc_vesicle.gro <tutorials/dppc_vesicle.gro>`
- **Trajectory file:** :download:`dppc_vesicle.xtc <tutorials/dppc_vesicle.xtc>`
- **Atoms selected as headgroups**: *PO4* (Phosphate bead)
- **Index file:** :download:`dppc_vesicle.ndx <tutorials/dppc_vesicle.ndx>`

.. note::

    Generating the (provided) index file using the atom selection is left as an exercise.

    Check the :ref:`above tutorial <tuto_make_ndx>` if needed.

Rather than running FATSLiM twice (once for membrane area, once for area per lipid), it is possible to
run the analysis once and store both results at the same time. Simply run the following command:

.. code-block:: bash

    fatslim apl -c dppc_vesicle.gro -n dppc_vesicle.ndx -t dppc_vesicle.xtc --plot-apl dppc_vesicle_apl.xvg --plot-area dppc_vesicle_area.xvg

This should give you an output similar to:

.. code-block:: bash

    FATSLiM - Fast Analysis Toolbox for Simulations of Lipid Membranes
    version 0.2.0
    Copyright (c) 2013-2016 Sébastien Buchoux

    Running command: 'apl'... This may take some time, be patient!
    Analysis will be performed using 8 threads.
    Analysing frame    11/   11 (time: 5000 ps)... done in 34 ms (Remaining: 0 s)
    Results:
    Average values over 11 processed frames:
    Area per lipid: Membrane=0.677±0.002nm^2 - Outer leaflet=0.793±0.002nm^2 - Inner leaflet=0.496±0.003nm^2
    'dppc_vesicle_apl.xvg' backed up to 'dppc_vesicle_apl.xvg'
    Area per lipid values saved to 'dppc_vesicle_apl.xvg'

    Area: Membrane=1026.241±3.145nm^2 - Outer leaflet=1468.016±3.455nm^2 - Inner leaflet=584.465±2.970nm^2
    'dppc_vesicle_area.xvg' backed up to 'dppc_vesicle_area.xvg'
    Area values saved to 'dppc_vesicle_area.xvg'

    'apl' command executed in 503.603 ms (CPU)
    Goodbye!

Two files, :download:`dppc_vesicle_apl.xvg <tutorials/dppc_vesicle_apl.xvg>` and
:download:`dppc_vesicle_area.xvg <tutorials/dppc_vesicle_area.xvg>` was created as a result.

.. _tuto_apl_protein:

Tutorial #4: Extracting raw APL data for further processing
===========================================================

- **Goal:** Map the area per lipid from a membrane with a protein
- **Configuration file:** :download:`bilayer_prot.gro <tutorials/bilayer_prot.gro>`
- **Atoms selected as headgroups**: *P8* (Phosphorus atom)
- **Protein atoms selected**: All but hydrogen atoms (*Protein-H* default group from ``gmx make_ndx``)
- **Index file:** :download:`bilayer_prot.ndx <tutorials/bilayer_prot.ndx>`

Index file
----------

In this tutorial, we need to take into account protein atoms for area per lipid calculation.
As for lipid headgroups, this involves ``gmx make_ndx``:

1. Launch ``gmx make_ndx`` using :download:`bilayer_prot.gro <tutorials/bilayer_prot.gro>` as input file:

.. code-block:: bash

    gmx make_ndx -f bilayer_prot.gro -o bilayer_prot.ndx

You should see the default groups GROMACS creates:

.. code-block:: bash

      0 System              : 26703 atoms
      1 Protein             :  1520 atoms
      2 Protein-H           :  1250 atoms
      3 C-alpha             :   160 atoms
      4 Backbone            :   480 atoms
      5 MainChain           :   645 atoms
      6 MainChain+Cb        :   805 atoms
      7 MainChain+H         :   805 atoms
      8 SideChain           :   715 atoms
      9 SideChain-H         :   605 atoms
     10 Prot-Masses         :  1520 atoms
     11 non-Protein         : 25183 atoms
     12 Other               :  5824 atoms
     13 POPC                :  5824 atoms
     14 Water               : 19359 atoms
     15 SOL                 : 19359 atoms
     16 non-Water           :  7344 atoms

Note that, as a protein is present, the number of groups is greater than :ref:`previously <tuto_make_ndx>`.

2. First, we will handle the lipid headgroups:

.. code-block:: bash

    r POPC & a P8

and renamed the newly created group to the default name so we do not need to specify if with :ref:`--hg-group <hg_group>`:

.. code-block:: bash

    name 17 headgroups

3. Then, we need to select the protein atoms we need to take into account we calculating area per lipid.
   In this tutorial, we will use only the heavy atoms but it is a good exercice to test other possibilities (all atoms, only the backbone, etc)
   to see how the selection affects the area per lipid calculation. In our case, the group of atoms is already defined (group 2)
   so we just need to rename it (once again this is just to avoid specifying :ref:`--interacting-group <parameter_interacting_group>`:

.. code-block:: bash

    name 2 protein

4. Optionally, it is also possible to "clean" our index file by deleting the groups we do not need before quitting ``gmx make_ndx``:

.. code-block:: bash

    del 0-1
    del 1-14

When done, we now have :download:`bilayer_prot.ndx <tutorials/bilayer_prot.ndx>` that contains only the two groups needed by FATSLiM.


Analysis
--------

Because, we use the default names (*headgroups* and *protein*) for the atom groups,
we do not need to specify them for FATSLiM to use the correct atom selection and
the following command will suffise:

.. code-block:: bash

    fatslim apl -c bilayer_prot.gro -n bilayer_prot.ndx --export-apl-raw bilayer_prot_apl.csv

This will create a file name :download:`bilayer_prot_apl_frame_00000.csv <tutorials/bilayer_prot_apl_frame_00000.csv>`.

.. note::

    As previously, the frame index (starting from 0) is appended to the filename as one csv file is created per frame.
    Hence *bilayer_prot_apl_0000.csv* instead of *bilayer_prot_apl.csv*

For illustration, here are the first lines of the file:

.. literalinclude:: tutorials/bilayer_prot_apl_frame_00000.csv
    :lines: 1-10

Further analysis
----------------

As you can see, :download:`bilayer_prot_apl_frame_00000.csv <tutorials/bilayer_prot_apl_frame_00000.csv>` contains
all the information we need to create a map with all the area per lipid values.
To actually create this map, we will use a Python script based on `Matplotlib <http://matplotlib.org/>`_.

Writing this script is not directly related to FATSLim and is thus beyond the scope of this tutorial,
this is why it will not be discussed here. Nonetheless, it is an example of post-processing that can be achieved
using FATSLiM raw results.

You download the script :download:`here <tutorials/show_apl_map.py>` or just take a look at its content:

.. literalinclude:: tutorials/show_apl_map.py
    :language: python
    :linenos:
    :caption: show_apl_map.py

and the generated figure:

.. figure:: images/bilayer_prot_apl.png
    :align: center

with the protein added:

.. figure:: images/bilayer_prot_apl2.png
    :align: center

.. |ang| unicode:: U+212B .. angstrom symbol
.. |--| unicode:: U+2013   .. en dash
