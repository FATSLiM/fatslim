.. _tutorials:

Tutorials & Usage examples
##########################

This chapter contains a few information and examples to get familiarized with FATSLiM.

Introduction
************

This section presents common things used later in the tutorials.

.. _tuto_test_systems:

Test systems
============

The following MD systems are used both in the tutorials and as tests to check FATSLiM accuracy (see
:ref:`thickness <thickness_accuracy>` and :ref:`area per lipid <apl_accuracy>` chapters)

.. _tuto_lipid_system:

Lipid-only system
"""""""""""""""""

.. figure:: images/bilayer_chol_nobox.png
    :align: center

This system is originally from Lukat et al.\ [1]_, in which it is called system M2.
It consists of a coarse-grained (Martini force field) planar bilayer made of
828 DPPC, 540 DLPC and 576 cholesterol molecules. Details about this system can be found in the
original paper and the corresponding files are freely available from
the APL@Voro `website <aplvoro_downloads>`_.

.. [1] TODO: Add formatted ref

- Related files used in tests or tutorials: `bilayer_chol.gro`_ and `bilayer_chol.ndx`_.

.. _aplvoro_downloads: http://www.aplvoro.org/index.php?section=downloads
.. _bilayer_chol.gro: tutorials/bilayer_chol.gro
.. _bilayer_chol.ndx: tutorials/bilayer_chol.ndx

.. _tuto_protein_system:

Protein system
""""""""""""""

.. figure:: images/bilayer_prot_nowater.png
    :align: center

This is originally from Krüger and Fischer\ [2]_ but is also mentioned in Lukat et al.\ [1]_ as system M3.
It is made of Vpu pore embedded in a membrane made of 112 DOPC. This pore results from the aggregation of
five VPU1−32 WT peptides and is then not a protein *per se*. Yet, the size and morphology of
this pentameric assembly is comparable to a real trans-membrane protein and the analysis
is quite similar. Details about this system can be found in the original paper and the corresponding
files are freely available from the APL@Voro `website <aplvoro_downloads>`_.

.. [2] TODO: Add DOI


.. _tuto_peptide_system:

Peptide system
""""""""""""""

.. figure:: images/bilayer_kalp.png
    :align: center

This system is based on work by Kandasamy and Larson [3]_ and consists of a single
trans-membrane KALP15 peptide (sequence: Ac-GKK(LA)\ :sub:`4` - LKKA-NH\ :sub:`2`\ ) embedded in
a bilayer made of 126 DPPC.

.. [3] Add DOI

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

.. _tuto_real_vesicle:

Real vesicle
""""""""""""

.. figure:: images/dppc_vesicle.png
    :align: center

This vesicle (3030 DPPC) was obtained from the self-aggregation of MARTINI lipids.

.. _tuto_generate_ndx:

Generating the index file
=========================

FATSLiM uses a `GROMACS index file <.ndx>`_ to identify the atoms corresponding to lipid headgroup so
every GROMACS users should be at ease and already know how to use the ``gmx make_ndx`` `utility <make_ndx>`_
providing by GROMACS.

.. _.ndx: http://manual.gromacs.org/current/online/ndx.html
.. _make_ndx: http://manual.gromacs.org/current/programs/gmx-make_ndx.html

In the following examples, the atom selection along with the index file will be provided to make things
as clear as possible.


Analysis examples
*****************

This section presents several examples of analysis one can do with FATSLiM.
All the following tutorial refer to a few example MD systems which are presented first.

.. _tuto_membrane_identification:

Leaflet and membrane identification
===================================

Tutorial #1: Simple membrane identification
"""""""""""""""""""""""""""""""""""""""""""

- **Goal:** Identify leaflets from the :ref:`lipid-only system <tuto_lipid_system>` and save them to
  an index file
- **Configuration file:** `bilayer_chol.gro`_

Index file
~~~~~~~~~~

- **Atoms selected as headgroups**: *PO4* for DUPC and DPPC residues  and *ROH* for CHOL residue.

When working with Martini lipids, the most common choice is to choose the phosphate moiety (*PO4* bead)
to describe phospholipid headgroups. For the cholesterol molecule, the alcool moiety (*ROH* bead) is usually used.
Before performing analysis, one must then create an index file containing these atoms so FATSLiM can use it to identify lipids from their headgroups:

1. Launch ``gmx make_ndx`` using `bilayer_chol.gro`_ as input file:

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



Analysis
~~~~~~~~




.. _tuto_apl:

Membrane area and Area per lipid calculation
============================================

.. |ang| unicode:: U+212B .. angstrom symbol
