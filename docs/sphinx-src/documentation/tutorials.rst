.. _tutorials:

Tutorials
#########

Common stuff
************

.. _tuto_generate_ndx:

Generating the index file
=========================


Analysis examples
*****************

This section presents several examples of analysis one can do with FATSLiM.
All the following tutorial refer to a few example MD systems which are presented first.

.. _tuto_example_systems:

Example systems
===============

.. _tuto_lipid_system:

Lipid-only system
"""""""""""""""""

.. figure:: images/bilayer_chol_nobox.png
    :align: center

This system is originally from Lukat et al.\ [1]_, in which it is called system M2.
It consists of a coarse-grained (Martini force field) planar bilayer made of
828 1,2-DiPalmitoyl-sn-glycero-3-PhosphoCholine (DPPC), 540 1,2-DiLinoleoyl-sn-glycero-
3-PhosphoCholine (DLPC) and 576 cholesterol molecules. All details regarding this system
and the related MD simulations can be found in the original paper and the corresponding
files are freely available from the APL@Voro `website <aplvoro_downloads>`_.

.. [1] TODO: Add formatted ref

.. _aplvoro_downloads: http://www.aplvoro.org/index.php?section=downloads

.. _tuto_protein_system:

Protein system
""""""""""""""

.. figure:: images/bilayer_prot_nowater.png
    :align: center

This is originally from Krüger and Fischer\ [2]_ but is also mentioned in Lukat et al.\ [1]_ as system M3.
It is made of Vpu pore embedded in a membrane made of 112 1,2-DiOleoyl-sn-glycero-3-PhosphoCholine (DOPC)
described using an "united-atom" force field. This pore results from the aggregation of
five VPU1−32 WT peptides and is then not a protein *per se*. Yet, the size and morphology of
this pentameric assembly is comparable to a real trans-membrane protein and the analysis
is quite similar. Details about this system are available in the original paper and MD files
can be retrieved for the APL@Voro `website <aplvoro_downloads>`_.

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
arbitrarily set to 5 nm and the area per lipid for both leaflets are arbitrarily set to 40 |ang|:sup:`2`
(inner) and 64 |ang|:sup:`2` (outer), respectively. Obviously, this vesicle is meant to be a test case for
which physical paramaters are clearly known.



.. _tuto_real_vesicle:

Real vesicle
""""""""""""

.. figure:: images/dppc_vesicle.png
    :align: center

This vesicle (3030 DPPC) was obtained from the self-aggregation of MARTINI lipids.


.. _tuto_membrane_identification:

Leaflet and membrane identification
===================================


.. _tuto_apl:

Membrane area and Area per lipid calculation
============================================

.. |ang| unicode:: U+212B .. angstrom symbol