Thickness calculation
#####################

How FATSLiM estimates membrane thickness
****************************************

.. _algo_thickness:

Step 0. Membrane identification
===============================

Obviously, before estimating membrane thickness (or inter-leaflet distance), the major prerequisite is to identify the actual membrane and its leaflets!
Please refer to the :ref:`corresponding section <algo_membrane_id>` in the :ref:`previous chapter <chapter_leaflet_membrane>` for details.

Once leaflets are clearly identified, membrane thickness is estimated for every single lipid:
successively, each lipid is taken as a reference and the inter-leaflet distance is calculated as described in the following.



Estimating inter-leaflet distance from reference lipid
======================================================

Because fluctuations of lipid positions are very likely to introduce noise when estimating membrane thickness,
FATSLiM does not rely on single lipid coordinates to calculate inter-leaflet distance but rather
uses neighborhood-averaged coordinates to smooth the individual fluctuations out.


Reference lipid position
""""""""""""""""""""""""

As described the position of the reference is modulated according to its neighborhood.
This procedure can be summarized as follows:

.. figure:: images/thickness_same.png
    :align: center

1. A neighbor search is performed to identify lipids (green) surrounding the reference (purple).
For consistency (as well as computational efficiency) purposes, the results from the neighbor search
done previously to :ref:`calculate normals <algo_local_normals>` are used here.

2. When its normal is almost parallel (a hard-coded 10 degree tolerance is used) to the reference normal,
a lipid is selected (yellow) to modulate the reference position.
This step may be unnecessary when dealing with flat membranes but it can become necessary to avoid bias
in the case of membranes with higher curvature.

3. An averaged position (orange) is calculated from the selected neighbors.
This position, with the addition of the local normal (yellow arrow and black line) will be used to estimate inter-leaflet distance.


Average position of the other leaflet
"""""""""""""""""""""""""""""""""""""

Similarly to what was described previously, an average position is also calculated for the other leaflet:

.. figure:: images/thickness_other.png
    :align: center

4. A neighbor search is performed to identify the lipids which (i) are close to the reference position
and (ii) belong the other leaflet. Obviously, the cutoff distance used for this neighbor search needs
to be big enough to reach the other leaflet (default value: 6.0 nm. See :ref:`--thickness-cutoff <parameter_thickness_cutoff>`).

5. If the distance vector between the neighbor and the reference position is almost parallel
(a hard-coded 10 degree tolerance is used) with the reference normal, the corresponding lipid is selected (yellow).

6. An averaged position (green) is calculated from the selected neighbors.


Thickness estimation
""""""""""""""""""""

The inter-leaflet distance is estimated as the projection of the distance vector **dx** between the two averaged positions:

.. figure:: images/thickness_results.png
    :align: center


Associated command and parameters
*********************************

Command
=======

If you want FATSLiM to estimate membrane thickness, use the following command:

.. code-block:: bash

    fatslim thickness

Parameters
==========

In addition to the common :ref:`analytical parameters <analytical_parameters>`,
Some parameters are specific to the ``thickness`` command.

Analytical parameters
"""""""""""""""""""""

.. _parameter_thickness_cutoff:

Cutoff distance for inter-leaflet neighbor search
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **Associated parameter:** ``--thickness-cutoff``

- **Purpose:** This option allows user to specify the cutoff distance (in nm) to be used when
  performing the inter-leaflet neighbor search needed by the thickness calculation algorithm.

- **Default value:** ``6.0``

