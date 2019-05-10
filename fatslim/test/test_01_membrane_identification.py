# -*- coding: utf-8; Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
# 
# This file is part of FATSLiM --- http://fatslim.github.io/
# 
# Copyright (c) 2013-2018, Sébastien Buchoux
# Copyright (c) 2019, by the FATSLiM development team (see AUTHORS file)
#
# FATSLiM is free software and is released under the GNU Public Licence,
# v3 or any higher version
#
# If you use FATSLiM in publised work, please cite:
#
# S. Buchoux.
# FATSLiM: a fast and robust software to analyze MD simulations of membranes.
# Bioinformatics 33(1) (2017), 133--134, doi:10.1093/bioinformatics/btw563
#

import numpy as np
from numpy.testing import assert_allclose, assert_equal
import pytest

from . import system_model_flat, system_model_vesicle, system_model_bulged, system_model_curved, system_model_bicelle
from . import system_big_deformed
from .data import MODELS_METADATA


def test_aggregate_model_flat(system_model_flat):
    system = system_model_flat

    aggregates = system.aggregates

    expected_aggregate_ids = [
        MODELS_METADATA["flat"]["upper_leaflet_ids"],
        MODELS_METADATA["flat"]["lower_leaflet_ids"]
    ]

    assert len(aggregates) == len(expected_aggregate_ids)

    for i, aggregate in enumerate(aggregates):
        assert_equal(aggregate.indices, expected_aggregate_ids[i], err_msg="Bad lipid ids for aggregate #{}".format(i))


def test_aggregate_model_vesicle(system_model_vesicle):
    system = system_model_vesicle

    aggregates = system.aggregates

    expected_aggregate_ids = [
        MODELS_METADATA["vesicle"]["outer_leaflet_ids"],
        MODELS_METADATA["vesicle"]["inner_leaflet_ids"]
    ]

    assert len(aggregates) == len(expected_aggregate_ids)

    for i, aggregate in enumerate(aggregates):
        assert_equal(aggregate.indices, expected_aggregate_ids[i], err_msg="Bad lipid ids for aggregate #{}".format(i))


def test_aggregate_model_bicelle(system_model_bicelle):
    system = system_model_bicelle

    aggregates = system.aggregates

    assert len(aggregates) == 1
    assert len(aggregates[0]) == len(system)


def test_aggregate_big_deformed(system_big_deformed):
    system = system_big_deformed
    expected_sizes = [12152, 11903, 1]

    aggregates = system.aggregates

    assert len(aggregates) == len(expected_sizes)

    for i, aggregate in enumerate(aggregates):
        assert aggregate.size == expected_sizes[i]


@pytest.mark.skip
def test_temptest(system_model_flat,
                  system_model_vesicle,
                  system_model_bicelle,
                  system_model_curved,
                  system_model_bulged,
                  system_big_deformed):
    old_dir = "/Users/sbuchoux/Hacking/fatslim-old/fatslimlib/test/data"
    from ..coreobjects import LipidSystem
    import MDAnalysis as mda
    systems = {
        #"Flat": system_model_flat,
        #"Vesicle": system_model_vesicle,
        "Bicelle": system_model_bicelle,
        #"Curved": system_model_curved,
        #"Bulged": system_model_bulged,
        #"Big deformed": system_big_deformed
    }

    for name, system in systems.items():

        print("Getting naive aggregates for {}:".format(name))

        current_aggregate_id = -1

        aggregate_ids = np.ones(system.nlipids, dtype=np.int) * -1

        stack = np.zeros(system.nlipids, dtype=np.int)
        stack_size = 0

        current_aggregate_lipid_ids = np.zeros(system.nlipids, dtype=np.int)
        current_aggregate_size = 0


        aggregates = []

        # First step find potential leaflets
        for seed_id in range(system.nlipids):

            if aggregate_ids[seed_id] > -1:  # This lipid already belong to an aggregate, we can skip it
                continue

            # Increment the aggregate id as we are dealing with another aggregate
            current_aggregate_id += 1

            # Add this seed bead to the new aggregate
            aggregate_ids[seed_id] = current_aggregate_id
            current_aggregate_size = 1
            current_aggregate_lipid_ids[0] = seed_id

            # Reset the stack
            stack_size = 1
            stack[0] = seed_id

            # Search for all the lipids that belong to the current aggregate
            # This is done by exhausting the stack that contains all the candidates
            while stack_size > 0:

                # Pop the stack
                ref_nid = stack[stack_size - 1]
                stack_size -= 1

                for nid in system.lipid_neighbours[ref_nid]:

                    if aggregate_ids[nid] > -1:  # This lipid already belong to an aggregate, we can skip it
                        continue

                    # If still here, add bead to current aggregate
                    aggregate_ids[nid] = current_aggregate_id
                    current_aggregate_lipid_ids[current_aggregate_size] = nid
                    current_aggregate_size += 1

                    # Append bead to stack
                    stack[stack_size] = nid
                    stack_size += 1

            # Store the aggregate
            ids = np.sort(current_aggregate_lipid_ids[:current_aggregate_size])
            aggregates.append(ids)

        print("Found {} aggregates:".format(len(aggregates)))
        for aggregate in aggregates:
            print("-> {} lipids".format(len(aggregate)))

        # Allocate buffer
        angles = np.empty(system.nlipids, dtype=np.float32)
        angles_bins = np.empty(19, dtype=np.int)

        # Z axis
        ref_axis = np.zeros(3, dtype=np.float32)
        ref_axis[2] = 1.0

        cleaned_aggregates = []
        actual_leftovers = []

        for i, aggregate in enumerate(aggregates):
            print("Checking aggregate #{}:".format(i+1))

            angles_bins[:] = 0
            aggregate_size = len(aggregate)

            n_angles_used = 0
            for seed_id in aggregate:

                dprod = np.sum(system.lipid_normals[seed_id] * ref_axis)
                if dprod > 1:
                    dprod = 1
                if dprod < -1:
                    dprod = -1

                angle = np.arccos(dprod) / np.pi * 180
                angles[seed_id] = angle
                bin_id = int(angle / 10)

                if angles_bins[bin_id] == 0:
                    n_angles_used += 1
                angles_bins[bin_id] += 1
            #
            #max_prob = 0
            # most_prob_angle = -1
            # best_bin_id = -1
            # for bin_id, count in enumerate(angles_bins):
            #     prob = count / aggregate_size
            #
            #     if bin_id == 18:
            #         angle_range = "angle = 180°"
            #     else:
            #         angle_range = "{}° <= angle < {}°".format(bin_id*10, (bin_id+1) * 10)
            #
            #     if prob > max_prob:
            #         max_prob = prob
            #         most_prob_angle = angle_range
            #         best_bin_id = bin_id

                #print("Probability that {}: {:.3f} ({} angles)".format(angle_range, prob, count))

            #print("Most probable angle: {}".format(most_prob_angle))

            if n_angles_used < 9:
                print("Aggregate is more likely planar")

            else:
                print("Aggregate needs checking")
                print(np.argwhere(angles_bins > angles_bins.mean()))


            # Trying to find edges
            edgers = []
            for seed_id in aggregate:
                min_angle = 180
                max_angle = 0
                for nid in system.lipid_neighbours[seed_id]:
                    angle = angles[nid]

                    if angle > max_angle:
                        max_angle = angle

                    if angle < min_angle:
                        min_angle = angle

                angle_extend = max_angle - min_angle

                # print("Angle for resid {}: {:.1f}°".format(system.lipids[seed_id].resid, angle_extend), end="")

                max_extend = 90

                if angle_extend > max_extend:
                    # print(" -> Edge!")
                    edgers.append(seed_id)
                # else:
                    # print()

            print("{} Edgers: resid {}".format(len(edgers), " ".join([str(system.lipids[val].resid) for val in edgers])))


            aggregate_ids[:] = -1

            stack_size = 0

            current_aggregate_size = 0

            splitted_aggregates = []

            current_aggregate_ref_angle = np.empty(3, dtype=np.float32)

            current_aggregate_id = -1

            for seed_id in aggregate:

                if aggregate_ids[seed_id] > -1:  # This lipid already belong to an aggregate, we can skip it
                    continue

                if seed_id in edgers:
                    continue

                # Increment the aggregate id as we are dealing with another aggregate
                current_aggregate_id += 1

                # Add this seed bead to the new aggregate
                aggregate_ids[seed_id] = current_aggregate_id
                current_aggregate_size = 1
                current_aggregate_lipid_ids[0] = seed_id
                current_aggregate_ref_angle = angles[seed_id]

                # Reset the stack
                stack_size = 1
                stack[0] = seed_id

                # Search for all the lipids that belong to the current aggregate
                # This is done by exhausting the stack that contains all the candidates

                print("Starting aggregate from resid {}".format(system.lipids[seed_id].resid))

                while stack_size > 0:

                    # Pop the stack
                    ref_nid = stack[stack_size - 1]
                    stack_size -= 1

                    #if np.dot(angles[ref_nid], angles[edge_id]) < np.cos(10 / 180 * np.pi):
                        #print("Updating Edge from resid #{} to resid #{}".format(system.lipids[edge_id].resid,
                        #                                                         system.lipids[ref_nid].resid))
                    #    edge_id = ref_nid
                    #    ref_angle = angles[edge_id]

                    for nid in system.lipid_neighbours[ref_nid]:

                        if aggregate_ids[nid] > -1:  # This lipid already belong to an aggregate, we can skip it
                            continue

                        if nid in edgers:
                            continue

                        dprod = np.dot(system.lipid_normals[ref_nid], system.lipid_normals[nid])

                        if dprod <= np.cos(45/180*np.pi):
                            continue

                        # if abs(ref_angle - angles[nid]) > 60:
                        #     dx = system.pbcbox.pbc_dx(system.lipid_positions[ref_nid],
                        #                               system.lipid_positions[nid])
                        #
                        #     d = np.sqrt(np.sum(dx**2))

                            #print(d, system.lipids[edge_id].resid, system.lipids[nid].resid,
                            #      ref_angle, angles[nid], abs(ref_angle - angles[nid]))

                            # print("Potential edger: resid {} (d={:.3f}, ref angle: {:.1f}°, angle: {:.1f}°, delta angle: {:.1f}°)".format(
                            #     system.lipids[nid].resid,
                            #     d,
                            #     current_aggregate_ref_angle,
                            #     angles[nid],
                            #     abs(ref_angle - angles[nid])
                            # ))
                            # if d < 60:
                            #     continue
                        #
                        # print("  -> Checking resid {} (angle: {:.1f}°, last angle: {:.1f}°, delta: {:.1f}°)".format(
                        #     system.lipids[nid].resid,
                        #     angles[nid],
                        #     last_angle,
                        #     np.abs(last_angle - angles[nid])), end="")
                        #
                        # if np.abs(last_angle - angles[nid]) > 30:
                        #     print(" -> REJECTED")
                        #     continue
                        # else:
                        #     print(" -> OK")

                        # If still here, add bead to current aggregate
                        aggregate_ids[nid] = current_aggregate_id
                        current_aggregate_lipid_ids[current_aggregate_size] = nid
                        current_aggregate_size += 1


                        # Append bead to stack
                        stack[stack_size] = nid
                        stack_size += 1


                # Store the aggregate
                ids = np.sort(current_aggregate_lipid_ids[:current_aggregate_size])
                splitted_aggregates.append(ids)

            n_added = 1
            while n_added > 0:
                n_added = 0

                leftovers = []
                for ref_nid in edgers:
                    if aggregate_ids[ref_nid] > -1:
                        continue

                    highest_dprod = -2
                    aggregate_id = -1

                    for nid in system.lipid_neighbours[ref_nid]:
                        dprod = np.dot(system.lipid_normals[ref_nid], system.lipid_normals[nid])
                        if dprod > highest_dprod and aggregate_ids[nid] > -1:
                            highest_dprod = dprod
                            aggregate_id = aggregate_ids[nid]

                        # print("  -> checking edger index {} with neighbor index {} (aggid: {}): dport={:.3f}".format(
                        #     system.lipids[ref_nid].hg_atoms[0].ix,
                        #     system.lipids[nid].hg_atoms[0].ix,
                        #     aggregate_ids[nid],
                        #     dprod
                        # ))

                    # print("Best dprod for edger index {}: {:.3f}\n".format(system.lipids[ref_nid].hg_atoms[0].ix, highest_dprod))

                    if highest_dprod > 0.70710:
                        #print("Adding edger {} to aggregate {}".format(system.lipids[ref_nid].resid, aggregate_id))
                        splitted_aggregates[aggregate_id] = np.concatenate((splitted_aggregates[aggregate_id], [ref_nid]))
                        aggregate_ids[ref_nid] = aggregate_id
                        n_added += 1
                    else:
                        leftovers.append(ref_nid)

            print("Aggregate was split into {} aggregates:".format(len(splitted_aggregates)))
            for aggregate in splitted_aggregates:
                #print("  -> {} lipids: resid {}".format(len(clean_aggregate), " ".join([str(val+1) for val in clean_aggregate])))
                print("  -> {} lipids".format(len(aggregate)))
            print("{} lipids are left over".format(len(leftovers)))

            cleaned_aggregates.extend(splitted_aggregates)
            actual_leftovers.extend(leftovers)

        print("Potential leaflets: ({})".format(len(cleaned_aggregates)))
        for i, aggregate in enumerate(cleaned_aggregates):
            print("-> Aggregate #{}:  {} lipids".format(i+1, len(aggregate)), end="")
            if len(aggregate) < 50:
                print(": index {}".format(" ".join([str(system.lipids[val].hg_atoms[0].ix) for val in sorted(aggregate)])))
            else:
                print("")

        print("{} lipids are left over".format(len(actual_leftovers)), end="")
        if len(actual_leftovers) > 0:
            print(": index {}".format(" ".join([str(system.lipids[val].hg_atoms[0].ix) for val in sorted(actual_leftovers)])))
        else:
            print("")
        print("\n\n\n")

    assert False


def test_membrane_model_flat(system_model_flat):
    system = system_model_flat

    expected_leaflets = [
        MODELS_METADATA["flat"]["upper_leaflet_ids"],
        MODELS_METADATA["flat"]["lower_leaflet_ids"]
    ]

    assert len(system.membranes) == 1

    membrane = system.membranes[0]
    for i, leaflet in enumerate(membrane):
        assert_allclose(leaflet.indices, expected_leaflets[i], err_msg="Bad lipid ids for leaflet #{}".format(i))


def test_membrane_model_vesicle(system_model_vesicle):
    system = system_model_vesicle

    expected_leaflets = [
        MODELS_METADATA["vesicle"]["outer_leaflet_ids"],
        MODELS_METADATA["vesicle"]["inner_leaflet_ids"]
    ]

    assert len(system.membranes) == 1

    membrane = system.membranes[0]
    for i, leaflet in enumerate(membrane):
        assert_allclose(leaflet.indices, expected_leaflets[i], err_msg="Bad lipid ids for leaflet #{}".format(i))


def test_membrane_model_bicelle(system_model_bicelle):
    system = system_model_bicelle

    expected_leaflets = [
        np.array([val for val in MODELS_METADATA["bicelle"]["lower_leaflet_ids"] if val not in MODELS_METADATA["bicelle"]["leaflet_pivot_ids"]]),
        np.array([val for val in MODELS_METADATA["bicelle"]["upper_leaflet_ids"] if val not in MODELS_METADATA["bicelle"]["leaflet_pivot_ids"]]),
    ]

    assert len(system.membranes) == 1

    membrane = system.membranes[0]

    for i, leaflet in enumerate(membrane):
        indices = np.array([val for val in leaflet.indices if val not in MODELS_METADATA["bicelle"]["leaflet_pivot_ids"]])
        assert_allclose(indices, expected_leaflets[i], err_msg="Bad lipid ids for leaflet #{}".format(i))


def test_membrane_big_deformed(system_big_deformed):
    system = system_big_deformed
    expected_sizes = [12152, 11861]

    assert len(system.membranes) == 1

    membrane = system.membranes[0]
    for i, leaflet in enumerate(membrane):
        assert leaflet.size == expected_sizes[i]
