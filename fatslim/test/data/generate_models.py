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

"""
This module contains some quick and dirty code to generate model system to test FATSLiM.

Do not expect anything polished here...
"""

import numpy as np
import scipy.integrate
import MDAnalysis as mda
import pickle


# Lipid definitions
lipids = {
    'DPPC': mda.Universe("DPPC-em.gro"),
    'DTPC': mda.Universe("DTPC-em.gro")
}

refs = {
    'DPPC': lipids['DPPC'].select_atoms("name PO4").positions[0],
    'DTPC': lipids['DTPC'].select_atoms("name PO4").positions[0],
}


# Useful functions
def create_dummy_universe(composition: dict):
    n_residues = 0
    atom_resindex = []
    atom_names = []
    resids = np.empty(0, dtype=np.int)
    resnames = []
    n_resids = 0

    for name, n in composition.items():
        lipid = lipids[name]

        size = len(lipid.atoms)

        for i in range(n):
            atom_resindex.extend([n_resids, ] * size)
            n_resids += 1

        atom_names.extend(list(lipid.atoms.names) * n)

        n_residues += 1
        resids = np.concatenate([resids, np.arange(n_resids - n, n_resids)+1])
        resnames.extend([name, ] * n)

    u = mda.Universe.empty(len(atom_resindex), n_resids, atom_resindex=atom_resindex, trajectory=True)
    u.add_TopologyAttr("resids", resids)
    u.add_TopologyAttr("resnames", resnames)
    u.add_TopologyAttr("names", atom_names)

    current_index = 0
    for name, n in composition.items():
        positions = lipids[name].atoms.positions
        size = len(positions)

        for i in range(n):
            u.trajectory.ts.positions[current_index: current_index + size] = positions.copy()
            current_index += size
    return u


def move_to_cartesian_position(positions, ref_position, new_position):
    vector = new_position - ref_position
    return positions + vector


def move_to_spherical_position(positions, ref_position, r, theta, phi):
    new_ref_position = np.empty_like(ref_position)
    new_ref_position[0] = r * np.sin(theta) * np.cos(phi)
    new_ref_position[1] = r * np.sin(theta) * np.sin(phi)
    new_ref_position[2] = r * np.cos(theta)

    new_positions = move_to_cartesian_position(positions, ref_position, new_ref_position)
    return new_positions


def rotate(positions, ref_position, axis, angle):
    new_positions = positions - ref_position

    ux = np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0]
    ])

    R = np.cos(angle) * np.identity(3, dtype=np.float) + np.sin(angle) * ux + (1 - np.cos(angle)) * np.outer(axis, axis)

    new_positions = np.dot(new_positions, R)

    new_positions += ref_position

    return new_positions


def get_sunflower_seeds_arragement(nseeds):
    b = np.sqrt(nseeds)
    phi = (np.sqrt(5) + 1) / 2 # golden ratio

    thetas = 2 * np.pi * np.arange(nseeds) / phi ** 2
    rs = np.sqrt(np.arange(nseeds)) / b

    return rs, thetas


def phyllotaxic_surface(x_func, z_func, t, radius=None):
    def find_nearest_idx(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

    golden_ratio = (np.sqrt(5) + 1) / 2
    golden_angle = 2 * np.pi * (1 - 1 / golden_ratio)

    # Calculate total area
    if not radius:
        x = x_func(t)
        z = z_func(t)
    else:
        x = x_func(t, radius)
        z = z_func(t, radius)

    x_prime = np.gradient(x, t)
    z_prime = np.gradient(z, t)
    y = 2 * np.pi * x * np.sqrt(x_prime ** 2 + z_prime ** 2)
    total_area = np.abs(scipy.integrate.trapz(y, t))
    areas = np.abs(scipy.integrate.cumtrapz(y, t, initial=0))

    n = 10 * t.size

    # Interpolate values with more points to have smooth values
    interpol_t = np.interp(np.linspace(0, t.size, n), np.arange(t.size), t)
    interpol_areas = np.interp(interpol_t, t, areas)

    points = np.zeros((t.size, 3), dtype=np.float)
    for i in range(t.size):
        area = total_area * i/t.size
        t_val = interpol_t[find_nearest_idx(interpol_areas, area)]
        theta = i * golden_angle

        if not radius:
            points[i] = [
                x_func(t_val) * np.cos(theta),
                x_func(t_val) * np.sin(theta),
                z_func(t_val)]
        else:
            points[i] = [
                x_func(t_val, radius) * np.cos(theta),
                x_func(t_val, radius) * np.sin(theta),
                z_func(t_val, radius)]

    return points


def planar_surface(nx, ny, apl, z_func=None):
    if not z_func:
        def z_func(x, y):
            return 0.0

    positions = np.empty((nx * ny, 3), dtype=np.float32)

    coord_inc = np.sqrt(apl)

    for i in range(nx):
        x = i * coord_inc
        for j in range(ny):
            y = j * coord_inc
            if i % 2 == 1:
                y += 0.5 * coord_inc

            z = z_func(x, y)

            positions[i * nx + j][0] = x
            positions[i * nx + j][1] = y
            positions[i * nx + j][2] = z

    return positions

# Common parameters
APL = 64
THICKNESS = 35

# Exported data

metadata = {}

#################
# Build bicelle #
#################
bicelle_metadata = {}

nlipids = 1000
box_z = 150
ratio = 0.75

n_dppc = int(np.round(nlipids * ratio))
n_dtpc = nlipids - n_dppc

bicelle_radius = np.sqrt((n_dppc * APL) / np.pi)

rs, phis = get_sunflower_seeds_arragement(int(n_dppc/2) + 1)
rs *= bicelle_radius

bicelle_metadata["size"] = nlipids
bicelle_metadata["n_dppc"] = n_dppc
bicelle_metadata["n_dtpc"] = n_dtpc
bicelle_metadata["radius"] = bicelle_radius
bicelle_metadata["apl"] = APL
bicelle_metadata["thickness"] = THICKNESS

bicelle_u = create_dummy_universe({"DPPC": n_dppc, "DTPC": n_dtpc})

# Build bicelle plane
bicelle_plane_positions = np.empty((n_dppc, 3), dtype=np.float32)
for i in range(n_dppc):
    resname = bicelle_u.atoms.resnames[i]

    lipid = lipids[resname]
    ref_position = refs[resname]
    positions = bicelle_u.residues[i].atoms.positions

    if i >= n_dppc / 2:
        z = -THICKNESS / 2
        positions = rotate(positions, ref_position, np.array([1., 0., 0]), np.pi)
        r = rs[i - n_dppc]
        phi = phis[i - n_dppc]
    else:
        z = THICKNESS / 2
        r = rs[i]
        phi = phis[i]

    bicelle_u.residues[i].atoms.positions = move_to_spherical_position(positions, ref_position, r, np.pi / 2, phi)
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    bicelle_plane_positions[i] = [x, y, z]

    bicelle_u.residues[i].atoms.positions += np.array([0, 0, z])

# Build bicelle edges


def x_func(t):
    return bicelle_radius + THICKNESS * 0.5 * np.sin(t)


def z_func(t):
    return THICKNESS * 0.5 * np.cos(t)


t = np.linspace(0, np.pi, n_dtpc)

edge_positions = phyllotaxic_surface(x_func, z_func, t)

for i in range(n_dtpc):
    resindex = i+n_dppc
    resname = bicelle_u.atoms.residues[resindex].resname

    lipid = lipids[resname]
    ref_position = refs[resname]
    positions = bicelle_u.residues[resindex].atoms.positions

    target_position = edge_positions[i]

    # get rotation around z (i.e. we work in the XY plane)
    position = target_position.copy()
    position[2] = 0  # We don't care about the Z coordinate
    position /= np.linalg.norm(position)
    phi = np.arccos(position[0])
    if position[1] > 0:
        phi *= -1

    # get rotation around y, easy just the distribution use to build the bicelle edge
    theta = -t[i]

    positions = rotate(positions, ref_position, np.array([0., 1., 0.]), theta)
    positions = rotate(positions, ref_position, np.array([0, 0., 1.]), phi)

    bicelle_u.residues[resindex].atoms.positions = move_to_cartesian_position(positions, ref_position, target_position)

bicelle_u.atoms.positions += np.array([2 * bicelle_radius, 2 * bicelle_radius, 0.5 * box_z])
bicelle_u.dimensions = np.array([4 * bicelle_radius, 4 * bicelle_radius, box_z, 90, 90, 90], dtype=np.float32)
bicelle_u.atoms.write("model_bicelle.gro")

bicelle_metadata["dimensions"] = np.array([4 * bicelle_radius, 4 * bicelle_radius, box_z, 90, 90, 90], dtype=np.float32)
positions = np.concatenate((bicelle_plane_positions, edge_positions)) + \
            np.array([2 * bicelle_radius, 2 * bicelle_radius, 0.5 * box_z])
bicelle_metadata["positions"] = positions
bicelle_metadata["upper_leaflet_ids"] = np.argwhere(positions[:, 2] >= 0.5 * box_z).flatten()
bicelle_metadata["lower_leaflet_ids"] = np.argwhere(positions[:, 2] < 0.5 * box_z).flatten()
bicelle_metadata["leaflet_pivot_ids"] = np.argwhere(np.abs(positions[:, 2] - 0.5 * box_z) / box_z < 0.01).flatten()

metadata["bicelle"] = bicelle_metadata

# Build vesicle

outer_radius = 60
inner_radius = outer_radius - THICKNESS

n_lipids_outer = int(np.round(4 * np.pi * outer_radius**2 / APL))
n_lipids_inner = int(np.round(4 * np.pi * inner_radius**2 / (APL * 0.75)))

vesicle_metadata = {
    "size": n_lipids_inner + n_lipids_outer,
    "n_outer": n_lipids_outer,
    "n_inner": n_lipids_inner,
    "outer_radius": outer_radius,
    "inner_radius": inner_radius,
    "outer_apl": APL,
    "inner_apl": APL * 0.75,
    "thickness": THICKNESS
}

vesicle_u = create_dummy_universe({"DPPC": n_lipids_inner + n_lipids_outer})


def x_func(t, radius):
    return radius * np.sin(t)

def z_func(t, radius):
    return radius * np.cos(t)


t_outer = np.linspace(0, np.pi, n_lipids_outer)
sphere_positions_outer = phyllotaxic_surface(x_func, z_func, t_outer, radius=outer_radius)

t_inner = np.linspace(0, np.pi, n_lipids_inner)
sphere_positions_inner = phyllotaxic_surface(x_func, z_func, t_inner, radius=inner_radius)

for resindex in range(n_lipids_outer + n_lipids_inner):
    if resindex < n_lipids_outer:
        i = resindex
        t = t_outer
        sphere_positions = sphere_positions_outer
        theta_offset = 0
    else:
        i = resindex - n_lipids_outer
        t = t_inner
        sphere_positions = sphere_positions_inner
        theta_offset = np.pi

    resname = vesicle_u.atoms.residues[resindex].resname
    lipid = lipids[resname]
    ref_position = refs[resname]
    positions = vesicle_u.residues[resindex].atoms.positions

    target_position = sphere_positions[i]

    # get rotation around z (i.e. we work in the XY plane)
    position = target_position.copy()
    position[2] = 0  # We don't care about the Z coordinate

    if np.linalg.norm(position) < 1e-3:
        phi = 0.
        theta = 0.
    else:
        position /= np.linalg.norm(position)
        if abs(position[0]) < 1e-3:
            phi = 0
        else:
            phi = np.arccos(position[0])
    if position[1] > 0:
        phi *= -1

    # get rotation around y, easy just the distribution use to build the bicelle edge
    theta = -t[i] + theta_offset

    positions = rotate(positions, ref_position, np.array([0., 1., 0.]), theta)
    positions = rotate(positions, ref_position, np.array([0, 0., 1.]), phi)

    vesicle_u.residues[resindex].atoms.positions = move_to_cartesian_position(positions, ref_position, sphere_positions[i])

vesicle_u.atoms.positions += np.array([1.5 * outer_radius, 1.5 * outer_radius, 1.5 * outer_radius])
vesicle_u.dimensions = np.array([3 * outer_radius, 3 * outer_radius, 3 * outer_radius, 90, 90, 90], dtype=np.float32)
vesicle_u.atoms.write("model_vesicle.gro")

vesicle_metadata["dimensions"] = vesicle_u.dimensions
vesicle_metadata["positions"] = np.concatenate((sphere_positions_outer, sphere_positions_inner)) + \
                                np.array([1.5 * outer_radius, 1.5 * outer_radius, 1.5 * outer_radius])
vesicle_metadata["outer_leaflet_ids"] = np.arange(n_lipids_outer)
vesicle_metadata["inner_leaflet_ids"] = np.arange(n_lipids_inner) + n_lipids_outer

metadata["vesicle"] = vesicle_metadata


# Build flat membrane
nx = 8
ny = 8

box_x = nx * np.sqrt(APL)
box_y = ny * np.sqrt(APL)

nlipids = 2 * nx * ny

flat_metadata = {
    "size": nlipids,
    "apl": APL,
    "thickness": THICKNESS
}
flat_u = create_dummy_universe({"DPPC": nlipids})

flat_positions = planar_surface(nx, ny, APL)
for resindex in range(nlipids):
    if resindex < nlipids / 2:
        i = resindex
        theta = 0
        z_offset = THICKNESS / 2
    else:
        i = resindex - nx * ny
        theta = np.pi
        z_offset = - THICKNESS / 2

    resname = flat_u.atoms.residues[resindex].resname
    lipid = lipids[resname]
    ref_position = refs[resname]
    positions = flat_u.residues[resindex].atoms.positions

    positions = rotate(positions, ref_position, np.array([0., 1., 0.]), theta)
    positions[:, 2] += z_offset

    flat_u.residues[resindex].atoms.positions = move_to_cartesian_position(positions, ref_position, flat_positions[i])

flat_u.atoms.positions += np.array([0, 0, box_z / 2])
flat_u.dimensions = np.array([box_x, box_y, box_z, 90, 90, 90], dtype=np.float32)
flat_u.atoms.write("model_flat.gro")

flat_metadata["dimensions"] = flat_u.dimensions
positions = flat_positions + np.array([0, 0, box_z / 2])
upper = positions + np.array([0, 0, THICKNESS/2])
lower = positions - np.array([0, 0, THICKNESS/2])
flat_metadata["positions"] = np.concatenate((upper, lower))
flat_metadata["upper_leaflet_ids"] = np.arange(nx * ny)
flat_metadata["lower_leaflet_ids"] = np.arange(nx * ny) + nx * ny

metadata["flat"] = flat_metadata

# Build curved membrane
nx = 20
ny = 20

box_x = nx * np.sqrt(APL)
box_y = ny * np.sqrt(APL)

nlipids = 2 * nx * ny

curved_metadata = {
    "size": nlipids,
    "apl": APL,
    "thickness": THICKNESS
}

curved_u = create_dummy_universe({"DPPC": nlipids})
current_u = curved_u


def z_func(x, y):
    sigma = 0.2 * min(box_x, box_y)
    x0 = box_x * 0.5
    y0 = box_y * 0.5

    A = box_z / 3

    x_value = (x - x0)**2 / (2 * sigma**2)
    y_value = (y - y0)**2 / (2 * sigma**2)

    return A * np.exp(-(x_value + y_value))

curved_positions = planar_surface(nx, ny, APL, z_func=z_func)
current_positions = curved_positions

for resindex in range(nlipids):
    if resindex < nlipids / 2:
        i = resindex
        theta = 0
        z_offset = THICKNESS / 2
    else:
        i = resindex - nx * ny
        theta = np.pi
        z_offset = - THICKNESS / 2

    resname = current_u.atoms.residues[resindex].resname
    lipid = lipids[resname]
    ref_position = refs[resname]
    positions = current_u.residues[resindex].atoms.positions

    positions = rotate(positions, ref_position, np.array([0., 1., 0.]), theta)
    positions[:, 2] += z_offset

    current_u.residues[resindex].atoms.positions = move_to_cartesian_position(positions, ref_position, current_positions[i])

current_u.atoms.positions += np.array([0, 0, box_z / 2])
current_u.dimensions = np.array([box_x, box_y, box_z, 90, 90, 90], dtype=np.float32)

curved_metadata["dimensions"] = current_u.dimensions
positions = current_positions + np.array([0, 0, box_z / 2])
upper = positions + np.array([0, 0, THICKNESS/2])
lower = positions - np.array([0, 0, THICKNESS/2])
curved_metadata["positions"] = np.concatenate((upper, lower))
curved_metadata["upper_leaflet_ids"] = np.arange(nx * ny)
curved_metadata["lower_leaflet_ids"] = np.arange(nx * ny) + nx * ny
metadata["curved"] = curved_metadata

current_u.atoms.write("model_curved.gro")


# Build curved membrane
nx = 20
ny = 20

box_x = nx * np.sqrt(APL)
box_y = ny * np.sqrt(APL)

bulge_radius = 0.35 * min(box_x, box_y)
bulge_height = 0.75 * THICKNESS
radius = (2 * bulge_radius)**2 / (8 * bulge_height) + bulge_height / 2
bulge_offset = bulge_height - radius

nlipids = 2 * nx * ny

bulged_metadata = {
    "size": nlipids,
    "apl": APL,
    "thickness": THICKNESS,
    "bulge_radius": bulge_radius,
    "bulge_height": bulge_height
}
bulged_u = create_dummy_universe({"DPPC": nlipids})
current_u = bulged_u


def z_func(x, y):
    x_centered = x - box_x / 2
    y_centered = y - box_y / 2

    if np.sqrt(x_centered**2 + y_centered**2) > bulge_radius:
        return 0
    else:
        return np.sqrt(radius**2 - x_centered**2 - y_centered**2) + bulge_offset


bulged_positions = planar_surface(nx, ny, APL, z_func=z_func)
current_positions = bulged_positions

for resindex in range(nlipids):
    if resindex < nlipids / 2:
        i = resindex
        theta = 0
        z_offset = THICKNESS / 2
        z_switch = 1
    else:
        i = resindex - nx * ny
        theta = np.pi
        z_offset = - THICKNESS / 2
        z_switch = -1

    resname = current_u.atoms.residues[resindex].resname
    lipid = lipids[resname]
    ref_position = refs[resname]
    positions = current_u.residues[resindex].atoms.positions

    positions = rotate(positions, ref_position, np.array([0., 1., 0.]), theta)
    positions[:, 2] += z_offset

    target_position = current_positions[i].copy()
    target_position[2] *= z_switch

    current_u.residues[resindex].atoms.positions = move_to_cartesian_position(positions, ref_position, target_position)

current_u.atoms.positions += np.array([0, 0, box_z / 2])
current_u.dimensions = np.array([box_x, box_y, box_z, 90, 90, 90], dtype=np.float32)
current_u.atoms.write("model_bulged.gro")

bulged_metadata["dimensions"] = current_u.dimensions
positions = current_positions
upper = positions + np.array([0, 0, THICKNESS/2])
lower = positions + np.array([0, 0, THICKNESS/2])
lower[:, 2] *= -1
bulged_metadata["positions"] = np.concatenate((upper, lower)) + np.array([0, 0, box_z / 2])
bulged_metadata["upper_leaflet_ids"] = np.arange(nx * ny)
bulged_metadata["lower_leaflet_ids"] = np.arange(nx * ny) + nx * ny

metadata["bulged"] = bulged_metadata


# Export metadata to help writing tests
with open("models_metadata.pickle", "wb") as fp:
    pickle.dump(metadata, fp)
