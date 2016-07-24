# -*- coding: utf-8 -*-
"""
Created on Sun Mar  2 11:33:04 2014

@author: sebastien
"""
import sys
import os
import numpy as np

try:
    from pymol import cmd
    from pymol import cgo
    from pymol.cgo import CONE, CYLINDER, COLOR, SPHERE, VERTEX, BEGIN, LINEWIDTH, LINES, END, \
        ALPHA
except ImportError:
    print("Not inside PyMOL... Exiting!")
    sys.exit(1)

FATSLIM_DIR = os.path.expanduser("~/Hacking/fatslim")
sys.path.insert(0, FATSLIM_DIR)

try:
    from fatslimlib.datareading import load_trajectory
    from fatslimlib.core_ns import neighbor_search
except ImportError:
    print("Could not find FATSLiM!")
    sys.exit(1)
else:
    import fatslimlib

    print("Fatslim v.%s found here: %s" %
          (fatslimlib.__version__, os.path.dirname(fatslimlib.__file__)))

EPSILON = 1e-6
BILAYER = "bilayer"
VESICLE = "vesicle"
BILAYER_PROT = "bilayer_prot"
# BILAYER_REF = 54 - 1
BILAYER_REF = 336 - 1
VESICLE_REF = 1115 - 1
BILAYER_PROT_REF = 205 - 161
REF_BEAD = 0
NS_RADIUS = 2
THICKNESS_DEFAULT_CUTOFF = 6.0
THICKNESS_MAXIMUM_MINMAX_RATIO = 1.5
THICKNESS_MIN_COS_DX = 0.93969262078590843 # Cosine 20 deg
THICKNESS_MIN_COS_NORMAL = 0.93969262078590843 # Cosine 20 deg
THICKNESS_DEBUG_BEADID = 160
XX = 0
YY = 1
ZZ = 2
COS_45 = 0.70710678
APL_DEFAULT_CUTOFF = 3.0


def setup():
    cmd.delete("all")
    print("Everything removed!")

    cmd.set("bg_rgb", "white")
    cmd.set("ray_opaque_background", "off")
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_trace_gain", 0.5)
    cmd.set("depth_cue", 0)
    cmd.set("orthoscopic", 1)
    print("Display settings OK!")


def get_box(obj):
    dimensions = cmd.get_symmetry(obj)[:6]

    B = np.zeros((3, 3), dtype=np.float32)
    x, y, z, a, b, c = dimensions

    if np.all(dimensions[:3] == 0):
        return B

    B[0][0] = x
    if a == 90. and b == 90. and c == 90.:
        B[1][1] = y
        B[2][2] = z
    else:
        a = np.deg2rad(a)
        b = np.deg2rad(b)
        c = np.deg2rad(c)
        B[1][0] = y * np.cos(c)
        B[1][1] = y * np.sin(c)
        B[2][0] = z * np.cos(b)
        B[2][1] = z * (np.cos(a) - np.cos(b) * np.cos(c)) / np.sin(c)
        B[2][2] = np.sqrt(z * z - B[2][0] ** 2 - B[2][1] ** 2)
    return B


def draw_pbc_box(obj, linewidth=1.0, colorRGB=(0.1, 0.3, 0.15)):
    box = get_box(obj)

    r, g, b = colorRGB
    minX = minY = minZ = 0
    maxX = box[XX][XX]
    maxY = box[YY][YY]
    maxZ = box[ZZ][ZZ]

    boxCGO = [
        LINEWIDTH, float(linewidth),

        BEGIN, LINES,
        COLOR, float(r), float(g), float(b),

        VERTEX, minX, minY, minZ,  # 1
        VERTEX, minX, minY, maxZ,  # 2

        VERTEX, minX, maxY, minZ,  # 3
        VERTEX, minX, maxY, maxZ,  # 4

        VERTEX, maxX, minY, minZ,  # 5
        VERTEX, maxX, minY, maxZ,  # 6

        VERTEX, maxX, maxY, minZ,  # 7
        VERTEX, maxX, maxY, maxZ,  # 8

        VERTEX, minX, minY, minZ,  # 1
        VERTEX, maxX, minY, minZ,  # 5

        VERTEX, minX, maxY, minZ,  # 3
        VERTEX, maxX, maxY, minZ,  # 7

        VERTEX, minX, maxY, maxZ,  # 4
        VERTEX, maxX, maxY, maxZ,  # 8

        VERTEX, minX, minY, maxZ,  # 2
        VERTEX, maxX, minY, maxZ,  # 6


        VERTEX, minX, minY, minZ,  # 1
        VERTEX, minX, maxY, minZ,  # 3

        VERTEX, maxX, minY, minZ,  # 5
        VERTEX, maxX, maxY, minZ,  # 7

        VERTEX, minX, minY, maxZ,  # 2
        VERTEX, minX, maxY, maxZ,  # 4

        VERTEX, maxX, minY, maxZ,  # 6
        VERTEX, maxX, maxY, maxZ,  # 8

        END
    ]

    cmd.load_cgo(boxCGO, "pbcbox")
    cmd.set("cgo_transparency", 0.5, "pbcbox")


cmd.extend("draw_pbc_box", draw_pbc_box)


def draw_vector(x_begin, direction, cgo_obj=None, color=(0.7, 0.7, 0.2), alpha=1.0, factor=1.0,
                length=10.0):
    if cgo_obj is None:
        cgo_obj = []

    head_length = 0.4 * length * factor
    radius = 1.0 * factor
    head_radius = 1.6 * radius

    red, green, blue = color

    tail = ALPHA, alpha, \
        CYLINDER, \
        x_begin[XX], x_begin[YY], x_begin[ZZ], \
        x_begin[XX] + length * direction[XX], x_begin[YY] + length * direction[YY], x_begin[
                                                                                        ZZ] + \
                                                                                    length * \
                                                                                              direction[
                                                                                                  ZZ], \
        radius, red, green, blue, red, green, blue  # Radius and RGB for each cylinder tail
    head = ALPHA, alpha, \
        CONE, \
        x_begin[XX] + length * direction[XX], x_begin[YY] + length * direction[YY], x_begin[
                                                                                        ZZ] + \
                                                                                    length * \
                                                                                              direction[
                                                                                                  ZZ], \
        x_begin[XX] + (length + head_length) * direction[XX], x_begin[YY] + (
                                                                                length +
                                                                                head_length) * \
                                                                            direction[YY], x_begin[
                                                                                               ZZ] + (
                                                                                                     length + head_length) * \
                                                                                                     direction[
                                                                                                         ZZ], \
        head_radius, \
        0.0, \
        red, green, blue, red, green, blue, \
        1.0, 1.0
    cgo_obj.extend(tail)
    cgo_obj.extend(head)

    return cgo_obj


def dprod(a, b):
    val = 0
    for i in range(3):
        val += a[i] * b[i]
    return val


def norm(a):
    val = 0
    for i in range(3):
        val += a[i] * a[i]
    return np.sqrt(val)

def complete_basis(z_axis):
    z_axis = z_axis / norm(z_axis)
    x_axis = np.zeros(3)
    y_axis = np.zeros(3)

    if z_axis[ZZ] < -0.9999999:
        x_axis[XX] = 0.0
        x_axis[YY] = -1.0
        x_axis[ZZ] = 0.0

        y_axis[XX] = -1.0
        y_axis[YY] = 0.0
        y_axis[ZZ] = 0.0

    else:
        a = 1.0  / (1.0 + z_axis[ZZ])
        b = -z_axis[XX] * z_axis[YY] * a

        x_axis[XX] = 1.0 - z_axis[XX] * z_axis[XX] * a
        x_axis[YY] = b
        x_axis[ZZ] = -z_axis[XX]

        y_axis[XX] = b
        y_axis[YY] = 1.0 - z_axis[YY] * z_axis[YY] * a
        y_axis[ZZ] = -z_axis[YY]

    return x_axis, y_axis, z_axis


def to_basis(v, conv_mat):
    v1 = conv_mat[XX][XX] * v[XX] + conv_mat[XX][YY] * v[YY] + conv_mat[XX][ZZ] * v[ZZ]
    v2 = conv_mat[YY][XX] * v[XX] + conv_mat[YY][YY] * v[YY] + conv_mat[YY][ZZ] * v[ZZ]
    v3 = conv_mat[ZZ][XX] * v[XX] + conv_mat[ZZ][YY] * v[YY] + conv_mat[ZZ][ZZ] * v[ZZ]

    return np.array([v1, v2, v3])


def put_atoms_on_plane(ref_position, ref_normal,  neighbors, coords, box):

    # Build plane basis
    plane_x, plane_y, plane_z = complete_basis(ref_normal)

    # Get conversion matrix
    conv_mat = np.linalg.inv(np.stack((plane_x, plane_y, plane_z)))

    coords_2d = np.zeros((len(neighbors), 2))
    for i, nid in enumerate(neighbors):
        # Get the closest image of the neighbor
        dx = box.dx(ref_position,
                    coords[nid])

        # Get the coordinates into the new basis (plane_x, plane_y, plane_z)
        proj = to_basis(dx, conv_mat)

        # Store the projection onto the (plane_x, plane_y) plane
        coords_2d[i][XX] = proj[XX]
        coords_2d[i][YY] = proj[YY]

    return np.array(coords_2d)


def show_positions(frame):
    positions = frame.bead_coords * 10.0
    pivot = positions.mean(axis=0)[2]
    obj_name = "positions"
    cgo_upper = []
    cgo_lower = []
    for position in positions:
        x, y, z = position
        if z < pivot:
            cgo_lower.extend([COLOR, 1.0, 0.6, 0.1, SPHERE, x, y, z, 2.0])
        else:
            cgo_upper.extend([COLOR, 1.0, 0.6, 0.1, SPHERE, x, y, z, 2.0])

    cmd.load_cgo(cgo_lower, "%s_lower" % obj_name)
    cmd.load_cgo(cgo_upper, "%s_upper" % obj_name)
    print("Fatslim positions drawn")


def show_directions(frame):
    directions = -frame.directions
    positions = frame.bead_coords * 10.0
    pivot = positions.mean(axis=0)[2]
    obj_name = "directions"
    cgo_upper = []
    cgo_lower = []
    for i, position in enumerate(positions):
        if position[2] < pivot:
            cgo_lower = draw_vector(position,
                                    directions[i],
                                    cgo_obj=cgo_lower,
                                    color=(1.0, 1.0, 0.22),
                                    alpha=1.0)
        else:
            cgo_upper = draw_vector(position,
                                    directions[i],
                                    cgo_obj=cgo_upper,
                                    color=(1.0, 1.0, 0.22),
                                    alpha=1.0)
    cmd.load_cgo(cgo_lower, "%s_lower" % obj_name)
    cmd.load_cgo(cgo_upper, "%s_upper" % obj_name)
    print("Fatslim directions drawn")


def show_ref(frame):
    cmd.create("Ref_full", "resid %i and resname DMPC" % (REF_BEAD + 1))
    cmd.hide("lines", "Ref_full")
    cmd.show("spheres", "Ref_full")

    cmd.create("Ref_Headgroup", "resid %i and name P" % (REF_BEAD + 1))
    cmd.hide("lines", "Ref_Headgroup")
    cmd.show("spheres", "Ref_Headgroup")
    cmd.set("sphere_scale", 2.0, "Ref_Headgroup")
    cmd.color("orange", "Ref_Headgroup")
    cmd.set("sphere_transparency", 0.5, "Ref_Headgroup")

    position = frame.bead_coords[REF_BEAD] * 10
    x, y, z = position
    cgo_ref = [COLOR, 1.0, 0.6, 0.1, SPHERE, x, y, z, 3.0]
    cgo_ref = draw_vector(position,
                          -frame.directions[REF_BEAD],
                          cgo_obj=cgo_ref,
                          color=(1.0, 1.0, 0.22),
                          alpha=1.0)
    cmd.load_cgo(cgo_ref, "Ref_simplified")

    position = frame.bead_coords[REF_BEAD] * 10
    x, y, z = position
    cgo_ref = [COLOR, 1.0, 0.6, 0.1, SPHERE, x, y, z, 3.0]
    cgo_ref = draw_vector(position,
                          -frame.directions[REF_BEAD],
                          cgo_obj=cgo_ref,
                          color=(1.0, 1.0, 0.22),
                          alpha=1.0)
    cmd.load_cgo(cgo_ref, "Ref_simplified")


def show_ref_ns(frame):
    beadid = REF_BEAD
    position = frame.bead_coords[beadid] * 10
    normal = frame.normals[beadid]
    x, y, z = position

    cmd.load_cgo([COLOR, 0.80, 0.22, 0.22, SPHERE, x, y, z, 3.0], "NS_ref")

    neighbors = neighbor_search(frame.box, frame.bead_coords, cutoff=NS_RADIUS)
    cmd.load_cgo([COLOR, 0.8, 1.0, 0.8, SPHERE, x, y, z, NS_RADIUS * 10.0], "NS_cutoff")
    cmd.set("cgo_transparency", 0.6, "NS_cutoff")

    cgo_neighbors = []
    for nid in neighbors[beadid]:
        x, y, z = frame.bead_coords[nid] * 10
        cgo_neighbors.extend([COLOR, 0.18, 0.53, 0.18, SPHERE, x, y, z, 2.1])
    cmd.load_cgo(cgo_neighbors, "NS_neighbors")

    cgo_pca = []
    cgo_pca = draw_vector(position,
                          normal,
                          cgo_obj=cgo_pca,
                          color=(0.33, 0.67, 0.33),
                          alpha=1.0,
                          factor=1.5,
                          length=20)
    cmd.load_cgo(cgo_pca, "NS_normal")

    print("Fatslim reference bead drawn")


def show_normals(frame):
    normals = frame.normals
    positions = frame.bead_coords * 10.0
    pivot = positions.mean(axis=0)[2]
    obj_name = "normals"
    cgo_upper = []
    cgo_lower = []
    for i, position in enumerate(positions):
        if position[2] < pivot:
            cgo_lower = draw_vector(position,
                                    normals[i],
                                    cgo_obj=cgo_lower,
                                    color=(0.33, 0.67, 0.33),
                                    alpha=1.0,
                                    factor=0.9)
        else:
            cgo_upper = draw_vector(position,
                                    normals[i],
                                    cgo_obj=cgo_upper,
                                    color=(0.33, 0.67, 0.33),
                                    alpha=1.0,
                                    factor=0.9)
    cmd.load_cgo(cgo_lower, "%s_lower" % obj_name)
    cmd.load_cgo(cgo_upper, "%s_upper" % obj_name)
    print("Fatslim normals drawn")


def show_leaflets(frame, z_limit=1e9):
    membrane = frame.get_membranes()[0]
    colors = [(1.0, 0.2, 0.2), (0.2, 0.2, 1.0)]
    for i, leaflet in enumerate(membrane):

        obj_name = "leaflets_hg%i" % (i + 1)
        cgo_positions = []
        positions = leaflet.coords * 10.0
        for position in positions:
            x, y, z = position
            red, green, blue = colors[i % 2]
            if z < z_limit:
                cgo_positions.extend([COLOR, red, green, blue, SPHERE, x, y, z, 2.0])

        cmd.load_cgo(cgo_positions, obj_name)
    print("Fatslim leaflets drawn")


def show_leaflet_normals(frame, z_limit=1e9):
    membrane = frame.get_membranes()[0]
    colors = [(1.0, 0.6, 0.6), (0.6, 0.6, 1.0)]

    for i, leaflet in enumerate(membrane):
        obj_name = "leaflets_normals%i" % (i + 1)
        cgo_normals = []

        coords = leaflet.coords * 10.0

        for j in range(len(leaflet.normals)):
            if coords[j][2] < z_limit:
                cgo_normals = draw_vector(coords[j],
                                          leaflet.normals[j],
                                          cgo_obj=cgo_normals,
                                          color=colors[i % 2])

        if len(cgo_normals) > 0:
            cmd.load_cgo(cgo_normals, obj_name)
    print("Fatslim leaflet normals drawn")


def ClosestPointOnLine(a, b, p):
    ap = p - a
    ab = b - a
    result = a + np.dot(ap, ab) / np.dot(ab, ab) * ab
    return result


def show_thickness(frame):

    beadid = REF_BEAD
    ref_position = frame.bead_coords[beadid]
    x, y, z = ref_position * 10

    membrane = frame.get_membranes()[0]

    cmd.load_cgo([COLOR, 1.0, 0.2, 1.0, SPHERE, x, y, z, 2.5],
                 "THICKNESS_ref")

    cmd.load_cgo([COLOR, 0.8, 1.0, 0.8, SPHERE, x, y, z, THICKNESS_DEFAULT_CUTOFF * 10.0],
                 "THICKNESS_cutoff")
    cmd.set("cgo_transparency", 0.6, "THICKNESS_cutoff")

    if beadid in membrane.leaflet1.beadids:
        other_leaflet = membrane.leaflet2
        same_leaflet = membrane.leaflet1
    else:
        other_leaflet = membrane.leaflet1
        same_leaflet = membrane.leaflet2

    ref_leaflet_beadid = list(same_leaflet.beadids).index(beadid)
    ref_normal = same_leaflet.normals[ref_leaflet_beadid]

    cgo_ref_normal = draw_vector(ref_position * 10,
                                 ref_normal,
                                 cgo_obj=[],
                                 color=(1.0, 1.0, 0.22),
                                 alpha=1.0)

    cmd.load_cgo(cgo_ref_normal, "THICKNESS_ref_normal")

    # Step 1: get XCM
    same_neighbors = neighbor_search(frame.box,
                                     np.array([ref_position, ]),
                                     neighbor_coords=same_leaflet.coords,
                                     cutoff=NS_RADIUS
                                     )

    total_weight = 0.0
    ref_xcm = np.zeros(3)
    cgo_neighbors = []
    cgo_useful = []
    for nid in same_neighbors[0]:
        dprod_normal = dprod(same_leaflet.normals[nid], ref_normal)

        if dprod_normal < THICKNESS_MIN_COS_NORMAL:
            continue

        dx = frame.box.dx(ref_position,
                          same_leaflet.coords[nid])

        x, y, z = same_leaflet.coords[nid] * 10
        cgo_neighbors.extend([COLOR, 0.18, 0.53, 0.18, SPHERE, x, y, z, 2.1])

        dx_norm = norm(dx)

        if dx_norm < EPSILON:
            continue

        weight = 1 - dx_norm / NS_RADIUS
        weight = 1
        total_weight += weight

        dx *= weight

        ref_xcm += dx
        cgo_useful.extend([COLOR, 0.90, 0.80, 0.18, SPHERE, x, y, z, 2.5])

    if total_weight > EPSILON:
        ref_xcm /= total_weight

    ref_xcm += ref_position

    x, y, z = ref_xcm * 10.0
    cmd.load_cgo([COLOR, 0.0, 0.3, 1.0, SPHERE, x, y, z, 3.5],
                 "THICKNESS_REF_XCM")
    cmd.load_cgo(cgo_neighbors, "THICKNESS_REF_neighbors")
    cmd.load_cgo(cgo_useful, "THICKNESS_REF_used")


    neighbors = neighbor_search(frame.box,
                                np.array([ref_position, ]),
                                neighbor_coords=other_leaflet.coords,
                                cutoff=THICKNESS_DEFAULT_CUTOFF)
    cgo_directions = []
    cgo_normals = []
    cgo_neighbors = []
    cgo_useful = []
    normals = other_leaflet.normals

    proj_lines = []
    avg_dx = np.zeros(3)
    total_weight = 0
    for nid in neighbors[0]:
        other_coord = other_leaflet.coords[nid]
        other_normal = other_leaflet.normals[nid]

        dprod_normal = dprod(other_normal, ref_normal)

        x, y, z = other_coord * 10
        cgo_neighbors.extend([COLOR, 0.18, 0.53, 0.18, SPHERE, x, y, z, 2.1])

        cgo_normals = draw_vector(other_coord * 10,
                                  other_normal,
                                  cgo_obj=cgo_normals,
                                  color=(1.0, 1.0, 0.22),
                                  alpha=1.0)

        # Make sure that both beads belong to different leaflet
        if dprod_normal > -COS_45:
            continue

        # Get the distance (through the bilayer) between the ref bead and its twin
        dx = frame.box.dx_leaflet(ref_xcm, other_coord, ref_normal)
        dx_norm = norm(dx)

        # we check dx because the image which is on the right side of the bilayer may not be a good twin (too far)
        if dx_norm > THICKNESS_DEFAULT_CUTOFF:
            continue

        dprod_dx = np.abs(dprod(dx, ref_normal))
        cos_trial = dprod_dx / dx_norm

        if cos_trial < THICKNESS_MIN_COS_DX:
            continue

        weight = (np.abs(dprod_normal) - THICKNESS_MIN_COS_NORMAL) / \
                 (1.0 - THICKNESS_MIN_COS_NORMAL)
        # weight = 1
        if weight > 0:
            dx *= weight
            total_weight += weight

            avg_dx += dx

            cgo_useful.extend([COLOR, 0.90, 0.80, 0.18, SPHERE, x, y, z, 2.5])


        dx = frame.box.dx_leaflet(ref_position,
                                  other_leaflet.coords[nid],
                                  ref_normal)

        coords = ref_position + dx

        cgo_normals = draw_vector(coords * 10,
                                  normals[nid],
                                  cgo_obj=cgo_normals,
                                  color=(1.0, 1.0, 0.22),
                                  alpha=1.0)
    if total_weight > EPSILON:
        avg_dx /= total_weight

    other_xcm = ref_xcm + avg_dx
    x, y, z = other_xcm * 10
    cmd.load_cgo([COLOR, 0.0, 1.0, 0.2, SPHERE, x, y, z, 3.5],
                 "THICKNESS_OTHER_XCM")

    cmd.load_cgo(cgo_neighbors, "THICKNESS_neighbors")
    cmd.load_cgo(cgo_useful, "THICKNESS_used")
    cmd.load_cgo(cgo_normals, "THICKNESS_normals")

    x, y, z = ref_xcm * 10
    dx, dy, dz = 100 * ref_normal

    ref_line = [
        LINEWIDTH, 1.0,

        BEGIN, LINES,
        COLOR, 0.5, 0.5, 0.5,

        VERTEX, x - dx, y - dy, z - dz,  # 1
        VERTEX, x + dx, y + dy, z + dz,  # 2
        END
    ]
    ref_line = draw_vector(ref_xcm * 10,
                ref_normal,
                cgo_obj=ref_line,
                color=(1.0, 1.0, 0.22),
                alpha=1.0)
    cmd.load_cgo(ref_line, "THICKNESS_normal_xcm")

    print("Fatslim thickness OK")


def show_apl(frame):
    from fatslimlib.core_geometry import get_zoi, is_inside_polygon, clip_zoi
    from scipy.spatial import Voronoi

    beadid = REF_BEAD
    ref_position = frame.bead_coords[beadid]
    x, y, z = ref_position * 10

    membrane = frame.get_membranes()[0]

    cmd.load_cgo([COLOR, 0.8, 0.2, 1.0, SPHERE, x, y, z, 2.5],
                 "APL_ref")
    cmd.load_cgo([COLOR, 0.8, 0.2, 1.0, SPHERE, 0, 0, 0, 2.5],
                 "APL_ref_2D")

    cmd.load_cgo([COLOR, 0.8, 1.0, 0.8, SPHERE, x, y, z, APL_DEFAULT_CUTOFF * 10.0],
                 "APL_cutoff")

    cmd.set("cgo_transparency", 0.6, "APL_cutoff")

    if beadid in membrane.leaflet1.beadids:
        same_leaflet = membrane.leaflet1
    else:
        same_leaflet = membrane.leaflet2

    ref_leaflet_beadid = list(same_leaflet.beadids).index(beadid)
    ref_normal = same_leaflet.normals[ref_leaflet_beadid]

    cgo_ref_normal = draw_vector(ref_position * 10,
                                 ref_normal,
                                 cgo_obj=[],
                                 color=(1.0, 1.0, 0.22),
                                 alpha=1.0)

    cmd.load_cgo(cgo_ref_normal, "APL_ref_normal")

    lipid_neighbors = neighbor_search(frame.box,
                                      np.array([ref_position, ]),
                                      neighbor_coords=same_leaflet.coords,
                                      cutoff=APL_DEFAULT_CUTOFF
                                      )
    lipid_coords_2d = put_atoms_on_plane(ref_position,
                                         ref_normal,
                                         lipid_neighbors[0],
                                         same_leaflet.coords,
                                         frame.box)

    cgo_2d_plane_axes = []
    cgo_plane_axes = []
    for axis in complete_basis(ref_normal):
        cgo_plane_axes = draw_vector(ref_position * 10,
                                 axis,
                                 cgo_obj=cgo_plane_axes,
                                 color=(0.22, 1.0, 0.5),
                                 alpha=1.0,
                                        factor=0.5)
    for i in range(2):
        axis = np.zeros(3)
        axis[i] = 1.0
        cgo_2d_plane_axes = draw_vector([0.0, 0.0, 0.0],
                                        axis,
                                        cgo_obj=cgo_2d_plane_axes,
                                        color=(0.22, 1.0, 0.5),
                                        alpha=1.0,
                                        factor=0.5)

    cmd.load_cgo(cgo_2d_plane_axes, "APL_2d_plane")
    cmd.load_cgo(cgo_plane_axes, "APL_plane")

    cgo_lipid_neighbors = []
    cgo_lipid_neighbors_2d = []
    cgo_lipid_projection_lines = []
    cgo_lipid_neighbors_proj = []
    for i, nid in enumerate(lipid_neighbors[0]):
        cur_coords = same_leaflet.coords[nid]
        x, y, z = cur_coords * 10
        cgo_lipid_neighbors.extend([COLOR, 0.18, 0.53, 0.18, SPHERE, x, y, z, 2.1])

        x2d, y2d = lipid_coords_2d[i] * 10
        cgo_lipid_neighbors_2d.extend([COLOR, 0.4, 0.6, 1.0, SPHERE, x2d, y2d, 0.0, 1.9])

        dx_vector = frame.box.dx(ref_position, cur_coords)

        proj = cur_coords - np.dot(dx_vector, ref_normal) * ref_normal
        proj *= 10

        cgo_lipid_projection_lines.extend([
            LINEWIDTH, 1.0,

            BEGIN, LINES,
            COLOR, 0.5, 0.5, 0.5,

            VERTEX, x, y, z,  # 1
            VERTEX, proj[XX], proj[YY], proj[ZZ],  # 2
            END
        ])
        cgo_lipid_neighbors_proj.extend([COLOR, 0.4, 0.6, 1.0, SPHERE, proj[XX], proj[YY], proj[ZZ], 1.9])

    cmd.load_cgo(cgo_lipid_neighbors, "APL_lipid_neighbors")
    cmd.load_cgo(cgo_lipid_neighbors_2d, "APL_2D_lipid_neighbors")
    cmd.load_cgo(cgo_lipid_neighbors_proj, "APL_proj_lipid_neighbors")
    cmd.load_cgo(cgo_lipid_projection_lines, "APL_lipid_projection_lines")

    zoi = get_zoi(np.zeros(2), lipid_coords_2d) * 10
    if len(zoi) == 0:
        print("APL: Error no ZOI!")
    else:
        print("APL: ZOI polygon has %i vertices" % len(zoi))

    cgo_zoi = []
    for i in range(len(zoi)):
        start = zoi[i]
        end = zoi[(i + 1) % len(zoi)]

        cgo_zoi.extend([
            LINEWIDTH, 1.0,

            BEGIN, LINES,
            COLOR, 0.5, 0.5, 0.5,

            VERTEX, start[XX], start[YY], 0,  # 1
            VERTEX, end[XX], end[YY], 0,  # 2
            END
        ])
    cmd.load_cgo(cgo_zoi, "APL_ZOI")

    voro_pts = lipid_coords_2d.tolist()
    voro_pts.append([0.0, 0.0])
    voro_pts = np.array(voro_pts)
    voro_diagram = Voronoi(voro_pts)
    vor_vertices = voro_diagram.vertices * 10

    cgo_voro = []
    for region in voro_diagram.regions:
        if len(region) == 0 or -1 in region:
            continue

        for i in range(len(region)):
            start = vor_vertices[region[i]]
            end = vor_vertices[region[(i + 1) % len(region)]]

            cgo_voro.extend([
                LINEWIDTH, 1.0,

                BEGIN, LINES,
                COLOR, 0.5, 0.5, 0.5,

                VERTEX, start[XX], start[YY], 0,  # 1
                VERTEX, end[XX], end[YY], 0,  # 2
                END
            ])
    cmd.load_cgo(cgo_voro, "APL_VORO")

    if len(frame.interacting_group_coords) != 0:
        print("APL: %i interacting coords to take into account. " %
              len(frame.interacting_group_coords))

        interacting_neighbors = neighbor_search(frame.box,
                                          np.array([ref_position, ]),
                                          neighbor_coords=frame.interacting_group_coords,
                                          cutoff=APL_DEFAULT_CUTOFF
                                          )
        interacting_coords_2d = put_atoms_on_plane(ref_position,
                                             ref_normal,
                                                   interacting_neighbors[0],
                                                   frame.interacting_group_coords,
                                             frame.box)

        cgo_interacting_neighbors = []
        cgo_interacting_neighbors_2d = []
        cgo_interacting_projection_lines = []
        cgo_interacting_neighbors_proj = []
        cgo_interacting_xcm_used = []

        total_weight = 0
        xcm_interacting = np.zeros(2)

        zoi /= 10
        for i, nid in enumerate(interacting_neighbors[0]):
            cur_coords = frame.interacting_group_coords[nid]
            x, y, z = cur_coords * 10
            cgo_interacting_neighbors.extend([COLOR, 1.0, 1.0, 0.18, SPHERE, x, y, z, 1.0])

            x2d, y2d = interacting_coords_2d[i] * 10
            cgo_interacting_neighbors_2d.extend([COLOR, 1.0, 0.7, 0.2, SPHERE, x2d, y2d, 0.0, 0.9])

            dx_vector = frame.box.dx(ref_position, cur_coords)
            interacting_weight = np.abs(dprod(dx_vector, ref_normal)) / APL_DEFAULT_CUTOFF

            proj = cur_coords - np.dot(dx_vector, ref_normal) * ref_normal
            proj *= 10

            cgo_interacting_projection_lines.extend([
                LINEWIDTH, 1.0,

                BEGIN, LINES,
                COLOR, 0.5, 0.5, 0.5,

                VERTEX, x, y, z,  # 1
                VERTEX, proj[XX], proj[YY], proj[ZZ],  # 2
                END
            ])
            cgo_interacting_neighbors_proj.extend(
                [COLOR, 0.4, 0.6, 1.0, SPHERE, proj[XX], proj[YY], proj[ZZ], 1.9])

            if is_inside_polygon(zoi, interacting_coords_2d[i]):
                total_weight += interacting_weight

                xcm_interacting[XX] += interacting_weight * interacting_coords_2d[i][XX]
                xcm_interacting[YY] += interacting_weight * interacting_coords_2d[i][YY]

                cgo_interacting_xcm_used.extend(
                    [COLOR, 1.0, 0.2, 0.2, SPHERE,
                        interacting_coords_2d[i][XX] * 10,
                        interacting_coords_2d[i][YY] * 10, 0, 1.0])

        if total_weight > 0:
            xcm_interacting[XX] /= total_weight
            xcm_interacting[YY] /= total_weight

            clipped_zoi = clip_zoi(zoi, np.zeros(2), xcm_interacting) * 10

            print("DEBUG: %i" % len(clipped_zoi))

            cmd.load_cgo([COLOR, 0.8, 0.2, 1.0, SPHERE,
                             xcm_interacting[XX] * 10,
                             xcm_interacting[YY] * 10,
                             np.zeros(2), 2.5],
                         "APL_intera_XCM_2D")

            cgo_clipped_zoi = []
            for j, start in enumerate(clipped_zoi):
                end = clipped_zoi[(j + 1) % len(clipped_zoi)]

                cgo_clipped_zoi.extend([
                    LINEWIDTH, 1.0,

                    BEGIN, LINES,
                    COLOR, 0.5, 0.5, 0.5,

                    VERTEX, start[XX], start[YY], 0,  # 1
                    VERTEX, end[XX], end[YY], 0,  # 2
                    END
                ])
            cmd.load_cgo(cgo_clipped_zoi, "APL_CLIPPED_ZOI")

        cmd.load_cgo(cgo_interacting_neighbors, "APL_interacting_neighbors")
        cmd.load_cgo(cgo_interacting_neighbors_2d, "APL_2D_interacting_neighbors")
        cmd.load_cgo(cgo_interacting_neighbors_proj, "APL_proj_interacting_neighbors")
        cmd.load_cgo(cgo_interacting_projection_lines, "APL_interacting_projection_lines")
        cmd.load_cgo(cgo_interacting_xcm_used, "APL_interacting_xcm_used")




    print("FATSLiM APL OK")


def fatslim_bilayer():
    global REF_BEAD
    REF_BEAD = BILAYER_REF
    setup()

    # Load file
    cmd.load("%s.pdb" % BILAYER)
    main_obj = "bilayer"
    cmd.disable(main_obj)
    traj = load_trajectory("%s.gro" % BILAYER, "%s.ndx" % BILAYER)
    traj.initialize()
    frame = traj[0]
    draw_pbc_box(main_obj)
    print("Bilayer Loaded!")

    # Show lipids
    cmd.create("lipids", "resname DMPC")
    cmd.hide("lines", "lipids")
    cmd.show("spheres", "lipids")

    # Show water
    cmd.create("water", "resname SOL")
    cmd.hide("lines", "water")
    cmd.set("solvent_radius", 2)
    cmd.show("surface", "water")
    cmd.color("skyblue", "water")
    cmd.set("transparency", 0.5, "water")
    # cmd.rebuild()
    # cmd.refresh()

    # Show positions
    show_positions(frame)

    # Show directions
    show_directions(frame)

    # Show ref bead and its neighbors
    show_ref(frame)

    show_ref_ns(frame)

    # Show normals
    show_normals(frame)

    # Identify leaflets
    show_leaflets(frame)

    # Calculate and show normals
    show_leaflet_normals(frame)

    # Show stuff related to thickness
    show_thickness(frame)

    # Show stuff related to APL
    show_apl(frame)

    # Zoom on leaflets
    cmd.zoom("all", 5)


cmd.extend("fatslim_bilayer", fatslim_bilayer)


def fatslim_vesicle():
    setup()
    global REF_BEAD
    REF_BEAD = VESICLE_REF

    # Load file
    cmd.load("%s.pdb" % VESICLE)
    main_obj = "vesicle"
    cmd.disable(main_obj)
    traj = load_trajectory("%s.gro" % VESICLE, "%s.ndx" % VESICLE)
    traj.initialize()
    frame = traj[0]
    draw_pbc_box(main_obj)
    print("Vesicle Loaded!")

    # Show positions
    show_positions(frame)

    # Show directions
    show_directions(frame)

    # Show ref bead and its neighbors
    show_ref(frame)

    show_ref_ns(frame)

    # Show normals
    show_normals(frame)

    # Identify leaflets
    show_leaflets(frame)

    # Calculate and show normals
    show_leaflet_normals(frame)

    # Show stuff related to thickness
    show_thickness(frame)

    # Show stuff related to APL
    show_apl(frame)


cmd.extend("fatslim_vesicle", fatslim_vesicle)


def fatslim_apl_prot():
    global REF_BEAD
    setup()

    FILENAME = BILAYER_PROT
    REF_BEAD = BILAYER_PROT_REF

    # Load file
    cmd.load("%s.pdb" % FILENAME)
    main_obj = "%s" % FILENAME
    cmd.disable(main_obj)
    traj = load_trajectory("%s.gro" % FILENAME, "%s.ndx" % FILENAME)
    traj.initialize()
    frame = traj[0]
    draw_pbc_box(main_obj)

    cmd.create("protein", "resi 1-160")
    cmd.hide("lines", "protein")
    cmd.color("yellow", "protein")
    cmd.show("cartoon", "protein")
    cmd.show("surface", "protein")
    cmd.set("transparency", 0.5, "protein")
    cmd.set("surface_color", "yelloworange", "protein")

    # Show stuff related to APL
    show_apl(frame)

    print("Bilayer with protein loaded!")

cmd.extend("fatslim_apl_prot", fatslim_apl_prot)



def render_to_file(fname, top=True, side=True, vesicle=False):
    fname = os.path.splitext(fname)[0]

    cmd.reset()
    if vesicle:
        cmd.zoom("NS_ref", 120)
        cmd.turn("z", -90)
        cmd.move("y", -50)
    if top:
        cmd.ray("2048")
        cmd.png("%s_top.png" % fname)

    if side:
        cmd.turn("x", -90)
        if vesicle:
            cmd.move("y", 150)
        cmd.ray("2048")
        cmd.png("%s_side.png" % fname)
    cmd.reset()
    if vesicle:
        cmd.zoom("NS_ref", 120)
        cmd.turn("z", -90)
        cmd.move("y", -50)


cmd.extend("render_to_file", render_to_file)
