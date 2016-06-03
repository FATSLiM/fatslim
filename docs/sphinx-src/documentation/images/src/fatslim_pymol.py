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
    from pymol.cgo import CONE, CYLINDER, COLOR, SPHERE, VERTEX, BEGIN, LINEWIDTH, LINES, END, ALPHA
except ImportError:
    print ("Not inside PyMOL... Exiting!")
    sys.exit(1)

FATSLIM_DIR =  os.path.expanduser("/home/sebastien/Hacking/fatslim-rewrite")
sys.path += [FATSLIM_DIR, ]

try:
    from fatslimlib.datareading import load_trajectory
    from fatslimlib.core_ns import neighbor_search
except ImportError:
    print ("Could not find FATSLiM!")
    sys.exit(1)





BILAYER = "bilayer"
#BILAYER_REF = 54 - 1
BILAYER_REF = 336 - 1
NS_RADIUS = 2
XX = 0
YY = 1
ZZ = 2

def setup():
    cmd.delete("all")
    print ("Everything removed!")
    
    cmd.set("bg_rgb", "white")
    cmd.set("ray_opaque_background", "off")
    cmd.set("ray_trace_mode",  1)
    cmd.set("ray_trace_gain", 0.5)
    cmd.set("depth_cue", 0)
    cmd.set("orthoscopic", 1)
    print("Display settings OK!")

def get_box(obj):
    dimensions = cmd.get_symmetry(obj)[:6]
    
    B = np.zeros((3,3), dtype=np.float32)
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
        B[1][0] = y*np.cos(c)
        B[1][1] = y*np.sin(c)
        B[2][0] = z*np.cos(b)
        B[2][1] = z*(np.cos(a)-np.cos(b)*np.cos(c))/np.sin(c)
        B[2][2] = np.sqrt(z*z-B[2][0]**2-B[2][1]**2)
    return B
    
def draw_pbc_box(obj, linewidth=1.0, colorRGB=[0.1, 0.3, 0.15]):
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

        VERTEX, minX, minY, minZ, #1
        VERTEX, minX, minY, maxZ, #2

        VERTEX, minX, maxY, minZ, #3
        VERTEX, minX, maxY, maxZ, #4

        VERTEX, maxX, minY, minZ, #5
        VERTEX, maxX, minY, maxZ, #6

        VERTEX, maxX, maxY, minZ, #7
        VERTEX, maxX, maxY, maxZ, #8

        VERTEX, minX, minY, minZ, #1
        VERTEX, maxX, minY, minZ, #5

        VERTEX, minX, maxY, minZ, #3
        VERTEX, maxX, maxY, minZ, #7

        VERTEX, minX, maxY, maxZ, #4
        VERTEX, maxX, maxY, maxZ, #8

        VERTEX, minX, minY, maxZ, #2
        VERTEX, maxX, minY, maxZ, #6


        VERTEX, minX, minY, minZ, #1
        VERTEX, minX, maxY, minZ, #3

        VERTEX, maxX, minY, minZ, #5
        VERTEX, maxX, maxY, minZ, #7

        VERTEX, minX, minY, maxZ, #2
        VERTEX, minX, maxY, maxZ, #4

        VERTEX, maxX, minY, maxZ, #6
        VERTEX, maxX, maxY, maxZ, #8

        END
    ]

    cmd.load_cgo(boxCGO, "pbcbox")
    cmd.set("cgo_transparency",0.5, "pbcbox")


cmd.extend("draw_pbc_box", draw_pbc_box)

def draw_vector(x_begin, direction, cgo_obj=None, color=(0.7,0.7,0.2), alpha=1.0, factor=1.0,
                length=10.0):
    if cgo_obj is None:
        cgo_obj = []
    
    head_length = 0.4 * length * factor
    radius = 1.0 * factor
    head_radius = 1.6 * radius
    
    red, green, blue = color

    tail = ALPHA, alpha,\
           CYLINDER,\
           x_begin[XX], x_begin[YY], x_begin[ZZ],\
           x_begin[XX] + length * direction[XX], x_begin[YY] + length * direction[YY], x_begin[ZZ] + length * direction[ZZ],\
		radius, red,green,blue, red, green, blue # Radius and RGB for each cylinder tail
    head = ALPHA, alpha,\
           CONE,\
           x_begin[XX] + length * direction[XX], x_begin[YY] + length * direction[YY], x_begin[ZZ] + length * direction[ZZ],\
           x_begin[XX] + (length + head_length) * direction[XX], x_begin[YY] + (length + head_length) * direction[YY], x_begin[ZZ] + (length + head_length) * direction[ZZ],\
           head_radius, \
           0.0,\
           red,green,blue, red, green, blue,\
           1.0, 1.0
    cgo_obj.extend(tail)
    cgo_obj.extend(head)
    
    return cgo_obj

def show_positions(frame, z_limit=1e9):
    positions = frame.bead_coords*10.0
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
    positions = frame.bead_coords*10.0
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
    cmd.create("Ref_full", "resid %i and resname DMPC" % (BILAYER_REF + 1))
    cmd.hide("lines", "Ref_full")
    cmd.show("spheres", "Ref_full")
    
    cmd.create("Ref_Headgroup", "resid %i and name P" % (BILAYER_REF + 1))
    cmd.hide("lines", "Ref_Headgroup")
    cmd.show("spheres", "Ref_Headgroup")
    cmd.set("sphere_scale", 2.0, "Ref_Headgroup")
    cmd.color("orange", "Ref_Headgroup")
    cmd.set("sphere_transparency", 0.5, "Ref_Headgroup")
    
    position = frame.bead_coords[BILAYER_REF] * 10
    x, y, z = position
    cgo_ref = [COLOR, 1.0, 0.6, 0.1, SPHERE, x, y, z, 3.0]
    cgo_ref = draw_vector(position,
                          -frame.directions[BILAYER_REF],
                          cgo_obj=cgo_ref,
                          color=(1.0, 1.0, 0.22),
                          alpha=1.0)
    cmd.load_cgo(cgo_ref, "Ref_simplified")
    
    position = frame.bead_coords[BILAYER_REF] * 10
    x, y, z = position
    cgo_ref = [COLOR, 1.0, 0.6, 0.1, SPHERE, x, y, z, 3.0]
    cgo_ref = draw_vector(position,
                          -frame.directions[BILAYER_REF],
                          cgo_obj=cgo_ref,
                          color=(1.0, 1.0, 0.22),
                          alpha=1.0)
    cmd.load_cgo(cgo_ref, "Ref_simplified")
    
def show_ref_ns(frame):
    beadid = BILAYER_REF
    position = frame.bead_coords[beadid] * 10
    normal = frame.normals[beadid]
    x, y, z = position
    
    cmd.load_cgo([COLOR, 0.80, 0.22, 0.22, SPHERE, x, y, z, 3.0], "NS_ref")
    
    
    neighbors = neighbor_search(frame.box, frame.bead_coords, cutoff=NS_RADIUS)
    cmd.load_cgo([COLOR, 0.8, 1.0, 0.8, SPHERE, x, y, z, NS_RADIUS*10.0], "NS_cutoff")
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
    positions = frame.bead_coords*10.0
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
    colors = [(1.0, 0.2, 0.2),(0.2, 0.2, 1.0)]
    for i, leaflet in enumerate(membrane):
        
        obj_name = "leaflets_hg%i" % (i+1)
        cgo_positions = []
        positions = leaflet.coords * 10.0
        for position in positions:
            x, y, z = position
            red, green, blue = colors[i%2]
            if z < z_limit:
                cgo_positions.extend([COLOR, red, green, blue, SPHERE, x, y, z, 2.0])
            
        cmd.load_cgo(cgo_positions, obj_name)
    print("Fatslim leaflets drawn")

def show_leaflet_normals(frame, z_limit=1e9):
    membrane = frame.get_membranes()[0]
    colors = [(1.0, 0.6, 0.6),(0.6, 0.6, 1.0)]
    
    for i, leaflet in enumerate(membrane):
        obj_name = "leaflets_normals%i" % (i+1)
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
    
def fatslim_bilayer():
    setup()
    
    # Load file
    cmd.load("%s.pdb" % BILAYER)
    main_obj = "bilayer"
    cmd.disable(main_obj)
    traj = load_trajectory("%s.gro" % BILAYER, "%s.ndx" % BILAYER)
    traj.initialize()
    frame = traj[0]
    draw_pbc_box(main_obj)
    print ("Bilayer Loaded!")
    
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
    cmd.rebuild()
    cmd.refresh()
    
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
    
    # Zoom on leaflets
    cmd.zoom("all", 5)
    
cmd.extend("fatslim_bilayer", fatslim_bilayer)
    
def render_to_file(fname, top=True, side=True):
    fname = os.path.splitext(fname)[0]
    
    cmd.reset()
    if top:
        cmd.ray("2048")
        cmd.png("%s_top.png" % fname)
    
    if side:
        cmd.turn("x", -90)
        cmd.ray("2048")
        cmd.png("%s_side.png" % fname)
    cmd.reset()
    
    

cmd.extend("render_to_file", render_to_file)