#!/usr/bin/env python
# -*- coding: utf8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

CSV_FILENAME = "bilayer_prot_apl_frame_00000.csv"
GRO_FILENAME = "bilayer_prot.gro"
PNG_FILENAME = "bilayer_prot_apl_frame_00000.png"

# Get Box vectors
last_line = ""
with open(GRO_FILENAME) as fp:
    for line in fp:
        line = line.strip()

        if len(line) == 0:
            continue

        last_line = line

box_x, box_y = [float(val) for val in line.split()[:2]]

# Get values
membrane_property = "Area per lipid"
x_values = []
y_values = []
z_values = []
property_values = []
with open(CSV_FILENAME) as fp:
    for lino, line in enumerate(fp):
        if lino == 0:
            membrane_property = line.split(",")[-1].strip()

        else:
            line = line.strip()

            if len(line) == 0:
                continue

            resid, leaflet, x, y, z, value = line.split(",")

            x_values.append(float(x))
            y_values.append(float(y))
            property_values.append(float(value))

# Building data from plotting
grid_x, grid_y = np.mgrid[0:box_x:50j, 0:box_y:50j]
points = np.stack((np.array(x_values).T, np.array(y_values).T), axis=-1)
values = np.array(property_values)
grid = griddata(points, values, (grid_x, grid_y), method='cubic')

# Plot map
plt.contourf(grid_x, grid_y, grid)
cbar = plt.colorbar()

plt.title(membrane_property)
plt.xlabel("Box X (nm)")
plt.ylabel("Box Y (nm)")

plt.tight_layout()
plt.savefig(PNG_FILENAME)

