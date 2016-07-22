run fatslim_pymol.py

fatslim_bilayer

disable all

enable pbcbox
enable positions_lower
enable directions_lower
enable normals_lower
enable normals_upper
zoom all, 5
render_to_file("04_Normals")

disable positions_lower
disable positions_upper
disable normals_lower
disable normals_upper
enable leaflets_hg1
enable leaflets_hg2
render_to_file("05_Leaflets")

disable directions_upper
disable directions_lower
render_to_file("05_Leaflets_positions")
