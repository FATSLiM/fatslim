run fatslim_pymol.py

fatslim_bilayer

disable all

enable pbcbox
enable positions_lower
enable directions_lower
enable positions_upper
enable directions_upper

enable NS_ref
enable NS_cutoff
enable NS_neighbors
zoom all, 5
render_to_file("03_NS_all")

disable NS_cutoff
render_to_file("03_NS_nocutoff")

disable positions_lower
disable directions_lower
render_to_file("03_NS_no_lower")

disable positions_upper
disable directions_upper
render_to_file("03_NS_only_neighbors")

enable NS_normal
render_to_file("03_NS_only_neighbors_normal")

enable positions_upper
enable directions_upper
enable positions_lower
enable directions_lower
render_to_file("03_NS_normal")
