run fatslim_pymol.py

fatslim_bilayer

disable all

enable lipids
#enable water
refresh
#render_to_file("00_raw_system", top=False)


disable lipids
disable water
enable Ref_full
enable Ref_Headgroup
zoom Ref_full, 20
#render_to_file("01_simplification_full", top=False)

enable Ref_simplified
zoom Ref_full, 20
#render_to_file("01_simplification_both", top=False)

disable Ref_full
disable Ref_Headgroup
#render_to_file("01_simplification_simplified", top=False)


disable Ref_Headgroup
disable Ref_simplified
enable pbcbox
enable positions_lower
enable directions_lower
enable positions_upper
enable directions_upper
zoom all, 5
#render_to_file("02_simplified_system", top=False)


enable NS_ref
enable NS_cutoff
enable NS_neighbors
zoom all, 5
#render_to_file("03_NS_all")

disable NS_cutoff
#render_to_file("03_NS_nocutoff")

disable positions_lower
disable directions_lower
#render_to_file("03_NS_no_lower")

disable positions_upper
disable directions_upper
#render_to_file("03_NS_only_neighbors")

enable NS_normal
#render_to_file("03_NS_only_neighbors_normal")

enable positions_upper
enable directions_upper
enable positions_lower
enable directions_lower
#render_to_file("03_NS_normal")


disable NS_ref
disable NS_cutoff
disable NS_neighbors
disable NS_normal
enable positions_lower
enable directions_lower
enable normals_lower
enable normals_upper
zoom all, 5
#render_to_file("04_Normals")


disable positions_lower
disable positions_upper
disable normals_lower
disable normals_upper
enable leaflets_hg1
enable leaflets_hg2
#render_to_file("05_Leaflets")


disable directions_upper
disable directions_lower
enable NS_ref
enable THICKNESS_normal
enable THICKNESS_cutoff
#render_to_file("06_Thickness_all", top=False)

disable THICKNESS_cutoff
#render_to_file("06_Thickness_nocutoff", top=False)

enable THICKNESS_neighbors
#render_to_file("06_Thickness_neighbors", top=False)

disable leaflets_hg1
disable leaflets_hg2
#render_to_file("06_Thickness_neighbors_only", top=False)

enable THICKNESS_used
enable THICKNESS_projections
#render_to_file("06_Thickness_used", top=False)

disable THICKNESS_neighbors
#render_to_file("06_Thickness_used_only")



