run fatslim_pymol.py

fatslim_bilayer

disable all

enable pbcbox
zoom pbcbox, 10
enable leaflets_hg2
render_to_file("07_APL_01-leaflets_upper", side=False)

disable leaflets_hg2
enable leaflets_hg1
render_to_file("07_APL_01-leaflets_lower", side=False)

enable leaflets_hg2
render_to_file("07_APL_02-leaflets_both")

disable leaflets_hg1
disable leaflets_hg2
enable APL_ref
render_to_file("07_APL_03-ref")

enable APL_cutoff
enable APL_lipid_neighbors
render_to_file("07_APL_04-NS")

disable APL_cutoff
render_to_file("07_APL_05-NS_nocutoff")

disable pbcbox
enable APL_lipid_neighbors
disable leaflets_hg1
disable leaflets_hg2
render_to_file("07_APL_06-NS_projection_neighbors", zoom_obj="APL_plane", zoom_distance=30)

disable APL_lipid_neighbors
enable APL_ref
enable APL_plane
enable APL_lipid_projection_lines
enable APL_proj_lipid_neighbors
render_to_file("07_APL_06-NS_projection", zoom_obj="APL_plane", zoom_distance=30)

disable all
enable APL_ref_2D
enable APL_2d_plane
enable APL_2D_lipid_neighbors
zoom APL_2D_lipid_neighbors, 10
render_to_file("07_APL_07-2D_neighbors", side=False, zoom_obj="APL_2D_lipid_neighbors")

enable APL_VORO
render_to_file("07_APL_08-voro_complete", side=False, zoom_obj="APL_2D_lipid_neighbors")

disable APL_VORO
enable APL_ZOI
render_to_file("07_APL_09-zoi", side=False, zoom_obj="APL_2D_lipid_neighbors")








