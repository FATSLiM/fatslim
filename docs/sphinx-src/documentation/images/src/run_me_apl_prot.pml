run fatslim_pymol.py

fatslim_apl_prot

disable all

enable pbcbox
zoom pbcbox, 10
enable leaflets_hg2
enable leaflets_hg1
enable protein
enable APL_ref
render_to_file("07_APLProt_01-system_start", side=False)

disable all
enable APL_ref
enable APL_lipid_neighbors
enable APL_interacting_neighbors
render_to_file("07_APLProt_01-system_start_NS", side=False)

disable all
enable APL_ref_2D
enable APL_2d_plane
enable APL_2D_lipid_neighbors
zoom APL_2D_lipid_neighbors, 10
enable APL_ZOI
render_to_file("07_APLProt_02-zoi_lipid", side=False, zoom_obj="APL_2D_lipid_neighbors")

enable APL_2D_interacting_neighbors
render_to_file("07_APLProt_03-2D_interacting_neighbors", side=False, zoom_obj="APL_2D_lipid_neighbors")

disable APL_2D_interacting_neighbors
enable APL_interacting_xcm_used
render_to_file("07_APLProt_04-interacting_used", side=False, zoom_obj="APL_2D_lipid_neighbors")

disable APL_interacting_xcm_used
enable APL_intera_XCM_2D
render_to_file("07_APLProt_05-interacting_xcm", side=False, zoom_obj="APL_2D_lipid_neighbors")

enable APL_CLIPPED_ZOI
disable APL_ZOI
render_to_file("07_APLProt_06-clipped_zoi", side=False, zoom_obj="APL_2D_lipid_neighbors")











