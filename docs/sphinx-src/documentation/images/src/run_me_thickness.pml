run fatslim_pymol.py

fatslim_bilayer

disable all

enable pbcbox
enable leaflets_hg1
enable leaflets_hg2
zoom all, 5

enable THICKNESS_ref
enable THICKNESS_ref_normal
enable THICKNESS_REF_neighbors
render_to_file("06_Thickness_NS_same")

disable leaflets_hg1
disable leaflets_hg2
render_to_file("06_Thickness_NS_same_neighbors_only")

enable THICKNESS_REF_used
render_to_file("06_Thickness_NS_same_used")

disable THICKNESS_ref
enable THICKNESS_REF_XCM
render_to_file("06_Thickness_NS_same_xcm")

disable THICKNESS_REF_neighbors
disable THICKNESS_REF_used
disable THICKNESS_ref_normal
enable THICKNESS_normal_xcm
render_to_file("06_Thickness_NS_same_xcm_only")

enable THICKNESS_cutoff
render_to_file("06_Thickness_other_cutoff", top=False)

disable THICKNESS_cutoff
enable THICKNESS_neighbors
render_to_file("06_Thickness_other_neighbors", top=False)

disable THICKNESS_neighbors
enable THICKNESS_used
render_to_file("06_Thickness_other_used", top=False)

enable THICKNESS_OTHER_XCM
render_to_file("06_Thickness_other_xcm", top=False)

disable THICKNESS_used
disable THICKNESS_neighbors
render_to_file("06_Thickness_other_xcm_only", top=False)
