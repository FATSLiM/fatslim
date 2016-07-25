run fatslim_pymol.py

fatslim_vesicle

disable all
zoom NS_ref, 120
turn z,-90
turn x,-90

enable leaflets_hg1
enable leaflets_hg2
enable THICKNESS_ref
enable THICKNESS_ref_normal
enable THICKNESS_REF_neighbors
render_to_file("06_Vesicle_NS_same", vesicle=True, zoom_obj="NS_ref", zoom_distance=120)

disable leaflets_hg1
disable leaflets_hg2
render_to_file("06_Vesicle_NS_same_neighbors_only", vesicle=True, zoom_obj="NS_ref", zoom_distance=120)

enable THICKNESS_REF_used
render_to_file("06_Vesicle_NS_same_used", vesicle=True, zoom_obj="NS_ref", zoom_distance=120)

disable THICKNESS_ref
enable THICKNESS_REF_XCM
render_to_file("06_Vesicle_NS_same_xcm", vesicle=True, zoom_obj="NS_ref", zoom_distance=120)

disable THICKNESS_REF_neighbors
disable THICKNESS_REF_used
disable THICKNESS_ref_normal
enable THICKNESS_normal_xcm
render_to_file("06_Vesicle_NS_same_xcm_only", top=False, vesicle=True, zoom_obj="NS_ref", zoom_distance=120)

enable THICKNESS_cutoff
render_to_file("06_Vesicle_other_cutoff", top=False, vesicle=True, zoom_obj="NS_ref", zoom_distance=120)

disable THICKNESS_cutoff
enable THICKNESS_neighbors
render_to_file("06_Vesicle_other_neighbors", top=False, vesicle=True, zoom_obj="NS_ref", zoom_distance=120)

disable THICKNESS_neighbors
enable THICKNESS_used
render_to_file("06_Vesicle_other_used", top=False, vesicle=True, zoom_obj="NS_ref", zoom_distance=120)

enable THICKNESS_OTHER_XCM
render_to_file("06_Vesicle_other_xcm", top=False, vesicle=True, zoom_obj="NS_ref", zoom_distance=120)

disable THICKNESS_used
disable THICKNESS_neighbors
render_to_file("06_Vesicle_other_xcm_only", top=False, vesicle=True, zoom_obj="NS_ref", zoom_distance=120)
