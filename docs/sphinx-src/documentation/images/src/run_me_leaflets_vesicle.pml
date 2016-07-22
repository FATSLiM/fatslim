run fatslim_pymol.py

fatslim_vesicle

disable all
zoom NS_ref, 120
turn z,-90
turn x,-90

enable leaflets_hg1
enable leaflets_hg2
render_to_file("05_Vesicle_leaflets", vesicle=True)
