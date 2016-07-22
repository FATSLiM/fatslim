run fatslim_pymol.py

fatslim_bilayer

disable all

enable lipids
enable water
refresh
render_to_file("00_raw_system", top=False)


disable lipids
disable water
enable Ref_full
enable Ref_Headgroup
zoom Ref_full, 20
render_to_file("01_simplification_full", top=False)

enable Ref_simplified
zoom Ref_full, 20
render_to_file("01_simplification_both", top=False)

disable Ref_full
disable Ref_Headgroup
render_to_file("01_simplification_simplified", top=False)


disable Ref_Headgroup
disable Ref_simplified
enable pbcbox
enable positions_lower
enable directions_lower
enable positions_upper
enable directions_upper
zoom all, 5
render_to_file("02_simplified_system", top=False)
