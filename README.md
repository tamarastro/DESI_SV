# DESI_SV
Repository for code relating to Dark Energy Spectroscopic Instrument (DESI) survey validation and visual inspection.  
To merge VI you:
1. Run VI_merge_sets_v1.py (changing the filenames in the code to the ones you want to merge)
2. Run VI_merge_merged_files.py (which merges all of the output of VI_merge_sets_v1.py and generates a truth table)
3. Run VI_merge_plots.py (to get some diagnostic plots)

To compare VI outputs:
1. Run VI_compare.py (chaging the filenames in the code to the ones you want to compare, usually one VI vs the truth table if truth exists)

Also:
1. VI_convert_v1pt1_to_v1pt2.py converts the old v1.1 from Rita to our common v1.2 format. 

