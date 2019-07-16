#!/usr/bin/env python
'''
Script that uses functions in annotationCoordsReorganization_functions.py
For each cell line:
    - get counts of each combo of 5(r) with 3(r)
        combos: 53, 5r3, 5r3r, 53r, 35, 3r5, 3r5r, 35r
    - exports these counts to csv so I can use R ggplot to visualize them.
'''

import annotationCoordsReorganization_functions as af
import numpy as np

#cells = ['2_B01', '3_C01', '4_D01']
# TODO RESTORE CELLS TO ALL CELL LINES
# FOR NOW, JUST 2_B01 BECAUSE THAT'S THE ONLY ONE WITH COMPLETE COORD ANNOTATION
cells = ['2_B01','test']
path_stub = '/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/profile_hmm/annotated_ccs/'
for i in np.arange(len(cells)):
    cell = cells[i]
    #0: create zmw dict
    initial_zmw_dict = af.gen_ZMW_dict(cell)
#    print("class of initial zmw dict is:", type(initial_zmw_dict))
#    print("initial zmw dict is: ", initial_zmw_dict )

    #1: remove 5(r) coords that overlap with 3(r) regions (and make csv for future use)
    output_filename = path_stub + cell + '_annotation_coords_noOverlap.csv'
    af.remove_overlapping_fiveCoords(initial_zmw_dict,output_filename, True )
#    print("in coordinating script")
#    print("after remove overlap function, zmw_dict is", initial_zmw_dict)
    #2: select only ZMWs with {5 or 5r} AND {3 or 3r}
    pair_zmw_dict = af.select_fiveThreePairs(initial_zmw_dict)
#    print("in coordinating script")
#    print("pair_zmw_dict is", pair_zmw_dict)
    #3: get combo counts
    combo_counts = af.count_combos(pair_zmw_dict)
    print()
    print("combo counts is",combo_counts)

    #4: export combo counts to csv
    output_filename = path_stub + cell + '_comboCounts.csv'
    af.combo_counts_to_csv(combo_counts, output_filename)
