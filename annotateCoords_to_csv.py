#!/usr/bin/env python
'''
For each cell line, runs Annotate_CCS annotateseq function on extracted ccs
csv file to get start and stop coords for each region
(3' barcode, 5' barcode, adapter).
Writes these to a csv in following format:
first line is headers. then:
ZMW#, region character, start, stop, yval
    'C' for full CCS
yval is height at which we'll place our line segment
    based on what region it is
C gets height 100
5' barcode or ref complement: 75
3' barcode or rev complement: 50
adapter or rev complement: 25


'''
import Annotate_CCS as ac
import numpy as np

cells = ['2_B01', '3_C01', '4_D01']
barcode_chars = ['2','3','4']
barcodes3p = [ac.barcode2, ac.barcode3, ac.barcode4]
adapters = [ac.two_adapter, ac.three_adapter, ac.four_adapter]
path_stub = '/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/'
input_path_stub = path_stub + 'demultiplexed_full_bams/no_passes/'
# input_path_stub = path_stub + 'profile_hmm/'
# Resolved RESTORE INPUT_PATH_STUB TO THAT OF ALL EXTRACTED CCS
output_path_stub = path_stub + 'profile_hmm/annotated_ccs/'

for i in np.arange(len(cells)):
    cell = cells[i]
    ref_list = [[ac.fivePBarcode, '5'], [barcodes3p[i], barcode_chars[i]], [adapters[i], 'A']]
    input_filename = input_path_stub + cell + '/extracted_ccs.csv'
    # Resolved: RESTORE INPUT FILENAME TO THAT OF ALL EXTRACTED CCS
    # input_filename = input_path_stub + cell + '_first10lines.csv'
    output_filename = output_path_stub + cell + '_annotation_coords.csv'
    ordered_filename = output_path_stub + cell + '_annotationChars_by_zmw.fasta'
    #input file format: zmw#, ccs
    ccs_list = [] #store tuples of (zmw number, ccs strings)
    with open(input_filename) as input:
        zmw_list = []
        for line in input:
            if ',' in line:
                strip_line = line.strip()
                split_list = strip_line.split(',')
                ccs_list.append((split_list[0],split_list[1]))

    with open(output_filename, 'w') as output:
        with open(ordered_filename,'w') as ordered_output:
            #first write headers line. order:
            # ZMW#, region character, start, stop, yval
            output.write('ZMW,region,start,stop,y\n')
            i = 0 #lets us know what spot we're in in the ccs list
                    #so we know when we hit last line
            for pair in ccs_list: #iterate through cell line's CCSes
                zmw_num = pair[0]
                ccs = pair[1]

                #first: write whole CCS
                to_write_str = ",".join([str(zmw_num), 'CCS', str(0),
                                        str(len(ccs)), str(100)])
                to_write_str = to_write_str + '\n'
    #            print("string to write is: ", to_write_str)
    #            print('\n\n')
                output.write(to_write_str)

                #now do regions for this ccs
                annotated_ccs = ac.annotate_seq(ccs, ref_list, justCoords=True)
                #returns final_coord_list: [(start,stop,char), ...]
                for grouping in annotated_ccs:
                    start = grouping[0]
                    stop = grouping[1]
                    region_char = grouping[2]
                    if region_char == 'A' or region_char == 'Ar':
                        y = 25
                        if region_char == 'A':
                            region_char = 'Adapter'
                        else:
                            region_char = 'Adapter_Reverse'
                    elif region_char == '5' or region_char == '5r':
                        y = 75
                        if region_char == '5':
                            region_char = 'Five_Barcode'
                        else:
                            region_char = 'Five_Barcode_Reverse'
                    else: #3' barcode or its rev complement
                        y = 50
                        if len(region_char) == 1:
                            region_char = 'Three_Barcode'
                        else:
                            region_char = 'Three_Barcode_Reverse'
                    # ZMW#, region character, start, stop, yval
                    to_write_str = ",".join([str(zmw_num), region_char, str(start),
                                             str(stop), str(y)])
                    to_write_str = to_write_str + '\n'
                    output.write(to_write_str)

                #let's sort the coord tuples so they're in order by start position
                ordered_coord_tuples = annotated_ccs.copy()
                ordered_coord_tuples.sort(key=lambda group:group[0])
                #now write it to a file
                zmw_line = '>' + str(zmw_num)
                ordered_output.write(zmw_line)
                ordered_refChars = [group[2] for group in ordered_coord_tuples]
                refChars_string = " ".join(ordered_refChars)
                if i < len(ccs_list)-1: #not at end of list of zmws yet
                    refChars_string +=1
                ordered_output.write(refChars_string)


                #increment i
                i += 1
