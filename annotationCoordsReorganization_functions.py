#!/usr/bin/env python
'''
Set of functions to reorganize info from the {cell}_annotation_coords.csv file
for each cell line (2_B01, 3_C01, 4_D01)
NOTE: stop coord is inclusive
'''

import numpy as np
import Annotate_CCS as ac
cells = ['2_B01', '3_C01', '4_D01']
barcode_chars = ['2','3','4']
fiveBar = ac.fivePBarcode
fiveBarRC = ac.fivePBarcodeRC
threeBar_list = [ac.barcode2, ac.barcode3, ac.barcode4]
threeBarRC_list = [ac.barcode2RC, ac.barcode3RC, ac.barcode4RC]

five_name = 'Five_Barcode'
fiveR_name = 'Five_Barcode_Reverse'
three_name = 'Three_Barcode'
threeR_name = 'Three_Barcode_Reverse'

path_stub = '/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/profile_hmm/annotated_ccs/'

'''
function to organize {cell}_annotation_coords.csv by ZMW
    Dictionary where key is ZMW
    value is another dictionary: 
        {region: [(start, stop), ...], diffRegion: [(start,stop)],...}
        key are region
        value: list of (start,stop) pairs of coords. This allows for a region 
        to occur multiple times in a zmw
PARAMETERS: cell name
RETURNS: dictionary for that cell's ZMWs
'''
def gen_ZMW_dict(cell):
    input_filename = path_stub + cell + '_annotation_coords.csv'
    return_dict = {}
    with open(input_filename) as input:
        for line in input:
            if ',' in line:
                strip_line = line.strip()
                split_list = strip_line.split(',')
                zmw = split_list[0]
                region = split_list[1]
                start = split_list[2]
                stop = split_list[3]
                #don't need to record entire CCS coords (ie start/stop coords)
                if region != 'CCS':
                    if zmw in return_dict:
                        if region in return_dict[zmw]:
                            return_dict[zmw][region].append((start,stop))
                        else: #first occurence of this region for this zmw
                            return_dict[zmw][region] = [(start,stop)]
                    else: #first occurence of this zmw's listings
                        #so also first occurence of this region for this zmw
                        return_dict[zmw] = {region:[(start,stop)]}
    return return_dict

'''
function to remove 5' barcode or its rev complement if its coords overlap with 
coords of 3' barcode or its rev complement.
Updates the dictionary given and if desired, writes a new csv file (so the results are reusable)
PARAMETERS: zmw dictionary, output filename, writeCSV
    writeCSV is true if you want it to make a CSV
RETURNS: null
'''
def remove_overlapping_fiveCoords(zmw_dict, output_filename = '', writeCSV=False):
    print("in remove_overlapping_fiveCoords")
    for zmw, zmw_region_dict in zmw_dict.items():
        print("made it past first zmw_dict.items")
        new_five_regions = [] #holds kept fiveBar coords
        new_fiveRC_regions = [] # holds kept fiveBarRC coords
        # list of sets of coords for fiveBar and fiveRC regions
        five_regions = zmw_region_dict[five_name]
        fiveRC_regions = zmw_region_dict[fiveR_name]
        # list of sets of coords for threeBar and threeBarRC regions
        three_regions = zmw_region_dict[threeR_name] + zmw_region_dict[three_name]

        #filter out undesireable five coords

        for five_coords in five_regions:
            for three_coords in three_regions:
                three_range = [n for n in np.arange(three_coords[0], three_coords[1] + 1)]
                five_range  = [n for n in np.arange( five_coords[0],  five_coords[1] + 1)]
                if not((five_coords[0] in three_range) or (five_coords[1] in three_range) or\
                        (three_coords[0] in five_range) or (three_coords[1] in five_range)):
                    #then we do want to keep the five_coords
                    new_five_regions.append(five_coords)
        for five_coords in fiveRC_regions:
            for three_coords in three_regions:
                three_range = [n for n in np.arange(three_coords[0], three_coords[1] + 1)]
                five_range  = [n for n in np.arange( five_coords[0],  five_coords[1] + 1)]
                if not((five_coords[0] in three_range) or (five_coords[1] in three_range) or\
                        (three_coords[0] in five_range) or (three_coords[1] in five_range)):
                    #then we do want to keep the five_coords
                    new_fiveRC_regions.append(five_coords) #appends a tuple
        print("new five regions:", new_five_regions)
        print("new fiveRC regions:", new_fiveRC_regions)
        #now update dict
        #forward
        if len(five_regions) != len(new_five_regions): #need to update
            del zmw_region_dict[five_name]
            zmw_region_dict[five_name] = new_five_regions

        #reverse
        if len(fiveRC_regions) != len(new_fiveRC_regions): #need to update
            del zmw_region_dict[fiveR_name]
            zmw_region_dict[fiveR_name] = new_fiveRC_regions
        print("after deletions and insertions, zmw_region_dict:", zmw_region_dict)
    print("after all edits, final zmw_dict", zmw_dict)
    #############
    #############
    #now write new csv, if desired
    if writeCSV:
        with open(output_filename, 'w') as output:
            #header line
            header_line = 'ZMW,region,start,stop\n'
            output.write(header_line)
        #have zmw, need to iterate through region dict to write lines of format:
        # zmw, region, start, stop
        for zmw, zmw_region_dict in zmw_dict.items():
            for region, coords_list in zmw_region_dict.items():
                for coord_pair in coords_list:
                    write_line = ",".join([str(zmw), region, str(coord_pair[0]), str(coord_pair[1])])
                    write_line += '\n'
                    output.write(write_line)



'''
Function to select only zmws with both 5' barcode (or rev complement) and 
3' barcode (or its rev complement).
Given zmw dictionary, will shallow copy entries of zmws that make the cut. 
ZMW dictionary format:
    Dictionary where key is ZMW
    value is another dictionary: 
        {region: [(start, stop), ...], diffRegion: [(start,stop)],...}
        key are region
        value: list of (start,stop) pairs of coords.

PARAMETERS: zmw_dict
RETURN: new dict containing only zmws that make cut
'''
def select_fiveThreePairs(zmw_dict):
    reject_zmws_list = []
    for zmw, zmw_regions_dict in zmw_dict.items():
        isFive = False #mark true if 5' barcode or rev complement present
        isThree = False #mark true if 3' barcode or rev complement present
        for region in zmw_regions_dict.keys():
            if region == five_name or region== fiveR_name:
                isFive = True
            elif region == three_name or region == threeR_name:
                isThree = True
        if isFive == False or isThree == False: #don't keep this zmw
            reject_zmws_list.append(zmw)

    #now make new zmw dict
    return_zmw_dict = zmw_dict.copy()
    for zmw in reject_zmws_list:
        del return_zmw_dict[zmw]

    return return_zmw_dict

'''
function to iterate through a zmw dict and tally instances of each possible pairing 
of 5 or 5r with 3 or 3r
PARAMETERS: zmw dictionary with all overlapping 5(r) coords removed adn only containing
            zmws with 5(r) - 3(r) pairs
                (see this file's functions:
                    remove_overlapping_fiveCoords and select_fiveTreePairs, respectively
                    
RETURNS:new dict cataloguing counts of each pair combo
        key: combo name
        value: count
'''
def count_combos(zmw_dict):
    combo_counts = {'five_three':0, 'fiveR_three':0, 'fiveR_threeR':0, 'five_threeR':0,
                    'three_five':0, 'threeR_five':0, 'threeR_fiveR':0, 'three_fiveR':0}
    for zmw, zmw_region_dict in zmw_dict:
        five_type = '' #will be Five_Barcode or Five_Barcode_Reverse
        three_type = ''# will be Three_Barcode or Three_Barcode_Reverse

        #what regions present?
        for region in zmw_region_dict.keys():
            if region == five_name or fiveR_name:
                five_type = region
            elif region == three_name or threeR_name:
                three_type = region

        #now update combo_counts
        five_coords_list = zmw_region_dict[five_type]
        three_coords_list = zmw_region_dict[three_type]
        for five_coords in five_coords_list:
            for three_coords in three_coords_list:
                #five type before three type: five end < three start
                if five_coords[1] < three_coords[0]:
                    if five_type == five_name: #5
                        if three_type == three_name:#53
                            combo_counts['five_three'] += 1
                        else: #53r
                            combo_counts['five_threeR'] += 1
                    else: #5r
                        if three_type == three_name:  # 5r3
                            combo_counts['fiveR_three'] += 1
                        else:  # 5r3r
                            combo_counts['fiveR_threeR'] += 1

    return combo_counts


'''
function to export 5(R) barcode - 3(R) barcode combo counts to csv
format of csv: first line is column headers
                following lines: one per combo type
                    comboName, count
PARAMETERS combo_counts: dictionary of counts of each possible 5(r) 3(r) combo
                            keys are combos, vals are counts
            output_filename
RETURNS: null
'''
def combo_counts_to_csv(combo_counts, output_filename):
    with open(output_filename, 'w') as output:
        header_line = 'Combination,Count\n'
        output.write(header_line)
        for combo, count in combo_counts.items():
            write_string = ",".join([combo,str(count)])
            write_string += '\n'
            output.write(write_string)











