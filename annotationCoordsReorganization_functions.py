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
def remove_overlapping_fiveCoords(zmw_dict, output_filename, writeCSV):
    for zmw, zmw_region_dict in zmw_dict.items():
        new_five_regions = [] #holds kept fiveBar coords
        new_fiveRC_regions = [] # holds kept fiveBarRC coords
        # list of sets of coords for fiveBar and fiveRC regions
        five_regions = zmw_region_dict['Five_Barcode_Reverse']
        fiveRC_regions = zmw_region_dict['Five_Barcode']
        # list of sets of coords for threeBar and threeBarRC regions
        three_regions = zmw_region_dict['Three_Barcode_Reverse'] + zmw_region_dict['Three_Barcode']

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

        #now update dict
        #forward
        if len(five_regions) != len(new_five_regions): #need to update
            del zmw_region_dict['Five_Barcode']
            zmw_region_dict['Five_Barcode'] = new_five_regions

        #reverse
        if len(fiveRC_regions) != len(new_fiveRC_regions): #need to update
            del zmw_region_dict['Five_Barcode_Reverse']
            zmw_region_dict['Five_Barcode_Reverse'] = new_fiveRC_regions

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
            if region == 'Five_Barcode' or region== 'Five_Barcode_Reverse':
                isFive = True
            elif region == 'Three_Barcode' or region == 'Three_Barcode_Reverse':
                isThree = True
        if isFive == False or isThree == False: #don't keep this zmw
            reject_zmws_list.append(zmw)

    #now make new zmw dict
    return_zmw_dict = zmw_dict.copy()
    for zmw in reject_zmws_list:
        del return_zmw_dict[zmw]

    return return_zmw_dict








