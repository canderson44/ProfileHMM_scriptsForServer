#!/usr/bin/env python
#for each cell line, creates dictionary of emission probabilities for the following regions of the CCS sequence: 

#sequence
    #5p barcode
    #RNA insert
    #3p barcode
#insertion state

#FOR NOW: just emission probs for forward 
#TODO:    emission probs for rev complement branch

import Annotate_CCS as ac
import numpy as np

barcode2 = ac.barcode2RC #as found analyzing annotation coords
barcode2rc = ac.barcode2 #as found analyzing annotation coords
barcode3 = ac.barcode3
barcode3rc = ac.barcode3RC
barcode4 = ac.barcode4
barcode4rc = ac.barcode4RC
fiveBarcode = ac.fivePBarcode
fiveBarcodeRC = ac.fivePBarcodeRC
adapter2 = ac.two_adapter
adapter2rc = ac.two_adapter_RC
adapter3 = ac.three_adapter
adapter3rc = ac.three_adapter_RC
adapter4 = ac.four_adapter
adapter4rc = ac.four_adapter_RC

threeP_barcodes = [barcode2, barcode3, barcode4]
adapters = [adapter2, adapter3, adapter4]

RCthreeP_barcodes = [barcode2rc, barcode3rc, barcode4rc]
RCadapters = [adapter2rc, adapter3rc, adapter4rc]

cell_emission_list = [] #entry 0 for 2B01, 1 for 3C01, 2 for 4D01
#key: index in sequence (0 for first character)
#value: emission probs dict for that position


#RNAinsert and Insert states
random_dict = {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25}
###########################
###########################
###########################
#FORWARD SEQUENCE COMPONENTS
###########################
###########################
###########################


#5p barcode
fiveP_emission_list = [0] *len(fiveBarcode)
minor_percentage = 0.02
major_percentage = 0.94
for i in np.arange(len(fiveBarcode)):
    fiveP_emission_list[i] = {'A':minor_percentage, 'C':minor_percentage, 'G':minor_percentage, 'T':minor_percentage}
    bp = fiveBarcode[i]
    fiveP_emission_list[i][bp] = major_percentage

#5p barcode rc
RCfiveP_emission_list = [0] *len(fiveBarcode)
for i in np.arange(len(fiveBarcode)):
    RCfiveP_emission_list[i] = {'A':minor_percentage, 'C':minor_percentage, 'G':minor_percentage, 'T':minor_percentage}
    bp = fiveBarcodeRC[i]
    RCfiveP_emission_list[i][bp] = major_percentage


cells = ['2_B01', '3_C01', '4_D01']
for i in np.arange(3):
    cell = cells[i]
    barcode3p = threeP_barcodes[i]
    adapter = adapters[i]
    #adapter_dict = {}
    #bar_emission_dict = {}
    #3p barcode:
    this_bar_emission_list = [0]*len(barcode3p)
    for j in np.arange(len(barcode3p)):
        this_bar_emission_list[j] = {'A':minor_percentage, 'C':minor_percentage, 'G':minor_percentage, 'T':minor_percentage}
        bp = barcode3p[j]
        this_bar_emission_list[j][bp] = major_percentage
    #adapter
    this_adapter_emission_list = [0]*len(adapter)
    for adapter_index in np.arange(len(adapter)):
        this_adapter_emission_list[adapter_index] = {'A':minor_percentage, 'C':minor_percentage, 'G':minor_percentage, 'T':minor_percentage}
        bp = adapter[adapter_index]
        this_adapter_emission_list[adapter_index][bp]=major_percentage

    #combine into this cell's list of emission dicts

    this_dict = {}
    # 5p barcode
    for index in np.arange(25): # 0 to 24; 5p barcode
        match_str = "M" + str(index)
        insert_str = "I" + str(index)
        this_dict[match_str] = fiveP_emission_list[index]
        this_dict[insert_str] = random_dict

    #RNAinsert
    this_dict["RNA"] = random_dict
    this_dict["IRNA"] = random_dict

    #3p barcode
    for index in np.arange(25,96): #25-95; len 71
        match_str = "M" + str(index)
        insert_str = "I" + str(index)
        this_dict[match_str] = this_bar_emission_list[index-25]
        this_dict[insert_str] = random_dict

    #initial junk insert state
    this_dict["IS"] = random_dict

    ###########################
    ###########################
    ###########################
    # BACKWARD SEQUENCE COMPONENTS
    ###########################
    ###########################
    ###########################
    RCbarcode3p = RCthreeP_barcodes[i]
    RCadapter = RCadapters[i]
    # 3p barcode:
    this_bar_emission_list = [0] * len(RCbarcode3p)
    for j in np.arange(len(RCbarcode3p)):
        this_bar_emission_list[j] = {'A':minor_percentage, 'C':minor_percentage, 'G':minor_percentage, 'T':minor_percentage}
        bp = RCbarcode3p[j]
        this_bar_emission_list[j][bp] = major_percentage
    # adapter
    this_adapter_emission_list = [0] * len(RCadapter)
    for rc_index in np.arange(len(RCadapter)):
        this_adapter_emission_list[rc_index] = {'A':minor_percentage, 'C':minor_percentage, 'G':minor_percentage, 'T':minor_percentage}
        bp = RCadapter[rc_index]
        this_adapter_emission_list[rc_index][bp] = major_percentage

    # combine into this cell's list of emission dicts
    # 3p barcode
    for index in np.arange(71):  # 0 to 70; 3p barcode
        match_str = "Mr" + str(index)
        insert_str = "Ir" + str(index)
        this_dict[match_str] = this_bar_emission_list[index]
        this_dict[insert_str] = random_dict

    # RNAinsert
    this_dict["RNAr"] = random_dict
    this_dict["IRNAr"] = random_dict

    # 5p barcode
    for index in np.arange(71, 96):  # 71-96
        match_str = "Mr" + str(index)
        insert_str = "Ir" + str(index)
        this_dict[match_str] = RCfiveP_emission_list[index - 71]
        this_dict[insert_str] = random_dict

    #initial junk insert state
    this_dict["ISr"] = random_dict


    #finally, add dictionary to list
    cell_emission_list.append(this_dict)

#Getters.
#getter
#returns list of emissions dicts.
## One dict each for 2B01, 3C01, 4D01 cell lines, in that order
def get_all_cells_emissions():
    return cell_emission_list
#getter
#returns 2B01's emission dict
def get_2B01_emissions():
    return cell_emission_list[0]
# getter
#returns 3C01's emission dict
def get_3C01_emissions():
    return cell_emission_list[1]
# getter
#returns 4D01's emission dict
def get_4D01_emissions():
    return cell_emission_list[2]


