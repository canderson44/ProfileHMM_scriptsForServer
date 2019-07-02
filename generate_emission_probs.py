#!/usr/bin/env python
#for each cell line, creates dictionary of emission probabilities for the following regions of the CCS sequence: 

#sequence
    #adapter (TODO: WHERE IS IT? BEGINNING OR END?)
    #5p barcode
    #RNA insert
    #3p barcode
#insertion state

#FOR NOW: just emission probs for forward 
#TODO:    emission probs for rev complement branch

import Annotate_CCS as ac
import numpy as np

barcode2 = ac.barcode2
barcode2rc = ac.barcode2RC
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

bar2emission_dict = {}
bar3emission_dict = {}
bar4emission_dict = {}
threeP_emission_lists = []
adapter_emission_lists = []
cell_emission_list = [] #entry 0 for 2B01, 1 for 3C01, 2 for 4D01
###########################
###########################
###########################
#FORWARD SEQUENCE COMPONENTS
###########################
###########################
###########################


#RNAinsert and Insert states
random_dict = {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25}




#key: index in sequence (0 for first character)
#value: emission probs dict for that position  

#5p barcode
fiveP_emission_list = [0] *len(fiveBarcode)
for i in np.arange(len(fiveBarcode)):
    fiveP_emission_list[i] = {'A':0.05, 'C':0.05, 'G':0.05, 'T':0.05}
    bp = fiveBarcode[i]
    fiveP_emission_list[i][bp] = 0.95

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
        this_bar_emission_list[j] = {'A':0.05, 'C':0.05, 'G':0.05, 'T':0.05}
        bp = barcode3p[j]
        this_bar_emission_list[j][bp] = 0.95
    #adapter
    this_adapter_emission_list = [0]*len(adapter)
    for i in np.arange(len(adapter)):
        this_adapter_emission_list[i] = {'A':0.05, 'C':0.05, 'G':0.05, 'T':0.05}
        bp = adapter[i]
        this_adapter_emission_list[i][bp]=0.95

    #combine into this cell's list of emission dicts

    this_dict = {}
    # 5p barcode
    for index in np.arange(25): # 0 to 24; 5p barcode
        match_str = "M" + str(index)
        insert_str = "I" + str(index)
        this_dict[match_str] = fiveP_emission_list[index]
        this_dict[insert_str] = random_dict

    #RNAinsert
    this_dict["R"] = random_dict
    this_dict["IR"] = random_dict

    #adapter
    for index in np.arange(25,70): # 25-69
        match_str = "M" + str(index)
        insert_str = "I" + str(index)
        this_dict[match_str] = this_adapter_emission_list[index - 25]
        this_dict[insert_str] = random_dict

    #3p barcode
    for index in np.arange(70,141): #70-140
        match_str = "M" + str(index)
        insert_str = "I" + str(index)
        this_dict[match_str] = this_bar_emission_list[index-70]
        this_dict[insert_str] = random_dict
    cell_emission_list.append(this_dict)

print(cell_emission_list[2])
#
# #BACKWARD SEQUENCE COMPONENTS
# RCthreeP_barcodes = [barcode2rc, barcode3rc, barcode4rc]
# RCadapters = [adapter2rc, adapter3rc, adapter4rc]
# RCbar2emission_dict = {}
# RCbar3emission_dict = {}
# RCbar4emission_dict = {}
# RCthreeP_emission_lists = [0,0,0]
# RCadapter_emission_lists = [0,0,0]
#
# #RNA insert
# RCrna_insert_emission_dict = {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25}
#
# #key: index in sequence (0 for first character)
# #value: emission probs dict for that position
#
# #5p barcode rc
# starting_keys = [n for n in np.arange(len(fiveBarcodeRC))]
# RCfiveP_dict = dict.fromkeys(starting_keys, '!')
# for i in np.arange(len(fiveBarcodeRC)):
#     RCfiveP_dict[i] = {'A':0.05, 'C':0.05, 'G':0.05, 'T':0.05}
#     bp = fiveBarcodeRC[i]
#     RCfiveP_dict[i][bp] = 0.95
#
#
# cells = ['2_B01', '3_C01', '4_D01']
# for i in np.arange(3):
#     cell = cells[i]
#     RCbarcode3p = RCthreeP_barcodes[i]
#     RCadapter = RCadapters[i]
#     RCadapter_dict = RCadapter_emission_lists[i]
#     RCbar_emission_dict = RCthreeP_emission_lists[i]
#     #3p barcode:
#     starting_keys =[n for n in np.arange(len(RCbarcode3p))]
#     RCbar_emission_dict = dict.fromkeys(starting_keys, '!')
#     for i in np.arange(len(RCbarcode3p)):
#         RCbar_emission_dict[i] = {'A':0.05, 'C':0.05, 'G':0.05, 'T':0.05}
#         bp = RCbarcode3p[i]
#         RCbar_emission_dict[i][bp] = 0.95
#     #adapter
#     starting_keys = [n for n in np.arange(len(RCadapter))]
#     adapter_dict = dict.fromkeys(starting_keys, '!')
#     for i in np.arange(len(adapter)):
#         RCadapter_dict[i] = {'A':0.05, 'C':0.05, 'G':0.05, 'T':0.05}
#         bp = RCadapter[i]
#         RCadapter_dict[i][bp]=0.95
#
# #######
# #######
# #STATES
# ######
# ######
#
# #INSERT STATE
# insertState_emission_dict = {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25}
#
# #SEQUENCE STATE: FORWARD
# forwardSeqState2B01 = {}
# forwardSeqState3C01 = {}
# forwardSeqState4D01 = {}
#
# #SEQUENCE STATE: BACKWARD
# backwardSeqState2B01 = {}
# backwardSeqState3C01 = {}
# backwardSeqState4D01 = {}
