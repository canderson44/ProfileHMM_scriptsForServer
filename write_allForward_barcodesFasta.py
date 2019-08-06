#!/usr/bin/env python

'''
Writes new barcodes fasta file, this one with all barcodes in the empirically 
determined forward direction. 
FOR NOW: only 2_B01 has been analyzed. the other two cell lines are same as in barcodes_for_profileHMM.fasta
'''
import Annotate_CCS as ac
with open('/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/barcodes_allForward_forHMM.fasta','w') as output:
    #5' barcode: already in forward direction
    five_str = ">primer_5p\n" + ac.fivePBarcode + '\n'
    output.write(five_str)
    
    #4_D01 3' barcode
    #TODO update after annotation analysis complete
    four_str = ">CD1_4D01_3p\n" + ac.barcode4 + '\n'
    output.write(four_str)
    
    #3_C01 3' barcode
    three_str = ">CD2_3C01_3p\n" + ac.barcode3 + '\n'
    output.write(three_str)

    #2_B01 3' barcode: reverse complement of ac's barcode2
    two_str = ">CD3_2B01_3p\n" + ac.barcode2RC + '\n'
    output.write(two_str)
