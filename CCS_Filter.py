#!/usr/bin/env python
# coding: utf-8

# Program for separating useable ZMW CCSs from unusable. Useable if don't begin with 3' barcode and do end with 3' barcode. 
# 
# kept sequence must have 1 and 2 as follows:
#     1. 
#         3' barcode in end region of sample (sample is forward transcript)
#         OR
#         3' barcode reverse complement in beginning region of sample (sample is rev compliment)
#     AND
#     2. 
#         no 3' barcode in beginning region of sample (forward transcript)
#         no 3' barcode reverse complement in end region of sample (rc transcript)

# In[12]:


import numpy as np
import random as rd


# In[13]:


#given DNA nucleotide, returns its complement
def give_nuc_complement (nucleotide):
    if nucleotide == "A":
        return 'T'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    elif nucleotide == 'T': 
        return 'A'
    else: #
        return '?'


# In[14]:


#generates reverse compliment. Returns rc in standard 5' to 3' direction
def gen_rev_complement(sequence):
    complement_list = []
    for nuc in sequence:
        complement_list.append(give_nuc_complement(nuc))
    #have reverse complement, but it's backwards
    complement_list.reverse()
    return "".join(complement_list)

'''#test this
test_one = 'ACTG' #expected reverse complement: CAGT
test_two = 'AAATG' #expect rc CATTT
assert gen_rev_complement(test_one) == 'CAGT', "expected CAGT but got: " + gen_rev_complement(test_one)
assert gen_rev_complement(test_two) == 'CATTT', "expected CATTT but got: " + gen_rev_complement(test_two)'''
a=1


# In[15]:


#gets your beginning and ending sample from sequence. Either one is 1.5 * passed length. 
def get_beg_and_end(sample, barcode_length):
    segment_length = int(1.5*barcode_length)
    beg = sample[:segment_length]
    end = sample[len(sample)-segment_length:]
    return [beg,end]


# In[16]:


#barcodes_list elements: [5', CD1, CD2, CD3] 
                     #ie [5', 4_D01, 3_C01, 2_B01]
barcodes_list = []
with open('barcodes_for_profileHMM.fasta') as f: 
    for line in f:
        line = line.strip()
        if line[0]!='>':
            barcodes_list.append(line)
barcode_4 = barcodes_list[1]
barcode_3 = barcodes_list[2]
barcode_2 = barcodes_list[3]
barcode_2_rc = gen_rev_complement(barcode_2)
print('len barcode 2',len(barcode_2))
print('len barcode 5', len(barcodes_list[0]))

#adapter
two_adapter = ''
with open('two_B01.adapters.fasta') as f:
    for line in f: 
        line = line.strip()
#        print('line',line)
        if line[0]!='>':
            two_adapter = line
print("len two adapter: ",len(two_adapter))
            


# In[17]:


#test "CCSs"
nucleotides = ['A','C','G','T']


#generates a list of test strings given 3p barcode and its rev complement
#list is: 
#exact matches to barcode:
##endsIn3RC, begIn3RC, endIn3, begIn3, rand, endIn3
#substituted, inserted, and deleted versions of barcode:
#endIn3RC * 3, begIn3RC * 3, endIn3 * 3, begIn3 * 3, endIn3 * 3
def gen_test_strings(barcode, barcode_rc):
    test_strings = []

    #First: exact matches
    #exclude
    test_end_in_3pBarcodeRC = "".join(rd.choices(nucleotides, k=1000)) + barcode_rc
    test_strings.append(test_end_in_3pBarcodeRC)
    #keep
    test_beg_3pBarcodeRC = barcode_rc + "".join(rd.choices(nucleotides, k=1000))
    test_strings.append(test_beg_3pBarcodeRC)
    #keep
    test_5pBarcodeRandom3pBarcode = barcodes_list[0] + "".join(rd.choices(nucleotides, k=1000)) + barcode
    test_strings.append(test_5pBarcodeRandom3pBarcode)
    #exclude
    test_3pBarcodeThenRandom= barcode + "".join(rd.choices(nucleotides, k=(1000)))
    test_strings.append(test_3pBarcodeThenRandom)
    #exclude
    test_random="".join(rd.choices(nucleotides, k=(1000)))
    test_strings.append(test_random)
    #keep
    test_endIn3pBarcode= "".join(rd.choices(nucleotides, k=1000)) + barcode
    test_strings.append(test_endIn3pBarcode)


    #Next: mutated matches
    bar2_list = list(barcode)
    rand_subs_list = rd.choices(np.arange(len(barcode)), k=int(len(barcode)/10)) #pacbio MUCH more accurate
                                                            #here, error in sequencing is 10%
                                                            #pacbio has error of about 0.001%
                                                #source: https://www.pacb.com/uncategorized/a-closer-look-at-accuracy-in-pacbio/
    sub_bar2_list = bar2_list.copy()
    rand_nuc_list = rd.choices(nucleotides, k=len(rand_subs_list))
    for i in np.arange(len(rand_nuc_list)):
        sub_bar2_list[i] = rand_nuc_list[i]
    sub_bar2_str = "".join(sub_bar2_list)
    sub_bar2_str_rc = gen_rev_complement(sub_bar2_str)

    insert_bar2_list = bar2_list.copy()
    for i in np.arange(len(rand_nuc_list)):
        insert_bar2_list.insert(i, rand_nuc_list[i])
    insert_bar2_str = "".join(insert_bar2_list)
    insert_bar2_str_rc = gen_rev_complement(insert_bar2_str)

    delet_bar2_list = bar2_list.copy()
    for i in np.arange(int(len(rand_nuc_list))):     
        del bar2_list[i]
    delet_bar2_str = "".join(delet_bar2_list)
    delet_bar2_str_rc = gen_rev_complement(delet_bar2_str)

    #exclude
    test_end3RC_sub = "".join(rd.choices(nucleotides, k=1000)) + sub_bar2_str_rc
    test_end3RC_insert = "".join(rd.choices(nucleotides, k=1000)) + insert_bar2_str_rc
    test_end3RC_del = "".join(rd.choices(nucleotides, k=1000)) + delet_bar2_str_rc
    test_strings.append(test_end3RC_sub)
    test_strings.append(test_end3RC_insert)
    test_strings.append(test_end3RC_del)
    #keep
    test_beg3RC_sub = sub_bar2_str_rc + "".join(rd.choices(nucleotides, k=1000))
    test_beg3RC_insert = insert_bar2_str_rc + "".join(rd.choices(nucleotides, k=1000))
    test_beg3RC_del = delet_bar2_str_rc + "".join(rd.choices(nucleotides, k=1000))
    test_strings.append(test_beg3RC_sub)
    test_strings.append(test_beg3RC_insert)
    test_strings.append(test_beg3RC_del)

    #keep
    test_5pRand3p_sub = barcodes_list[0] + "".join(rd.choices(nucleotides, k=1000)) + sub_bar2_str
    test_5pRand3p_insert = barcodes_list[0] + "".join(rd.choices(nucleotides, k=1000)) + insert_bar2_str
    test_5pRand3p_del = barcodes_list[0] + "".join(rd.choices(nucleotides, k=1000)) + delet_bar2_str
    test_strings.append(test_5pRand3p_sub)
    test_strings.append(test_5pRand3p_insert)
    test_strings.append(test_5pRand3p_del)

    #exclude
    test_beg3p_sub = sub_bar2_str + "".join(rd.choices(nucleotides, k=1000))
    test_beg3p_insert = insert_bar2_str + "".join(rd.choices(nucleotides, k=(1000)))
    test_beg3p_del = delet_bar2_str + "".join(rd.choices(nucleotides, k=(1000)))
    test_strings.append(test_beg3p_sub)
    test_strings.append(test_beg3p_insert)
    test_strings.append(test_beg3p_del)

    #keep
    test_end3p_sub = "".join(rd.choices(nucleotides, k=1000)) + sub_bar2_str
    test_end3p_insert = "".join(rd.choices(nucleotides, k=1000)) + insert_bar2_str
    test_end3p_del = "".join(rd.choices(nucleotides, k=1000)) + delet_bar2_str
    test_strings.append(test_end3p_sub) 
    test_strings.append(test_end3p_insert) 
    test_strings.append(test_end3p_del)
    
    return test_strings

test_2B01_strings = gen_test_strings(barcode_2, barcode_2_rc)
test_3C01_strings = gen_test_strings(barcode_3, gen_rev_complement(barcode_3))
test_4D01_strings = gen_test_strings(barcode_4, gen_rev_complement(barcode_4))

expected_keep_results = [False, True, True, False, False, True, False, False, False, True, True, True, 
                             True, True, True, False, False, False, 
                             True, True, True]
#expected_keep_results


# In[18]:



'''#test csv file 
f = open('test_strings.csv','w+')
i = 0
for string in test_strings:
    f.write("" + str(i)+',' + string + '\n')
    i += 1
f.close()'''


# In[19]:


score_dict = {}
mismatch_penalty = -2
score_dict[('A','A')] = 1
score_dict[('A','C')] = mismatch_penalty
score_dict[('A','G')] = mismatch_penalty
score_dict[('A','T')] = mismatch_penalty
score_dict[('C','A')] = mismatch_penalty
score_dict[('C','C')] = 1
score_dict[('C','G')] = mismatch_penalty
score_dict[('C','T')] = mismatch_penalty
score_dict[('G','A')] = mismatch_penalty
score_dict[('G','C')] = mismatch_penalty
score_dict[('G','G')] = 1
score_dict[('G','T')] = mismatch_penalty
score_dict[('T','A')] = mismatch_penalty
score_dict[('T','C')] = mismatch_penalty
score_dict[('T','G')] = mismatch_penalty
score_dict[('T','T')] = 1


#test matrix

'''test = np.zeros((4,3))
print(test)
test[3,0] = 9
print(test)'''


# In[ ]:





# In[ ]:





# In[ ]:





# In[20]:


#performs SmithWaterman adaptation (glocal alignment). Returns score; no traceback
    #FOR NOW: linear gap penalty d
    #parameters: gap penalty, seq x (vertical), barcode reference (horizontal),  
    #initialization: 0s in first row; (i,0) = -i*d
def score_alignment(gap_penalty, sequence, barcode, score_dict):
    gap_penalty = np.abs(gap_penalty)
    #matrix
    num_rows = len(sequence) +1
    num_cols = len(barcode) +1
    F = np.zeros((num_rows, num_cols))
    #initialization\n",
    #row 0 already initialized to zeros
    for i in np.arange(1,num_rows):
        F[i][0] = -1*i*gap_penalty

    for i in np.arange(1,num_rows): # i represents the row
        for j in np.arange(1,num_cols): #j represents the col
#            print('i is', i)
#            print('j is', j)
#            print(\"looking at sequence char\", sequence[i-1])
#            print(\"looking at barcode char\", barcode[j-1])
            m = F[i-1][j-1] + score_dict[(sequence[i-1],barcode[j-1])]
            ix = F[i-1][j] - gap_penalty
            iy = F[i][j-1] - gap_penalty
#            print('m is', m)
#            print('ix is',ix)
#            print('iy is', iy)

            F[i][j] = max(m,ix,iy)
    #now matrix is filled\n",
    score = max(F[num_rows-1])
    return score


# In[21]:


#compute s prime = score + log(pm/pr)
#score: alignment score
#pm: prior prob of m
#pr: prior prob of r
def compute_s_prime(score, pm, pr):
    return score + np.log(pm/pr)


# In[22]:


#computes sigma function e^x / (1+e^x)
def compute_sigma(x):
    return np.exp(x)/(1.0+np.exp(x))


# In[23]:


#calculates score significance with Bayesian approach: P(M|x,y)
#if P(M|x,y)>0.9, x and y related
#possibilities for sequence of interest: 3' adapter, 3' barcode, 
                        #RNA insert, 5' barcode, 5' adapter. Only one is M
#P(M) = 1/5 Match model
#P(R) = 4/5 Random model
def calc_score_sig(score):
    pr = 4.0/5
    pm = 1.0/5
    s_prime = compute_s_prime(score, pr, pm)  
    return compute_sigma(s_prime)


# In[24]:


#given sequence, of CCS, 
#determines whether it should be kept
#True means keep, false means exclude

#kept sequence must have 1 and 2 as follows:
    #1. 
        #3' barcode in end region of sample (if sample is forward transcript)
        #OR
        #3' barcode reverse compliment in beginning region of sample (if sample is rev compliment)
    #AND
    #2. 
        #no 3' barcode in beginning region of sample (forward transcript)
        #no 3' barcode reverse complement in end region of sample (rc transcript)

#simplified. 
def should_we_keep(sample, barcode):
    segments = get_beg_and_end(sample, len(barcode))
    barcode_rc = gen_rev_complement(barcode)
    #beg
    #score_alignment(gap_penalty, sequence, barcode, score_dict):
    barcode_beg_score_sig = calc_score_sig(score_alignment(1,segments[0],barcode,score_dict))
    barcodeRC_beg_score_sig = calc_score_sig(score_alignment(1,segments[0],barcode_rc,score_dict))
    barcode_end_score_sig = calc_score_sig(score_alignment(1,segments[1],barcode,score_dict))
#    print("barcode beg score_sig:",barcode_beg_score_sig)
#    print("barcodeRC_beg_score_sig: ", barcodeRC_beg_score_sig)
#    print("barcode_end_score_sig:",barcode_end_score_sig)
#    print()
    
    if barcode_beg_score_sig >0.99999999:
#        print("beg: barcode in beg. Returning False.")
#        print()
        return False            
    elif  barcodeRC_beg_score_sig>0.99999999:
#        print('beg: barcodeRC in beg. Returning True')
#        print()
        return True
    #end                        
    elif barcode_end_score_sig > 0.99999999:
#        print("end: barcode in end. Returning True")
#        print()
        return True
    else: #if we've gottent to this point, either the barcode isn't in teh 
        #sample and/or the barcode rev comp is at the end, which we don't want
#        print('end: Barcode not in end and/or barcode rc in end; Returning False')
#        print()
        return False
        

#test it
#tests should_we_keep function given barcode and test strings
def test_should_we_keep(test_strings, barcode): 
    keep_list = []
    for entry in test_strings:
        keep_list.append(should_we_keep(entry, barcode))
    print('results ',keep_list)
    print('expected',expected_keep_results)
    
#test_should_we_keep(test_4D01_strings, barcode_4)    
a=1
#endsIn3RC, begIn3RC, endIn3, begIn3, rand, endIn3
#endIn3RC, begIn3RC, endIn3, begIn3, rand, endIn3


# In[25]:


#given csv file and desired file name for kept CCSs, filters CCS
def filter_ccs(all_ccs_filename, filtered_filename, barcode):
    retained_csv_list = []
    with open(all_ccs_filename) as f:
        for line in f:
            line = line.strip()
            line_list = line.split(',')
            ccs = line_list[1]
            if should_we_keep(ccs, barcode):
                retained_csv_list.append(line_list)
    f = open(filtered_filename, 'w+')
    for entry in retained_csv_list:
        f.write(entry[0] + ',' + entry[1])
        f.write('\n')
    f.close()  


# In[26]:


#test the function
'''
filter_ccs('test_strings.csv', 'filtered_test_strings.csv', barcode_2)
'''


# In[ ]:





# In[ ]:




