#!/usr/bin/env Rscript
#### Description ####
# Counts combinations of regions identified per zmw. Outputs a csv of format RegionCombo, count
# Goal: determine where Adapter, 3' barcode, 5' barcode, and their reverse complements
# lie within a given CCS
# The pipe taken from script by Colin Dewey

#### load packages and data ####
library(tidyverse)

install.packages("gridExtra")
library(gridExtra)

#FOR SERVER 
coords_2B01 <- read.csv("/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/profile_hmm/annotated_ccs/2_B01_annotation_coords.csv")
coords_3C01 <- read.csv("/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/profile_hmm/annotated_ccs/3_C01_annotation_coords.csv")
coords_4D01 <- read.csv("/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/profile_hmm/annotated_ccs/4_D01_annotation_coords.csv")

#FOR CATHERINE MAC
#coords_2B01 <- read.csv("../data/2_B01_annotation_coords_2019-08-5.csv")
# coords_3C01 <- read.csv("../data/toy_coords_B.csv")
# coords_4D01 <- read.csv("../data/toy_coords_C.csv")

#### pipe by Colin Dewey ####
#2_B01
pattern_counts_2 <- coords_2B01 %>% 
  filter(region != "CCS") %>%
  arrange(ZMW, start) %>%
  group_by(ZMW) %>%
  summarise(pattern = paste(region, collapse = " ")) %>%
  count(pattern) %>%
  arrange(desc(n))

#3_C01
pattern_counts_3 <- coords_3C01 %>% 
  filter(region != "CCS") %>%
  arrange(ZMW, start) %>%
  group_by(ZMW) %>%
  summarise(pattern = paste(region, collapse = " ")) %>%
  count(pattern) %>%
  arrange(desc(n))

#4_D01
pattern_counts_4 <- coords_4D01 %>% 
  filter(region != "CCS") %>%
  arrange(ZMW, start) %>%
  group_by(ZMW) %>%
  summarise(pattern = paste(region, collapse = " ")) %>%
  count(pattern) %>%
  arrange(desc(n))

#### save result ####
# for server
write_csv(pattern_counts_2, "/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/profile_hmm/annotated_ccs/2_B01_annotationAllRegionCombos.csv")
write_csv(pattern_counts_3, "/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/profile_hmm/annotated_ccs/3_C01_annotationAllRegionCombos.csv")
write_csv(pattern_counts_4, "/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/profile_hmm/annotated_ccs/4_D01_annotationAllRegionCombos.csv")
#for Catherine Mac
#write_csv(pattern_counts_2,"../2_B01_annotationRegionCombos.csv")


