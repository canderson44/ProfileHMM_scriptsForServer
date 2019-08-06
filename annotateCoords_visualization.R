#!/usr/bin/env Rscript
#### Description ####
# Counts combinations of regions identified per zmw. Outputs a csv of format RegionCombo, count
# Goal: determine where Adapter, 3' barcode, 5' barcode, and their reverse complements
# lie within a given CCS

library(tidyverse)

install.packages("gridExtra")
library(gridExtra)
#### Load data ####
#FOR SERVER 
#coords_2B01 <- read.csv("/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/profile_hmm/annotated_ccs/2_B01_annotation_coords.csv")
#FOR MY MAC
coords_2B01 <- read.csv("../data/2_B01_annotation_coords_2019-08-5.csv")
# coords_3C01 <- read.csv("../data/toy_coords_B.csv")
# coords_4D01 <- read.csv("../data/toy_coords_C.csv")


test <- coords_2B01 %>% filter(ZMW <4260000) %>% filter(region != 'CCS') %>% select (ZMW,region) %>% 
  unite("region_combos",region, sep=" ")

for (x in seq_along(coords_2B01)){
  print(coords_2B01[[x]])
}
