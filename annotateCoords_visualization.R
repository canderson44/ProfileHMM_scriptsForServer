#!/usr/bin/env Rscript
#### Description ####
# plots coord data. For now, it's toy data. 
# Goal: determine where Adapter, 3' barcode, 5' barcode, and their reverse complements
# lie within a given CCS

library(tidyverse)

install.packages("gridExtra")
library(gridExtra)
#### Load data ####
#FOR SERVER 
#coords_2B01 <- read.csv("/tier2/deweylab/scratch/ipsc_pacbio/demultiplexing/profile_hmm/annotated_ccs/2_B01_annotation_coords.csv")
#FOR MY MAC
coords_2B01 <- read.csv("../data/2_B01_annotation_coords.csv")
# coords_3C01 <- read.csv("../data/toy_coords_B.csv")
# coords_4D01 <- read.csv("../data/toy_coords_C.csv")


#### Plots ####

#FOR NOW: only 2_B01 and rest toy data
#2_B01
ggplot(coords_2B01, aes(x=start, y=y)) + 
  facet_wrap(~ZMW) + 
  geom_segment(aes(xend=stop, yend=y, color=region, size=10)) +
  labs(title="2_B01 Annotations by CCS", 
       x = "Index within CCS", y = "") +
  scale_color_brewer(breaks = c("CCS", "Adapter", "Adapter_Reverse","Five_Barcode", 
                                "Five_Barcode_Reverse", "Three_Barcode", "Three_Barcode_Reverse"),
                     palette = "Paired") +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks=element_blank())
# #3_C01
# ggplot(coords_3C01, aes(x=start, y=y)) + 
#   facet_wrap(~ZMW) + 
#   geom_segment(aes(xend=stop, yend=y, color=region, size=10)) +
#   labs(title="3_C01 Annotations by CCS", 
#        x = "Index within CCS", y = "") +
#   scale_color_brewer(breaks = c("CCS", "Adapter", "Adapter_Reverse","Five_Barcode", 
#                                 "Five_Barcode_Reverse", "Three_Barcode", "Three_Barcode_Reverse"),
#                      palette = "Paired") +
#   theme_bw() +
#   theme(axis.text.y = element_blank(), axis.ticks=element_blank())
# #4_D01
# ggplot(coords_4D01, aes(x=start, y=y)) + 
#   facet_wrap(~ZMW) +
#   geom_segment(aes(xend=stop, yend=y, color=region, size=10)) +
#   labs(title="4_D01 Annotations by CCS", 
#        x = "Index within CCS", y = "") +
#   scale_color_brewer(breaks = c("CCS", "Adapter", "Adapter_Reverse","Five_Barcode", 
#                                 "Five_Barcode_Reverse", "Three_Barcode", "Three_Barcode_Reverse"),
#                     palette = "Paired") +
#   theme_bw() +
#   theme(axis.text.y = element_blank(), axis.ticks=element_blank())





