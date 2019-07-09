#!/usr/bin/env Rscript
#### Description ####
# plots coord data. For now, it's toy data. 
# Goal: determine where Adapter, 3' barcode, 5' barcode, and their reverse complements
# lie within a given CCS

library(tidyverse)

install.packages("gridExtra")
library(gridExtra)
#### Load data ####
coords_2B01 <- read.csv("../data/toy_coords_A.csv")
coords_3C01 <- read.csv("../data/toy_coords_B.csv")
coords_4D01 <- read.csv("../data/toy_coords_C.csv")

#### Plots ####

#FOR NOW: toy plots
#2_B01
ggplot(coords_2B01, aes(x=start, y=y)) + 
  facet_wrap(~ZMW#)








