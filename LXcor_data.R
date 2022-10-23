
setwd("D:/R-lin study/R packages/LXcor")
library(openxlsx)

gene_data_example <- read.xlsx("gene_data.xlsx")
meta_data_example <- read.xlsx("meta_data.xlsx")


usethis::use_data(gene_data_example,overwrite = T)
usethis::use_data(meta_data_example,overwrite = T)

rm(list=ls())

data(gene_data_example)
data(meta_data_example)

