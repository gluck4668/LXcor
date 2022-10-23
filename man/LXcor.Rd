\name{LXcor}
\alias{LXcor}
\title{LXcor::LXcor}
\usage{
LXcor(gene_file,meta_file,group1,group2,row_n)
}
\description{
LXcor can be used for the correlation analysis of transctiptome and metabolome
}
\examples{

if(!requireNamespace("devtools"))
  install.packages("devtools")

library(devtools)

install_github("gluck4668/LXcor")

library(LXcor)

??LXcor

#------------------------------------
data(gene_data_example)
data(meta_data_example)

#------------------------------------

setwd("D:/LXcor_example")

gene_file="genes_data.xlsx"
meta_file="meta_data.xlsx"

group1="Model"
group2="Normal"

row_n=1

LXcor(gene_file,meta_file,group1,group2,row_n)


}
