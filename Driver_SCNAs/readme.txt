## Driver-PCGs
To identify driver SCNAs combining network propagation and linear regression model

You should do the following things first
##########1#########################loading functions,data,and R packages###########################################
library(clusterProfiler)
library (tidyverse)
library(GSVA)
load("STRING.RData") ####loading PPI data
load("H1.RData") #####loading hallmark genesets
load("m_hallmark.RData")####loading hallmark geneset list

source("gene_deal.r")
source("cna_deal.r")
source("deal_data.r")
source("rw_weight.r")
source("weight_ppi.r")
source("dys_hallmark_profile.r")
source("last_program.r")

#########2#######################################
Then,run the first "last_program.r" first to invoke other programs
last_driver_net <- last_program(STRING=STRING,H1=H1,path_input,path_output,m_hallmark)

The results are stored in "output" folder,The genes in the active gene-hallmark network were identified as driver genes.
###########
Functions
##deal_data.r
To obtain copy number genetic alteration profile and gene expression profile
###cna_deal.r
To obtain copy number genetic alteration profile
###gene_deal.r
To obtain expression profile
####weight_ppi.r
Building a weighted PPI network
####rw_weight.r
To run RWR algorithm

####dys_hallmark_profile.r
Constructing hallmark dysregulation profile

####candidate_gene_hallmark.r
Building the candidate gene-hallmark Network using random walk with restart(RWR)

####plrs.r
running plrs regression to select the important regulatory factors for activity of a hallmark





