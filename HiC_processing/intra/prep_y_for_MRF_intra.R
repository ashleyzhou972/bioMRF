#!/usr/bin/env Rscript
######################################
# For general intra HiC networks
########################################
rm(list=ls())
setwd('/home/nzhou/hic/rao2014')
source('../get_edge_list.R')
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  print(args[1])
}
homefolder = '/home/nzhou/hic/rao2014/'
genename_folder = '/home/nzhou/hic/rao2014/gene_names/'
##update 20190517
##remove '1','7','17','X' from list of chromosomes to loop
##because already done
#new_chrms = chrms[-c(1,7,17,23)]

#args should be in the format of "$chrm"_"$quant"_"$method"_TRUE_"$del"

chrm_long = args[1]
#key = 'chrm22_0.9_mean_TRUE_TRUE'
cat(chrm_long, "\n")
chrm = substr(chrm_long, 5,6)

rnaseq_files = paste0("/home/nzhou/hic/IMR90/rnaseq/rnaseq_", c('RIW.txt', 'SBP.txt'))
output_folder = './IMR90_10kb/intra/by_gene/data/y/'

y <-output_y(genename_folder, chrm, rnaseq_files)
name = paste0("chrm", chrm)
save_y(output_folder, y, name)
