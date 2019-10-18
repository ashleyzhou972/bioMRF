#!/usr/bin/env Rscript
######################################
# For general inter HiC networks
########################################
rm(list=ls())
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

chrm1_long = args[1]
chrm2_long = args[2]
cat(chrm1_long, "\n")
cat(chrm2_long, "\n")
chrm1 = substr(chrm1_long, 5,6)
chrm2 = substr(chrm2_long, 5,6)

rnaseq_files = paste0("/home/nzhou/hic/IMR90/rnaseq/rnaseq_", c('RIW.txt', 'SBP.txt'))
output_folder = './IMR90_10kb/inter/by_gene/data/y/'

if (compare_chrms(chrm1, chrm2)){
	chrms = c(chrm1, chrm2)
	y <-output_y(genename_folder, chrms, rnaseq_files)
	name = paste0("chrm", chrm1, "_chrm", chrm2)
	save_y(output_folder, y, name)
}else {
	cat("Invalid chromosome pair,", chrm1, chrm2, "\n")
}
