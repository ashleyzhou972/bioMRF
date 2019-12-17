rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("Three arguments must be supplied\n", call.=FALSE)
} 
source('./__subset_mat.R')

chrms = as.character(seq(1,22,1))
chrms = c(chrms, 'X')

TADfolder = args[1]
by_chrm_data_folder = args[2] # data path for full chromosomes gene networks
outfolder = args[3]

#nonTADs
for (chrm in chrms){
  chrm_long = paste0("chrm", chrm)
  nonTADgenes = read.table(paste0(TADfolder, '/', chrm_long, '_nonTAD.txt'), header = F)[,1]
  ret = subset_mat(by_chrm_data_folder, paste0(chrm_long, '_0.9_mean_TRUE_FALSE'), chrm_long, TRUE, nonTADgenes, 
                   outfolder, paste0("nonTADs_", chrm_long) , paste0("nonTADs_", chrm_long))
  
}