#!/usr/bin/env Rscript
######################################
# For general inter HiC networks
########################################
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
if (length(args)<5) {
  stop("Five arguments must be supplied\n", call.=FALSE)
} 

compare_chrms<-function(chrm1, chrm2){
  #return true if chrm1 < chrm2
  if (chrm1=='X' || chrm1=='chrmX'){
    return(FALSE)
  }
  else {
    if (chrm2=='X' || chrm2=='chrmX'){
      return(TRUE)
    }
  }
  #check if input is long form or short form
  if (nchar(chrm1)>2 & nchar(chrm2)>2) { #long form
    num1 = as.numeric(substr(chrm1, 5,6))
    num2 = as.numeric(substr(chrm2, 5,6))
  }
  else if (!grepl("\\D", chrm1) & !grepl("\\D", chrm2)){# short form only numbers
    num1 = as.numeric(chrm1)
    num2 = as.numeric(chrm2)
  }
  else {
    stop("Input correct form of chromsome names, such as chrm1 or 1\n")
  }
  if (num1<num2) return(TRUE)
  else return(FALSE);
}

read_genes<-function(genename_folder, chrms){
  all = c()
  for (chrm in chrms){
    genes = readRDS(paste0(genename_folder, '/chrm', chrm,'_genes.rds'))
    all = c(all, genes)
  }
  return(all)
}

output_y<-function(genename_folder, chrms, rnaseq_files){
  y_output = c()
  genes = read_genes(genename_folder, chrms)
  #rnaseq_files should be a vector of filepaths
  for (f in rnaseq_files){
    rnaseq = read.table(f,header=T)
    #match y to have the same order as the gene list in genes (neighborhood)
    y<-as.numeric(rnaseq[match(genes, rnaseq[,1]),2])
    y_output = cbind(y_output, y)
  }
  return(y_output)
}

save_y<-function(folder, y, name){
  saveRDS(y, paste0(folder,'/', name, '_y.rds'))
}


homedir = args[3]
genename_folder = args[4]
output_folder = paste0(homedir, '/data/y/')
rnaseq_files = args[5:length(args)]

chrm_long1 = args[1]
chrm_long2 = args[2]
chrm1 = substr(chrm_long1, 5,6)
chrm2 = substr(chrm_long2, 5,6)

if (compare_chrms(chrm1, chrm2)){
  chrms = c(chrm1, chrm2)
  y <-output_y(genename_folder, chrms, rnaseq_files)
  name = paste0("chrm", chrm1, "_chrm", chrm2)
  save_y(output_folder, y, name)
}else {
  cat("Invalid chromosome pair,", chrm1, chrm2, "\n")
}
