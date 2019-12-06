#!/usr/bin/env Rscript
######################################
# For general intra HiC networks
########################################
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("Four arguments must be supplied\n", call.=FALSE)
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

homedir = args[2]
genename_folder = args[3]
output_folder = paste0(homedir, '/data/y/')
rnaseq_files = args[4]
print(rnaseq_files)

chrm_long = args[1]
cat(chrm_long, "\n")
chrm = substr(chrm_long, 5,6)

y <-output_y(genename_folder, chrm, rnaseq_files)
name = paste0("chrm", chrm)
save_y(output_folder, y, name)
