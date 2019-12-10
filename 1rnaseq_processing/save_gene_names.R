#!/usr/bin/env Rscript
rm(list=ls())

args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Two argument must be supplied\n", call.=FALSE)
}
chrms = as.character(seq(1,22,1))
chrms = c(chrms, 'X')

rnaseq_file = args[1]
outputfolder = args[2]
rnaseq = read.table(rnaseq_file, header = T, sep = '\t')
all_names = c()
for (chrm in chrms){
  genes = as.character(rnaseq[rnaseq$chromosome==chrm,"Ensembl_ID"])  
  all_names = c(all_names, genes)
  output_name = paste0(outputfolder, '/chrm',chrm, '_genes.rds')
  saveRDS(genes, output_name)
}
write.table(all_names, paste0(outputfolder, "/all_gene_names.txt"), row.names = F, col.names = F, quote = F, sep = '\n')

