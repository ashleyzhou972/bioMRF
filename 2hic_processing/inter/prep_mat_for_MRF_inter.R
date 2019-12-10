#!/usr/bin/env Rscript
######################################
# For general intra HiC networks
########################################
rm(list=ls())
source('./__get_edge_list.R')
args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("Three arguments must be supplied\n", call.=FALSE)
} 

homedir = args[2]
genename_folder = args[3]
output_folder = paste0(homedir, '/data/')

key = args[1]
#key should be in the format of "$chrm1"_"$chrm2"_"$quant"_"$method"_TRUE_"$del"
cat(key, "\n")
keys = strsplit(key, '_')[[1]]
chrm_long1 = keys[1] #chrm_long contains the string "chrm"
chrm_long2 = keys[2]
chrm1 = substr(chrm_long1, 5,6) #chrm only contains the number of the chromosome
chrm2 = substr(chrm_long2, 5,6)
chrms = c(chrm1, chrm2)
quant = as.numeric(keys[3])
method = keys[4]
add_singles = as.logical(keys[5])
del_neighbors = as.logical(keys[6])
  
if ( compare_chrms(chrm1, chrm2)) {
  edgelist_list = read_edge_list_from_genepairs_inter(paste0(homedir, '/genepairs/'), chrm1, chrm2, quant, method)
  chrms = c(chrm1, chrm2)
  if (method=='all'){
    method_names = c("mean", "median", "max", "min")
    i = 1
    for (el in edgelist_list$el){
      net = output_graph(homefolder, genename_folder, el, add_singles, del_neighbors, chrms)
      mat = output_adjacency_matrix(genename_folder, net, chrms, add_singles)
      name = paste0('chrm', chrm1, '_', 'chrm', chrm2, '_', quant, '_', method_names[i],'_', add_singles, '_', del_neighbors)
      #print(name)
      save_mat(output_folder,mat, name)
      i = i + 1
    }  
    } else {
	  el = edgelist_list$el
	  net = output_graph(homefolder,genename_folder, el, add_singles, del_neighbors, chrms)
	  mat = output_adjacency_matrix(genename_folder, net, chrms, add_singles)
	  name = paste0('chrm', chrm1, '_', 'chrm', chrm2, '_', quant, '_', method,'_', add_singles, '_', del_neighbors)
	  save_mat(output_folder, mat, name)
    }
} else {
	cat("Invalid chromosome pair,", chrm1, chrm2, "\n")
}

