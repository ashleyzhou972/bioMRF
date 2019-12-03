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
#key should be in the format of "$chrm"_"$quant"_"$method"_TRUE_"$del"
cat(key, "\n")
keys = strsplit(key, '_')[[1]]
chrm = substr(keys[1], 5,6)
quant = as.numeric(keys[2])
method = keys[3]
add_singles = as.logical(keys[4])
del_neighbors = as.logical(keys[5])

  
edgelist_list = read_edge_list_from_genepairs(paste0(homedir, '/genepairs/'), chrm, quant, method)
if (method=='all'){
  method_names = c("mean", "median", "max", "min")
  i = 1
  for (el in edgelist_list$el){
    net = output_graph(homefolder, genename_folder, el, add_singles, del_neighbors, chrm)
    mat = output_adjacency_matrix(genename_folder, net, chrm, add_singles)
    name = paste0('chrm', chrm,  '_', quant, '_', method_names[i],'_', add_singles, '_', del_neighbors)
    #print(name)
    save_mat(output_folder,mat, name)
    i = i + 1
  }  
} else {
  el = edgelist_list$el
  net = output_graph(homefolder,genename_folder, el, add_singles, del_neighbors, chrm)
  mat = output_adjacency_matrix(genename_folder, net, chrm, add_singles)
  name = paste0('chrm', chrm,  '_', quant, '_', method,'_', add_singles, '_', del_neighbors)
  save_mat(output_folder, mat, name)
}

