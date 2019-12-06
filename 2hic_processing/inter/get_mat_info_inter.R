#!/usr/bin/env Rscript
######################################
# For general inter HiC networks
########################################
rm(list=ls())
source('./__get_edge_list.R')
args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("Three arguments must be supplied\n", call.=FALSE)
} 


homedir = args[2]
genename_folder = args[3]

key = args[1]
#key should be in the format of "$chrm"_"$quant"_"$method"_TRUE_"$del"
#key = 'chrm1_0.9_all_TRUE_FALSE'
cat(key, "\n")
keys = strsplit(key, '_')[[1]]
chrm_long = keys[1]
chrm = substr(chrm_long, 5,6)
quant = as.numeric(keys[2])
method = keys[3]
add_singles = as.logical(keys[4])
del_neighbors = as.logical(keys[5])

info.ecount = c()
info.degrees_median = c()
info.degrees_mean = c()
info.degrees_sd = c()
info.percent_linear = c()

if (method=='all'){
  #read in unique neighbor network
  linear_neighbor_adj = readRDS(paste0(homedir, '/linear/',chrm_long,'_neighbors_trans.rds'))
  linear_net<-graph_from_adjacency_matrix(linear_neighbor_adj, weighted = TRUE, mode = "undirected")

  edgelist_list = read_edge_list_from_genepairs_inter(paste0(homedir,'/genepairs/'), chrm, quant, method)
  cat("read edgelist\n")
  method_names = c("mean", "median", "max", "min")
  i = 1
  for (el in edgelist_list$el){
    net = output_graph(homefolder, genename_folder, el, add_singles, del_neighbors, chrm)
    info.ecount[i] = ecount(net)
    degrees = degree(net)
    info.degrees_median[i] = median(degrees)
    info.degrees_mean[i] = mean(degrees)
    info.degrees_sd[i] = sd(degrees)
    inter_with_linear = intersection(linear_net, net, keep.all.vertices = F)
    cat("\n ratio of edges from linear neighbors in edges from HiC\n")
    r3 = ecount(inter_with_linear)/ ecount(net)
    cat(r3)
    cat("\n ratio of edges from HiC in edges from linear neighbors\n")
    r4 = ecount(inter_with_linear)/ ecount(linear_net)
    cat(r4)
    cat("\n")
    info.percent_linear[i] = r4
    i = i + 1
  }  
} else {
  stop("Use 'all' as method in input key\n")
}

table = cbind(info.ecount, info.degrees_median, info.degrees_mean, info.degrees_sd, info.percent_linear)
name = paste0('chrm', chrm, '_', quant, '_', method,'_', add_singles, '_', del_neighbors)
saveRDS(table, paste0(homedir, '/info/', name, '_info.rds'))
