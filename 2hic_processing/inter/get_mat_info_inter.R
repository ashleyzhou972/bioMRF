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
#key should be in the format of "$chrm1"_"chrm2"_"$quant"_"$method"_TRUE_"$del"
#key = 'chrm1_chrm2_0.9_all_TRUE_FALSE'
cat(key, "\n")
keys = strsplit(key, '_')[[1]]
chrm_long1 = keys[1]
chrm_long2 = keys[2]
chrm1 = substr(chrm_long1, 5,6)
chrm2 = substr(chrm_long2, 5,6)
quant = as.numeric(keys[3])
method = keys[4]
add_singles = as.logical(keys[5])
del_neighbors = as.logical(keys[6])

info.ecount = c()
info.degrees_median = c()
info.degrees_mean = c()
info.degrees_sd = c()

if (compare_chrms(chrm1, chrm2)) {
  chrms = c(chrm1, chrm2)
  if (method=='all'){
    edgelist_list = read_edge_list_from_genepairs_inter(paste0(homedir,'/genepairs/'), chrm1, chrm2, quant, method)
    cat("read edgelist\n")
    method_names = c("mean", "median", "max", "min")
    i = 1
    for (el in edgelist_list$el){
      net = output_graph(homefolder, genename_folder, el, add_singles, del_neighbors, chrms)
      info.ecount[i] = ecount(net)
      degrees = degree(net)
      info.degrees_median[i] = median(degrees)
      info.degrees_mean[i] = mean(degrees)
      info.degrees_sd[i] = sd(degrees)
      cat("\n")
      i = i + 1
    }  
  } else {
    stop("Use 'all' as method in input key\n")
  }
}

table = cbind(info.ecount, info.degrees_median, info.degrees_mean, info.degrees_sd)
name = paste0('chrm', chrm1, '_', 'chrm', chrm2, '_', quant, '_', method,'_', add_singles, '_', del_neighbors)
saveRDS(table, paste0(homedir, '/info/', name, '_info.rds'))
