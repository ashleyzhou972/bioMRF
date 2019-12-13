rm(list=ls())
source('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/analysis/change_neighbor_mat.R')
datafolder = '/home/nzhou/hic/rao2014/IMR90_10kb/intra/by_gene/data/'
chrms = as.character(seq(1,22,1))
chrms = c(chrms, 'X')
#chrms = c("1")

#TADs: inter-only
# TAD genes with interTAD edges only
# the adjacency matrix is simply the element-wise difference 
# between the original TADs matrix and the intraTAD matrix
# need to align the order of the genes!!!
outputfolder = "/home/nzhou/hic/rao2014/IMR90_10kb/intra/by_gene/TAD_inter_data/"
for (chrm in chrms){
  chrm_long = paste0("chrm", chrm)
  allTAD = readRDS(paste0(datafolder, "../TAD_data/TADs_", chrm_long, "_neighbors_trans.rds"))
  gene_names_all = colnames(allTAD)
  intraTAD = readRDS(paste0(datafolder, "../TAD_intra_data/TADs_intra_", chrm_long, "_neighbors_trans.rds"))
  gene_names_intra = colnames(intraTAD)
  #### re-order the intraTAD matrix according to the order of the allTAD matrix
  order = match(gene_names_all, gene_names_intra)
  new_intraTAD = intraTAD[order, order]
  # check if all the row indices of the new_intraTAD sparse matrix is a subset 
  # of all the row indices of the allTAD sparse matrix
  if (!all(attributes(new_intraTAD)$i%in%attributes(allTAD)$i)){
    stop("Edges in intraTAD is not a subset of edges in allTAD\n")
  }
  ##un-transform back to basic adj matrix with nonzero values = 1
  attributes(new_intraTAD)$x<-rep(1, length(attributes(new_intraTAD)$x))
  attributes(allTAD)$x<-rep(1, length(attributes(allTAD)$x))
  interTAD_notrans = allTAD - new_intraTAD
  #coerse the resulting zeros out of the nonzero values
  #@TODO this is expensive
  interTAD_notrans_sparse = as(as.matrix(interTAD_notrans), "dgCMatrix")
  attributes(interTAD_notrans_sparse)$x<-rep(1, length(attributes(interTAD_notrans)$x))
  interTAD_trans = change_mat(interTAD_notrans_sparse)
  ## sum of nnz of interTAD and intraTAD should be nnz of allTAD
  if (length(attributes(allTAD)$x)!=length(attributes(interTAD_trans)$x)+length(attributes(new_intraTAD)$x)){
    stop("Sum of nnz in intraTAD and interTAD not nnz of allTAD\n")
  }
  
  # Read in y from intraTAD and check for colnames before save as interTAD
  y = readRDS(file = paste0(outputfolder, "../TAD_data/y/", paste0("TADs_", chrm_long), "_y.rds"))
  if (sum(rownames(y)!=colnames(interTAD_trans))>0){
    stop("Gene names for y and interTAD do not agree\n")
  }
  saveRDS(interTAD_trans, file = paste0(outputfolder, "/", paste0("TADs_inter_", chrm_long), "_neighbors_trans.rds"))
  saveRDS(y, file = paste0(outputfolder, "/y/", paste0("TADs_inter_", chrm_long), "_y.rds"))
}
