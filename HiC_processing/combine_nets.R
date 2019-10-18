#############################################
##### Read in transformed adjacency matrices from 
##### both intra and inter HiC data
##### produce a combined NON-transformed adjacency matrix 
##### saved in sparse format 
##### Transformation needs to be done AFTER the whole matrix is constructed
#############################################
rm(list=ls())
library(Matrix)
source('/home/nzhou/hic/rao2014/script/get_edge_list.R')
datafolder_intra = '/home/nzhou/hic/rao2014/IMR90_10kb/intra/by_gene/data/'
datafolder_inter = '/home/nzhou/hic/rao2014/IMR90_10kb/inter/by_gene/data/'
chrms = as.character(seq(1,22,1))
chrms = c(chrms,'X')
#chrms = c(1,7,14,17,'X')
long_chrms = paste0("chrm", chrms)

##toggle
singleton = "TRUE"
del="FALSE"
method = "mean"
quan = "0.9"
##
get_block_inter<-function(inter_mat, index1, index2){
  # for each inter matrix, separate into four blocks. 
  # The two diagonal blocks should be empty(all zero)
  # The two off-diagonal blocks should be the symmetric (transposed of each other)
  # Return the off-diagonal block in the LOWER triangle
  dim1 = dims[index1]
  dim2 = dims[index2]
  off_lower = inter_mat[(dim1+1):(dim1+dim2), 1:dim1]
  off_upper = inter_mat[1:dim1,(dim1+1):(dim1+dim2)] 
  on_1 = inter_mat[1:dim1, 1:dim1]
  on_2 = inter_mat[(dim1+1):(dim1+dim2),(dim1+1):(dim1+dim2)]
  if (check_mat_zero(on_1) & check_mat_zero(on_2)) {
    if (sum(off_lower!=t(off_upper))==0) return(off_lower)  
    else stop("Lower triagle not the transpose of upper triangle\n");
  } else {
    stop("On diagonal not zero")
  }
}

check_mat_zero<-function(mat){
  if (sum(mat)==0) return(TRUE)
  else return(FALSE)
}
check_if_names_agree<-function(big, small, i, j){
  ## check if the colnames and rownames of the part of the (big) matrix 
  ## to be replaced
  ## is the same as the (small) matrix replacing it
  part = big[(blocks[j,1]:blocks[j,2]), (blocks[i,1]:blocks[i,2])]
  big_1 = colnames(part)
  big_2 = rownames(part)
  small_1 = colnames(small)
  small_2 = rownames(small)
  if (length(big_1)==length(small_1) & length(big_2)==length(small_2)){
    if (sum(big_1!=small_1)==0 & sum(big_2!=small_2)==0)
      return(TRUE)
    else return(FALSE)
  }else return(FALSE)
}
read_intra<-function(long_chrms){
  #read in all inTRAchromosomal matrices
  #combine them into a block diagonal matrix
  #size of this block diagonal matrix is the size of the final desired matrix
  #use read_inter() to add the off diagonal parts
  dims = rep(NA, length(long_chrms))
  blocks = matrix(NA, ncol = 2, nrow = length(long_chrms))
  intras = list()
  names = c()
  total_nnz = 0
  for (i in 1:length(long_chrms)){
    chrm1 = long_chrms[i]
    cat(chrm1, "\n")
    key1 = paste(chrm1, quan, method, singleton, del, sep = "_")
    intra = readRDS(paste0(datafolder_intra, key1,"_neighbors_trans.rds"))
    nnz = length(attributes(intra)$x)
    total_nnz = total_nnz + nnz
    cat("Number of nonzero values for intra", chrm1, ":", nnz, "\n")
    print(length(colnames(intra)))
    dims[i] = dim(intra)[1]
    if (i==1) blocks[i,]<-c(1, dims[i])
    else blocks[i,]<-c((blocks[i-1,2]+1), (blocks[i-1,2]+dims[i]))
    intras = c(intras, intra)
    names = c(names, colnames(intra))
  }
  big = bdiag(intras)
  colnames(big)<-names
  rownames(big)<-names
  return(list(mat = big, dims = dims, blocks = blocks, nnz = total_nnz))
}



read_inter<-function(long_chrms, big, dims, blocks){
  nnz_vec = c()
  total_nnz = 0
  for (i in 1:length(long_chrms)){
  chrm1 = long_chrms[i]
    for (j in 1:length(long_chrms)){
      chrm2 = long_chrms[j]
      if (compare_chrms(chrm1, chrm2)){
        cat(chrm1, chrm2, "\n")
        key12 = paste(chrm1, chrm2, quan, method, singleton, del, sep = "_")
        inter = readRDS(paste0(datafolder_inter, key12,"_neighbors_trans.rds"))
        nnz = length(attributes(inter)$x)
        nnz_vec = c(nnz_vec, nnz)
        total_nnz = total_nnz + nnz
        cat("Number of nonzero values for inter", chrm1, chrm2, ":",nnz , "\n")
        off_lower = get_block_inter(inter, i, j)
        if (check_if_names_agree(big, off_lower, i, j)){
          #replace
          big[(blocks[j,1]:blocks[j,2]), (blocks[i,1]:blocks[i,2])]<-off_lower
          big[(blocks[i,1]:blocks[i,2]), (blocks[j,1]:blocks[j,2])]<-t(off_lower)
        }
      }
    }
  }
  return(list(mat = big, nnz = total_nnz, nnz_vec = nnz_vec))
}


combine_y<-function(long_chrms){
  #y should follow the order in gene_names
  #which is preserved in the matrices
  y_all = c()
  for (i in 1:length(long_chrms)){
    chrm1 = long_chrms[i]
    y = readRDS(paste0(datafolder_intra, '/y/', chrm1,"_y.rds"))
    y_all = rbind(y_all, y)
  }
  return(y_all)
}

#####################################################################
##main
ret1 = read_intra(long_chrms)
big = ret1$mat
dims = ret1$dims
blocks = ret1$blocks
ret2 = read_inter(long_chrms, big, dims, blocks)
newbig = ret2$mat

#count total nonzero values as sanity check
length(attributes(newbig)$x)
ret1$nnz+ret2$nnz

#make all nonzero values one instead of ratio (for preprocessing)
att = attributes(newbig)
attributes(newbig)$x<-rep(1, length(att$x))

key_all = paste("ALL", quan, method, singleton, del, sep = "_")
datafolder_all = "/home/nzhou/hic/rao2014/IMR90_10kb/inter/by_gene/all_data/"
saveRDS(newbig, paste0(datafolder_all,key_all, "_neighbors.rds"))

y_all <-combine_y(long_chrms)
##check if length agrees
dim(y_all)
#dim(newbig)
saveRDS(y_all, paste0(datafolder_all, "/y/", 'ALL', "_y.rds"))



