##########################################
# For general HiC intrachromosomal networks
##########################################
#rm(list=ls())
library(igraph, warn.conflicts = F)
library(data.table)
library(Matrix)
source('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/analysis/change_neighbor_mat.R')
#chrms = as.character(seq(1,23,1))
#chrms = c(chrms,'X')
#chromosome 1
#chrms = c('1')


read_genes<-function(genename_folder, chrms){
  all = c()
  for (chrm in chrms){
    genes = readRDS(paste0(genename_folder, '/chrm', chrm,'_genes.rds'))
    all = c(all, genes)
  }
  return(all)
}

delete_neighbor_edges<-function(homefolder, net, chrm){
  # read neighborhood network directly from hard link
  print(homefolder)
  linear_neighbor_cutoff = readRDS(paste0(homefolder, '/linear_neighbors/data/cutoff_10kb/chrm',chrm,'_neighbors_trans.rds'))
  linear_net_cutoff<-graph_from_adjacency_matrix(linear_neighbor_cutoff,weighted = TRUE, 
                                                 mode = "undirected")
  intersection = intersection(net, linear_net_cutoff, keep.all.vertices = F)
  if (vcount(intersection)==vcount(net)){
    #make sure we are not deleting nodes
    newnet = difference(net, intersection)  
    return(newnet)
  }
  else{
   stop("Error creating intersection graph in deleting neighboring edges from HiC edges\n")
  }
}

add_singletons<-function(homefolder, genename_folder, net, chrms){
  # read gene names of each chrmosome directly from hard link
  genes = read_genes(genename_folder, chrms)
  singletons = genes[!(genes%in%V(net)$name)]
  newnet = add_vertices(net, length(singletons), name=singletons)
  newnet<-simplify(newnet, remove.multiple = T, remove.loops = T)
  return(newnet)
}

read_edge_list_from_hiccups<-function(folder, chrms){
  edge_list = data.frame(t(c(NA,NA)))
  colnames(edge_list)<-c('Gene1','Gene2')
  for (chrm in chrms){
    edge_list_chrm = tryCatch(read.table(paste0(folder,'edge_list_chrm',chrm,'.txt')),error=function(e) NULL)
    tryCatch(colnames(edge_list_chrm)<-c('Gene1','Gene2'), error=function(e) NULL)
    edge_list = rbind(edge_list, edge_list_chrm)
  }
  edge_list = edge_list[-1,]
  return(as.matrix(edge_list))
}

read_edge_list_from_genepairs<-function(infolder, chrm, quant, method){
  #three methods: according to p-value, q-value and hard threshold on mean
  #for those gene_pair that do not have repetition, 
    #add them if the mean(from one value) is above the 
  result = as.data.table(read.table(paste0(infolder, '/','/intra_chr', chrm,'_10kb_quantile', 
                                           quant,'.genepairs'), header = T, sep='\t'))
  edgelist1 = as.matrix(t(data.frame(strsplit(as.character(result[result[,direct_mean],pair]), split='_'))))
  edgelist2 = as.matrix(t(data.frame(strsplit(as.character(result[result[,direct_median],pair]), split='_'))))
  edgelist3 = as.matrix(t(data.frame(strsplit(as.character(result[result[,direct_max],pair]), split='_'))))
  edgelist4 = as.matrix(t(data.frame(strsplit(as.character(result[result[,direct_min],pair]), split='_'))))

  if (method=="mean"){
    return(list(el = edgelist1))
  } else if (method=="median"){
    return(list(el = edgelist2))
  } else if (method=="max"){
    return(list(el = edgelist3))
  } else if (method =="min"){
    return(list(el = edgelist4))
  } else if (method =="all"){
    return(list(el = list(edgelist1, edgelist2, edgelist3, edgelist4)))
  }
}


read_edge_list_from_genepairs_inter<-function(infolder, chrm1, chrm2, quant, method){
  #three methods: according to p-value, q-value and hard threshold on mean
  #for those gene_pair that do not have repetition, 
    #add them if the mean(from one value) is above the 
  result = as.data.table(read.table(paste0(infolder, '/','/inter_chr', chrm1, '_chr', chrm2, '_10kb_quantile', 
                                           quant,'.genepairs'), header = T, sep='\t'))
  edgelist1 = as.matrix(t(data.frame(strsplit(as.character(result[result[,direct_mean],pair]), split='_'))))
  edgelist2 = as.matrix(t(data.frame(strsplit(as.character(result[result[,direct_median],pair]), split='_'))))
  edgelist3 = as.matrix(t(data.frame(strsplit(as.character(result[result[,direct_max],pair]), split='_'))))
  edgelist4 = as.matrix(t(data.frame(strsplit(as.character(result[result[,direct_min],pair]), split='_'))))

  if (method=="mean"){
    return(list(el = edgelist1))
  } else if (method=="median"){
    return(list(el = edgelist2))
  } else if (method=="max"){
    return(list(el = edgelist3))
  } else if (method =="min"){
    return(list(el = edgelist4))
  } else if (method =="all"){
    return(list(el = list(edgelist1, edgelist2, edgelist3, edgelist4)))
  }
}


output_graph<-function(homefolder, genename_folder, edge_list, add_singles, delete_neighbors, chrms){
  #print(edge_list)
  if (dim(edge_list)[1]==0) return(NULL)
  else {
    net<-graph_from_edgelist(edge_list,directed = F)
    net0<-simplify(net, remove.multiple = TRUE, remove.loops = TRUE,
                 edge.attr.comb = igraph_opt("edge.attr.comb"))
    cat("Number of nodes", vcount(net0),"\n")
    cat("Number of edges", ecount(net0),"\n")
    if (add_singles){
      cat("Adding singletons...\n")
      net1 = add_singletons(homefolder, genename_folder, net0, chrms)
      cat("New number of nodes", vcount(net1),"\n")
      cat("New number of edges", ecount(net1),"\n")
      if (delete_neighbors){
	if (length(chrms)>1){
	  stop("Cannot delete neighbors with INTER-chromosomal networks\n")
	}
        cat("Deleting neighboring edges...\n")
        net2 = delete_neighbor_edges(homefolder, net1, chrms)
        cat("New number of nodes", vcount(net2),"\n")
        cat("New number of edges", ecount(net2),"\n")
      }else {
        net2 = net1
      }
    }else {
      if (delete_neighbors){
        stop("Cannot delete neighbors if singletons not added\n")
      }
      net2 = net0
    }
    return(net2)
  }
}
 
output_adjacency_matrix<-function(genename_folder, net, chrms, added_singles){
  genes = read_genes(genename_folder, chrms)
  #check if the order of genes in the network (mat) agrees with the saved one
  if (!added_singles){
    cat("Singletons were not added to the network\n")
  }
  else {
    if (length(genes)!=vcount(net)){
      stop("Adjacency matrix does not contain all genes in the chromosome\n")
    }
    else {
      matched = match(V(net)$name, genes)
      if (sum(is.na(matched))>0){
        stop("Some nodes in the network not matched to genes\n")
      }
      else {
        net1 = permute(net, matched) #order of genes in the adjacency matrix is the same as in the genename_folder 
        matched1 = match(V(net1)$name, genes)
      }
    }
  }
  mat0 = as_adjacency_matrix(net1,type='both')
  mat1 = change_mat(mat0)
  return(mat1)
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

#output_y_inter<-function(genename_folder, chrm1, chrm2, rnaseq_files){
#  y_output = c()
#  genes1 = read_genes(genename_folder, chrm1)
#  genes2 = read_genes(genename_folder, chrm2)
#  #rnaseq_files should be a vector of filepaths
#  for (f in rnaseq_files){
#    rnaseq = read.table(f,header=T)
#    #match y to have the same order as the gene list in genes (neighborhood)
#    y1<-as.numeric(rnaseq[match(genes1, rnaseq[,1]),2])
#    y2<-as.numeric(rnaseq[match(genes2, rnaseq[,1]),2])
#    y = c(y1, y2)
#    y_output = cbind(y_output, y)
#  }
#  return(y_output)
#}

save_mat<-function(folder, mat, name){
  #mat need to be already converted
  saveRDS(mat, paste0(folder,'/', name, '_neighbors_trans.rds'))
}

save_y<-function(folder, y, name){
  saveRDS(y, paste0(folder,'/', name, '_y.rds'))
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
