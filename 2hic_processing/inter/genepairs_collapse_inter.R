#!/usr/bin/env Rscript

#########################################################
#Updated 20190620: varying threshold for each chromosome
#based on quantile
#########################################################
rm(list=ls())
chrs = as.character(seq(1,22,1))
chrs = c(chrs, 'X')
args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("At least three arguments must be supplied: homedir, quant, resolution. Optional argument: chrs\n", call.=FALSE)
} else if (length(args)==3) {
  chrms = chrs
} else if (length(args)>3){
  chrms = c(args[4])
}
print(args)
###
setwd(args[1])
quant = as.numeric(args[2])
resolution = args[3]
###

library(data.table)
compare_chromosomes<-function(chrm1, chrm2){
  if (chrm1=='X') return(FALSE)
  else if (chrm2=='X') return(TRUE)
  else if (as.numeric(chrm1)< as.numeric(chrm2)) return(TRUE)
  else return(FALSE)
}
join_gene_pair<-function(genepair){
  return(paste(as.character(sort(unlist(genepair))), collapse='_'))
}

join_gene_pair2<-function(gene1, gene2){
  return(paste(sort(c(as.character(gene1), as.character(gene2))), collapse='_'))
}
compute_values = function(values, thres){
  #print(values)
  #num = as.numeric(levels(values))[values]
  num = as.numeric(values)
  #print('####')
  #print(num)
  mean = mean(num)
  sd = sd(num)
  ##updated 20190621: do not use p-values
  #p=-1
  #if (length(num)>5) {
  #  tryCatch({p = t.test(num, alternative = 'greater', mu = thres)$p.value}, 
  #           error = function(e) {print(num); print(e)})
  #print(p)
  #}
  #return(list(mean=mean, sd=sd, pval = p, len = length(num)))
  max = max(num)
  min = min(num)
  median = median(num)
  len = length(num)
  return(list(mean= mean, sd = sd, max = max, min = min, median = median, len = len))
}
compute_threshold<-function(values, quantile){
  #values should be all pairs of (mapped) interactions in this chromosome
  num = as.numeric(values)
  thres = quantile(num, probs=c(quantile))
  return(thres)
}
cat("quantile", quant, "\n")
for (chrm1 in chrms){
  for (chrm2 in chrms){
    if (compare_chromosomes(chrm1, chrm2)){
      cat("Chromosome1", chrm1, "\n")
      cat("Chromosome2", chrm2, "\n")
      data.read<-read.table(file = paste0('./norms/inter_chr', chrm1,'_chr',chrm2,'_',resolution,'.norm'), header=T)
      dt = as.data.table(data.read)
      dt[,pair:=join_gene_pair2(gene1, gene2), by = seq_len(nrow(dt))]
      dt<-dt[,c(4,3)]
      #newrows=apply(dat[,1:2],1, join_gene_pair)
      #newdat =cbind(newrows, dat$inter_norm)
      colnames(dt)<-c('pair', 'val')
      thres = compute_threshold(dt$val, quant)
      cat("begin collapsing...\n")
      result = dt[, compute_values(val,thres), by = pair]
      cat("\n###number of total gene pairs\n")
      cat(nrow(result))
      #cat("\n###number of total interactions with more than one bin\n")
      #cat(sum(!is.na(result$pval)))
      cat("\n###number of total interactions with le five\n")
      cat(sum(result$len<=5))
      #cat("\n###number of interactions with significant pvalue\n")
      #cat(sum(result$pval<0.05, na.rm=T))
      #cat("\n###number of interactions with significant qvalue\n")
      #cat(sum(q_value<0.05, na.rm=T))
      cat("\n###number of interactions with mean interaction greater than treshold\n")
      direct_mean = result$mean>thres
      cat(sum(direct_mean))
      cat("\n###number of interactions with median interaction greater than treshold\n")
      direct_median = result$median>thres
      cat(sum(direct_median))
      cat("\n###number of interactions with max interaction greater than treshold\n")
      direct_max = result$max>thres
      cat(sum(direct_max))
      cat("\n###number of interactions with min interaction greater than treshold\n")
      direct_min = result$min>thres
      cat(sum(direct_min))
      cat("\n")
      newresult = cbind(result, direct_mean, direct_median, direct_max, direct_min)
      write.table(newresult, file = paste0('./genepairs/inter_chr', chrm1,'_chr', chrm2,'_',resolution,'_quantile', quant,'.genepairs'), quote=F, row.names = F, sep='\t')
	  }
	}
}
