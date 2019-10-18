#########################################################
#Updated 20190620: varying threshold for each chromosome
#based on quantile
#########################################################
rm(list=ls())
library(data.table)
setwd('/home/nzhou/hic/rao2014/IMR90_10kb/intra/by_gene')
####toggle
chrms = as.character(seq(1,22,1))
chrms = c(chrms, 'X')
quant = 0.90
#chrms = c("1")
resolution = '10kb'
###
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
for (chrm in chrms){
  cat("Chromosome", chrm, "\n")
  data.read<-read.table(file = paste0('./newnorms/intra_chr', chrm,'_',resolution,'.norm'), header=T)
  dt = as.data.table(data.read)
  dt[,pair:=join_gene_pair2(gene1, gene2), by = seq_len(nrow(dt))]
  dt<-dt[,c(4,3)]
  #newrows=apply(dat[,1:2],1, join_gene_pair)
  #newdat =cbind(newrows, dat$inter_norm)
  colnames(dt)<-c('pair', 'val')
  thres = compute_threshold(dt$val, quant)
  print(thres)
  cat("begin collapsing...\n")
  result = dt[, compute_values(val,thres), by = pair]
  #result[result$pval==-1,'pval']<-NA
  #q_value = p.adjust(result$pval, 'fdr')
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
  write.table(newresult, file = paste0('./genepairs/intra_chr', chrm,'_',resolution,'_quantile', quant,'.genepairs'), quote=F, row.names = F, sep='\t')
}




# compute_values = function(values, thres){
#   #print(values)
#   num = as.numeric(levels(values))[values]
#   #print('####')
#   #print(num)
#   mean = mean(num)
#   sd = sd(num)
#   if (length(num)==1){
#     p = NA
#   }
#   else {
#     tryCatch({p = t.test(num, alternative = 'greater', mu = thres)$p.value}, 
#              error = function(e) {print(num);print(e)}, finally={p = NA})
#   }
#   return(list(mean=mean, sd=sd,pval = p))
# }
# 
# result = aggregate(val~pair, data = newdat[1:100,], FUN = compute_values, thres=70)
# colnames(result)<-c('genepair', 'mean', 'sd', 'p_value')

#############################

