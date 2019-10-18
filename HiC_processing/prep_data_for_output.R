#####################################
# prepare presentation of data for intra model
#Two outputs: 
  #1. side by side boxplot of four methods, linear and loop
  #2. Table showing the eta intervals for above runs
#####################################
rm(list = ls())
source("/home/nzhou/hic/IMR90/work/MRF_HIC_GE/analysis/post_analysis.R")
chrms = c('1','7', '14', '17', 'X')
methods = c("mean", "median", "max", "min")
outfolder = "/home/nzhou/hic/IMR90/work/MRF_HIC_GE/real_by_chrm/out/new/"
retfolder = "/home/nzhou/hic/IMR90/work/MRF_HIC_GE/real_by_chrm/returns/new/"
infofolder = "/home/nzhou/hic/rao2014/IMR90_10kb/intra/by_gene/info/"
resultsfolder = "/home/nzhou/hic/rao2014/IMR90_10kb/results/"
linearfolder = "/home/nzhou/hic/IMR90/work/MRF_HIC_GE/real_by_chrm/returns/linear/cutoff_10kb/"
loopfolder = "/home/nzhou/hic/IMR90/work/MRF_HIC_GE/real_by_chrm/returns/loop/"
B = 400

#output tables
for (chrm in chrms){
  tab = c()
  for (del in c(TRUE, FALSE)){
    key0 = paste(paste0('chrm', chrm), 0.9, 'all', TRUE, del,sep = '_')
    info= readRDS(paste0(infofolder, '/', key0, '_info.rds'))
    eta_4 = c()
    for (method in methods){
      key = paste(paste0('chrm', chrm), 0.9, method, TRUE, del,sep = '_')
      ret = readRDS(paste0(retfolder, '/', 'ret_', key, '.rds'))
      newret = delete_burn_in(ret, B, NULL)
      eta_est = save_eta_estimates(newret)
      eta_4 = c(eta_4, eta_est)
    }
    tab = cbind(tab, info, eta_4)
  }
  newtab = tab[, c(8,1,9,2,10,3,11,4,12,5,13,6,14,7)]
  colnames(newtab)<-c("ecount.w", "ecount.wo", "deg_median.w", "deg_median.wo", "deg_mean.w", "deg_mean.wo",
                      "deg_sd.w", "deg_sd.wo", "percet_loop.w", "percent_loop.wo", "percent_linear.w",
                      "percent_linear.wo","eta.w", "eta.wo")
  write.table(paste0("Chromosome ", chrm, "\n"), file = paste0(resultsfolder,"/results_20190628.csv"), sep = ',',
              append = T, quote = T, row.names = F, col.names = F)
  write.table(newtab, file = paste0(resultsfolder,"/results_20190628.csv"), sep = ',',
             append = T, quote = T, row.names = F, col.names = T)
  write.table(paste0( "\n"), file = paste0(resultsfolder,"/results_20190628.csv"), sep = ',',
              append = T, quote = T, row.names = F, col.names = F)
}

B1 = 500
#output boxplots
adaptive = TRUE
#if adaptive the ylims adapt to data
set.seed(12)
for (chrm in chrms){
  eta_all = c()
  linear = readRDS(paste0(linearfolder, "/", "linear_ret_chrm", chrm, ".rds"))
  new_linear = delete_burn_in(linear, B1, NULL)
  linear_eta = sample(new_linear$eta, size = 1602)
  
  loop = readRDS(paste0(loopfolder, "/", "loop_retchrm", chrm, ".rds"))
  new_loop = delete_burn_in(loop, B1, NULL)
  loop_eta = sample(new_loop$eta, size = 1602)
  
  for (method in methods){
    #key0 = paste(paste0('chrm', chrm), 0.9, 'all', TRUE, del,sep = '_')
    #info= readRDS(paste0(infofolder, '/', key0, '_info.rds'))
    eta_4 = c()
    for (del in c(FALSE, TRUE)){
      key = paste(paste0('chrm', chrm), 0.9, method, TRUE, del,sep = '_')
      ret = readRDS(paste0(retfolder, '/', 'ret_', key, '.rds'))
      newret = delete_burn_in(ret, B, NULL)
      eta = newret$eta
      eta_all = cbind(eta_all, eta)
    }
  }
  eta_all = cbind(eta_all, loop_eta, linear_eta)
  colnames(eta_all)<-c("mean.w", "mean.wo", "median.w", "median.wo", "max.w", "max.wo", 
                      "min.w", "min.wo", "loop", "linear")
  if (adaptive){
    png(paste0(resultsfolder, "/boxplot_adpt_chrm", chrm, ".png"), width = 900)
    boxplot(eta_all, outline = F, col = c(rep(0,8),2,4), ylab = 'estimated eta')
  }
  else{
    png(paste0(resultsfolder, "/boxplot_noadpt_chrm", chrm, ".png"), width = 900)
    boxplot(eta_all, outline = F, ylim = c (-2.5, 2.5), col = c(rep(0,8),2,4), ylab = 'estimated eta')
  }
  abline(h = 0, col = 2)
  dev.off()
}