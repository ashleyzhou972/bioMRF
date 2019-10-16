##################################
#Post MCMC analysis
#including test for eta
#get fitted values
##################################

get_jump_frequency<-function(ret, Ti, N){
  #ret is a list of returned values from c
  #Ti is total number of iterations
  #N is data size
  
  jumpcount = ret[[5]]
  rate_w = jumpcount[1]/(N*Ti)
  rate_alpha = jumpcount[2]/Ti
  rate_eta = jumpcount[3]/Ti
  rate_tau2 = jumpcount[4]/Ti
  
  return(list(w = rate_w, alpha = rate_alpha, eta = rate_eta, tau2 = rate_tau2))
}

get_jump_frequency_gaussian<-function(ret, Ti, N){
  #ret is a list of returned values from c
  #Ti is total number of iterations
  #N is data size
  
  jumpcount = ret[[4]]
  rate_alpha = jumpcount[1]/Ti
  rate_eta = jumpcount[2]/Ti
  rate_tau2 = jumpcount[3]/Ti
  
  return(list(alpha = rate_alpha, eta = rate_eta, tau2 = rate_tau2))
}


delete_burn_in <-function(ret, burn_in, neighbor){
  total_iter = length(ret$alpha)
  new_w = ret$w[,(burn_in:total_iter)]
  new_alpha = ret$alpha[burn_in:total_iter]
  new_eta = ret$eta[burn_in:total_iter]
  new_tau2 = ret$tau2[burn_in:total_iter]
  return(list(w=new_w, alpha=new_alpha, eta=new_eta, tau2=new_tau2, neighbor=neighbor))
}


print_param_estimates<-function(ret){
  #updated 20190228
  #Just take the 2.5% and 97.5% percentile!!!
  alpha_percentile = quantile(ret$alpha, probs=c(0.025, 0.975))
  cat(paste(round(mean(ret$alpha),5), " (", round(alpha_percentile[1],5), ", ", round(alpha_percentile[2],5), ") ", sep=""))
  cat("\n")
  eta_percentile = quantile(ret$eta, probs=c(0.025, 0.975))
  cat(paste(round(mean(ret$eta),5), " (", round(eta_percentile[1],5), ", ", round(eta_percentile[2],5), ") ", sep=""))
  cat("\n")
  tau2_percentile = quantile(ret$tau2, probs=c(0.025, 0.975))
  cat(paste(round(mean(ret$tau2),5), " (", round(tau2_percentile[1],5), ", ", round(tau2_percentile[2],5), ") ", sep=""))
  cat("\n")
  #options(digits=10)
}

print_param_estimates_latex<-function(ret, sim, vars, index){
  #updated 20190314
  alpha_percentile = quantile(ret$alpha, probs=c(0.025, 0.975))
  eta_percentile = quantile(ret$eta, probs=c(0.025, 0.975))
  tau2_percentile = quantile(ret$tau2, probs=c(0.025, 0.975))
  row1 = paste("Dataset", index, "& True &", sim$alpha, "&", sim$eta, "&", sim$tau2, "\\\\")
  alpha_est = paste0(round(mean(ret$alpha),5), " (", round(alpha_percentile[1],5), ",", round(alpha_percentile[2],5), ") ")
  eta_est = paste0(round(mean(ret$eta),5), " (", round(eta_percentile[1],5), ", ", round(eta_percentile[2],5), ") " )
  tau2_est = paste0(round(mean(ret$tau2),5), " (", round(tau2_percentile[1],5), ", ", round(tau2_percentile[2],5), ") ")
  row2 = paste(" & estimated &",alpha_est,"&",eta_est,"&", tau2_est,"\\\\", sep=" ")
  row3 = paste(" & variance &", vars[2], "&", vars[3], "&", vars[4], "\\\\")
  print=paste(row1, "\n", row2, "\n", row3, "\n")
  cat(print)
}

save_param_estimates_to_table<-function(ret, sim, vars, jump_count){
	results_table = matrix(NA, nrow=4, ncol=3)
	#colnames(results_table)<-c('alpha', 'eta', 'tau2')
	rownames(results_table)<-c('true', 'estimates', 'jumps', 'var')
	alpha_percentile = quantile(ret$alpha, probs=c(0.025, 0.975))
  	eta_percentile = quantile(ret$eta, probs=c(0.025, 0.975))
  	tau2_percentile = quantile(ret$tau2, probs=c(0.025, 0.975))
	alpha_est = paste0(round(mean(ret$alpha),5), " (", round(alpha_percentile[1],5), ",", round(alpha_percentile[2],5), ") ")
  	eta_est = paste0(round(mean(ret$eta),5), " (", round(eta_percentile[1],5), ", ", round(eta_percentile[2],5), ") " )
  	tau2_est = paste0(round(mean(ret$tau2),5), " (", round(tau2_percentile[1],5), ", ", round(tau2_percentile[2],5), ") ")
	results_table[1,]<-c(sim$alpha, sim$eta, sim$tau2)
	results_table[3,]<-c(jump_count$alpha, jump_count$eta, jump_count$tau2)
	#since we are not reporting jump_count for w, warning is printed if it is out of allowed range
	if (jump_count$w<0.2 || jump_count$w>0.6) {
		warnings(paste0("Jump count for w is out of bounds: ", jump_count$w, ".\n"))
	}
	results_table[4,]<-c(vars[2], vars[3], vars[4])
	results_table[2,]<-c(alpha_est, eta_est, tau2_est)
	return(results_table)
}


save_eta_estimates<-function(ret){
  eta_percentile = quantile(ret$eta, probs=c(0.025, 0.975))
  eta_est = paste0(round(mean(ret$eta),5), " (", round(eta_percentile[1],5), ", ", round(eta_percentile[2],5), ") " )
  return(eta_est)
}

save_param_estimates_to_table_3<-function(ret, sim, vars, jump_count){
	#3 digits instead of 5
	results_table = matrix(NA, nrow=4, ncol=3)
	#colnames(results_table)<-c('alpha', 'eta', 'tau2')
	rownames(results_table)<-c('true', 'estimates', 'jumps', 'var')
	alpha_percentile = quantile(ret$alpha, probs=c(0.025, 0.975))
  	eta_percentile = quantile(ret$eta, probs=c(0.025, 0.975))
  	tau2_percentile = quantile(ret$tau2, probs=c(0.025, 0.975))
	alpha_est = paste0(round(mean(ret$alpha),3), " (", round(alpha_percentile[1],3), ",", round(alpha_percentile[2],3), ") ")
  	eta_est = paste0(round(mean(ret$eta),3), " (", round(eta_percentile[1],3), ", ", round(eta_percentile[2],3), ") " )
  	tau2_est = paste0(round(mean(ret$tau2),3), " (", round(tau2_percentile[1],3), ", ", round(tau2_percentile[2],3), ") ")
	results_table[1,]<-c(sim$alpha, sim$eta, sim$tau2)
	results_table[3,]<-c(jump_count$alpha, jump_count$eta, jump_count$tau2)
	#since we are not reporting jump_count for w, warning is printed if it is out of allowed range
	if (jump_count$w<0.2 || jump_count$w>0.6) {
		warnings(paste0("Jump count for w is out of bounds: ", jump_count$w, ".\n"))
	}
	results_table[4,]<-c(vars[2], vars[3], vars[4])
	results_table[2,]<-c(alpha_est, eta_est, tau2_est)
	return(results_table)
}

get_fitted_y_full<-function(newret, neighbor, n){
  #each column is an iteration
  if (dim(neighbor)[1]!= n ){
    stop("Dimension n does not agree with neighborhood matrix\n")
  }
  total_iter = length(newret$alpha)
  alpha_mat = t(replicate(n, newret$alpha)) #n*T matrix
  tau2_mat = t(replicate(n, newret$tau2)) #n*T matrix
  mu =  alpha_mat + neighbor%*%(newret$w-alpha_mat)%*%diag(newret$eta)
  mean = mu + tau2_mat/2
  yhat = rowMeans(exp(mean))
  return(yhat)
}

posterior_predictive_y<-function(newret, n){
  total_iter = length(newret$alpha)
  #matrices are collapsed by column in R
  y<-matrix(rpois(n*total_iter, exp(newret$w)),byrow=F,ncol = total_iter)
  yhat = rowMeans(y)
  return(yhat)
}

RMSE<-function(y_ob, y_fitted) {
  naid = c(which(is.na(y_ob)), which(is.na(y_fitted)))
  if (length(naid)>0) print("NAs exist in data\n")
  rmse=sqrt(mean((y_ob-y_fitted)^2, na.rm=T))
  return(rmse)
}

PMSE<-function(y_ob, y_fitted) {
  naid = c(which(is.na(y_ob)), which(is.na(y_fitted)))
  if (length(naid)>0) print("NAs exist in data\n")
  mse=mean((y_ob-y_fitted)^2, na.rm=T)
  return(mse)
}

correlation<-function(y, sub_neighbor){
  neighbors = t(y%*%sub_neighbor)/rowSums(sub_neighbor)
  pearson = cor(y, neighbors, method = "pearson")
  return(pearson)
}


# RRSE<-function(y_ob, y_fitted, ret) {
#   #ret should be deleted of burn-in
#   fitted_mean = mean(exp(ret$alpha))
#   rrse = sqrt(sum((y_ob-y_fitted)^2,na.rm=T)/sum((y_fitted-fitted_mean)^2,na.rm=T))
#   return(rrse)
# }
# 
# RRSE1<-function(y_ob, y_fitted, alpha) {
#   #ret should be deleted of burn-in
#   fitted_mean = mean(exp(alpha))
#   rrse = sqrt(sum((y_ob-y_fitted)^2,na.rm=T)/sum((y_fitted-fitted_mean)^2,na.rm=T))
#   return(rrse)
# }
# 
# RRSE2<-function(y_ob, y_fitted) {
#   rrse = sqrt(sum((y_ob-y_fitted)^2,na.rm=T)/sum((y_ob-mean(y_ob, na.rm=T))^2,na.rm=T))
#   return(rrse)
# }

