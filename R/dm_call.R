############################################
# R wrapper for double metropolis
############################################

#' main pmrf
#'
#' @useDynLib bioMRF double_metropolis_cont
#' @export
dm_call_wrapper<-function(total_iter,N, y, mat_dgcmatrix, vars, bounds_a, bounds_e, bounds_t, inis, wInis){
	if (!is.numeric(total_iter) || !is.numeric(y) || !is.numeric(N) || !is.numeric(vars) || !is.numeric(bounds_a) || !is.numeric(bounds_e) || !is.numeric(bounds_t) || !is.numeric(inis) || !is.numeric(wInis)){
		stop("input data not numeric\n")
	}
	else if (class(mat_dgcmatrix)!='dgCMatrix') {
		stop("Please input neighborhood matrix in 'dgCMatrix' class.\n Use as(yourmatrix, 'dgCMatrix') in the Matrix package\n")
	}
  print("Hello")
  att = attributes(mat_dgcmatrix)
  dim = att$Dim[1]
  print(is.integer(dim))
  val = att$x
  row_ind = att$i
  col_ptr = att$p
	ret = .Call(double_metropolis_cont, T_in = as.integer(total_iter),
		    N_in = as.integer(N), y_in = as.double(y) ,
		    dim_in = as.integer(dim),
		    val_in = as.double(val), row_ind_in = as.integer(row_ind),
		    col_ptr_in = as.integer(col_ptr), vars_in =as.double(vars),
		    bounds_alpha = as.double(bounds_a),
		    bounds_eta = as.double(bounds_e),
		    bounds_tau2 = as.double(bounds_t),
		    initials = as.double(inis), wInitials = as.double(wInis))
	return(list(w = ret[[1]], alpha = ret[[2]], eta = ret[[3]],
		    tau2 = ret[[4]], jump_count = ret[[5]]))
}




