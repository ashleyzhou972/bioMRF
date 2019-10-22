library(PhiMRF)
adj_trans = transform_small(irregular_lattice)
n = dim(irregular_lattice)[1]
y = simulate_y(n = n, adj_mat = adj_trans, alpha = 2, eta = 2, tau2 = 3, m = 2, M = 2000)
bounds_e = get_eta_param_space_small(adj_trans)
ret = pmrf(total_iter = 1000, N = n , y, adj_trans, c(1, 0.5, NA, 2), bounds_e)
get_jump_frequency(ret, 1000, 1000)
new_ret = delete_burn_in(ret, 200)
print_param_estimates(ret)
