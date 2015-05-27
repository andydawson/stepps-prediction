library(Rcpp)
library(inline)
library(ggplot2)
library(rstan)
library(reshape)
library(fields)

source('r/utils/pred_plot_funs.r')
source('r/utils/pred_helper_funs.r')
source('r/read_stanbin.r')

# edit this file to process different runs
source('r/runs.r')

for (run in runs){
  suff_dat  = run$suff_dat
  suff_fit  = run$suff_fit
  suff_figs = run$suff_figs
  #   source('r/pred_process_full_test.r')
  
  # where to put the figures
  subDir <- paste("figures/", suff_fit, sep='')
  create_figure_path(subDir)
  
  # load the data and posterior draws
  load(paste0('r/dump/', suff_dat, '.rdata'))
  post_dat = load_stan_output(suff_fit)
  
  process_out = build_r(post_dat, T, K)
#   save(process_out, file=paste0(subDir, '/process_out.rdata'))
  
#   save(post_dat, rho, eta, T, K, d, d_inter, d_knot, od, mpp, mu0, file='build_mu_g_test.rdata')
  process_mean = build_mu_g(post_dat, rho, eta, T, K, d, d_inter, d_knots, od, mpp, mu0)  
#   save(process_mean, file=paste0(subDir, '/process_mean.rdata'))

  # for full model
  N_pars = 3*(K-1) + 1
  write_par_vals(post_dat, taxa, subDir, N_pars)
  
#   trace_plot_pars(post, N_knots, T, N_pars, taxa=taxa, suff=suff, save_plots=save_plots)
#   trace_plot_mut(post, N_knots, T, N_pars, mean_type='other', suff=suff, save_plots=save_plots)
#   
#   ## fix this!!!
#   #trace_plot_knots(fit, N_knots, T, K, N_pars=N_pars, suff=suff, save_plots=save_plots)

  # figure out nugget size
  adj = get_corrections(post_dat, rho, eta, T, K, d_inter, d_knots)
  adj_s = adj$adj_s
  adj_t = adj$adj_t

  t(apply(adj_s + adj_t, 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975))))

  diff_g_mug = process_out$g - process_mean$mu_g
  
  summary_diff_g_mug = matrix(nrow=K-1, ncol=3)
  for (k in 1:(K-1)){
    summary_diff_g_mug[k,] = quantile(abs(as.vector(diff_g_mug[,k,])), probs=c(0.025, 0.5, 0.975))
  }
  
    source('r/pred_plot.r')
  
}