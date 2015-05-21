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
  save(process_out, file=paste0(subDir, '/process_out.rdata'))
  
  process_mean = build_mu_g(post_dat, rho, eta, T, K, d, d_inter, d_knots, od, mpp, mu0)  
  save(process_mean, file=paste0(subDir, '/process_mean.rdata'))
  
#   source('r/pred_plot.r')

}