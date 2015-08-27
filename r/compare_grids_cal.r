library(Rcpp)
library(inline)
library(ggplot2)
library(rstan)
library(reshape)
library(fields)

source('r/utils/pred_plot_funs.r')
source('r/utils/pred_helper_funs.r')
source('r/utils/build_mu_g.r')
source('r/read_stanbin.r')
source('r/mugp.r')

# edit this file to process different runs
source('r/runs_cal_3by.r')
runs3by = runs
source('r/runs_cal_1by.r')
runs1by = runs

# where to put the figures
subDir <- paste('figures/compare_grids')
create_figure_path(subDir)

# source('data/comp_data_12taxa_mid_ALL_v0.3.rdata')

for (run1by in runs1by){
  
  suff_dat  = run$suff_dat
  suff_fit  = run$suff_fit
  suff_figs = run$suff_figs
  
  run3by = runs3by[[which(sapply(runs3by, function(x) x$suff_figs) %in% suff_figs)]]

  
  # load the data and posterior draws
  load(paste0('r/dump/', run1by$suff_dat, '.rdata'))
  post_dat = load_stan_output(suff_fit)
  
  # hack for now, fix later
  post_dat$par_names = post_dat$par_names[7:length(post_dat$par_names)]
  
  process_out = build_r(post_dat, N, T, K)

  r_pred = process_out$r
  r_mean = matrix(NA, nrow=N*T, ncol=K)
  niter  = dim(r_pred)[3]
  
  for (i in 1:(N*T)){
    r_mean[i,]       = rowSums(r_pred[i,,])/niter
  }
  
  limits = get_limits(centers_pls)
  
  ####################################################################################################
  # chunk: plots 
  ####################################################################################################
  # suff1=paste(suff_fit, '_props', sep='')
  suff = paste0('props_', suff_figs)
  
  # plot_pred_maps(r_mean, centers_veg, taxa=taxa, ages, N, K, T, thresh=0.5, 
                  #limits, type='prop', suff=suff1,  save_plots=save_plots)
  
  # suff1.1=paste(suff_fit, '_props_select', sep='')
  plot_pred_maps_select(r_mean, centers_veg, taxa=taxa, ages, N, K, T, thresh=0.5, 
                        limits, type='prop', suff=suff,  save_plots=save_plots)
  
  breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
  # p_binned <- plot_pred_maps_binned(r_mean, centers_veg, breaks, taxa, ages, N, K, T, limits, suff=suff1, save_plots=save_plots)
  p_binned <- plot_pred_maps_binned_select(r_mean, centers_veg, breaks, taxa, ages, N, K, T, 
                                           limits, suff=suff_figs, save_plots, fpath=subDir)
  
#   ####################################################################################################
#   # chunk: plot raw pls
#   ####################################################################################################
#   plot_data_maps(y_veg, centers=centers_pls, taxa=taxa, ages, N_pls, K, T, thresh=0.5, limits, suff=suff_figs, save_plots=save_plots)
#   plot_data_maps_binned(y_veg, centers=centers_pls, taxa=taxa, ages, N_pls, K, T, breaks, limits, suff=suff_figs, save_plots=save_plots)
#   
}