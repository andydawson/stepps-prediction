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
source('r/runs_cal_1by.r')
source('r/runs_cal_5by.r')

# source('data/comp_data_12taxa_mid_ALL_v0.3.rdata')

for (run in runs){
  
  # post_process_run(run)
  # gc()
  
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
  
  # hack for now, fix later
  post_dat$par_names = post_dat$par_names[7:length(post_dat$par_names)]
  
  # N       = nrow(d_inter)
  # N_knots = ncol(d_inter)
  
  process_out = build_r(post_dat, N, T, K)
  print('Built r')
  
  process_mean = build_mu_g_no_time(post_dat, rho, eta, T, K, d, d_inter, d_knots, od, mpp, mu0) 
  
  print('Built mu_g')
  
  # for full model
  W = K-1
  N_pars = 3*(K-1) + 1
  write_par_vals(post_dat, taxa, subDir, N_pars)
  
  #   trace_plot_pars(post, N_knots, T, N_pars, taxa=taxa, suff=suff, save_plots=save_plots)
  #   trace_plot_mut(post, N_knots, T, N_pars, mean_type='other', suff=suff, save_plots=save_plots)
  #   
  #   ## fix this!!!
  #   #trace_plot_knots(fit, N_knots, T, K, N_pars=N_pars, suff=suff, save_plots=save_plots)
  
  #   # figure out nugget size
  #   adj = get_corrections(post_dat, rho, eta, T, K, d_inter, d_knots)
  #   adj_s = adj$adj_s
  #   adj_t = adj$adj_t
  
  #   t(apply(adj_s + adj_t, 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975))))
  
  diff_g_mug = process_out$g - process_mean$mu_g
  
  summary_diff_g_mug = matrix(nrow=K-1, ncol=3)
  for (k in 1:(K-1)){
    summary_diff_g_mug[k,] = quantile(abs(as.vector(diff_g_mug[,k,])), probs=c(0.025, 0.5, 0.975))
  }
  
  summary_diff_g_mug
  
  ###############################################################################################################
  # chunk: load processed output
  ###############################################################################################################
  r_pred   = process_out$r
  g        = process_out$g
  mu_g     = process_mean$mu_g
  mu       = process_mean$mu
  Halpha_s = process_mean$Halpha_s
  
  r_mu_g = build_r_from_mu_g(mu_g, N, T, K)
  
  niter = dim(g)[3]
  
  mean_Halpha_s = array(NA, dim=c(W, niter))
  for (k in 1:W){
    mean_Halpha_s[k, ] = colSums(Halpha_s[,k,])/N
  }
  
  pdf(file=paste0(subDir, '/trace_mu.pdf'), width=8, height=12)
  par(mfrow=c(5,2))
  par(oma=c(0,0,2,0))
  for (k in 1:W){
    plot(mu[,k], type="l", ylab=paste0('mu[', k, ']'))
    abline(h=mean(mu[,k]), col="blue")
    abline(h=quantile(mu[,k],probs=0.025), col='blue', lty=2)
    abline(h=quantile(mu[,k],probs=0.975), col='blue', lty=2)
    #   lines(mu_t[t+1,1,], col="blue")
    #     plot(sum_Halpha_t[k,t,], type="l", ylab=paste0('sum_Halpha[', k, ',', t, ']'))
    #     title(main=taxa[k], outer=TRUE)
  }
  dev.off()
  
  pdf(file=paste0(subDir, '/trace_Halpha_s.pdf'), width=8, height=12)
  par(mfrow=c(5,2))
  par(oma=c(0,0,2,0))
  for (k in 1:W){
    plot(mean_Halpha_s[k,], type="l", ylab=paste0('mean_Halpha_s[', k, ']'))
    abline(h=mean(mean_Halpha_s[k,]), col="blue")
    abline(h=quantile(mean_Halpha_s[k,],probs=0.025), col='blue', lty=2)
    abline(h=quantile(mean_Halpha_s[k,],probs=0.975), col='blue', lty=2)
    #   lines(mu_t[t+1,1,], col="blue")
    #     plot(sum_Halpha_t[k,t,], type="l", ylab=paste0('sum_Halpha[', k, ',', t, ']'))
    #     title(main=taxa[k], outer=TRUE)
  }
  dev.off()
  
  
  # trace_plot_process(mu_g, suff='mu_g', save_plots=save_plots)
  trace_plot_process(r_pred, suff='r', save_plots=save_plots)
  trace_plot_process(g, suff='g', save_plots=save_plots)
  
  
  r_mean = matrix(NA, nrow=N*T, ncol=K)
  r_mu_g_mean = matrix(NA, nrow=N*T, ncol=K)
  g_mean = matrix(NA, nrow=N*T, ncol=W)
  mu_g_mean = matrix(NA, nrow=N*T, ncol=W)
  niter = dim(r_pred)[3]
  
  for (i in 1:(N*T)){
    r_mean[i,]       = rowSums(r_pred[i,,])/niter
    r_mu_g_mean[i,]  =  rowSums(r_mu_g[i,,])/niter
    g_mean[i,]       = rowSums(g[i,,])/niter
    mu_g_mean[i,]    = rowSums(mu_g[i,,])/niter
  }
  
  print('Computed mean pred vals')
  
  limits = get_limits(centers_pls)
  ####################################################################################################
  # chunk: plot predicted distributions
  ####################################################################################################
  # suff1=paste(suff_fit, '_props', sep='')
  suff = paste0('props_', suff_figs)
  
  # plot_pred_maps(r_mean, centers_veg, taxa=taxa, ages, N, K, T, thresh=0.5, limits, type='prop', suff=suff1,  save_plots=save_plots)
  
  # suff1.1=paste(suff_fit, '_props_select', sep='')
  plot_pred_maps_select(r_mean, centers_veg, taxa=taxa, ages, N, K, T, thresh=0.5, limits, type='prop', suff=suff,  save_plots=save_plots)
  
  breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
  # p_binned <- plot_pred_maps_binned(r_mean, centers_veg, breaks, taxa, ages, N, K, T, limits, suff=suff1, save_plots=save_plots)
  p_binned <- plot_pred_maps_binned_select(r_mean, centers_veg, breaks, taxa, ages, N, K, T, limits, suff=suff_figs, save_plots, fpath=subDir)
  
  print('Plotted predictions')
  
  ####################################################################################################
  # chunk: predicted process maps
  ####################################################################################################
  # suff2=paste(suff_fit, '_process', sep='')
  # 
  # plot_pred_maps(g_mean, centers_veg, taxa=taxa, ages, N, K-1, T, thresh=NA, limits, type='process', suff=suff2,  save_plots=save_plots)
  
  suff=paste0('process_', suff_figs)
  plot_pred_maps_select(g_mean, centers_veg, taxa=taxa, ages, N, K-1, T, thresh=NA, limits, type='process', suff=suff,  save_plots=save_plots)
  
  suff=paste0('mug_', suff_figs)
  plot_pred_maps_select(r_mu_g_mean, centers_veg, taxa=taxa, ages, N, K-1, T, thresh=NA, limits, type='prop', suff=suff,  save_plots=save_plots)
  
  print('Plotted process')
  
  ####################################################################################################
  # chunk: plot raw pls
  ####################################################################################################
  plot_data_maps(y_veg, centers=centers_pls, taxa=taxa, ages, N_pls, K, T, thresh=0.5, limits, suff=suff_figs, save_plots=save_plots)
  plot_data_maps_binned(y_veg, centers=centers_pls, taxa=taxa, ages, N_pls, K, T, breaks, limits, suff=suff_figs, save_plots=save_plots)
  
}