library(Rcpp)
library(inline)
library(ggplot2)
library(rstan)
library(reshape)
library(fields)

source('r/utils/pred_plot_funs.r')
source('r/utils/pred_helper_funs.r')
source('r/read_stanbin.r')

# where to put the figures
subDir <- paste("figures/", suff_fit, sep='')

create_figure_path(subDir)

load(paste('r/dump/', suff_dat, '.rdata', sep=''))

W = K-1

# for full model
N_pars = 3*W + 1

trace_plot_pars(post, N_knots, T, N_pars, taxa=taxa, suff=suff, save_plots=save_plots)
trace_plot_mut(post, N_knots, T, N_pars, mean_type='other', suff=suff, save_plots=save_plots)

## fix this!!!
#trace_plot_knots(fit, N_knots, T, K, N_pars=N_pars, suff=suff, save_plots=save_plots)

# write mean vals to file
# sink(sprintf('%s/%s/summary.txt', wd, path_figs1), type='output')
sink(sprintf('%s/summary.txt', subDir), type='output')
print('The taxa modelled are:')
taxa
cat('\n')
print('Summary of posterior parameter vals:')
get_quants(post, N_pars)
sink()

adj = get_corrections(post, rho, eta, T, K, d_inter, d_knots)
adj_s = adj$adj_s
adj_t = adj$adj_t

####################################################################################################
# chunk: compute and plot proportion chains
####################################################################################################
# process_out = build_props_full(post, rho, eta, T, K, d, d_inter, d_knots, od, mpp, mu0)
process_out = build_r(post, rho, eta, T, K, d, d_inter, d_knots, od, mpp, mu0)
save(process_out, file=paste0(subDir, '/process_out.rdata'))
# rm(post)

# r_pred = process_out$r
# g      = process_out$g
# rm(process_out)


process_mean = build_mu_g(post, rho, eta, T, K, d, d_inter, d_knots, od, mpp, mu0, res)
# build_mu_g <- function(post, rho, eta, T, K, d, d_inter, d_knots, od, mpp, mu0, res){

# mu_g     = process_mean$mu_g
# Halpha_t = process_mean$Halpha_t
# Halpha_s = process_mean$Halpha_s

save(process_mean, file=paste0(subDir, '/process_mean.rdata'))


process_out$g - process_mean$mu_g


# niter = dim(post[,1,])[1]
# mean_Halpha_t = array(NA, dim=c(W, T-1, niter))
# 
# for (k in 1:W){
#   for (t in 1:(T-1)){
#     mean_Halpha_t[k, t, ] = colSums(Halpha_t[,(k-1)*(T-1)+t,])/N
#   }
# }
# 
# mu_t = get_mut(post, N_pars, W)
# 
# pdf(file=paste0(subDir, '/trace_Halpha_t.pdf'), width=8, height=12)
# for (k in 1:W){
#     par(mfrow=c(5,2))
#     par(oma=c(0,0,2,0))
#     for (t in 1:(T-1)){
#       plot(mean_Halpha_t[k,t,], type="l", ylab=paste0('mean_Halpha_t[', k, ',', t, ']'))
#       abline(h=mean(mean_Halpha_t[k,t,]), col="blue")
#       abline(h=quantile(mean_Halpha_t[k,t,],probs=0.025), col='blue', lty=2)
#       abline(h=quantile(mean_Halpha_t[k,t,],probs=0.975), col='blue', lty=2)
#   #   lines(mu_t[t+1,1,], col="blue")
# #     plot(sum_Halpha_t[k,t,], type="l", ylab=paste0('sum_Halpha[', k, ',', t, ']'))
#   
#     }
# #     title(main=taxa[k], outer=TRUE)
# }
# dev.off()


# niter = dim(post[,1,])[1]
# mean_Halpha_s = array(NA, dim=c(W, niter))
# 
# for (k in 1:W){
#     mean_Halpha_s[k, ] = colSums(Halpha_s[,k,])/N
# }
# 
# mu = get_mu(post, W)
# pdf(file=paste0(subDir, '/trace_mu.pdf'), width=8, height=12)
# par(mfrow=c(5,2))
# par(oma=c(0,0,2,0))
# for (k in 1:W){
#   plot(mu[k,], type="l", ylab=paste0('mu[', k, ']'))
#   abline(h=mean(mu[k,]), col="blue")
#   abline(h=quantile(mu[k,],probs=0.025), col='blue', lty=2)
#   abline(h=quantile(mu[k,],probs=0.975), col='blue', lty=2)
#   #   lines(mu_t[t+1,1,], col="blue")
#   #     plot(sum_Halpha_t[k,t,], type="l", ylab=paste0('sum_Halpha[', k, ',', t, ']'))
#   #     title(main=taxa[k], outer=TRUE)
# }
# dev.off()
# 
# pdf(file=paste0(subDir, '/trace_Halpha_s.pdf'), width=8, height=12)
# par(mfrow=c(5,2))
# par(oma=c(0,0,2,0))
# for (k in 1:W){
#   plot(mean_Halpha_s[k,], type="l", ylab=paste0('mean_Halpha_s[', k, ']'))
#   abline(h=mean(mean_Halpha_s[k,]), col="blue")
#   abline(h=quantile(mean_Halpha_s[k,],probs=0.025), col='blue', lty=2)
#   abline(h=quantile(mean_Halpha_s[k,],probs=0.975), col='blue', lty=2)
#     #   lines(mu_t[t+1,1,], col="blue")
#     #     plot(sum_Halpha_t[k,t,], type="l", ylab=paste0('sum_Halpha[', k, ',', t, ']'))
#   #     title(main=taxa[k], outer=TRUE)
# }
# dev.off()


# ######################################
# # check if need further OD
# mu_t       = get_mut(post, N_pars=W)
# sum_Halpha = process_out$sumHalpha
# 
# pdf(file=paste0(subDir, '/compare_mu_Halpha.pdf'), width=8, height=8)
# for (k in 1:W){
#   par(mfrow=c(3,2))
#   par(oma=c(0,0,2,0))
#   for (t in 1:3){
#     plot(mu_t[t,k,], type="l", ylab=paste0('mu_t[', t, ',', k, ']'))
#   #   lines(mu_t[t+1,1,], col="blue")
#     plot(sum_Halpha[t,k,], type="l", ylab=paste0('sum_Halpha[', t, ',', k, ']'))
# 
#   }
#   title(main=taxa[k], outer=TRUE)
# }
# dev.off()
# 
# mu_t_int = colSums(mu_t[,1,])
# 
# mu = get_mu(post, N_pars=W)
# 
# mu_t_int - mu[,1]
# 
# pdf(file=paste0(subDir, '/compare_mus.pdf'), width=8, height=6)
# for (k in 1:W){
#   mu_t_int = colSums(mu_t[,k,])
#   par(mfrow=c(2,1))
#   plot(mu_t_int, type='l', ylab=paste0('int(mu_t[t, ', k, '])'))
#   plot(mu[,k], type='l',  ylab=paste0('mu[', k, ']'), col='black')
# }
# dev.off()
#######################################

# # trace_plot_process(mu_g, suff='mu_g', save_plots=save_plots)
# trace_plot_process(r_pred, suff='r', save_plots=save_plots)
# trace_plot_process(g, suff='g', save_plots=save_plots)
# 
# 
# r_mean = matrix(NA, nrow=N*T, ncol=K)
# g_mean = matrix(NA, nrow=N*T, ncol=W)
# niter = dim(r_pred)[3]
# 
# for (i in 1:(N*T)){
#   r_mean[i,] = rowSums(r_pred[i,,])/niter
#   g_mean[i,] = rowSums(g[i,,])/niter
# }
# 
# rm(g)
# rm(r_pred)
# 
# 
# limits = get_limits(centers_pls)
# ####################################################################################################
# # chunk: plot predicted distributions
# ####################################################################################################
# suff1=paste(suff_fit, '_props', sep='')
# 
# # plot_pred_maps(r_mean, centers_veg, taxa=taxa, ages, N, K, T, thresh=0.5, limits, type='prop', suff=suff1,  save_plots=save_plots)
# 
# suff1.1=paste(suff_fit, '_props_select', sep='')
# plot_pred_maps_select(r_mean, centers_veg, taxa=taxa, ages, N, K, T, thresh=0.5, limits, type='prop', suff=suff1.1,  save_plots=save_plots)
# 
# breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
# # p_binned <- plot_pred_maps_binned(r_mean, centers_veg, breaks, taxa, ages, N, K, T, limits, suff=suff1, save_plots=save_plots)
# p_binned <- pred_maps_binned_select(r_mean, centers_veg, breaks, taxa, ages, N, K, T, limits, suff1, save_plots, fpath=subDir)
# 
# # FIX THIS SO it is ONLY select locations...
# # ####################################################################################################
# # # chunk: predicted process maps
# # ####################################################################################################
# # suff2=paste(suff_fit, '_process', sep='')
# # 
# # plot_pred_maps(g_mean, centers_veg, taxa=taxa, ages, N, K-1, T, thresh=NA, limits, type='process', suff=suff2,  save_plots=save_plots)
# 
# ####################################################################################################
# # chunk: plot observed proportions
# ####################################################################################################
# suff3=paste(suff_fit, '_data', sep='')
# 
# plot_data_maps(y_veg, centers=centers_pls, taxa=taxa, ages, N_pls, K, T, thresh=0.5, limits, suff=suff3, save_plots=save_plots)
# 
# suff4=paste(suff_fit, '_data_binned', sep='')
# plot_data_maps_binned(y_veg, centers=centers_pls, taxa=taxa, ages, N_pls, K, T, breaks, limits, suff=suff4, save_plots=save_plots)
# 
# 
# # N_cores and centers_polU broken for split domain data....
# idx.keep  = c(1,2,length(ages)/2,T)
# # plot_core_locations(y, centers_polU, centers_pls, ages, limits)
# plot_core_locations_select(y, centers_pol, centers_pls, ages, idx.keep, limits, fpat=subDir)
# 
# # # ####################################################################################################
# # # # chunk: plot predicted and observed proportions in same frame
# # # ####################################################################################################
# # # suff4=paste(suff, '_compare', sep='')
# # # 
# # # plot_both_maps(r_mean, y_veg, centers=centers_pls, taxa=taxa, ages, N, K, T, thresh=1, suff=suff3, save_plots=save_plots)
