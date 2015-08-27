
# ## relative to working directory!
# path_veg = 'data/composition_v0.3.csv'                                        # composition model results
# veg = read.table(file=file.path(getwd(), path_veg),    sep=',', row.names=NULL, header=TRUE)


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

####################################################################################################
# chunk: plot composition estimates
####################################################################################################
suff=paste0('comp_', suff_figs)
plot_data_maps(r_comp, centers=centers_comp/rescale, taxa=taxa, ages, N_comp, K, T, thresh=0.5, limits, suff=suff, save_plots=save_plots)
plot_data_maps_binned(r_comp, centers=centers_comp/rescale, taxa=taxa, ages, N_comp, K, T, breaks, limits, suff=suff, save_plots=save_plots)

# ####################################################################################################
# # chunk: residuals
# ####################################################################################################
# 
# r_pls = matrix(nrow=nrow(y_veg), ncol=ncol(y_veg))
# 
# for (i in 1:nrow(y_veg)){
#   
#   if (sum(y_veg[i,]) > 0){
#     r_pls[i, ] = y_veg[i,]/sum(y_veg[i,])
#   } else {
#     r_pls[i, ] = rep(0, ncol(y_veg))
#   }
# }
# 
# # measure of distance; remove MI:LP first
# dist_mat=rdist(as.matrix(centers_pls), as.matrix(centers_comp)/rescale)
# dist_mat[dist_mat < 1e-6] = 0
# 
# foo=unlist(apply(dist_mat, 2, function(x) any(x == 0)))
# 
# bar = r_comp[foo, ]
# 
# 
# resids = sqrt(rowSums((r_pls - bar)^2))
# 
# dat = data.frame(resids=resids, x=centers_pls$x*rescale, y=centers_pls$y*rescale)
# 
# p <- ggplot() + geom_tile(data=dat, aes(x=x, y=y, fill=resids)) + 
#   scale_fill_gradientn(colours=tim.colors(), name="Residuals") + coord_fixed() #+
# #scale_x_continuous(limits$xlims*1000) + scale_y_continuous(limits$ylims*1000)
# p <- add_map_albers(p, map_data=us.fort, limits)
# p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.3)), strip.text.x = element_text(size = rel(1.0)))
# print(p)
# 
# # try binning?
# resids_taxon = (r_pls - bar)#^2
# colnames(resids_taxon) = taxa
# dat_taxon = data.frame(resids_taxon, x=centers_pls$x*rescale, y=centers_pls$y*rescale)
# dat_taxon_melt = melt(dat_taxon, id=c('x', 'y'))
# 
# 
# p <- ggplot() + geom_tile(data=dat_taxon_melt, aes(x=x, y=y, fill=value)) + 
#   scale_fill_gradientn(colours=tim.colors(), name="Residuals") + coord_fixed() #+
# #scale_x_continuous(limits$xlims*1000) + scale_y_continuous(limits$ylims*1000)
# p <- p + facet_grid(~variable)
# p <- add_map_albers(p, map_data=us.fort, limits)
# p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.3)), strip.text.x = element_text(size = rel(1.0)))
# print(p)
# 
# # probably should normalize somehow...
# dat_taxon_melt = dat_taxon_melt[which(abs(dat_taxon_melt$value)>0.5), ]
# p <- ggplot() + geom_tile(data=dat_taxon_melt, aes(x=x, y=y, fill=value)) + 
#   scale_fill_gradientn(colours=tim.colors(), name="Residuals") + coord_fixed() #+
# #scale_x_continuous(limits$xlims*1000) + scale_y_continuous(limits$ylims*1000)
# p <- p + facet_grid(~variable)
# p <- add_map_albers(p, map_data=us.fort, limits)
# p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.3)), strip.text.x = element_text(size = rel(1.0)))
# print(p)
# 
# ####################################################################################################
# # chunk: see if we can compare preds on 3by with comp
# ####################################################################################################
# # measure of distance; remove MI:LP first
# dist_mat=rdist(as.matrix(centers_veg), as.matrix(centers_comp)/rescale)
# dist_mat[dist_mat < 1e-6] = 0
# 
# foo=unlist(apply(dist_mat, 1, function(x) which.min(x)))
# 
# bar = r_comp[foo, ]
# centers_comp_sub = centers_comp[foo,]
# 
# 
# resids = sqrt(rowSums((r_mean - bar)^2))
# 
# dat = data.frame(resids=resids, x=centers_comp_sub[,1], y=centers_comp_sub[,2])
# 
# p <- ggplot() + geom_tile(data=dat, aes(x=x, y=y, fill=resids)) + 
#   scale_fill_gradientn(colours=tim.colors(), name="Residuals") + coord_fixed() #+
# #scale_x_continuous(limits$xlims*1000) + scale_y_continuous(limits$ylims*1000)
# p <- add_map_albers(p, map_data=us.fort, limits)
# p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.3)), strip.text.x = element_text(size = rel(1.0)))
# print(p)
# 
# ####################################################################################################
# # chunk: residuals
# ####################################################################################################
# 
# # N_cores and centers_polU broken for split domain data....
# idx.keep  = c(1,2,length(ages)/2,T)
# # plot_core_locations(y, centers_polU, centers_pls, ages, limits)
# plot_core_locations_select(y, centers_pol, centers_pls, ages, idx.keep, limits, suff=suff_figs, fpat=subDir)

# # ####################################################################################################
# # # chunk: plot predicted and observed proportions in same frame
# # ####################################################################################################
# # suff4=paste(suff, '_compare', sep='')
# # 
# # plot_both_maps(r_mean, y_veg, centers=centers_pls, taxa=taxa, ages, N, K, T, thresh=1, suff=suff3, save_plots=save_plots)



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