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


  compute_distance <- function(comp, preds){
  
    ntaxa  = ncol(comp)
    ncells = nrow(comp)
  
    score = rep(0, ncells)
    for (n in 1:ncells){
#     score[n] = score[n] + sum((preds[n,] - comp[n,])^2)#sqrt(sum((preds[n,] - comp[n,])^2))
#       score[n] = score[n] + sqrt(sum((preds[n,] - comp[n,])^2))
     score[n] = score[n] + sqrt((preds[n,10] - comp[n,10])^2)
    }
  
    return(score)
  }

pls.raw = data.frame(read.table(file='../stepps-data/data/composition/composition_v0.3.csv' , sep=",", row.names=NULL, header=TRUE))
# pls.raw = data.frame(read.table(file='../stepps-data/data/composition/composition_v0.4.csv' , sep=",", row.names=NULL, header=TRUE))
colnames(pls.raw) = tolower(colnames(pls.raw))

# figure out why chestnut and atlantic.white.cedar are all NA!
pls.raw = pls.raw[,!(colnames(pls.raw) %in% c('chestnut', 'atlantic.white.cedar'))]

# conversion table
convert = read.table('data/dict-comp2stepps.csv', sep=',', row.names=1, header=TRUE)

# pull the subset of proportions
taxa.start.col = min(match(tolower(rownames(convert)), colnames(pls.raw)), na.rm=TRUE)

pls_dat  = pls.raw[,taxa.start.col:ncol(pls.raw)]
colnames(pls_dat) = as.vector(convert[match(colnames(pls_dat), tolower(rownames(convert))),1])
pls_dat_collapse  = sapply(unique(colnames(pls_dat)), 
                           function(x) rowSums( pls_dat[ , grep(x, names(pls_dat)), drop=FALSE]) )
comp_counts = data.frame(pls_dat_collapse[,sort(colnames(pls_dat_collapse))])
comp_meta   = pls.raw[,1:(taxa.start.col-1)]

comp_counts = comp_counts[which(comp_meta$region %in% c('michigan_north', 'wisconsin', 'minnesota')),]
comp_meta = comp_meta[which(comp_meta$region %in% c('michigan_north', 'wisconsin', 'minnesota')),]
colnames(comp_meta)[colnames(comp_meta) == 'region'] = 'state'

comp_meta = split_mi(comp_meta)
comp_counts = comp_counts[which(comp_meta$state2 %in% c('michigan:north', 'wisconsin', 'minnesota')),]
comp_meta = comp_meta[which(comp_meta$state2 %in% c('michigan:north', 'wisconsin', 'minnesota')),]

# edit this file to process different runs
source('r/runs_cal_3by.r')
runs3by = runs
source('r/runs_cal_1by.r')
runs1by = runs

# where to put the figures
subDir <- paste('figures/compare_grids')
create_figure_path(subDir)
# source('data/comp_data_12taxa_mid_ALL_v0.3.rdata')

# score_pred = data.frame(score=numeric(0), x=numeric(0), y=numeric(0), model=character(0))
score_comp = data.frame(score=numeric(0), x=numeric(0), y=numeric(0), model=character(0), grid=character(0))

for (run1by in runs1by){
  
  suff_dat  = run1by$suff_dat
  suff_fit  = run1by$suff_fit
  suff_figs = run1by$suff_figs
  
  run3by = runs3by[[which(sapply(runs3by, function(x) x$suff_figs) %in% suff_figs)]]

  # load the data and posterior draws
  load(paste0('r/dump/', run1by$suff_dat, '.rdata'))
  centers_veg_1by = centers_veg
  N_1by = N

  load(paste0('r/dump/', run3by$suff_dat, '.rdata'))
  centers_veg_3by = centers_veg
  N_3by = N

  grids = c('1by', '3by')
  for (grid in grids){
    
    #grid = '1by'

    N = get(paste0('N_', grid))
    run = get(paste0('run', grid))$suff_fit

    post_dat = load_stan_output(run)
    post_dat$par_names = post_dat$par_names[7:length(post_dat$par_names)]
    r_pred   = build_r(post_dat, N, T, K)$r

    r_mean = matrix(NA, nrow=N*T, ncol=K)
    niter  = dim(r_pred)[3]
  
    for (i in 1:(N*T)){
      r_mean[i,]       = rowSums(r_pred[i,,])/niter
    }
    
    assign(paste0('r_mean_', grid), r_mean)

  }
  
#   # match each 1by center with a 3by estimate
#   r_mean_3by_res = matrix(NA, nrow=N_1by*T, ncol=K)
#   centers_res    = matrix(NA, nrow=N_1by*T, ncol=2)
#   colnames(centers_res) = c('x', 'y')  
# 
#   for (i in 1:N_1by){
#       d_across = rdist(as.matrix(centers_veg_1by[i,]), as.matrix(centers_veg_3by))
#       r_mean_3by_res[i * T,] = r_mean_3by[which.min(d_across),]
#       centers_res[i,]    = as.matrix(centers_veg_1by[i,])
#   }
# 
#   score      = compute_distance(r_mean_1by, r_mean_3by_res)
#   score_pred = rbind(score_pred, data.frame(score=score, centers_res, model = suff_figs))

  # match each estimate with PLS data
  r_mean_1by_comp = matrix(NA, nrow=N_1by*T, ncol=K)
  r_mean_3by_comp = matrix(NA, nrow=N_1by*T, ncol=K)

  for (i in 1:N_1by){
      d_across = rdist(as.matrix(comp_meta[i,1:2]), as.matrix(centers_veg_1by)*1e6)
      r_mean_1by_comp[i * T,] = r_mean_1by[which.min(d_across),]
      
      d_across = rdist(as.matrix(comp_meta[i,1:2]), as.matrix(centers_veg_3by)*1e6)
      r_mean_3by_comp[i * T,] = r_mean_3by[which.min(d_across),]
  }

  score = compute_distance(comp_counts, r_mean_1by_comp)
  score_comp = rbind(score_comp, data.frame(score=score, comp_meta[,1:2], model = suff_figs, grid = '1by'))
  
  score = compute_distance(comp_counts, r_mean_3by_comp)
  score_comp = rbind(score_comp, data.frame(score=score, comp_meta[,1:2], model = suff_figs, grid = '3by'))

}



limits = get_limits(centers_veg_1by)

p <- ggplot() + geom_tile(data=score_pred, aes(x=x, y=y, fill=score)) + 
  scale_fill_gradientn(colours=tim.colors(), name="Distance") + coord_fixed() 
#p <- add_map_albers(p, map_data=us.fort, limits)
#p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.4)), 
#                            strip.text.x = element_text(size = rel(1.4)),
#                            legend.title=element_text(size = rel(1.1)),
#                            legend.text=element_text(size = rel(1.0)))
#p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# p <- p + facet_grid(model~.)
p <- p + facet_wrap(~model)
print(p)
fname = paste0(subDir, '/compare_grids.pdf')
ggsave(p, file=fname, width=14)
sys_str = paste("pdfcrop", fname, fname, sep=' ')
system(sys_str)

aggregate(. ~ model + grid, data=score_comp, FUN=sum)

limits = get_limits(centers_veg_1by)

p <- ggplot() + geom_tile(data=score_comp, aes(x=x, y=y, fill=score)) + 
  scale_fill_gradientn(colours=tim.colors(), name="Distance") + coord_fixed() 
#p <- add_map_albers(p, map_data=us.fort, limits)
#p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.4)), 
#                            strip.text.x = element_text(size = rel(1.4)),
#                            legend.title=element_text(size = rel(1.1)),
#                            legend.text=element_text(size = rel(1.0)))
#p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# p <- p + facet_grid(model~.)
p <- p + facet_grid(grid~model)
print(p)
fname = paste0(subDir, '/distance_pred_PLS_pine.pdf')
ggsave(p, file=fname, width=14)
sys_str = paste("pdfcrop", fname, fname, sep=' ')
system(sys_str)

 

#   suff   = paste0('compare_', suff_figs)
#   plot_pred_maps_select(score, centers_veg_1by, taxa=taxa, ages, N_1by, K, T, thresh=NA, 
#                         limits, type='process', suff=suff,  save_plots=save_plots)
  
#   ####################################################################################################
#   # chunk: plots 
#   ####################################################################################################
#   # suff1=paste(suff_fit, '_props', sep='')
#   suff = paste0('props_', suff_figs)
  
#   # plot_pred_maps(r_mean, centers_veg, taxa=taxa, ages, N, K, T, thresh=0.5, 
#                   #limits, type='prop', suff=suff1,  save_plots=save_plots)
  
#   # suff1.1=paste(suff_fit, '_props_select', sep='')
#   plot_pred_maps_select(r_mean, centers_veg, taxa=taxa, ages, N, K, T, thresh=0.5, 
#                         limits, type='prop', suff=suff,  save_plots=save_plots)
  
#   breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
#   # p_binned <- plot_pred_maps_binned(r_mean, centers_veg, breaks, taxa, ages, N, K, T, limits, suff=suff1, save_plots=save_plots)
#   p_binned <- plot_pred_maps_binned_select(r_mean, centers_veg, breaks, taxa, ages, N, K, T, 
#                                            limits, suff=suff_figs, save_plots, fpath=subDir)
  
# #   ####################################################################################################
# #   # chunk: plot raw pls
# #   ####################################################################################################
# #   plot_data_maps(y_veg, centers=centers_pls, taxa=taxa, ages, N_pls, K, T, thresh=0.5, limits, suff=suff_figs, save_plots=save_plots)
# #   plot_data_maps_binned(y_veg, centers=centers_pls, taxa=taxa, ages, N_pls, K, T, breaks, limits, suff=suff_figs, save_plots=save_plots)
# #   
# }