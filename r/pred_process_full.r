library(ggplot2)
library(rstan)
library(reshape)
library(fields)

source('r/utils/pred_plot_funs.r')
source('r/utils/pred_helper_funs.r')

suff_dat = '12taxa_457cells_77knots_0to2000ypb_umwE_3by_v0.3'
suff_fit = '12taxa_457cells_77knots_0to2000ypb_umwE_3by_od_mpp_full'#_mut2'

# where to put the figures
subDir <- paste("figures/", suff_fit, sep='')
save_plots = TRUE
suff = ''

# od         = TRUE
# bt         = TRUE
# mpp        = FALSE
# mut        = FALSE
# save_plots = TRUE

od         = TRUE
bt         = TRUE
mpp        = TRUE
mut        = FALSE
save_plots = TRUE

create_figure_path(subDir)

load(paste('r/dump/pred_data_', suff_dat, '.rdata', sep=''))

if (!file.exists(paste0('output/', suff_fit,'.rdata'))){
  fname = sprintf('output/%s.csv', suff_fit)
  system(sprintf('r/fixup.pl %s', fname)) # is this broken now?
  fit = read_stan_csv(fname)
  post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
  save(post, file=paste0('output/', suff_fit,'.rdata'))
  rm(fit)
} else {
  load(paste0('output/', suff_fit,'.rdata'))
}

W = K-1

# summary(fit)$summary # take a while

# for full model
N_pars = 3*W + 1

trace_plot_pars(post, N_knots, T, N_pars, taxa=taxa, suff=suff, save_plots=save_plots)

# fix this!!!
# trace_plot_knots(fit, N_knots, T, K, N_pars=K, suff=suff, save_plots=save_plots)

#fix this!!!
#trace_plot_mut(post, N_knots, T, N_pars, mean_type='other', suff=suff, save_plots=save_plots)

# col_substr = substr(names(summary(fit)$summary[,'mean']),1,2)
# 
# tau   = summary(fit)$summary[,'mean'][which(col_substr == 'ta')]
# mu    = summary(fit)$summary[,'mean'][which(col_substr == 'mu')]
# alpha = summary(fit)$summary[,'mean'][which(col_substr == 'al')]
# 
# ksi     = summary(fit)$summary[,'mean'][which(col_substr == 'ks')]
# omega   = summary(fit)$summary[,'mean'][which(col_substr == 'om')]
# 
# summary(fit)$summary[,'mean'][1:(2+W-1)]

npars = 3*W+1

# write mean vals to file
# sink(sprintf('%s/%s/summary.txt', wd, path_figs1), type='output')
sink(sprintf('%s/summary.txt', subDir), type='output')
print('The taxa modelled are:')
taxa
cat('\n')
print('Summary of posterior parameter vals:')
get_quants(post, npars)
sink()

adj = get_corrections(post, rho, eta, T, K, d_inter, d_knots)
adj_s = adj$adj_s
adj_t = adj$adj_t

####################################################################################################
# chunk: compute and plot proportion chains
####################################################################################################
process_out = build_props_full(post, rho, eta, T, K, d, d_inter, d_knots, od, mpp, mut)
# save(process_out, file='r/pred/dump/process_out.rdata')
# rm(post)


# M = diag(N*T)
# M = M - P
# 
# foo = M %*% Halpha[,1,1]
# 
# # Halpha


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

# load(file='r/pred/dump/process_out.rdata')
r_pred = process_out$r
g = process_out$g
# trace_plots_props(r_pred, suff=suff, save_plots=save_plots)
rm(process_out)


r_mean = matrix(NA, nrow=N*T, ncol=K)
g_mean = matrix(NA, nrow=N*T, ncol=W)
niter = dim(r_pred)[3]

for (i in 1:(N*T)){
  r_mean[i,] = rowSums(r_pred[i,,])/niter
  g_mean[i,] = rowSums(g[i,,])/niter
}

rm(g)
rm(r_pred)


limits = get_limits(centers_pls)
####################################################################################################
# chunk: plot predicted distributions
####################################################################################################
suff1=paste(suff_fit, '_props', sep='')

# plot_pred_maps(r_mean, centers_veg, taxa=taxa, ages, N, K, T, thresh=0.5, limits, type='prop', suff=suff1,  save_plots=save_plots)

suff1.1=paste(suff_fit, '_props_select', sep='')
plot_pred_maps_select(r_mean, centers_veg, taxa=taxa, ages, N, K, T, thresh=0.5, limits, type='prop', suff=suff1.1,  save_plots=save_plots)

breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
p_binned <- plot_pred_maps_binned(r_mean, centers_veg, breaks, taxa, ages, N, K, T, limits, suff=suff1, save_plots=save_plots)
p_binned <- pred_maps_binned_select(r_mean, centers_veg, breaks, taxa, ages, N, K, T, limits, suff1, save_plots, fpath=subDir)
####################################################################################################
# chunk: predicted process maps
####################################################################################################
suff2=paste(suff_fit, '_process', sep='')

plot_pred_maps(g_mean, centers_veg, taxa=taxa, ages, N, K-1, T, thresh=NA, limits, type='process', suff=suff2,  save_plots=save_plots)

####################################################################################################
# chunk: plot observed proportions
####################################################################################################
suff3=paste(suff_fit, '_data', sep='')

plot_data_maps(y_veg, centers=centers_pls, taxa=taxa, ages, N_pls, K, T, thresh=0.5, limits, suff=suff3, save_plots=save_plots)

suff4=paste(suff_fit, '_data_binned', sep='')
plot_data_maps_binned(y_veg, centers=centers_pls, taxa=taxa, ages, N_pls, K, T, breaks, limits, suff=suff4, save_plots=save_plots)


# N_cores and centers_polU broken for split domain data....

# plot_core_locations(y, centers_polU, centers_pls, ages, limits)
plot_core_locations_select(y, centers_pol, centers_pls, ages, limits, fpath=subDir)


# plot5<-grid.arrange(arrangeGrob(c, p_binned, nrow=1), 
#                      widths=c(10))
# plot5<-grid.arrange(c, p_binned, heights=c(10,10), nrow=1)


# 
# # objects are c, d, p_binned
# # Get the widths
# g <- ggplotGrob(p_binned)
# keep <- !grepl("strip-top", g$layout$name)
# g$grobs <- g$grobs[keep]
# g$layout <- g$layout[keep, ]
# 
# # gA <- ggplot_gtable(ggplot_build(p_binned))
# # gB <- ggplot_gtable(ggplot_build(c))
# 
# c_grob <- ggplotGrob(c)
# c_grob$heights <- g$heights
# 
# # gtable_show_layout(g)
# gf <- gtable_add_cols(g, unit(5, "cm"), pos=0)
# gf <- gtable_add_cols(g, unit(0.25, "lines"), pos=0)
# gf <- gtable_add_grob(gf, c_grob,
#                      t = 3, l=1, b=nrow(g), r=1)
# g <- gtable_add_rows(g, unit(8, "cm"), pos=0)
# # g <- gtable_add_rows(g, g$heights[1], pos=0)
# g <- gtable_add_grob(g, ggplotGrob(d),
#                      t = 1, l=2, b=1, r=ncol(g)-2)
# grid.newpage()
# pdf('figures/pred/test_gtable.pdf', width=12, height=16)
# grid.draw(g)
# dev.off()
# # 
# 
# 
# ng = nullGrob()
# g <- gtable(widths = unit(c(5, 7), "in"),  # need lcm(3,2)=6 for the matrix rows
#             heights = unit(c(5, 7), "in"))
# gtable_show_layout(g)
# 
# # # Add a grob:
# # rect <- rectGrob(gp = gpar(fill = "black"))
# # a <- gtable_add_grob(a, rect, 1, 1)
# # a
# # plot(a)
# 
# g<- gtable_add_grob(g, c_grob, 2, 1)
# g<- gtable_add_grob(g, ggplotGrob(d), 1, 2)
# g<- gtable_add_grob(g, ggplotGrob(p_binned), 2, 2)
# 
# grid.newpage()
# grid.draw(g)


# x1 <- freeze_panels(d, draw=TRUE, width=unit(1,"in"),height=unit(1,"in"))
# x2 <- freeze_panels(c, draw=TRUE, width=unit(1,"in"),height=unit(1,"in"))
# x3 <- freeze_panels(p_binned, draw=TRUE, width=unit(1,"in"),height=unit(1,"in"))
# 
# 
# 
# 
# 
# grid.arrange(arrangeGrob(ng, x1, nrow=1, widths=c(0.5,3.5)),
#              arrangeGrob(x2, x3, nrow=1, widths=c(0.5,3.5)),
#              nrow=2)
# 
# grid.arrange(ng, x1, x2, x3, nrow=2, widths=c(0.5,3.5))
# 
# 
# library(grid)
# ht = 3
# grid.newpage()
# pushViewport(viewport(layout = grid.layout(2, 2, heights = unit(c(ht, length(ages*ht)), "null"))))
# grid.rect()
# # grid.text("title of this panel", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
# print(d, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
# print(c, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
# print(p_binned, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
# 
# w = 5
# grid.newpage()
# pushViewport(viewport(width=10, height=10, default.units="in")); 
# #gt = grobTree(rectGrob(), circleGrob()) 
# grid.arrange(c, p_binned, nrow=1, newpage=FALSE,
#              widths = unit(c(w, w*7), "null"))
# upViewport()
# 
# 
# top <- arrangeGrob(r, r, r, r, nrow=2)
# 
# grid.arrange(top, r, ncol=1, main = "Main Title", heights=c(2, 1))
# 
# grid.newpage()
# align_plots(c,p_binned)
# # ####################################################################################################
# # # chunk: plot predicted and observed proportions in same frame
# # ####################################################################################################
# # suff4=paste(suff, '_compare', sep='')
# # 
# # plot_both_maps(r_mean, y_veg, centers=centers_pls, taxa=taxa, ages, N, K, T, thresh=1, suff=suff3, save_plots=save_plots)
