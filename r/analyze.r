library(ggplot2)
library(fields)
library(reshape2)
library(RColorBrewer)

source('r/utils/pred_analyze_funs.r')
source('r/utils/pred_helper_funs.r')
source('r/utils/pred_plot_funs.r')

# PL_var = list(suff_run = '120knots_150to2150ybp_PL_umw_3by_v2.1_ar_varves_set1.2',
#               suff_figs = 'varves_fix')
# PL_var = list(suff_run = '120knots_0to2000ybp_PL_umw_3by_v2.1_ar_varves',
#               suff_figs = 'from0')
PL_var = list(suff_run = '120knots_150to2150ybp_PL_umw_3by_v2.3_ar_set1',
              suff_figs = 'v2.3')
run = PL_var
run_num  = '1'
subDir <- paste0("runs/", run$suff_run, '/figures')


## load the data
# comp data
load(sprintf('%s/cal_data_%s.rdata', "../stepps-calibration/r/dump", "12taxa_mid_comp_ALL_v0.3"))
centers_comp=centers_veg

# pred run data
load(paste0('runs/', run$suff_run, '/run', run_num, '/input.rdata'))

rIts = readRDS(file=paste0('runs/', run$suff_run, '/rIts.RDS'))
rIts = rIts[,,,seq(1, dim(rIts)[4], by=10)]

get_r_quants <- function(rIts){

  N = dim(rIts)[1]
  K = dim(rIts)[2]
  T = dim(rIts)[3]
  
  r_int = apply(rIts, c(1,2,3), function(x) quantile(x, probs=c(0.1, 0.5, 0.9)))
  r_sd = apply(rIts, c(1,2,3), function(x) sd(x))

  r_quants = array(NA, c(4, N*T, K))
  for (k in 1:K){
    for (n in 1:N){
      for (t in 1:T){
        r_quants[1:3,(n-1)*T+t,k] = r_int[,n,k,t]
        r_quants[4,(n-1)*T+t,k] = r_sd[n,k,t]
      }
    }
  }
  return(r_quants)
}

r_quants = get_r_quants(rIts)
# 
# 
# 
# r_quants = apply(rIts, c(1,2,3), function(x) quantile(x, probs=c(0.1, 0.5, 0.9)))


# get_r_pt_devs <- function(rIts){
#   
#   N = dim(rIts)[1]
#   K = dim(rIts)[2]
#   T = dim(rIts)[3]
#   niter = dim(rIts)[4]
#   
#   r_t = array(NA, c(N, K, niter))
#   for (k in 1:K){
#     for (n in 1:N){
#         r_t[n,k,] = apply(rIts[n,k,,], 2, mean)
#       }
#   }
#   
#   
#   r_int = apply(rIts, c(1,2,3), function(x) quantile(x, probs=c(0.1, 0.5, 0.9)))
#   r_sd = apply(rIts, c(1,2,3), function(x) sd(x))
#   
#   r_quants = array(NA, c(4, N*T, K))
#   for (k in 1:K){
#     for (n in 1:N){
#       for (t in 1:T){
#         r_quants[1:3,(n-1)*T+t,k] = r_int[,n,k,t]
#         r_quants[4,(n-1)*T+t,k] = r_sd[n,k,t]
#       }
#     }
#   }
#   return(r_quants)
# }
# 



rescale=1e6
limits = get_limits(centers_pls)
save_plots = TRUE
breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)

r_mean2 = r_quants[2,,]
r_sd2   = r_quants[4,,]
 
# # p_binned <- plot_pred_maps_binned(r_mean, centers_veg, breaks, taxa, ages, N, K, T, limits, suff=suff1, save_plots=save_plots)
# p_binned <- plot_pred_maps_binned_select(r_mean2, centers_veg, breaks, taxa, taxa_sub=taxa_sub, ages, N, K, T, limits, 
#                                          suff=run$suff_figs, save_plots=TRUE, fpath=subDir)

taxa_sub = taxa[c(2,5,7,10,11)]

if (!file.exists(subDir)){
  dir.create(file.path(subDir))
} 

# p_binned <- plot_pred_maps_binned(r_mean, centers_veg, breaks, taxa, ages, N, K, T, limits, suff=suff1, save_plots=save_plots)
p_binned <- plot_pred_maps_binned_select(r_mean2, centers_veg, breaks, taxa, taxa_sub=taxa_sub, ages, N, K, T, limits, 
                                         suff=run$suff_figs, save_plots=TRUE, fpath=subDir)

# plot_pred_maps_binned_gif
p_binned <- plot_pred_maps_binned_gif(r_mean2, centers_veg, breaks, taxa, taxa_sub=taxa, ages, N, K, T, limits, 
                                         suff=run$suff_figs, save_plots=TRUE, fpath=subDir)

p_binned <- plot_pred_maps_binned_gifun(r_mean2, r_sd2, centers_veg, breaks, taxa, taxa_sub=taxa, ages, N, K, T, limits, 
                                      suff=run$suff_figs, save_plots=TRUE, fpath=subDir)

suff = paste0('sd_', run$suff_figs)
p_binned <- plot_pred_maps_binned_select(r_sd2, centers_veg, breaks, taxa, taxa_sub=taxa, ages, N, K, T, limits, 
                                         suff=suff, save_plots, fpath=subDir)


## plot bw
bw_buff = 164000
bw_lims <- list(xlims=c(min(bw_fort$long)-bw_buff, max(bw_fort$long)+bw_buff)/1e6, 
                ylims=c(min(bw_fort$lat)-64000, max(bw_fort$lat)+bw_buff)/1e6)

r_bw = rowSums(r[,c(1,4)])
plot_data_maps_binned_bw(r_bw, centers_comp/1e6, taxa='bw', 
                         nrow(centers_comp), K=1, breaks, bw_lims, 
                         suff='pls-comp_bw', save_plots, fpath=subDir)

rp_bw = rowSums(r_mean2[seq(1,N*T, T),c(1,4)])
plot_data_maps_binned_bw(rp_bw, centers_veg, taxa='pbw', 
                         N, K=1, breaks, bw_lims, 
                         suff='pred-comp_bw', save_plots, fpath=subDir)

r2_bw = rbind(data.frame(props=r_bw, x=centers_comp[,1], y=centers_comp[,2], type='PLS composition'),
              data.frame(props=rp_bw, x=centers_veg[,1]*rescale, y=centers_veg[,2]*rescale, type='STEPPS'))

p <- ggplot() + geom_tile(data=r2_bw, aes(x=x, y=y, fill=factor(props))) + 
  scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportions') + 
  coord_fixed() + 
  scale_x_continuous(limits$xlims*rescale) + scale_y_continuous(limits$ylims*rescale)
p <- add_map_albers(p, us.shp, limits)
p <- p + geom_path(data=bw, aes(x=long, y=lat),  colour='grey', alpha=0.8,size=1.8)
p <- theme_clean(p) #+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))

p <- p + theme(strip.text.x = element_blank(),
               strip.text.y = element_blank())
p <- p + theme(strip.background = element_blank())

print(p)

p <- ggplot() + geom_tile(data=r2_bw, aes(x=x, y=y, fill=props)) + 
  #scale_fill_gradientn(colours=tim.colors()) +
  # scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportions') + 
  coord_fixed() + 
  scale_x_continuous(limits$xlims*rescale) + scale_y_continuous(limits$ylims*rescale)
p <- add_map_albers(p, us.shp, limits)
p <- p + geom_path(data=bw, aes(x=long, y=lat),  colour='grey', alpha=0.8,size=1.8)
p <- p + facet_grid(.~type)
p <- theme_clean(p) #+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))

p <- p + theme(strip.text.x = element_blank(),
               strip.text.y = element_blank())
p <- p + theme(strip.background = element_blank())
print(p)
# plot_data_maps_binned_bw <- function(y, centers, taxa, N, K, T, breaks, limits, suff, save_plots, fpath=subDir){
  
# 
# p_binned <- plot_pred_maps_binned_select(r_mean2, centers_veg, breaks, taxa, taxa_sub=taxa, ages, N, K, T, limits, 
#                                          suff=run$suff_figs, save_plots=TRUE, fpath=subDir)

##########################################################################################################################################
## posterior difference plots
##########################################################################################################################################
ages =ages

t_diffs = c(1, 2, 5, 10, 15, 20)
# t_diffs = c(1, 2, 7, 14, 20)
# t_diffs = c(1, 2, 7, 14, 20)
# t_diffs = c(1, 2, 10, 20)
es    = eco_sig(rIts, t_diffs, eco_cut=0.03)
diffs = post_diffs(rIts, t_diffs)

post_prob = apply(diffs, c(1,2,3,4), diff_q, prob=0.85)
post_mag  = apply(diffs, c(1,2,3,4), diff_mag, prob=0.85)

post_prob_es = post_prob
post_prob_es[which(!es)] = NA

post_mag_es = post_mag
post_mag_es[which(!es)] = NA

# plot_post_change <- function(post_change, t_diffs, ages, centers_veg, N, taxa, rescale, type, subDir){
# plot_post_change <- function(post_change, t_diffs, ages, centers_veg, N, taxa, rescale, type, subDir){

# function(post_change, t_diffs, ages, centers_veg, N, taxa, taxa_sub, limits, rescale, type, subDir)

plot_post_change(post_prob_es, t_diffs, ages, centers_veg, N, taxa, taxa_sub=taxa, limits, rescale, type='prob', subDir=subDir)
plot_post_change(post_mag_es, t_diffs, ages, centers_veg, N, taxa, taxa_sub=taxa, limits, rescale, type='mag', subDir=subDir)

# function(post_change, t_diffs, ages, centers_veg, N, taxa, rescale, limits, type, group, subDir)
plot_post_change_group(post_prob_es[,c(1,4),,], t_diffs, ages, centers_veg, N, taxa[c(1,4)], 
                       rescale, limits=bw_lims, type='prob', group='bw', subDir)

tsuga_lims = list(xlims=c(3*10^5/rescale, max(centers_veg$x)), ylims=c(min(centers_veg$y), 1200000/rescale))
# plot_post_change_bw(post_prob_es[,5,,], t_diffs, ages, centers_veg, N, taxa[5], 
#                     rescale, limits=tsuga_lims, type='prob', subDir)
plot_post_change(post_prob_es, t_diffs, ages, centers_veg, N, taxa, taxa_sub=taxa[5], tsuga_lims, rescale, type='prob', subDir=subDir)


# deciduous
plot_post_change_group(post_prob_es[,c(1,2,3,4,6,7,9),,], t_diffs, ages, centers_veg, N, taxa[c(1,2,3,4,6,7,9)], 
                       rescale, limits=limits, type='prob', group='deciduous', subDir=subDir)
# deciduous
plot_post_change_group(post_prob_es[,c(1,2,3,4,6,7,9),,], t_diffs, ages, centers_veg, N, taxa[c(1,2,3,4,6,7,9)], 
                       rescale, limits=limits, type='mag', group='deciduous', subDir=subDir)

# conifer
plot_post_change_group(post_prob_es[,c(5,8,10,11,12),,], t_diffs, ages, centers_veg, N, taxa[c(5,8,10,11,12)], 
                       rescale, limits=limits, type='prob', group='conifer', subDir)
plot_post_change_group(post_prob_es[,c(5,8,10,11,12),,], t_diffs, ages, centers_veg, N, taxa[c(5,8,10,11,12)], 
                       rescale, limits=limits, type='mag', group='conifer', subDir)
##########################################################################################################################################
## posterior distance plots
##########################################################################################################################################
t_diffs = c(1, 2, 5, 10, 15, 20)
n_diffs = length(t_diffs)
niter   = dim(rIts)[4]

# dist_sc = array(NA, c(N, n_diffs, n_diffs, niter))
# dist_can = array(NA, c(N, n_diffs, n_diffs, niter))
dist_bc = array(NA, c(N, n_diffs, n_diffs, niter))
for (t1 in 1:n_diffs){
  print(t1)
  for (t2 in 1:n_diffs){
    for (n in 1:N){
      # print(n)
      for (i in 1:niter){
        # dist_sc[n,t1,t2,i] = squared_chord(rIts[n,,t_diffs[t1],i], rIts[n,,t_diffs[t2],i]) 
        # dist_can[n,t1,t2,i] = canberra(rIts[n,,t_diffs[t1],i], rIts[n,,t_diffs[t2],i]) 
        dist_bc[n,t1,t2,i] = bray_curtis(rIts[n,,t_diffs[t1],i], rIts[n,,t_diffs[t2],i]) 
      }
    }
  }
}

# sc  = apply(dist_sc, c(1,2,3), mean)
# can = apply(dist_can, c(1,2,3), mean)
# bc     = apply(dist_bc, c(1,2,3), mean)
bc     = apply(dist_bc, c(1,2,3), median)
bc_all = apply(bc, c(2,3), sum)

bc_allm = melt(bc_all)
bc_allm$value[bc_allm$value==0] = NA
ggplot(data=bc_allm) + geom_tile(aes(x=Var1, y=Var2, fill=value))

# try to include uncertainty


# sc_melted = data.frame(d=numeric(0), x=integer(0), y=integer(0), t1=numeric(0), t2=numeric(0))
# for (t2 in 2:n_diffs){
#   for (t1 in (t2-1):1){
#     sc_melted  = rbind(sc_melted , data.frame(d = as.vector(sc[,t1,t2]), 
#                                               x     = centers_veg[,1]*rescale, 
#                                               y     = centers_veg[,2]*rescale, 
#                                               t1 = rep(t_diffs[t1], N),
#                                               t2 = rep(t_diffs[t2], N)))
#   }
# }
# 
# can_melted = data.frame(d=numeric(0), x=integer(0), y=integer(0), t1=numeric(0), t2=numeric(0))
# for (t2 in 2:n_diffs){
#   for (t1 in (t2-1):1){
#     can_melted  = rbind(can_melted , data.frame(d = as.vector(can[,t1,t2]), 
#                                                 x     = centers_veg[,1]*rescale, 
#                                                 y     = centers_veg[,2]*rescale, 
#                                                 t1 = rep(t_diffs[t1], N),
#                                                 t2 = rep(t_diffs[t2], N)))
#   }
# }

# start at 2?
bc_melted = data.frame(d=numeric(0), x=integer(0), y=integer(0), t1=numeric(0), t2=numeric(0))
for (t2 in 2:n_diffs){
  for (t1 in (t2-1):1){
    bc_melted  = rbind(bc_melted , data.frame(d = as.vector(bc[,t1,t2]), 
                                              x     = centers_veg[,1]*rescale, 
                                              y     = centers_veg[,2]*rescale, 
                                              t1 = rep(t_diffs[t1], N),
                                              t2 = rep(t_diffs[t2], N)))
  }
}


plot_dissim <- function(dm_dat, ages, type, r_type=r_type, subDir=subDir){
  
  # dm_dat = dm_dat[which((dm_dat$t1!=1)),]
  
  dm_dat$t1 = ages[dm_dat$t1]*100
  dm_dat$t2 = ages[dm_dat$t2]*100
  
  dm_dat$t1 = factor(dm_dat$t1, levels=rev(sort(unique(dm_dat$t1))))
  
  
  p <- ggplot(data=dm_dat) + geom_tile(aes(x=x, y=y, fill=d)) 
  #   p <- p + scale_fill_gradientn(colours=c("blue", "lightskyblue", "white", "pink", "red"),  
  #                                 values=values, limits=c(-0.6,0.6), na.value="white",
  #                                 rescaler = function(x, ...) x, oob = identity) 
  p <- p + scale_fill_gradientn(colours=tim.colors(8), name="Distance")
  # p <- p + scale_fill_gradientn(colours=brewer.pal(8,"Paired"), name="Distance")
  p <- p + coord_fixed()
  p <- add_map_albers(p, map_data=us.fort, limits)
  p <- p + facet_grid(t2~t1, switch='x')
  p <- theme_clean(p) + theme_bw() + theme(strip.text = element_text(size = 12),
                                           strip.background = element_rect(colour = 'grey'),
                                           axis.ticks = element_blank(), 
                                           axis.text = element_blank(),
                                           axis.title = element_blank(),
                                           legend.text = element_text(size = 12))#+ labs(color='NEW LEGEND TITLE') 
  print(p)
  
  ggsave(file=paste0(subDir, '/post_dist_', r_type, '_', type, '.pdf'))
  
}

r_type = 'r_mu_g'
# plot_dissim(can_melted, type='canberra_test', r_type=r_type, subDir=subDir)
# plot_dissim(sc_melted, type='squared_chord', r_type=r_type, subDir=subDir)
plot_dissim(bc_melted, ages, type='bray_curtis', r_type=r_type, subDir=subDir)

##########################################################################################################################################
## read in NH temp proxies
##########################################################################################################################################
# read in and organize proxies
mo05 = read.table(file='data/IPCC/moberg2005_orig.dat.txt', skip=96, fill=TRUE, header=TRUE)
mo05 = data.frame(year=mo05$Year, temp=mo05$T, line=rep('mo05'))
# mo05 = data.frame(year=mo05$Year, temp=mo05$LF, line=rep('mo05'))

cl12 = read.table(file='data/IPCC/christiansen_ljung2012_comb.dat.txt', skip=219, fill=TRUE, header=TRUE)
cl12 = data.frame(year=cl12$YearAD, temp=cl12$Recon, line=rep('cl12'))

lj10 = read.table(file='data/IPCC/ljungqvist2010decadal.csv', skip=0, fill=TRUE, header=TRUE, sep=',')
lj10 = data.frame(year=lj10$Year, temp=lj10$ReconT, line=rep('lj10'))

ma08 = read.table(file='data/IPCC/mannetal2008_cps_NHCRUTEM.dat.txt', skip=6, fill=TRUE, header=TRUE)
ma08 = data.frame(year=ma08$Year, temp=ma08$CPScomposite, line=rep('ma08'))

lm08 = read.table(file='data/IPCC/loehlemcculloch2008.dat.txt', skip=6, fill=TRUE, header=TRUE)
lm08 = data.frame(year=lm08$X.AD., temp=lm08$Anom., line=rep('lm08'))


plot(cl12$year, cl12$temp, type='l')
lines(mo05$year, mo05$temp, col='blue')
lines(lj10$year, lj10$temp, col='red')
lines(lm08$year, lm08$temp, col='green')
lines(ma08$year, ma08$temp, col='pink')


# recons = rbind(cl12, mo05, lj10, lm08, ma08)
recons = rbind(cl12, lj10)
# recons = rbind(cl12)
recons$year = 1950 - as.numeric(recons$year)
recons = recons[which(!is.nan(recons$temp)),]

ss = smooth.spline(recons$year, recons$temp, df=40)
# ss2 = ffcsaps(y=out$temp, x=out$year, nyrs = 100, f = 0.5)


plot(recons$year, recons$temp)
lines(ss$x, ss$y, col='blue', cex=5)

library(ggplot2)
library(splines)
ggplot(data=recons, aes(x=year, y=temp, colour=line)) + geom_line() +  stat_smooth(data = recons, aes(x = year, y = temp), method = 'lm', 
                                                                                   formula = y ~ ns(x, 50), col = 'red', size = 2, se= FALSE)

ggplot(data=recons, aes(x=year, y=temp)) + geom_line() +  stat_smooth(data = recons, aes(x = year, y = temp), method = 'lm', 
                                                                      formula = y ~ ns(x, 40), col = 'red', size = 2, se= FALSE)

##########################################################################################################################################
## posterior distance total consecutive
##########################################################################################################################################

t_diffs = seq(1, 20, by=2)
n_diffs = length(t_diffs)
niter   = dim(rIts)[4]

eco_cut   = 0.02
es        = eco_sig(rIts, t_diffs, eco_cut=eco_cut)
diffs     = post_diffs(rIts, t_diffs)
post_prob = apply(diffs, c(1,2,3,4), diff_q, prob=0.8)


dist_bc = array(NA, c(N, n_diffs-1, niter))
sig = array(NA, c(N, k, n_diffs-1))
for (t1 in 1:(n_diffs-1)){
  print(t1)
  for (n in 1:N){
    x = rIts[n,,t_diffs[t1],]
    y = rIts[n,,t_diffs[t1+1],]

    esig = rep(NA, K)
    for (k in 1:K){
      # eco sig; and or or?
      esig[k] = (mean(x[k,]) > eco_cut) | (mean(y[k,]) > eco_cut)
    }

    diffxy = x-y
    diffxy_sig = apply(diffxy, 1, diff_q, prob=0.7)

    sig[n,,t1] = (!is.na(esig)) & (!is.na(diffxy_sig))
    print(sum(sig[n,,t1]))
    if (sum(sig[n,,t1])<2){
      next
    } else {
      # xsig = x[sig[n,,t1],]
      # ysig = y[sig[n,,t1],]
      xsig = x#[sig[n,,t1],]
      ysig = y#[sig[n,,t1],]
      for (i in 1:niter){
        dist_bc[n,t1,i] = bray_curtis(xsig[,i]/sum(xsig[,i]), ysig[,i]/sum(ysig[,i]))
      }
    }
  }
}



# dist_bc = array(NA, c(N, n_diffs-1, niter))
# for (t1 in 1:(n_diffs-1)){
#   print(t1)
#   for (n in 1:N){
#     x = rIts[n,,t_diffs[t1],]
#     y = rIts[n,,t_diffs[t1+1],]
# 
#     for (i in 1:niter){
#       dist_bc[n,t1,i] = bray_curtis(x[,i]/sum(x[,i]), y[,i]/sum(y[,i]))
#     }
#   }
# }



# sc  = apply(dist_sc, c(1,2,3), mean)
# can = apply(dist_can, c(1,2,3), mean)
# bc  = apply(dist_bc, c(1,2), median)
# bc_all = colSums(bc, na.rm=TRUE)
# bc_all = bc_all[-1]
# bc  = apply(dist_bc, c(1,2), median)
# bc_all = colSums(bc, na.rm=TRUE)

# # only for locations with pollen cores
# dist_bc = dist_bc[idx_cores,,]

# try for all grid cells surrounding cores grid cells
d_close = rdist(centers_veg[idx_cores,], centers_veg)
idx_close=unique(unlist(apply(d_close, 1, function(x) which(x<=sqrt(0.024^2*2)))))


dist_bc_close = dist_bc[idx_close,,]
# bc_zero  = apply(dist_bc_close, c(2,3), function(x) length(which(x ==0)))
bc  = apply(dist_bc_close, c(2,3), sum, na.rm=TRUE)
bc_quants = apply(bc, 1, function(x) quantile(x, probs=c(0.05, 0.5, 0.95), na.rm=TRUE))

# bc_all = data.frame(bc=bc_all, time=(t_diffs[2:n_diffs]-t_diffs[1:(n_diffs-1)])/2 + t_diffs[1:(n_diffs-1)])
# ggplot(data=bc_all[-1,]) + geom_line(aes(x=time*100, y=bc))

dat = rbind(data.frame(value=bc_quants[2,], time=(((t_diffs[2:n_diffs]-t_diffs[1:(n_diffs-1)])/2 + t_diffs[1:(n_diffs-1)])*100), 
                       line=rep('mid'), type=rep('veg change'), lty=rep('solid'), colour=rep('veg'), size=rep(1)),
            data.frame(value=bc_quants[1,], time=(((t_diffs[2:n_diffs]-t_diffs[1:(n_diffs-1)])/2 + t_diffs[1:(n_diffs-1)])*100), 
                       line=rep('lb'), type=rep('veg change'), lty=rep('dashed'), colour=rep('veg'), size=rep(1)),
            data.frame(value=bc_quants[3,], time=(((t_diffs[2:n_diffs]-t_diffs[1:(n_diffs-1)])/2 + t_diffs[1:(n_diffs-1)])*100), 
                       line=rep('ub'), type=rep('veg change'), lty=rep('dashed'), colour=rep('veg'), size=rep(1)),
            data.frame(value=recons$temp, time=recons$year, line=recons$line, type=rep('northern hemispere temp'), lty=rep('solid'), 
                       colour=factor(recons$line), size=rep(1)),
            data.frame(value=ss$y, time=ss$x, line='smooth', type=rep('northern hemispere temp'), lty=rep('solid'), colour=rep('blue'), size=rep(2)))


blues <- brewer.pal(n = 4, name = "BrBG")

ggplot(data=dat) + geom_line(data=dat, aes(x=time, y=value, group=line, colour=colour, linetype=lty, size=factor(size))) +
  # scale_linetype_discrete(values=c("solid", "twodash", "twodash", "twodash")) +
  scale_linetype_manual(values=c("solid", "dashed"), guide=FALSE) +
  scale_size_manual(values=c(1,2), guide=FALSE)+
  scale_colour_manual(values=c('black', blues), guide=FALSE) +
  # scale_colour_brewer() +
  scale_x_reverse(limits=c(2000,0)) + 
  xlab('Years BP') + facet_grid(type~., scales="free_y")

ggsave('figures/total_change_veg_temp_time_no_sig.pdf')

# try to include uncertainty


# sc_melted = data.frame(d=numeric(0), x=integer(0), y=integer(0), t1=numeric(0), t2=numeric(0))
# for (t2 in 2:n_diffs){
#   for (t1 in (t2-1):1){
#     sc_melted  = rbind(sc_melted , data.frame(d = as.vector(sc[,t1,t2]), 
#                                               x     = centers_veg[,1]*rescale, 
#                                               y     = centers_veg[,2]*rescale, 
#                                               t1 = rep(t_diffs[t1], N),
#                                               t2 = rep(t_diffs[t2], N)))
#   }
# }
# 
# can_melted = data.frame(d=numeric(0), x=integer(0), y=integer(0), t1=numeric(0), t2=numeric(0))
# for (t2 in 2:n_diffs){
#   for (t1 in (t2-1):1){
#     can_melted  = rbind(can_melted , data.frame(d = as.vector(can[,t1,t2]), 
#                                                 x     = centers_veg[,1]*rescale, 
#                                                 y     = centers_veg[,2]*rescale, 
#                                                 t1 = rep(t_diffs[t1], N),
#                                                 t2 = rep(t_diffs[t2], N)))
#   }
# }

bc_melted = data.frame(d=numeric(0), x=integer(0), y=integer(0), t1=numeric(0), t2=numeric(0))
for (t2 in 2:n_diffs){
  for (t1 in (t2-1):1){
    bc_melted  = rbind(bc_melted , data.frame(d = as.vector(bc[,t1,t2]), 
                                              x     = centers_veg[,1]*rescale, 
                                              y     = centers_veg[,2]*rescale, 
                                              t1 = rep(t_diffs[t1], N),
                                              t2 = rep(t_diffs[t2], N)))
  }
}


plot_dissim <- function(dm_dat, type, r_type=r_type, subDir=subDir){
  
  dm_dat = dm_dat[which((dm_dat$t1!=1)),]
  
  p <- ggplot(data=dm_dat) + geom_tile(aes(x=x, y=y, fill=d)) 
  #   p <- p + scale_fill_gradientn(colours=c("blue", "lightskyblue", "white", "pink", "red"),  
  #                                 values=values, limits=c(-0.6,0.6), na.value="white",
  #                                 rescaler = function(x, ...) x, oob = identity) 
  p <- p + scale_fill_gradientn(colours=tim.colors(8))
  # p <- p + scale_fill_gradientn(colours=brewer.pal(8,"Accent"))
  p <- p + coord_fixed()
  p <- add_map_albers(p, map_data=us.fort, limits)
  p <- p + facet_grid(t2~t1)
  p <- theme_clean(p) # theme(strip.text.x = element_blank(), strip.text.y = element_text(size = rel(1.5)))
  print(p)
  
  ggsave(file=paste0(subDir, '/post_dist_', r_type, '_', type, '.pdf'))
  
}

r_type = 'r_mu_g'
# plot_dissim(can_melted, type='canberra_test', r_type=r_type, subDir=subDir)
# plot_dissim(sc_melted, type='squared_chord', r_type=r_type, subDir=subDir)
plot_dissim(bc_melted, type='bray_curtis', r_type=r_type, subDir=subDir)

##########################################################################################################################################
## posterior area totals
##########################################################################################################################################
# rIts = rIts[,,,1000:1800]

t_diffs = seq(2, 20, by=2)
n_diffs = length(t_diffs)
niter   = dim(rIts)[4]

years = ((ages[t_diffs[2:n_diffs]]-ages[t_diffs[1:(n_diffs-1)]])/2 + ages[t_diffs[1:(n_diffs-1)]])*100

eco_cut   = 0.03
es        = eco_sig(rIts, t_diffs, eco_cut=eco_cut)
diffs     = post_diffs(rIts, t_diffs)

probs = c(0.75, 0.8, 0.85)
area_con = array(NA, c(length(probs),n_diffs-1))
area_con_tax = array(NA, c(length(probs), K, n_diffs-1))

area_con_unique = array(NA, c(length(probs),n_diffs-1))

for (k in 1:length(probs)){
  print(k)
  
  prob = probs[k]
  post_prob = apply(diffs, c(1,2,3,4), diff_q, prob=prob)
  
  post_prob_es = post_prob
  post_prob_es[which(!es)] = NA

  area_unique = apply(post_prob_es, c(3,4), function(x) sum(apply(x,1,function(y) any(!is.na(y)))))

  area_tax = apply(post_prob_es, c(2,3,4), function(x) sum(!is.na(x)))
  area_tot = apply(area_tax, c(2,3), function(x) sum(x))

  for (i in 2:n_diffs){
    area_con[k, i-1] = area_tot[i-1,i]
    area_con_tax[k, ,i-1] = area_tax[,i-1,i]
    area_con_unique[k, i-1] = area_unique[i-1,i]
  }
}

area_con

dat = rbind(data.frame(value=area_con[2,], time=(((ages[t_diffs[2:n_diffs]]-ages[t_diffs[1:(n_diffs-1)]])/2 
                                                  + ages[t_diffs[1:(n_diffs-1)]])*100),
                       line=rep('mid'), type=rep('Significant change (n. cells)'), lty=rep('solid'), colour=rep('veg'), size=rep(1)),
            data.frame(value=recons$temp, time=recons$year, line=recons$line, type=rep('Temp. Anomaly (degrees C)'), lty=rep('solid'), 
                       colour=factor(recons$line), size=rep(1)))

mydf <- data.frame(ub=area_con[1,], mid=area_con[2,], lb=area_con[3,], type=rep('Significant change (n. cells)'), 
                   time=(((ages[t_diffs[2:n_diffs]]-ages[t_diffs[1:(n_diffs-1)]])/2 
                          + ages[t_diffs[1:(n_diffs-1)]])*100))



# cols <- brewer.pal(n = 4, name = "BrBG")
cols = c('black', 'grey', 'red')
p <- ggplot(data=dat) + geom_line(data=dat, aes(x=time, y=value, group=line, colour=colour, linetype=lty, size=factor(size))) +
  geom_ribbon(data=mydf, aes(x=time, ymin=lb, ymax=ub), colour='lightgrey', alpha=0.5)+
  # scale_linetype_discrete(values=c("solid", "twodash", "twodash", "twodash")) +
  scale_linetype_manual(values=c("solid", "twodash"), guide=FALSE) +
  scale_size_manual(values=c(1,2), guide=FALSE)+
  scale_colour_manual(values=c('black', cols), guide=FALSE) +
  # scale_colour_brewer() +
  # scale_x_continuous(breaks=seq(0, 2000, by=250))+
  scale_x_reverse(limits=c(2000,0),breaks=seq(0, 2000, by=250)) + 
  xlab('Years BP') + facet_grid(type~., scales="free_y")
p <- theme_clean(p) + theme_bw() + theme(strip.text = element_text(size = 12),
                                         axis.title.y = element_blank(),
                                         axis.text = element_text(size=14),
                                         axis.title = element_text(size=14))
print(p)

ggsave('figures/total_sig_change_area_ribbon.pdf')

##

sig_tax = t(area_con_tax[2,,])
colnames(sig_tax) = taxa

sig_tax = data.frame(sig_tax, time=(((ages[t_diffs[2:n_diffs]]-ages[t_diffs[1:(n_diffs-1)]])/2 
                                     + ages[t_diffs[1:(n_diffs-1)]])*100))
sig_taxm = melt(sig_tax, id.vars = 'time')  
  
p <- ggplot(data=sig_taxm) + geom_line(aes(x=time, y=value, colour=variable)) +
  # scale_linetype_manual(values=c("solid", "twodash"), guide=FALSE) +
  # scale_size_manual(values=c(1,2), guide=FALSE)+
  # scale_colour_manual(values=c('black', cols), guide=FALSE) +
  scale_x_reverse(limits=c(2000,0),breaks=seq(0, 2000, by=250)) + 
  xlab('Years BP') #+ facet_grid(type~., scales="free_y")
p <- theme_clean(p) + theme_bw() + theme(strip.text = element_text(size = 12),
                                         axis.title.y = element_blank(),
                                         axis.text = element_text(size=14),
                                         axis.title = element_text(size=14))
print(p)

ggsave('figures/total_change_area_taxon.pdf')


mid_ss = spline((((ages[t_diffs[2:n_diffs]]-ages[t_diffs[1:(n_diffs-1)]])/2 
         + ages[t_diffs[1:(n_diffs-1)]])*100), area_con_unique[2,], n=201)
ub_ss = spline((((ages[t_diffs[2:n_diffs]]-ages[t_diffs[1:(n_diffs-1)]])/2 
                  + ages[t_diffs[1:(n_diffs-1)]])*100), area_con_unique[1,], n=201)
lb_ss = spline((((ages[t_diffs[2:n_diffs]]-ages[t_diffs[1:(n_diffs-1)]])/2 
                 + ages[t_diffs[1:(n_diffs-1)]])*100), area_con_unique[3,], n=201)

dat = data.frame(value=area_con_unique[2,], 
                 time=(((ages[t_diffs[2:n_diffs]]-ages[t_diffs[1:(n_diffs-1)]])/2 
                                                  + ages[t_diffs[1:(n_diffs-1)]])*100),
                 line=rep('mid'), 
                 type=rep('Significant change (n. cells)'), 
                 lty=rep('solid'), 
                 colour=rep('veg'), 
                 size=rep(1))

dat2 = rbind(data.frame(value=mid_ss$y, 
                        time=mid_ss$x,
                        line=rep('mid'), 
                        type=rep('Significant change (n. cells)'), 
                        lty=rep('solid'), 
                        colour=rep('veg'), 
                        size=rep(1)),
            data.frame(value=recons$temp, 
                       time=recons$year, 
                       line=recons$line, 
                       type=rep('Temp. Anomaly (degrees C)'), 
                       lty=rep('solid'), 
                       colour=factor(recons$line), 
                       size=rep(1)))

# mydf <- data.frame(ub=area_con_unique[1,], 
#                    mid=area_con_unique[2,], 
#                    lb=area_con_unique[3,], 
#                    type=rep('Significant change (n. cells)'), 
#                    time=(((ages[t_diffs[2:n_diffs]]-ages[t_diffs[1:(n_diffs-1)]])/2 
#                           + ages[t_diffs[1:(n_diffs-1)]])*100))

mydf <- data.frame(ub=ub_ss$y, 
                   mid=mid_ss$y, 
                   lb=lb_ss$y, 
                   type=rep('Significant change (n. cells)'), 
                   time=mid_ss$x)


cols = c('black', 'grey', 'red')
p <- ggplot(data=dat) + geom_line(data=dat2, aes(x=time, y=value, group=line, colour=colour, linetype=lty, size=factor(size))) +
  geom_point(data=dat, aes(x=time, y=value, group=line, colour=colour, size=factor(size))) +
  geom_ribbon(data=mydf, aes(x=time, ymin=lb, ymax=ub), colour='lightgrey', alpha=0.5)+
  scale_linetype_manual(values=c("solid", "twodash"), guide=FALSE) +
  scale_size_manual(values=c(1,2), guide=FALSE)+
  scale_colour_manual(values=c('black', cols), guide=FALSE) +
  scale_x_reverse(limits=c(2000,0),breaks=seq(0, 2000, by=250)) + 
  xlab('Years BP') + facet_grid(type~., scales="free_y")
p <- theme_clean(p) + theme_bw() + theme(strip.text = element_text(size = 12),
                                         axis.title.y = element_blank(),
                                         axis.text = element_text(size=14),
                                         axis.title = element_text(size=14))
print(p)

ggsave('figures/community_change_area_unique.pdf')



dat = rbind(data.frame(value=area_con_unique[2,], 
                       time=years,
                       line=rep('mid'), 
                       type=rep('Total change (cells)'), 
                       lty=rep('solid'), 
                       colour=rep('veg'), 
                       size=rep(1)),
            data.frame(value=area_con_tax[2,5,], 
                       time=years, 
                       line=rep('hemlock'), 
                       type=rep('Taxon change (cells)'), 
                       lty=rep('solid'), 
                       colour=rep('Hemlock'), 
                       size=rep(1)),
            data.frame(value=area_con_tax[2,10,], 
                       time=years, 
                       line=rep('pine'), 
                       type=rep('Taxon change (cells)'), 
                       lty=rep('solid'), 
                       colour=rep('Pine'), 
                       size=rep(1)),
            data.frame(value=area_con_tax[2,4,], 
                       time=years, 
                       line=rep('elm'), 
                       type=rep('Taxon change (cells)'), 
                       lty=rep('solid'), 
                       colour=rep('Elm'), 
                       size=rep(1)),
            data.frame(value=area_con_tax[2,7,], 
                       time=years, 
                       line=rep('oak'), 
                       type=rep('Taxon change (cells)'), 
                       lty=rep('solid'), 
                       colour=rep('Oak'), 
                       size=rep(1)),
            data.frame(value=area_con_tax[2,9,], 
                       time=years, 
                       line=rep('oh'), 
                       type=rep('Taxon change (cells)'), 
                       lty=rep('solid'), 
                       colour=rep('OH'), 
                       size=rep(1)),
            data.frame(value=recons$temp, 
                       time=recons$year, 
                       line=recons$line, 
                       type=rep('Temp. anomaly (C)'), 
                       lty=rep('solid'), 
                       colour=factor(recons$line), 
                       size=rep(1)))

mydf <- data.frame(ub=area_con_unique[1,], 
                   mid=area_con_unique[2,], 
                   lb=area_con_unique[3,], 
                   type=rep('Total change (cells)'), 
                   time=years)


# dat=subset(dat, time>300)

cols <- brewer.pal(n = 5, name = "Set1")
# cols <- tim.colors(n = 6)
# cols = c('black', 'grey', 'red')
cols = c(cols, 'black', 'grey')
p <- ggplot(data=dat) + 
  geom_line(data=dat, aes(x=time, y=value, group=line, colour=colour, linetype=lty, size=factor(size))) +
  geom_ribbon(data=mydf, aes(x=time, ymin=lb, ymax=ub), colour='lightgrey', alpha=0.5)+
  # scale_linetype_discrete(values=c("solid", "twodash", "twodash", "twodash")) +
  scale_linetype_manual(values=c("solid", "twodash"), guide=FALSE) +
  scale_size_manual(values=c(1,2), guide=FALSE)+
  scale_colour_manual(values=c('black', cols), breaks=c('Elm', 'Hemlock', 'Oak', 'OH', 'Pine'), name='Taxon') +
  # scale_colour_brewer() +
  # scale_x_continuous(breaks=seq(0, 2000, by=250))+
  # scale_x_reverse(limits=c(2000,0),breaks=seq(0, 2000, by=250)) + 
  scale_x_reverse(limits=c(2000,400),breaks=seq(400, 2000, by=250)) + 
  xlab('Years BP') + facet_grid(type~., scales="free_y")
p <- theme_clean(p) + theme_bw() + theme(strip.text = element_text(size = 10),
                                         axis.title.y = element_blank(),
                                         axis.text = element_text(size=12),
                                         axis.title = element_text(size=12))
print(p)

ggsave('figures/community_change_area_unique.pdf')

#########################################################################################################################################
## contour plots
#########################################################################################################################################

prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0),time=integer(0), taxon=character())
for (k in 1:K){
  prop_dat = rbind(prop_dat, data.frame(props = r_mean2[,k], 
                                        x     = rep(centers_veg[,1], each=T)*rescale, 
                                        y     = rep(centers_veg[,2], each=T)*rescale, 
                                        time  = rep(ages,times=N), 
                                        taxon = rep(taxa[k], N*T)))
}

# find proportion transitions to define contours
prop_cuts = c(0.01, 0.1, 0.2, 0.4)
prop_dat_cuts = data.frame(prop_dat, cut=rep(prop_cuts[1], nrow(prop_dat)))
for (i in 2:length(prop_cuts)){
  prop_dat_shift = prop_dat
  prop_dat_shift$props =  prop_dat_shift$props - prop_cuts[i] + 0.01
  prop_dat_cuts = rbind(prop_dat_cuts, data.frame(prop_dat_shift, cut=rep(prop_cuts[i], nrow(prop_dat))))
}
prop_dat_cuts = data.frame(prop_dat_cuts)

pdf(file=paste0(subDir, '/contour_change_all.pdf'), width=16, height=12)
for (i in 1:length(taxa)){
  
  taxon = taxa[i]
  prop_sub = prop_dat_cuts[prop_dat_cuts$taxon == taxon,]
  p1 <- ggplot(prop_sub,aes(x=x,y=y))+ geom_contour(aes(z=props, group=time, color=time), breaks=0.01) + 
    facet_wrap(~cut, nrow=2) #+ geom_tile(aes(fill=props, alpha=0.2))+ 
  p1 <- p1 + scale_colour_gradientn(colours = tim.colors(10))
  p1 <- add_map_albers(p1, map_data=us.fort, limits)
  p1 <- theme_clean(p1) #+ theme(strip.text.x = element_blank(), strip.text.y = element_text(size = rel(1.5)))
  p1 <- p1 + annotate("text", y= max(prop_sub$y), x =max(prop_sub$x), label=taxon, hjust=1)
  print(p1)
  
  #ggsave(p, file=paste0(subDir, '/contour_change_', taxon, '.pdf'))
}
dev.off()

# pdf(file=paste0(subDir, '/contour_change_hemlock.pdf'), width=16, height=12)
taxon = taxa[5]
prop_sub = prop_dat_cuts[prop_dat_cuts$taxon == taxon,]
p1 <- ggplot(prop_sub,aes(x=x,y=y))+ geom_contour(aes(z=props, group=time, color=time), breaks=0.01) + facet_wrap(~cut, nrow=2) #+ geom_tile(aes(fill=props, alpha=0.2))+ 
p1 <- p1 + scale_colour_gradientn(colours = tim.colors(10), 
                                  name='Time', 
                                  breaks=c(0.5, 10.5, 19.5),
                                  labels=c(0.5, 10.5, 19.5))
p1 <- add_map_albers(p1, map_data=us.fort, limits)
p1 <- theme_clean(p1) + theme(strip.text.x = element_text(size = 14), 
                              legend.text=element_text(size=13), 
                              legend.title=element_text(size=16),
                              legend.title.align=0)
print(p1)
ggsave(p1,file=paste0(subDir, '/contour_change_hemlock.pdf'))

##########################################################################################################################################
## diagonal cut
##########################################################################################################################################
# diagonal cut
# p1 = c(-0.059, 0.862)
p1 = c(-0.059, 0.862) + c(0,0.024)
coord_old = p1
in_domain = TRUE
coords_diag = p1
while (in_domain) {
  # coord_new = coord_old + c(0.024,0.024)
  coord_new = coord_old + c(0.024,2*0.024)
  if (rdist(matrix(coord_new, nrow=1), matrix(c(0.205, 1.126), nrow=1)) < 1e-10) {
    in_domain == TRUE
  } else {
    in_domain = any(rdist(matrix(coord_new, nrow=1), centers_veg) < 1e-10)
  }
  
  if (in_domain) {
    coords_diag = rbind(coords_diag, coord_new)
    coord_old = coord_new
  }
}

coords_diag

pdf(file=paste0(subDir, '/ecotone_transect.pdf'))
plot(centers_veg$x, centers_veg$y)
points(coords_diag[,1], coords_diag[,2], pch=19)
dev.off()

diag_transect <- function(dat, coords_diag){
  dat = dat[with(dat, order(x)), ]
  d_diag = rdist(coords_diag[,1:2]*1e6, matrix(cbind(dat$x, dat$y), ncol=2))
  idx_keep = apply(d_diag, 2, function(x) any(x < 1e-8))
  # id_keep  = apply(d_diag, 2, function(x){if(any(x < 1e-10)){ which(x < 1e-10)}})
  dat_tran = dat[idx_keep, ]
  
  return(dat_tran)
}

# 
# can_tran = diag_transect(can_melted, coords_diag)
# can_tran = can_tran[can_tran$t1 != 1,]
# can_tran = can_tran[can_tran$t1 == 2,]
# 
# p <- ggplot() + 
#   geom_line(data=can_tran, aes(y = d, x=x, colour=factor(t2)))+
#   scale_colour_manual(values=tim.colors(4)) 
# print(p)
# 
# ggsave(file=paste0(subDir, '/diss_transect_canberra_200.pdf'))


bc_tran = diag_transect(bc_melted, coords_diag)
bc_tran = bc_tran[bc_tran$t1 != 1,]
bc_tran = bc_tran[bc_tran$t1 == 2,]

p <- ggplot() + 
  geom_line(data=bc_tran, aes(y = d, x=x, colour=factor(t2)))+
  scale_colour_manual(values=tim.colors(4)) 
print(p)

ggsave(file=paste0(subDir, '/diss_transect_bray_curtis_200.pdf'))

## ecotone proportions

prop_dat = data.frame(props=numeric(0), x=integer(0), y=integer(0),time=integer(0), taxon=character())
for (k in 1:K){
  prop_dat = rbind(prop_dat, data.frame(props = r_mean2[,k], 
                                        x     = rep(centers_veg[,1]*1000000, each=T), 
                                        y     = rep(centers_veg[,2]*1000000, each=T), 
                                        time  = rep(ages, times=N), 
                                        taxon = rep(taxa[k], N*T)))
}

prop_dat = prop_dat[with(prop_dat, order(x)), ]


# # diagonal cut
# p1 = c(-0.059, 0.862)
# coord_old = p1
# in_domain = TRUE
# coords_diag = p1
# while (in_domain) {
#   coord_new = coord_old + c(0.024,0.024)
#   if (rdist(matrix(coord_new, nrow=1), matrix(c(0.205, 1.126), nrow=1)) < 1e-10) {
#     in_domain == TRUE
#   } else {
#     in_domain = any(rdist(matrix(coord_new, nrow=1), centers_veg) < 1e-10)
#   }
#   
#   if (in_domain) {
#     coords_diag = rbind(coords_diag, coord_new)
#     coord_old = coord_new
#   }
# }
# 
# coords_diag

coords_diag = data.frame(coords_diag, id = seq(1, nrow(coords_diag)))

pdf(file=paste0(subDir, '/ecotone_transect.pdf'))
plot(centers_veg$x, centers_veg$y)
points(coords_diag[,1], coords_diag[,2], pch=19)
dev.off()


prop_dat = prop_dat[with(prop_dat, order(x)), ]
d_diag = rdist(coords_diag[,1:2]*1e6, matrix(cbind(prop_dat$x, prop_dat$y), ncol=2))
idx_keep = apply(d_diag, 2, function(x) any(x < 1e-8))
# id_keep  = apply(d_diag, 2, function(x){if(any(x < 1e-10)){ which(x < 1e-10)}})
prop_dat1 = prop_dat[idx_keep, ]



prop_dat2 = prop_dat1[which(prop_dat1$taxon == 'OAK'), ]
prop_dat3 = prop_dat2[order(prop_dat2$x),]

foo = rdist(coords_diag[,1:2]*1e6, matrix(cbind(prop_dat3$x, prop_dat3$y), ncol=2))
id_keep  = apply(foo, 2, function(x){if (any(x < 1e-8)){ which(x < 1e-8)}})
prop_dat3 = data.frame(prop_dat3, id=id_keep)

p <- ggplot() + geom_line(data=prop_dat3, aes(x=id, y=props, colour=factor(time))) + 
  # stat_smooth(data=prop_dat3, aes(y = props, x=id, colour=factor(time)), method="gam", formula = y ~ s(x), se=FALSE)+# level=0.5) +
  scale_colour_manual(values=tim.colors(20)) 
print(p)

prop_dat2 = prop_dat1[which(prop_dat1$taxon == 'OTHER.HARDWOOD'), ]
prop_dat3 = prop_dat2[order(prop_dat2$x),]

foo = rdist(coords_diag[,1:2]*1e6, matrix(cbind(prop_dat3$x, prop_dat3$y), ncol=2))
id_keep  = apply(foo, 2, function(x){if (any(x < 1e-8)){ which(x < 1e-8)}})
prop_dat3 = data.frame(prop_dat3, id=id_keep)

p <- ggplot() + geom_line(data=prop_dat3, aes(x=id, y=props, colour=factor(time))) + 
  # stat_smooth(data=prop_dat3, aes(y = props, x=id, colour=factor(time)), method="gam", formula = y ~ s(x), se=FALSE)+# level=0.5) +
  scale_colour_manual(values=tim.colors(20)) 
print(p)

prop_dat2 = prop_dat1[which(prop_dat1$taxon == 'PINE'), ]
prop_dat3 = prop_dat2[order(prop_dat2$x),]

foo = rdist(coords_diag[,1:2]*1e6, matrix(cbind(prop_dat3$x, prop_dat3$y), ncol=2))
id_keep  = apply(foo, 2, function(x){if (any(x < 1e-8)){ which(x < 1e-8)}})
prop_dat3 = data.frame(prop_dat3, id=id_keep)

p <- ggplot() + 
  # stat_smooth(data=prop_dat3, aes(y = props, x=id, colour=factor(time)), method="gam", formula = y ~ s(x), se=FALSE)+# level=0.5) +
  geom_line(data=prop_dat3, aes(x=id, y=props, colour=factor(time))) + 
  scale_colour_manual(values=tim.colors(20)) 
print(p)

prop_dat2 = prop_dat1[which(prop_dat1$taxon %in% c('PINE', 'OAK', 'OTHER.HARDWOOD')), ]
prop_dat3 = prop_dat2[order(prop_dat2$x),]
prop_dat3 = prop_dat3[prop_dat3$time %in% ages[c(2, 3, 5, 10, 15, 20)],]

foo = rdist(coords_diag[,1:2]*1e6, matrix(cbind(prop_dat3$x, prop_dat3$y), ncol=2))
id_keep  = apply(foo, 2, function(x){if (any(x < 1e-8)){ which(x < 1e-8)}})
prop_dat3 = data.frame(prop_dat3, id=id_keep)



# p <- ggplot() + 
#   stat_smooth(data=prop_dat3, aes(y = props, x=id, colour=factor(time), linetype=taxon), method="gam", formula = y ~ s(x), se=FALSE)+# level=0.5) +
#   # geom_line(data=prop_dat3, aes(x=id, y=props, colour=factor(time), linetype=taxon)) + 
#   xlim(10,19) +
#   scale_colour_manual(values=tim.colors(20)) 
# print(p)


p <- ggplot() + 
  # stat_smooth(data=prop_dat3, aes(y = props, x=id, colour=factor(time), linetype=taxon), method="gam", formula = y ~ s(x), se=FALSE)+# level=0.5) +
  # geom_line(data=prop_dat3, aes(x=id, y=props, colour=factor(taxon), linetype=factor(time))) + 
  geom_line(data=prop_dat3, aes(x=id, y=props, colour=factor(time), linetype=factor(taxon))) + 
  #xlim(10,19) +
  scale_colour_manual(values=tim.colors(6)) 
print(p)
ggsave(file=paste0(subDir, '/diag_transect_oop.pdf'))

p <- ggplot() + 
  # stat_smooth(data=prop_dat3, aes(y = props, x=id, colour=factor(time), linetype=taxon), method="gam", formula = y ~ s(x), se=FALSE)+# level=0.5) +
  # geom_line(data=prop_dat3, aes(x=id, y=props, colour=factor(taxon), linetype=factor(time))) + 
  geom_line(data=prop_dat3, aes(x=id, y=props, linetype=factor(taxon))) + 
  # xlim(9,15) +
  # scale_colour_manual(values=c("red", "blue", "black")) 
  facet_grid(time~.)
print(p)

ggsave(file=paste0(subDir, '/diag_transect_key.pdf'))

# prop_dat2 = prop_dat1[which(prop_dat1$taxon %in% c('PINE', 'OAK', 'OTHER.HARDWOOD')), ]
prop_dat2 = prop_dat1[which(prop_dat1$taxon %in% c('PINE', 'OAK', 'OTHER.HARDWOOD', 'TAMARACK', 'BIRCH', 'OTHER.CONIFER')), ]
prop_dat3 = prop_dat2[order(prop_dat2$x),]
prop_dat3 = prop_dat3[prop_dat3$time %in% ages[c(2, 5, 10, 15, 20)],]

foo = rdist(coords_diag[,1:2]*1e6, matrix(cbind(prop_dat3$x, prop_dat3$y), ncol=2))
id_keep  = apply(foo, 2, function(x){if (any(x < 1e-8)){ which(x < 1e-8)}})
prop_dat3 = data.frame(prop_dat3, id=id_keep)

# dat_rect = data.frame(x=numeric(0), y=numeric(0), time=numeric(0), taxon=character(0))
# for (t in 1:T) {
#   for (id in 1:nrow(coords_diag)) {
#     dat_sub = prop_dat3[which((prop_dat3$time == ages[t]) & (prop_dat3$id == id)), ]  
#     dat_rect = rbind(dat_rect, dat_sub[which.max(dat_sub$props),])
#   }
# } 
# 
# dat_rect = rbind(dat_rect, data.frame(0, 0, 0, 3, 'BIRCH', 3))

# droplevels(dat_rect$taxon, except=3)
# jDat <- droplevels(subset(dat_rect, taxon != c(unique(dat_rect$taxon), 'BIRCH')))

p <- ggplot() + 
  # stat_smooth(data=prop_dat3, aes(y = props, x=id, shape=factor(taxon), colour=factor(taxon)), method="gam", formula = y ~ s(x), se=FALSE)+# level=0.5) +
  geom_line(data=prop_dat3, aes(x=id, y=props, colour=factor(taxon)), size=1.6) + 
  geom_point(data=prop_dat3, aes(x=id, y=props, shape=factor(taxon), colour=factor(taxon)), size=2.5) + 
  # geom_rect(data=dat_rect, aes(xmin=id-0.5, xmax=id+0.5, ymin=0, ymax=1, fill=taxon), alpha=0.2) + 
  # scale_colour_discrete(drop=FALSE) + 
  # xlim(8,16) +
    scale_colour_manual(values=tim.colors(6)) +
  #   scale_fill_manual(values=tim.colors(6)) +
#   scale_colour_brewer(palette="Dark2") +
#   scale_fill_brewer(palette="Set1") +
  facet_grid(time~.)
print(p)

ggsave(file=paste0(subDir, '/ecotone_transect_comp.pdf'))

p <- ggplot() + 
  # stat_smooth(data=prop_dat3, aes(y = props, x=id, shape=factor(taxon), colour=factor(taxon)), method="gam", formula = y ~ s(x), se=FALSE)+# level=0.5) +
  geom_line(data=prop_dat3, aes(x=id, y=props, colour=factor(time)), size=1.6) + 
  geom_point(data=prop_dat3, aes(x=id, y=props, shape=factor(time), colour=factor(time)), size=2.5) + 
  # geom_rect(data=dat_rect, aes(xmin=id-0.5, xmax=id+0.5, ymin=0, ymax=1, fill=taxon), alpha=0.2) + 
  # scale_colour_discrete(drop=FALSE) + 
  # xlim(8,16) +
  scale_colour_manual(values=tim.colors(6)) +
  #   scale_fill_manual(values=tim.colors(6)) +
  #   scale_colour_brewer(palette="Dark2") +
  #   scale_fill_brewer(palette="Set1") +
  facet_grid(taxon~.)#, scales="free")
print(p)
ggsave(file=paste0(subDir, '/ecotone_transect_comp_by_taxon.pdf'))








# ndiffs = c(1, 2, 5, 10, 15, 20)
ages
ndiffs = seq(1,20) 

niter=dim(rIts)[4]
diffs_po = array(NA, c(N, length(ndiffs), niter))
for (i in 1:length(ndiffs)){
  diffs_po[,i,] = rIts[,7,ndiffs[i],] - rIts[,10,ndiffs[i],] 
}

avg_po = apply(diffs_po, c(1,2), mean)

po_melted = data.frame(diff=numeric(0), x=integer(0), y=integer(0), t=numeric(0))
for (t in 1:length(ndiffs)){
  po_melted  = rbind(po_melted , data.frame(diff = as.vector(avg_po[,t]), 
                                            x     = centers_veg[,1]*rescale, 
                                            y     = centers_veg[,2]*rescale, 
                                            t     = rep(ages[ndiffs[t]], N)))
}

# values=c(-1, -0.9, -0.899, 0.899, 0.9, 1)
values=c(-1, -0.5, -0.001, 0.001, 0.5, 1)

pdf(file=paste0(subDir, '/post_diff_', r_type, '_PO.pdf'))

p <- ggplot(data=po_melted) + geom_tile(data=po_melted, aes(x=x, y=y, fill=diff)) 
p <- p + scale_fill_gradientn(colours=c("blue", "lightskyblue", "white", "pink", "red"),  
                              values=values, limits=c(-1,1), na.value="white",
                              rescaler = function(x, ...) x, oob = identity) 
# p <- p + scale_fill_gradientn(colours=brewer.pal(11, 'BrBG'),  limits=c(-1,1), na.value="white")#,
#rescaler = function(x, ...) x, oob = identity) 
p <- p + coord_fixed()
p <- add_map_albers(p, map_data=us.fort, limits)
p <- p + facet_wrap(~t, ncol=4)
p <- theme_clean(p)
print(p)

dev.off()


lim=0
p <- ggplot(po_melted, aes(x=x,y=y))+ geom_contour(aes(z=diff, group=t, color=t), breaks=lim) #+ facet_wrap(~time) + geom_tile(aes(fill=props, alpha=0.2))+ 
p <- p + scale_colour_gradientn(colours = tim.colors(10))
p <- add_map_albers(p, map_data=us.fort, limits)
p <- theme_clean(p) + theme(strip.text.x = element_blank(), strip.text.y = element_text(size = rel(1.5)))
# p <- p + annotate("text", y= max(prop_sub$y), x =max(prop_sub$x), label=taxon, hjust=1)
print(p)
ggsave(p, file=paste0(subDir, '/post_ecotone_po.pdf'))

# foo = rdist(coords_diag[,1:2]*1e6, matrix(cbind(po_melted$x, po_melted$y), ncol=2))
# id_keep  = unlist(apply(foo, 2, function(x){if (any(x < 1e-8)){ which(x < 1e-8)} else {NA}}))
# po = data.frame(po_melted, id=id_keep)
# 
# po = po[which(!is.na(po$id)),]
# 
# p <- ggplot() + 
#   stat_smooth(data=po, aes(y = diff, x=id, colour=factor(t)), method="gam", formula = y ~ s(x, bs="cs"), se=FALSE)+# level=0.5) +
#   # geom_line(data=po, aes(x=id, y=diff, colour=factor(t))) + 
#   scale_colour_manual(values=tim.colors(20)) 
# print(p)


##########################################################################################################################################
## horizontal transect
##########################################################################################################################################
# horizontal transect
y_tran = rev(c(0.934, 0.982, 1.030, 1.078))
coords_hz = centers_veg[centers_veg$y %in% y_tran,]

tran = data.frame(x=rep(970000), y_tran = y_tran*1e6, label= seq(1,4))

# make a ggplot map of hemlock range
# then add horizontal transect lines with numbers

# pdf(file=paste0(subDir, '/post_diff_', r_type, '_PO.pdf'))
p <- ggplot(data=subset(prop_dat, taxon == 'HEMLOCK')) + geom_tile(aes(x=x, y=y, fill=props)) 
p <- p + scale_fill_gradientn(colours=c("white", "grey80", "grey58", "grey22"),  
                              values=c(0,0.05,0.15, 1), limits=c(0,1), na.value="white",
                              rescaler = function(x, ...) x, oob = identity, name='Proportion', guide=FALSE) 
p <- p + geom_hline(data=tran, aes(yintercept=y_tran), colour="royalblue", size=1.2, linetype='longdash')
p <- p + geom_text(data=tran, aes(x=x+1e3, y=y_tran+20e3, label=label), size=9) 
p <- p + coord_fixed()
p <- add_map_albers(p, map_data=us.fort, limits)
p <- theme_clean(p)
print(p)
ggsave(p, file=paste0(subDir,'/hz_tran_map.pdf'))

# dev.off()

## hz transect plots of proportions
prop_dat1 = prop_dat[which((prop_dat$y %in% c(y_tran*1e6)) & (prop_dat$taxon == 'HEMLOCK')),]
prop_dat2 = prop_dat1[order(prop_dat1$x),]
prop_dat2$y <- factor(prop_dat2$y)
prop_dat2$y <- factor(prop_dat2$y, levels(prop_dat2$y)[c(4,3,2,1)])
levels(prop_dat2$y) <- seq(1,4)


q <- ggplot() + geom_line(data=prop_dat2, aes(x=x, y=props, colour=factor(time))) + 
  scale_colour_manual(values=tim.colors(20), name='Time') + xlim(limits$xlims*1e6) + ylab('Proportion') + facet_grid(y~.) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        strip.text.y=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))#,
        # legend.position=c(0,0),legend.direction="horizontal")
print(q)
ggsave(q, file=paste0(subDir,'/hz_tran_hemlock.pdf'))

##########################################################################################################################################
## horizontal transect paper
##########################################################################################################################################
# horizontal transect
y_tran = rev(c(0.982))
coords_hz = centers_veg[centers_veg$y %in% y_tran,]

tran = data.frame(x=rep(970000), y_tran = y_tran*1e6)

# make a ggplot map of hemlock range
# then add horizontal transect lines with numbers

# pdf(file=paste0(subDir, '/post_diff_', r_type, '_PO.pdf'))
p <- ggplot(data=subset(prop_dat, taxon == 'HEMLOCK')) + geom_tile(aes(x=x, y=y, fill=props)) 
p <- p + scale_fill_gradientn(colours=c("white", "grey80", "grey58", "grey22"),  
                              values=c(0,0.05,0.15, 1), limits=c(0,1), na.value="white",
                              rescaler = function(x, ...) x, oob = identity, name='Proportion', guide=FALSE) 
p <- p + geom_hline(data=tran, aes(yintercept=y_tran), colour="black", size=1.2, linetype='longdash')
# p <- p + geom_text(data=tran, aes(x=x+1e3, y=y_tran+20e3, label=label), size=9) 
p <- p + coord_fixed()
p <- add_map_albers(p, map_data=us.fort, limits)
p <- p + theme_bw()
p <- theme_clean(p)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  
               panel.border = element_rect(colour = "black", fill=NA, size=0.8))#, panel.border=element_blank())
print(p)
ggsave(p, file=paste0(subDir,'/hz_tran_map_paper.pdf'))

# dev.off()

## hz transect plots of proportions
prop_dat1 = prop_dat[which((prop_dat$y %in% c(y_tran*1e6)) & (prop_dat$taxon == 'HEMLOCK')),]
prop_dat2 = prop_dat1[order(prop_dat1$x),]
prop_dat2$y <- factor(prop_dat2$y)

# prop_dat2$y <- factor(prop_dat2$y, levels(prop_dat2$y)[c(4,3,2,1)])
# levels(prop_dat2$y) <- seq(1,4)

prop_dat3 = data.frame(props=numeric(0), x=numeric(0), y=numeric(0), time=numeric(0), taxon=character(0))
for (t in 1:T){
  prop_sub = prop_dat2[which(prop_dat2$time == ages[t]),]
  prop_sub_ss = spline(prop_sub$x, prop_sub$props, n=201)
  
  prop_dat3 = rbind(prop_dat3, data.frame(props=prop_sub_ss$y, x=prop_sub_ss$x, y=rep(prop_sub$y[1]), time=rep(ages[t]), taxon=rep(prop_sub$taxon[1])))
}

prop_dat2$time = prop_dat2$time * 100
prop_dat3$time = prop_dat3$time * 100

prop_dat2$x = prop_dat2$x / 1000
prop_dat3$x = prop_dat3$x / 1000

q <- ggplot() + geom_point(data=prop_dat2, aes(x=x, y=props, colour=factor(time))) + 
  geom_line(data=prop_dat3, aes(x=x, y=props, colour=factor(time))) + 
  scale_colour_manual(values=tim.colors(20), name='Cal Yr BP') +
  # xlim(limits$xlims*1e6/1000) + 
  ylab('Proportion') + 
  xlab('Coordinate (km from projection origin)') #+
  # theme(axis.text=element_text(size=12),
  #       axis.title=element_text(size=14),
  #       strip.text.y=element_text(size=16),
  #       axis.text.x=element_blank(),
  #       axis.title.x=element_blank(),
  #       axis.ticks.x=element_blank(),
  #       legend.text=element_text(size=16),
  #       legend.title=element_text(size=16))#,
q <- q + theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        # strip.text.y=element_text(size=16),
        # axis.text.x=element_blank(),
        # axis.title.x=element_blank(),
        # axis.ticks=element_text(size=14),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18))#,
# legend.position=c(0,0),legend.direction="horizontal")
print(q)
ggsave(q, file=paste0(subDir,'/hz_tran_hemlock.pdf'))

#A viewport taking up a fraction of the plot area
vp <- viewport(width = 0.43, height = 0.43, x = 0.278, y = 0.779)

#Just draw the plot twice
pdf(paste0(subDir, "/hemlock_transect_panel.pdf"), width=8, height=6)
# png(paste0(subDir, "/hemlock_transect_panel.png"))
print(q)
print(p, vp = vp)
dev.off()

# # arrange p and q
# library(gridExtra)
# 
# multiplot(p,q, nrow=2, widths=c(2,2))
# grid.arrange(p, q, ncol = 1)
# 
# 
# fname = paste0(subDir, '/transect_hz_hemlock.pdf')	
# ggsave(file=fname, scale=1, width=12, height=12)
# sys_str = paste("pdfcrop", fname, fname, sep=' ')
# system(sys_str)

######################################################################################################################################


horizon_transect <- function(dat, coords_hz){
  dat = dat[with(dat, order(x)), ]
  d_diag = rdist(coords_hz[,1:2]*1e6, matrix(cbind(dat$x, dat$y), ncol=2))
  idx_keep = apply(d_diag, 2, function(x) any(x < 1e-8))
  # id_keep  = apply(d_diag, 2, function(x){if(any(x < 1e-10)){ which(x < 1e-10)}})
  dat_tran = dat[idx_keep, ]
  
  return(dat_tran)
}

# 
# can_tran = horizon_transect(can_melted, coords_hz)
# can_tran = can_tran[can_tran$t1 != 1,]
# can_tran = can_tran[can_tran$t1 == 2,]
# can_tran$y <- factor(can_tran$y, levels=sort(unique(can_tran$y), decreasing=TRUE))
# 
# p <- ggplot(data=can_tran) + 
#   geom_line(data=can_tran, aes(y = d, x=x, colour=factor(t2)))+
#   scale_colour_manual(values=tim.colors(4)) +
#   facet_grid(y~.)
# print(p)
# 
# ggsave(file=paste0(subDir, '/post_diff_maps/diss_hz_transect_canberra.pdf'))

bc_tran = horizon_transect(bc_melted, coords_hz)
bc_tran = bc_tran[bc_tran$t1 != 1,]
bc_tran = bc_tran[bc_tran$t1 == 2,]
bc_tran$y <- factor(bc_tran$y, levels=sort(unique(bc_tran$y), decreasing=TRUE))

p <- ggplot(data=bc_tran) + 
  geom_line(data=bc_tran, aes(y = d, x=x, colour=factor(t2)))+
  scale_colour_manual(values=tim.colors(4)) +
  facet_grid(y~.)
print(p)

ggsave(file=paste0(subDir, '/diss_hz_transect_bc.pdf'))

##################################################33

prop_dat = data.frame(props=numeric(0), props_lb=numeric(0), props_lb=numeric(0), x=integer(0), y=integer(0),time=integer(0), taxon=character())
for (k in 1:K){
  prop_dat = rbind(prop_dat, data.frame(props = r_quants[2,,k], 
                                        props_lb = r_quants[1,,k], 
                                        props_ub = r_quants[3,,k], 
                                        x     = rep(centers_veg[,1]*1000000, each=T), 
                                        y     = rep(centers_veg[,2]*1000000, each=T), 
                                        time  = rep(ages, times=N), 
                                        taxon = rep(taxa[k], N*T)))
}

prop_dat = prop_dat[with(prop_dat, order(x)), ]
prop_dat1 = prop_dat[which((prop_dat$y %in% c(y_tran*1e6)) & (prop_dat$taxon == 'HEMLOCK')),]
prop_dat2 = prop_dat1[order(prop_dat1$x),]
prop_dat2$y <- factor(prop_dat2$y)
prop_dat2$y <- factor(prop_dat2$y, levels(prop_dat2$y)[c(4,3,2,1)])
# levels(prop_dat2$y) <- sort( as.numeric(levels(prop_dat2$y)), decreasing=TRUE) 

prop_dat2 = prop_dat2[which(prop_dat2$time %in% ages[c(1, 5, 10, 20)]),]
prop_dat2 = prop_dat2[which(prop_dat2$y %in% c(1078000, 982000)),]

p <- ggplot() + geom_line(data=prop_dat2, aes(x=x, y=props, colour=factor(time))) + 
  geom_ribbon(data=prop_dat2, aes(x=x, ymin=props_lb, ymax=props_ub, fill=factor(time)), alpha=0.3) +
  scale_colour_manual(values=tim.colors(5)) + scale_fill_manual(values=tim.colors(5)) + facet_grid(y~.)
print(p)

#   p <- ggplot() + geom_bar(data=prop_dat2, aes(x=x, y=props, colour=factor(time))) + 
#     scale_colour_manual(values=tim.colors(20)) + facet_grid(y~.)
#   print(p)

#   p <- theme_clean(p) + theme(strip.text.x = element_blank(), strip.text.y = element_text(size = rel(1.5)))
# p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.3)), strip.text.x = element_text(size = rel(1.0)))
#theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
#   print(p)

fname = paste0(subDir, '/transect_hz_hemlock_ribbon.pdf')	
ggsave(file=fname, scale=1, width=12, height=12)
sys_str = paste("pdfcrop", fname, fname, sep=' ')
system(sys_str)

##########################################################################################################################################
## pollen/veg diagrams by pollen site
##########################################################################################################################################
hips_coords = read.table('data/HIPS_coords.csv', sep=',', header=TRUE)[3:5,]

#XXX: add pollen to albers function
hips_coordsa = data.frame(cbind(hips_coords$long, hips_coords$lat))
colnames(hips_coordsa) = c('x', 'y')

coordinates(hips_coordsa) <- ~ x + y
proj4string(hips_coordsa) <- CRS('+proj=longlat +ellps=WGS84')

hips_coordsa <- spTransform(hips_coordsa, CRS('+init=epsg:3175'))
hips_coordsa <- as.matrix(data.frame(hips_coordsa))/1e6

dmat=rdist(hips_coordsa, centers_veg)
idx_hips = apply(dmat, 1, which.min)


foo = rdist(hips_coordsa, centers_pol)
idx_hips_pol = apply(foo, 1, which.min)

# cell_id = 415

make_props <- function(x){
  if (sum(x) == 0){
    return(x)
  } else {
    return(x/sum(x))
  }
}


library(gridExtra)
library(reshape2)

# add the site id / site name



pdf(file=paste0(subDir, '/site_pollen_diagrams.pdf'), width=12)

for (i in 1:N_cores){
  
  print(i)
#   cell_id = idx_hips[i]
#   site = hips_coords$site[i]
  
  cell_id = idx_cores[i]
  site = i
  
  # idx = seq((cell_id-1)*T, (cell_id-1)*T + T-1)
  
  r_s = rIts[cell_id,,,]
  
  r_site = array(NA, dim=c(T, K, 4))
  for (t in 1:T){
    # r_site[t,,] = t(apply(r_s[,t,], 1, function(x) quantile(x, probs=c(0.1, 0.5, 0.9)))) 
    r_site[t,,] = t(apply(r_s[,t,], 1, function(x) c(mean(x), quantile(x, probs=c(0.1, 0.5, 0.9))))) 
  }
  
  r_melt = melt(r_site, varnames=c("time", "taxon", "quant"))
  r_melt = r_melt[order(r_melt$time),]
  
  #r_melt = r_melt[which(r_melt$quant == 2),]
  
  r_melt$taxon <- as.factor(r_melt$taxon)
  levels(r_melt$taxon) = taxa
  r_melt = r_melt[r_melt$quant %in% c(1, 2, 4),]
  
  # get the pollen data
  idx_pol = seq((i-1)*T+1, i*T)
  y_p = apply(y[idx_pol,], 1, make_props)
  
  t_pol = which(colSums(y_p) != 0)
  # y_p = y_p[,t_pol]
  
  y_melt = melt(y_p, varnames=c("taxon", "time"))
  y_melt = y_melt[y_melt$time %in% t_pol,]
  y_melt$time = ages[y_melt$time]
  
  y_melt$taxon <- as.factor(y_melt$taxon)
  levels(y_melt$taxon) = taxa
  
  
  # pdf(file=sprintf('runs/%s/run1/figures/cell%s_change.pdf', run$suff_run, cell_id))
  p <- ggplot() + geom_line(data=r_melt, aes(x=time, y=value, 
                                             linetype=as.factor(quant)), colour="black", size=0.7) + 
    scale_linetype_manual(name="preds", values=c( "solid", "dashed","dashed"), guide=FALSE)
  
  if (length(t_pol) == 1) {
    p <- p + geom_point(data=y_melt, aes(x=time, y=value), colour='blue', size=2)
  } else if (length(t_pol) > 1) {
    p <- p + geom_line(data=y_melt, aes(x=time, y=value), colour='blue', size=0.7) 
  }
  
  p <- p + scale_colour_manual(name='col', values=c("black", "blue"), guide=FALSE)+
    ylim(0,1) + scale_x_reverse() + coord_flip() 
  p <- p + facet_wrap(~taxon, ncol=6)
  p <- p + ggtitle(paste0(meta_pol[idx_pol[1],'sitename'], '; stat_id = ', site)) + ylab("Proportion") + xlab("Time")
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  #print(p)
  # ggsave(p, file=sprintf('runs/%s/run1/cell%s_change.pdf', run$suff_run, cell_id))
  
  loc = data.frame(centers_veg[cell_id,]*1e6)
  limits = get_limits(centers_pls)
  
  # where is this point
  q <- ggplot() + geom_point(data=loc, aes(x=x, y=y), colour='#FF6600', shape=19, size=5) + 
    coord_fixed() #+ scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
  q <- add_map_albers(q, map_data=us.fort, limits)#, xlims=c(xlo,xhi), ylims=c(ylo, yhi))
  q <- q + theme(strip.text.x = element_blank(),
                 strip.text.y = element_blank())
  q <- q + theme(strip.background = element_blank())
  q <- theme_clean(q) #+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  #print(q)

  grid.arrange(p,q, ncol=2, widths=c(2.5,1.5))
  
  # dev.off()
}

dev.off()

# hips = data.frame(hips_coordsa*1e6, site=hips_coords$site)
# 
# # where is this point
# p <- ggplot() + geom_point(data=hips, aes(x=x, y=y, colour=site, label=site), shape=19, size=7) + 
#   coord_fixed() #+ geom_text(data=hips, aes(x=x+8000, y=y+8000, colour=site, label=site), hjust=0, vjust=0)#+ scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
# p <- add_map_albers(p, map_data=us.fort, limits)#, xlims=c(xlo,xhi), ylims=c(ylo, yhi))
# p <- p + theme(strip.text.x = element_blank(),
#                strip.text.y = element_blank())
# p <- p + theme(strip.background = element_blank())
# p <- theme_clean(p) #+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
# print(p)
# ggsave(file=sprintf('runs/%s/run1/figures/hips_site_map.pdf', run$suff_run))
# 



# path_pollen   = '../stepps-data/data/bacon_ages/pollen_ts_bacon_meta_v4.csv'
# pollen_ts = read.table(path_pollen, header=TRUE, sep=',', stringsAsFactors=FALSE)


# add all pollen samples from 10 runs
y_all = list(length(10))
core_idx = list(length(10))
for (r in 1:10){
  load(paste0('runs/', run$suff_run, '/run', r, '/input.rdata')) 
  y_all[[r]] = y
  
  foo=rdist(as.matrix(centers_pol), as.matrix(meta_pol[,c(4,3)]))
  core_idx[[r]] = unlist(apply(foo, 1, function(x) which(x < 1e-8)))
  
  print(N_cores)
  print(dim(pollen_check))
}

# y_p = sapply(y_all, make_props)

pdf(file=paste0(subDir, '/site_pollen_var.pdf'), width=12)

for (i in 1:10){#nrow(pollen_check)){
  
  print(i)
  #   cell_id = idx_hips[i]
  #   site = hips_coords$site[i]
  
  y_melt = data.frame(time=numeric(0), taxon=numeric(0), value=numeric(0), run=numeric(0))
  for (r in 1:nruns){
  
    load(paste0('runs/', run$suff_run, '/run', r, '/input.rdata')) 
    
    if (!any(core_idx[[r]] == i)){
      next
    } else { 
      run_idx = which(core_idx[[r]] == i)
    }
    
    cell_id = idx_cores[run_idx]
    
    # get the pollen data
    idx_pol = seq((run_idx-1)*T+1, run_idx*T)
    
      
    y_sub = melt(apply(y_all[[r]][idx_pol,], 1, make_props))
    
    # y_sub = melt(y_p[[r]][idx_pol,])
    colnames(y_sub) = c('taxon', 'time', 'value')
    
    t_pol = which(rowSums(y_all[[r]][idx_pol,]) != 0)
    y_sub = y_sub[y_sub$time %in% t_pol,]
    
    y_melt = rbind(y_melt, data.frame(y_sub, run=rep(r, nrow(y_sub))))
  }
  
  y_melt$time = ages[y_melt$time]
  
  y_melt$taxon <- as.factor(y_melt$taxon)
  levels(y_melt$taxon) = taxa
  
  
  
  # pdf(file=sprintf('runs/%s/run1/figures/cell%s_change.pdf', run$suff_run, cell_id))
  p <- ggplot() 
  if (length(unique(y_melt$time)) == 1) {
    p <- p + geom_point(data=y_melt, aes(x=time, y=value, colour=factor(run)), size=2)
  } else if (length(unique(y_melt$time)) > 1) {
    p <- p + geom_line(data=y_melt, aes(x=time, y=value, colour=factor(run)), size=0.7) 
  }
  
  # p <- p + scale_colour_manual(name='col', values=c("black", "blue"), guide=FALSE)+
  p <- p +  ylim(0,1) + scale_x_reverse() + coord_flip() 
  p <- p + facet_wrap(~taxon, ncol=6)
  p <- p + ggtitle(site) + ylab("Proportion") + xlab("Time")
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p)
  # ggsave(p, file=sprintf('runs/%s/run1/cell%s_change.pdf', run$suff_run, cell_id))
  
#   loc = data.frame(centers_veg[cell_id,]*1e6)
#   limits = get_limits(centers_pls)
#   
#   # where is this point
#   q <- ggplot() + geom_point(data=loc, aes(x=x, y=y), colour='#FF6600', shape=19, size=5) + 
#     coord_fixed() #+ scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
#   q <- add_map_albers(q, map_data=us.fort, limits)#, xlims=c(xlo,xhi), ylims=c(ylo, yhi))
#   q <- q + theme(strip.text.x = element_blank(),
#                  strip.text.y = element_blank())
#   q <- q + theme(strip.background = element_blank())
#   q <- theme_clean(q) #+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
#   #print(q)
#   
#   grid.arrange(p,q, ncol=2, widths=c(2.5,1.5))
  
  # dev.off()
}

dev.off()


#######################################################################################################################################
# y_p = sapply(y_all, make_props)

pp = pollen_preds_sp(phi, gamma, w, sum_w_pot, d, idx_cores, r_mean)

pdf(file=paste0(subDir, '/site_pollen_diagrams_pp.pdf'), width=12)

for (i in 1:N_cores){
  
  #   cell_id = idx_hips[i]
  #   site = hips_coords$site[i]
  
  cell_id = idx_cores[i]
  site = i
  
  # idx = seq((cell_id-1)*T, (cell_id-1)*T + T-1)
  
  r_s = rIts[cell_id,,,]
  
  r_site = array(NA, dim=c(T, K, 4))
  for (t in 1:T){
    # r_site[t,,] = t(apply(r_s[,t,], 1, function(x) quantile(x, probs=c(0.1, 0.5, 0.9)))) 
    r_site[t,,] = t(apply(r_s[,t,], 1, function(x) c(mean(x), quantile(x, probs=c(0.1, 0.5, 0.9))))) 
  }
  
  r_melt = melt(r_site, varnames=c("time", "taxon", "quant"))
  r_melt = r_melt[order(r_melt$time),]
  
  #r_melt = r_melt[which(r_melt$quant == 2),]
  
  r_melt$taxon <- as.factor(r_melt$taxon)
  levels(r_melt$taxon) = taxa
  r_melt = r_melt[r_melt$quant %in% c(1, 2, 4),]
  
  # get the pollen data
  idx_pol = seq((i-1)*T+1, i*T)
  y_p = apply(y[idx_pol,], 1, make_props)
  
  t_pol = which(colSums(y_p) != 0)
  # y_p = y_p[,t_pol]
  
  y_melt = melt(y_p, varnames=c("taxon", "time"))
  y_melt = y_melt[y_melt$time %in% t_pol,]
  y_melt$time = ages[y_melt$time]
  
  y_melt$taxon <- as.factor(y_melt$taxon)
  levels(y_melt$taxon) = taxa
  
  yp_melt = melt(pp[,i,], varnames=c("taxon", "time"))
  yp_melt$time = ages[yp_melt$time]
  yp_melt$taxon <- as.factor(yp_melt$taxon)
  levels(yp_melt$taxon) = taxa
  
  # pdf(file=sprintf('runs/%s/run1/figures/cell%s_change.pdf', run$suff_run, cell_id))
  p <- ggplot() + geom_line(data=r_melt, aes(x=time, y=value, 
                                             linetype=as.factor(quant)), colour="black", size=0.7) + 
    scale_linetype_manual(name="preds", values=c("solid", "dashed", "dashed"), guide=FALSE)
  
  if (length(t_pol) == 1) {
    p <- p + geom_point(data=y_melt, aes(x=time, y=value), colour='blue', size=2)
  } else if (length(t_pol) > 1) {
    p <- p + geom_line(data=y_melt, aes(x=time, y=value), colour='blue', size=0.7) 
  }
  
  p <- p + geom_line(data=yp_melt, aes(x=time, y=value), colour='red', size=0.7) 
  
  p <- p + scale_colour_manual(name='col', values=c("black", "blue"), guide=FALSE)+
    ylim(0,1) + scale_x_reverse() + coord_flip() 
  p <- p + facet_wrap(~taxon, ncol=6)
  p <- p + ggtitle(paste0(meta_pol[idx_pol[1],'sitename'], '; stat_id = ', site)) + ylab("Proportion") + xlab("Time")
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  #print(p)
  # ggsave(p, file=sprintf('runs/%s/run1/cell%s_change.pdf', run$suff_run, cell_id))
  
  loc = data.frame(centers_veg[cell_id,]*1e6)
  limits = get_limits(centers_pls)
  
  # where is this point
  q <- ggplot() + geom_point(data=loc, aes(x=x, y=y), colour='#FF6600', shape=19, size=5) + 
    coord_fixed() #+ scale_x_continuous(limits$xlims) + scale_y_continuous(limits$ylims)
  q <- add_map_albers(q, map_data=us.fort, limits)#, xlims=c(xlo,xhi), ylims=c(ylo, yhi))
  q <- q + theme(strip.text.x = element_blank(),
                 strip.text.y = element_blank())
  q <- q + theme(strip.background = element_blank())
  q <- theme_clean(q) #+ theme(strip.text.y = element_text(size = rel(1.5)), strip.text.x = element_text(size = rel(1.5)))
  #print(q)
  
  grid.arrange(p,q, ncol=2, widths=c(2.5,1.5))
  
  # dev.off()
}

dev.off()



#######################################################################################################################################
## core locations and intervals
#######################################################################################################################################
# meta_pol = readRDS(file=paste0(subDir, '/pollen_meta.RDS'))
        
# plot_core_intervals(y, centers_pol, centers_pls, ages, limits, fpath=subDir)
# plot_core_intervals2(y, centers_pol, centers_pls, ages, limits, fpath=subDir)


pollen_ts = readRDS(file='data/pollen_ts.RDS')
# pollen_ts = pollen_ts[pollen_ts$id %in% meta_pol$id,]
# pollen_ts$stat_id = meta_pol$stat_id[match(pollen_ts$id, meta_pol$id)]
# # only plot cores used in model
# plot_core_intervals3(pollen_ts, fpath=subDir)
# plot_core_intervals3(subset(meta_pol, age_bacon>=150), fpath=subDir)


states_pol = c('minnesota', 'wisconsin', 'michigan:north')
pollen_ts_new = pollen_ts[which(pollen_ts$state %in% states_pol),]
plot_core_intervals3(pollen_ts_new, fpath=subDir)


plot_core_locations_select(y, centers_pol, centers_pls, ages, idx.keep=c(1,2,5,10,15,20), limits, suff='', fpath=subDir)

sj <- data.frame(read.csv(file='data/lake_coordinates_sj.csv', header=TRUE))
ms_dd <- function(d, m, s, dd = as.numeric(d), mm = as.numeric(m), 
              ss = as.numeric(s)) sign(dd) * (abs(dd) + mm / 60 + ss / 3600) 
for (i in 1:nrow(sj)){
    lat = sj[i,3:5]
    sj$lat[i] = ms_dd(lat[1], lat[2], lat[3])
}
for (i in 1:nrow(sj)){
  long = sj[i,6:8]
  sj$long[i] = ms_dd(long[1], long[2], long[3])
}

centers_sja = data.frame(cbind(x=-sj$long, y=sj$lat))

coordinates(centers_sja) <- ~ x + y
proj4string(centers_sja) <- CRS('+proj=longlat +ellps=WGS84')

centers_sja <- spTransform(centers_sja, CRS('+init=epsg:3175'))
centers_sja <- as.matrix(data.frame(centers_sja))/1000000

centers_sja = data.frame(centers_sja, 
                         names=sj[,1], 
                         handle=c('AL', 'LA', 'PL', 'SK', 'TL', 'YL', 
                                             'BL', 'HK', 'WL', 'BL', 'GL'))



# #XXX: add pollen to albers function
# hips_coordsa = data.frame(cbind(hips_coords$long, hips_coords$lat))
# colnames(hips_coordsa) = c('x', 'y')
# 
# coordinates(hips_coordsa) <- ~ x + y
# proj4string(hips_coordsa) <- CRS('+proj=longlat +ellps=WGS84')
# 
# hips_coordsa <- spTransform(hips_coordsa, CRS('+init=epsg:3175'))
# hips_coordsa <- as.matrix(data.frame(hips_coordsa))/1e6

plot_core_locations_select_sj(y, centers_pol, centers_sja, centers_pls, ages, 
                              idx.keep  = c(1,2,length(ages)/2,T), limits, fpath=subDir, 
                              save_plots=TRUE)

#######################################################################################################################################
## core locations and intervals
#######################################################################################################################################

covs <- read.csv(paste0('data/StatsGO2_soil_covariates.csv'))
colnames(covs)
coords_cov = cbind(covs$x, covs$y)

covs_in = data.frame(matrix(NA, nrow=N, ncol=ncol(covs)))
colnames(covs_in) = colnames(covs)
for (i in 1:N){
  print(i)
  d = rdist(as.matrix(centers_veg[i,]*1e6), coords_cov)
  idx = which(d<=sqrt(8000^2*2))
  if (length(idx)>1){
    covs_in[i,] = colMeans(covs[idx,]) 
  } else if (length(idx)==1) {
    covs_in[i,] = covs[idx,] 
  } else{
    next
  }
  # idx = which.min(d)
  # covs_in[i,] = covs[idx,] 
}

r_mean_s = r_mean
r_mean_s[r_mean_s < 0.05] = NA

for (i in 1:12){
  par(oma=c(0,0,2,0))
  par(mfrow=c(2,2))
  plot(covs_in$sand, r_mean_s[,i,1])
  plot(covs_in$silt, r_mean_s[,i,1])
  plot(covs_in$ksat, r_mean_s[,i,1])
  plot(covs_in$clay, r_mean_s[,i,1])
  title(main=taxa[i],outer=T)
}


p <- ggplot(data=covs) + geom_raster(aes(x=x, y=y, fill=sand)) 
#   p <- p + scale_fill_gradientn(colours=c("blue", "lightskyblue", "white", "pink", "red"),  
#                                 values=values, limits=c(-0.6,0.6), na.value="white",
#                                 rescaler = function(x, ...) x, oob = identity) 
p <- p + scale_fill_gradientn(colours=c("blue", "lightskyblue", "white", "pink", "red"),  
                              values = c(0,20,80,100),
                              limits=c(0,100), na.value="white",
                              rescaler = function(x, ...) x, oob = identity) 
# p <- p + scale_fill_gradientn(colours=tim.colors(8))
# p <- p + scale_fill_gradientn(colours=brewer.pal(8,"YlOrBr"))
p <- p + coord_fixed()
p <- add_map_albers(p, map_data=us.fort, limits)
p <- theme_clean(p) # theme(strip.text.x = element_blank(), strip.text.y = element_text(size = rel(1.5)))
print(p)