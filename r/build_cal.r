library(fields)
library(rstan)
library(sp)
library(rgdal)
library(ggplot2)
library(mvtnorm)
library(maptools)
library(maps)
library(plyr)

source('r/utils/pred_helper_funs.r')

######################################################################################################################################
# user defs
######################################################################################################################################

# stat model flags
decomp     = TRUE
bt         = TRUE
mpp        = TRUE
save_plots = TRUE

# us albers shape file
# us.shp <- readShapeLines('attic/r/data/map_data/us_alb.shp',
#                          proj4string=CRS('+init=epsg:3175'))
us.shp <- readOGR('data/map_data/us_alb.shp', 'us_alb')

# date of pollen_ts pull
mydate = '2014-07-22'

# use the veg knots? set to false for smaller test data sets
veg_knots = TRUE

n_cells = NA
# n_cells = seq(1,100)

# grid 
res  = res
side = side
# side = '' # 'E', 'W', or ''
# grid = 'MISP' 
grid = 'umw'
gridname = paste0(grid, side, '_', as.character(res), 'by')
#gridname = 'umwE_3by'

# reconstruction limits and bin-width
tmin = 50
tmax = 150
int  = 100

# rescale
rescale = 1e6

# knots
# nclust = 75
# clust.ratio = 6# approx clust.ratio knots per cell
clust.ratio = 7# approx clust.ratio knots per cell
clust.ratio = 15# approx clust.ratio knots per cell

# suff=''
version="v0.3"
suff    = paste(gridname, '_', version, sep='') 
# suff = '3by_v0.3_test'

# states_pol = c('minnesota', 'wisconsin', 'michigan:north', 'michigan:south')
# states_pls = c('minnesota', 'wisconsin', 'michigan:north', 'michigan:south')

states_pol = c('minnesota', 'wisconsin', 'michigan:north')
states_pls = c('minnesota', 'wisconsin', 'michigan:north')


# specify the taxa to use
# must be from the list: taxa_sub = c('oak', 'pine', 'maple', 'birch', 'tamarack', 'beeh', 'elm', 'spruce', 'ash', 'hemlock')
# always have: 'other.hardwood' and 'other.conifer'

taxa_all = toupper(c('oak', 'pine', 'maple', 'birch', 'tamarack', 'beech', 'elm', 'spruce', 'ash', 'hemlock'))
taxa_sub = toupper(c('oak', 'pine', 'maple', 'birch', 'tamarack', 'beech', 'elm', 'spruce', 'ash', 'hemlock'))
#taxa_sub = toupper(c('oak', 'pine', 'maple', 'birch', 'tamarack', 'beech'))#, 'elm', 'spruce', 'ash', 'hemlock'))
# taxa_sub = toupper(c('oak', 'pine', 'maple', 'birch'))

K = as.integer(length(taxa_sub) + 1)
W = K-1

##########################################################################################################################
## read in tables and data
##########################################################################################################################

# conversion tables
tree_type = read.table('data/assign_HW_CON.csv', sep=',', row.names=1, header=TRUE)
convert   = read.table('data/dict-comp2stepps.csv', sep=',', row.names=1, header=TRUE)

pls.raw = data.frame(read.table(file='data/western_comp_stepps_v0.3-1.csv', sep=",", row.names=NULL, header=TRUE))

# read in grid
load(file=paste('data/grid/', gridname, '.rdata', sep=''))

# read in calibration output
# load('calibration/r/dump/cal_data_12taxa_mid_comp_all.rdata')
# cal_fit = read_stan_csv('data/calibration_output/12taxa_mid_comp_long.csv')

cal_fit = rstan::read_stan_csv(paste0('data/calibration_output/', run$suff_fit,'.csv'))

# read in veg data and output
# veg data specifies which knots to use
suff_data = '12taxa_6341cells_120knots_v0.3'
load(paste('data/veg_data_', suff_data, '.rdata', sep=''))
load(file="data/12taxa_6341cells_53knots_cont_od_mpp.rdata")
y_veg = y

pollen_ts = read.table(paste('data/pollen_ts_', mydate, '.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)
load('data/calibration/cal_data_12taxa_mid_comp_ALL_v0.3.rdata')

##########################################################################################################################
## read in and organize pls data
##########################################################################################################################

colnames(pls.raw) = tolower(colnames(pls.raw))

# pull the subset of proportions
taxa.start.col = min(match(tolower(rownames(convert)), colnames(pls.raw)), na.rm=TRUE)

pls_dat  = pls.raw[,taxa.start.col:ncol(pls.raw)]
colnames(pls_dat) = as.vector(convert[match(colnames(pls_dat), tolower(rownames(convert))),1])
pls_dat_collapse  = sapply(unique(colnames(pls_dat)), 
                           function(x) rowSums( pls_dat[ , grep(x, names(pls_dat)), drop=FALSE]) )
counts = data.frame(pls_dat_collapse[,sort(colnames(pls_dat_collapse))])
meta   = pls.raw[,1:(taxa.start.col-1)]
# kilometers
# pls$X = pls$X/1000
# pls$Y = pls$Y/1000

meta        = split_mi(meta)
counts      = counts[which(meta$state2 %in% states_pls),]
meta        = meta[which(meta$state2 %in% states_pls),]
# if (length(num_cells) > 1){
#   counts   = counts[cells,]
#   meta = meta[cells,]
# }

centers_pls = data.frame(x=meta$x, y=meta$y)/rescale # megameters!
plot(centers_pls[,1]*rescale, centers_pls[,2]*rescale, asp=1, axes=F,  col='antiquewhite4', xlab='',ylab='', pch=19, cex=0.2)
plot(us.shp, add=T)

y_veg = convert_counts(counts, tree_type, taxa_sub)

taxa = colnames(y_veg)
y_veg = as.matrix(round(unname(y_veg)))
rownames(y_veg) = NULL
y_veg = unname(y_veg)
# y = y_build(counts, taxa_sub) # fix this if we want to use a subset of taxa

K = as.integer(ncol(y_veg))
W = K-1
N_pls = nrow(y_veg)

# make sure columns are in order! 
# y_veg = y_veg[,taxa]

##########################################################################################################################
## chunk: read in coarse grid and pollen data
##########################################################################################################################

# FIXME: ADD STATE TO GRID
# coarse_domain  = coarse_domain[coarse_domain$state %in% states_pls,]
coarse_centers = domain[,1:2]
if (length(n_cells) > 1){
  coarse_centers = coarse_centers[n_cells,]
}

plot(coarse_centers[,1]*rescale, coarse_centers[,2]*rescale, col='blue')
plot(us.shp, add=TRUE)

# assign grid to centers_veg
centers_veg = coarse_centers
N = nrow(centers_veg)

# subdomain boundaries
xlo = min(centers_veg$x)
xhi = max(centers_veg$x)
ylo = min(centers_veg$y)
yhi = max(centers_veg$y)

##########################################################################################################################
## chunk: reorganize pollen data
##########################################################################################################################

# ## pollen data!
# pollen_ts1 = pollen_ts[which((pollen_ts$ages <= 2500) & (pollen_ts$state %in% states_pol)),]
# 
# # reproject pollen coords from lat long to Albers
# pollen_ts2 <- pollen_to_albers(pollen_ts1)
# 
# # pollen_ts = pollen_ts[which((pollen_ts[,'x'] <= xhi) & (pollen_ts[,'x'] >= xlo) & 
# #                             (pollen_ts[,'y'] <= yhi) & (pollen_ts[,'y'] >= ylo)),]
# 
# pollen_locs = cbind(pollen_ts2$x, pollen_ts2$y)
# # pollen_int  = knots_in_domain4(unique(pollen_locs), centers_veg, cell_width = res*8000/rescale)
# # 
# # idx_pollen_int = apply(pollen_locs, 1, function(x) if (any(rdist(x, pollen_int) < 1e-8)) {return(TRUE)} else {return(FALSE)})
# # pollen_ts = pollen_ts[idx_pollen_int, ]
# 
pollen_int  = cores_near_domain(centers_polA/rescale, centers_veg, cell_width = res*8000/rescale)
# 
idx_pollen_int = apply(centers_polA/rescale, 1, 
                       function(x) if (any(rdist(x, pollen_int) < 1e-8)) {return(TRUE)} else {return(FALSE)})
centers_pol = centers_polA[idx_pollen_int, ]/rescale
colnames(centers_pol) = c('x', 'y')

y = y[idx_pollen_int,taxa_all %in% taxa_sub]

# 
# # check how does splitting affects weights... 
# pollen_check = pollen_ts2[,1:7]
# pollen_check$int = rep(FALSE, nrow(pollen_check))
# pollen_check$int[which(idx_pollen_int == TRUE)] = TRUE
# pollen_check=pollen_check[!duplicated(pollen_check),]

# plot domain and core locations 
par(mfrow=c(1,1))
plot(centers_veg$x*rescale, centers_veg$y*rescale)
points(centers_pol[,1]*rescale, centers_pol[,2]*rescale, col='blue', pch=19)
plot(us.shp, add=T, lwd=2)

# plot domain and core locations 
par(mfrow=c(1,1))
plot(centers_pol$x*rescale, centers_pol$y*rescale, col='blue', pch=19)
points(centers_veg$x*rescale, centers_veg$y*rescale)
plot(us.shp, add=T, lwd=2)

##########################################################################################################################
## chunk: prepare pollen data; aggregate over time intervals
##########################################################################################################################
# 
# # sum counts over int length intervals
# pollen_agg = build_pollen_counts(tmin=tmin, tmax=tmax, int=int, pollen_ts=pollen_ts3, taxa_all, taxa_sub)
# #pollen_agg = build_pollen_counts_fast_core(tmin=tmin, tmax=tmax, int=int, pollen_ts=pollen_ts)

meta_pol   = centers_pol
counts     = y

ages    = 150
T       = length(ages) 
lag     = unname(as.matrix(dist(matrix(ages), upper=TRUE)))
N_cores = nrow(meta_pol)

y = unname(y)

# indices for which cells the cores fall in
idx_cores <- build_idx_cores(centers_pol, centers_veg, N_cores)

plot(centers_veg$x*rescale, centers_veg$y*rescale, col='lightgrey')
points(centers_veg[idx_cores,'x']*rescale, centers_veg[idx_cores,'y']*rescale, col='red', pch=19)
points(centers_pol$x*rescale, centers_pol$y*rescale, col='blue', pch=4, cex=1.4)
plot(us.shp, add=TRUE)

##########################################################################################################################
## chunk 3: build distance matrices
##########################################################################################################################

if (!veg_knots){
  nclust = ceiling(N/clust.ratio)
  d_out = build_domain_objects(centers_veg, dx=20, cell_width=8, nclust=nclust)
  
  d = d_out$d
  # d_knots = d_out$d_knots
  # d_inter = d_out$d_inter
  # 
  knot_coords = d_out$knot_coords
} else {

  if (side == '') {
    # don't touch knot_coords
  } else {
    if (side == 'W'){
      if (res == 1) 
        cutlines = list(list(c(0.42, 1.0), c(0.0, 1.0)), list(c(0.397,1.15), c(0.168,0.119)))
      if (res == 5) 
        cutlines = list(list(c(0.386, 1.0), c(0.0, 1.0)), list(c(0.397,1.15), c(0.168,0.119)))
      if (res == 3) 
        cutlines = list(list(c(0.405, 1.0), c(0.0, 1.0)), list(c(0.397,1.15), c(0.168,0.119)))
    } else if (side == 'E'){ 
      if (res %in% c(1, 5)) cutlines = list(list(c(0.253, 1.0), c(0.0, -1.0)))
      if (res == 3) cutlines = list(list(c(0.27, 1.0), c(0.0, -1.0)))
    } 
    idx = choppy(knot_coords[,1], knot_coords[,2], cutlines)
    knot_coords = knot_coords[idx,]
#     knot_coords2 = knot_coords[idx,]
  } 
}

# 
# # knot_coords3 = knots_in_domain4(knot_coords, centers_veg, cell_width = res*8000/rescale)

# plot(domain[,1], domain[,2], asp=1)
plot(centers_veg[,1], centers_veg[,2], asp=1)
points(knot_coords[,1], knot_coords[,2], col='blue', pch=19)
# points(knot_coords2[,1], knot_coords2[,2], col='green', pch=19)

d = rdist(centers_veg, centers_veg)
diag(d) <- 0

d_knots = rdist(knot_coords, knot_coords)
diag(d_knots) <- 0

d_inter = rdist(centers_veg, knot_coords)
d_inter[which(d_inter<1e-8)]=0

d_pol = rdist(centers_pol, centers_veg)
d_pol[which(d_pol<1e-8)]=0

N_knots     = nrow(knot_coords)

##########################################################################################################################
## chunk: qr decompose X
##########################################################################################################################
KW     = FALSE
KGAMMA = FALSE

kernel    = run$kernel
post      = rstan::extract(cal_fit, permuted=FALSE, inc_warmup=FALSE)
col_names = colnames(post[,1,])
par_names  = unlist(lapply(col_names, function(x) strsplit(x, "\\[")[[1]][1]))

phi    = unname(colMeans(post[,1,which(par_names == 'phi')])[1:K])

one_gamma = run$one_gamma
if (one_gamma){
#   gamma = rep(mean(post[,1,which(par_names == 'gamma')]), K)
  gamma = unname(mean(post[,1,which(par_names == 'gamma')]))
} else {
  KGAMMA = TRUE
  gamma = unname(colMeans(post[,1,which(par_names == 'gamma')])[1:K])
}

if (kernel=='gaussian'){
  one_psi = run$one_psi
  if (one_psi){
#     psi   = rep(mean(post[,1,which(par_names == 'psi')]), K)
    psi   = unname(mean(post[,1,which(par_names == 'psi')]))
  } else {
    KW = TRUE
    psi   = unname(colMeans(post[,1,which(par_names == 'psi')])[1:K])
  }
} else if (kernel=='pl'){
  one_a = run$one_a
  if (one_a){
#     a = rep(mean(post[,1,which(par_names == 'a')]), K)
    a = unname(mean(post[,1,which(par_names == 'a')]))
  } else {
    KW = TRUE
    a = unname(colMeans(post[,1,which(par_names == 'a')])[1:K])
  }
  
  one_b = run$one_b
  if (one_b){
#     b = rep(mean(post[,1,which(par_names == 'b')]), K)
    b = unname(mean(post[,1,which(par_names == 'b')]))
  } else {
    KW = TRUE
    b = unname(colMeans(post[,1,which(par_names == 'b')])[1:K])
  }
}

w <- build_weight_matrix(post, d_pol, idx_cores, N, N_cores, run)
# head(apply(w, 1, rowSums))
#####################################################################################
# calculate potential d
# used to determine C normalizing constant in the non-local contribution term
#####################################################################################

library(plyr)

# xlo = min(centers_veg[,1])
# xhi = max(centers_veg[,1])
# ylo = min(centers_veg[,2])
# yhi = max(centers_veg[,2])

# x_pot =  
# 
# x_pot = seq(-528000, 528000, by=8000)
# y_pot = seq(-416000, 416000, by=8000)
# coord_pot = expand.grid(x_pot, y_pot)
# 
# d_pot = t(rdist(matrix(c(0,0), ncol=2), as.matrix(coord_pot, ncol=2))/dist.scale)
# d_pot = unname(as.matrix(count(d_pot)))
# 
# N_pot = nrow(d_pot)

coord_pot = seq(-700000, 700000, by=8000)
coord_pot = expand.grid(coord_pot, coord_pot)

d_pot = t(rdist(matrix(c(0,0), ncol=2), as.matrix(coord_pot, ncol=2))/rescale)
d_pot = unname(as.matrix(count(data.frame(d_pot))))

N_pot     = nrow(d_pot)
sum_w_pot = build_sumw_pot(post, K, N_pot, d_pot, run)

#####################################################################################
# recompute gamma
#####################################################################################
w_coarse  = build_sumw_pot(post, K, length(d_hood), cbind(t(d_hood), rep(1, length(d_hood))), run)

gamma_new = recompute_gamma(w_coarse, sum_w_pot, gamma)

# #####################################################################################
# # domain splitting check
# #####################################################################################
# w_all <- build_weight_matrix(post, d, idx_cores_all, N, length(idx_cores_all), run)
# 
# foo=apply(w_all, 1, rowSums)
# 
# pollen_check$sum_w = foo

#####################################################################################
# veg run pars
#####################################################################################
names_substr = substr(names(mean_pars),1,3)

eta = unname(mean_pars[which(names_substr == 'eta')])[1:W]
rho = unname(mean_pars[which(names_substr == 'rho')])[1:W]

# ##########################################################################################################################
# ## chunk: qr decompose X
# ##########################################################################################################################
# 
# x = matrix(1, nrow=(N*T), ncol=1)
# N_p = N*T
# 
# temp = qr(x)
# Q = qr.Q(temp)
# R = qr.R(temp)
# 
# P = Q %*% t(Q)
# # M = diag(N_p) - P
# 
# if (all(P-P[1,1]<1.0e-12)){
#   P = P[1,1]
#   N_p = 1
# }

##########################################################################################################################
## save the data; rdata more efficient, use for processing
##########################################################################################################################
if (kernel == 'gaussian'){ suff = paste0('G_', suff) } else if (kernel == 'pl'){suff = paste0('PL_', suff)}
if (KGAMMA) suff = paste0('KGAMMA_', suff)
if (KW) suff = paste0('KW_', suff)

# note that w is column-major 
save(K, N, T, N_cores, N_knots, res,
     gamma, psi, phi, rho, eta,
     y, 
     idx_cores, 
     d, d_knots, d_inter, w,
     #lag,
     #        P, N_p, sum_w_pot,
     sum_w_pot,
     knot_coords,
     centers_pls, centers_veg, centers_pol, taxa, ages, y_veg, N_pls,
     file=paste('r/dump/', K, 'taxa_', N, 'cells_', N_knots, 'knots_', 'cal_', suff, '.rdata',sep=""))

# convert to row-major
if (KW){
  w_new = vector(length=0)
  for (k in 1:K)
    w_new = c(w_new, as.vector(w[k,,]))
  w = array(w_new, c(K, N_cores, N))  
}

dump(c('K', 'N', 'T', 'N_cores', 'N_knots', 'res',
       'gamma', 'psi', 'phi', 'rho', 'eta',
       'y', 
       'idx_cores', 
       'd', 'd_knots', 'd_inter', 'w',
       #'lag',
       #        'P', 'N_p', 'sum_w_pot'),
       'sum_w_pot'),#, 'pollen_check'),
     #        'knot_coords',
     #        'centers_pls', 'centers_veg', 'centers_polU', 'taxa', 'ages', 'y_veg', 'N_pls'), 
     file=paste('r/dump/', K, 'taxa_', N, 'cells_', N_knots, 'knots_', 'cal_', suff, '.dump',sep=""))

# ##########################################################################################################################
# ## if base model copy gamma and w
# ##########################################################################################################################
# if ( (!KGAMMA) & (!KW) ) { 
#   suff = paste0('COPY_', suff)
#   
#   gamma = rep(gamma, K)
#   sum_w_pot = rep(sum_w_pot, K)
#   
# #   w_new = vector(length=K*N_cores*N)
# #   for (j in 1:N)
# #     for (i in 1:N_cores)
# #       for (k in 1:K)
# #         w_new[(k-1)*N*N_cores + (i-1)*N + j] = w[i, j]
# # #   
# #   w_new = array(0, c(K, N_cores, N))
# #   for (k in 1:K){
# #     w_new[k,,] = w
# #   }
#   
#   w = array(rep(as.vector(w), K), c(K, N_cores, N))
#   
#   dump(c('K', 'N', 'T', 'N_cores', 'N_knots', 'res',
#          'gamma', 'psi', 'phi', 'rho', 'eta',
#          'y', 
#          'idx_cores', 
#          'd', 'd_knots', 'd_inter', 'w',
#          'lag',
#          #        'P', 'N_p', 'sum_w_pot'),
#          'sum_w_pot'),#, 'pollen_check'),
#        #        'knot_coords',
#        #        'centers_pls', 'centers_veg', 'centers_polU', 'taxa', 'ages', 'y_veg', 'N_pls'), 
#        file=paste('r/dump/', K, 'taxa_', N, 'cells_', N_knots, 'knots_', tmin, 'to', tmax, 'ypb_', suff, '.dump',sep=""))
#   
#   save(K, N, T, N_cores, N_knots, res,
#        gamma, psi, phi, rho, eta,
#        y, 
#        idx_cores, 
#        d, d_knots, d_inter, w,
#        lag,
#        #        P, N_p, sum_w_pot,
#        sum_w_pot, #pollen_check,
#        knot_coords,
#        centers_pls, centers_veg, centers_pol, taxa, ages, y_veg, N_pls,
#        file=paste('r/dump/', K, 'taxa_', N, 'cells_', N_knots, 'knots_', tmin, 'to', tmax, 'ypb_', suff, '.rdata',sep=""))
#   
# }

# dump(c('K', 'N', 'T', 'N_cores', 'N_knots', 'res',
#        'gamma', 'psi', 'phi', 'rho', 'eta',
#        'y', 
#        'idx_cores', 
#        'd', 'd_knots', 'd_inter', 'w',
#        'lag',
# #        'P', 'N_p', 'sum_w_pot'),
#         'sum_w_pot', 'pollen_check'),
# #        'knot_coords',
# #        'centers_pls', 'centers_veg', 'centers_polU', 'taxa', 'ages', 'y_veg', 'N_pls'), 
#      file=paste('r/dump/pred_data_', K, 'taxa_', N, 'cells_', N_knots, 'knots_', tmin, 'to', tmax, 'ypb_', suff, '.dump',sep=""))
# 
# save(K, N, T, N_cores, N_knots, res,
#        gamma, psi, phi, rho, eta,
#        y, 
#        idx_cores, 
#        d, d_knots, d_inter, w,
#        lag,
# #        P, N_p, sum_w_pot,
#        sum_w_pot, pollen_check,
#        knot_coords,
#        centers_pls, centers_veg, centers_pol, taxa, ages, y_veg, N_pls,
#        file=paste('r/dump/pred_data_', K, 'taxa_', N, 'cells_', N_knots, 'knots_', tmin, 'to', tmax, 'ypb_', suff, '.rdata',sep=""))

# ##########################################################################################################################
# ## build initial conditions
# ##########################################################################################################################
# 
# tau   = 150
# mu    = rep(0,W) 
# omega = 0.5
# ksi   = 0.1
# 
# mu_t = array(0, dim=c(W, T))
# for (k in 1:W){
#   mu_t[k,1] = rnorm(1, mean= mu[k], sd= sqrt((ksi * ksi)/(1 - omega * omega)))
#   for (i in 2:T){
#     mu_t[k,i] = rnorm(1, mean = mu[k] + omega * mu_t[k, i-1], sd = ksi)     
#   }
# }
# 
# # alpha = build_alpha_init(W, N_knots, T, rep(rho,W), tau, rep(eta,W), d_knots, lag)
# inits = pred_build_inits(K, N, N_knots, eta, rho, mu, tau, d_knots, d_inter, lag)
# alpha = inits$alpha_init
# g     = inits$g_init
# 
# # eta   = 1#rep(1, W) 
# 
# dump(c('tau', 'mu', 'omega', 'ksi', 'mu_t', 'alpha', 'g'), 
#      file=paste('r/dump/pred_data_', K, 'taxa_', N, 'cells_', N_knots, 'knots_', tmin, 'to', tmax, 
#                 'ypb_', suff, '_inits.dump',sep=""))
# 
# ##########################################################################################################################
# ## build initial conditions for full
# ##########################################################################################################################
# 
# #tau   = 150
# mu    = rep(0,W)
# ksi   = 0.1
# 
# sigma  = rep(0.6, W)
# lambda = rep(0.2, W)
# 
# mu_t = array(0, dim=c(W, T))
# for (k in 1:W){
#   mu_t[k,1] = rnorm(1, mean= 0, sd= 20)
#   print(rnorm(1, mean= 0, sd= 20))
#   for (i in 2:T){
#     mu_t[k,i] = rnorm(1, mean = mu_t[k, i-1], sd = ksi)     
#   }
# }
# 
# inits = pred_build_inits_full(K, N, N_knots, eta, rho, mu, mu_t, tau, d_knots, d_inter, lag)
# alpha_s = inits$alpha_s_init
# alpha_t = inits$alpha_t_init
# g       = inits$g_init
# 
# # eta   = 1#rep(1, W) 
# 
# dump(c('tau', 'ksi', 'sigma', 'lambda', 'mu', 'mu_t', 'alpha_s', 'alpha_t', 'g'), 
#      file=paste('r/dump/pred_data_', K, 'taxa_', N, 'cells_', N_knots, 'knots_', tmin, 'to', tmax, 
#                 'ypb_', suff, '_inits_full.dump',sep=""))