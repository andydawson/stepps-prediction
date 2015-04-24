library(sp)
library(rgdal)
library(fields)
library(rstan)
library(maptools)

source('r/utils/pred_helper_funs.r')

us.shp <- readShapeLines('data/map_data/us_alb.shp',
                         proj4string=CRS('+init=epsg:3175'))

gridname = 'umw_3by'
load(file=paste('data/grid/', gridname, '.rdata',sep=''))

states_pol = c('wisconsin', 'michigan:north', 'minnesota')
# states_pol = c('wisconsin','michigan:north')
#states_pol = c('minnesota')

# read in pollen data
# pollen = read.csv('data/pollen_meta_2014-07-22.csv', header=TRUE)
# date of pollen_ts pull
mydate = '2014-07-22'
pollen = read.table(paste('data/pollen_ts_', mydate, '.csv', sep=''), header=TRUE)
fit    = read_stan_csv('data/12taxa_mid_comp_long.csv')


pollen = pollen[pollen$state %in% states_pol,]

# reproject pollen coords from lat long to Albers
centers_pol = data.frame(cbind(pollen$long, pollen$lat))
colnames(centers_pol) = c('x', 'y')
coordinates(centers_pol) <- ~ x + y
proj4string(centers_pol) <- CRS('+proj=longlat +ellps=WGS84')

centers_polA <- spTransform(centers_pol, CRS('+init=epsg:3175'))
centers_polA <- as.matrix(data.frame(centers_polA))#/1000

centers_polA <- unique(centers_polA)

plot(centers_polA[,1], centers_polA[,2])
plot(us.shp, add=T, lwd=2)

N_cores = nrow(centers_polA)

centers_polA =centers_polA/1000000

# find idx of the cell that the core is in
idx_cores = vector(length=N_cores)

for (i in 1:nrow(centers_polA)){
  core_site = centers_polA[i,]
  d_core = rdist(matrix(core_site, ncol=2), as.matrix(centers_pls))
  idx_cores[i] = which.min(d_core)
}

# visual checks
par(mfrow=c(1,1))
plot(centers_pls[idx_cores,1], centers_pls[idx_cores,2])
points(centers_polA[,1], centers_polA[,2], col='blue', pch=8)

gamma = unname(summary(fit)$summary[,'mean']['gamma'])
psi   = unname(summary(fit)$summary[,'mean']['psi'])

w <- build_weight_matrix(d, idx_cores, N, N_cores, psi)

k = 21

# indices of fine cells in a bigger coarse cell
idx_fine = which(d[idx_cores[k],] %in% d_hood)

weights = w[k,idx_fine]

C = sum(w[k,])

gamma_coarse = gamma + sum(weights)/C

save(gamma_coarse, file='prediction/gamma_coarse_5by.rdata')



