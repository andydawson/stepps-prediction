library(sp)
library(rgdal)
library(fields)
# library(gpclib)
library(maptools)

source('r/utils/pred_helper_funs.r')

# us.shp <- readShapeLines('data/map_data/us_alb.shp',
#                          proj4string=CRS('+init=epsg:3175'))
us.shp <- readOGR('data/map_data/us_alb.shp', 'us_alb')
is.projected(us.shp)
# busted!
# proj4string(us.shp)=CRS('+init=epsg:3175')

# gridname = 'wi_5by'
# gridname = 'mn_3by'

rescale    = 1e6
nside      = 1
cell_width = 8000/rescale

states_pls = c('wisconsin', 'michigan:north', 'minnesota')
states_pol = c('wisconsin', 'michigan:north', 'minnesota')

pls = data.frame(read.table(file='data/western_comp_stepps_v0.3-1.csv', sep=",", row.names=NULL, header=TRUE))

# get pls centers and counts
meta           = pls[,1:4]
centers_pls    = data.frame(x=pls$X, y=pls$Y)
colnames(meta) = tolower(colnames(meta))
meta = split_mi(meta) # this is a total hack...

plot(meta$x, meta$y) 

# get rid of mi:southern peninsula
meta = meta[meta$state2 %in% states_pls,]
meta$x = meta$x/rescale
meta$y = meta$y/rescale

plot(meta$x, meta$y) # but look, the southern peninsula is gone

centers_pls = cbind(meta$x, meta$y)
N           = nrow(centers_pls)

# distance matrix
d = rdist(as.matrix(centers_pls), as.matrix(centers_pls))
d[which(d<1e-8)]=0
# 
# gridname = paste0('umw)
# domain = domain_coarse
# save(domain, N, d, d_hood, centers_pls, file=paste('data/grid/', gridname, 'by.rdata',sep=''))

domain_fine = centers_pls
xlo = min(domain_fine[,1])
xhi = max(domain_fine[,1])
ylo = min(domain_fine[,2])
yhi = max(domain_fine[,2])

# pad out the domain to be divisible by nside
get_padding <- function(xlo, xhi, cell_width, nside){
  ncells     = (xhi-xlo)/cell_width
  ncells_new = round(ncells/nside)*nside

  pad = ncells_new-ncells

  return(pad)
}

x_pad = get_padding(xlo, xhi, cell_width, nside)
y_pad = get_padding(ylo, yhi, cell_width, nside)

# build the coarse rectangle that covers domain
x_coarse    = seq(xlo, xhi + x_pad*cell_width, by=cell_width*nside)
y_coarse    = seq(ylo, yhi + y_pad*cell_width, by=cell_width*nside)
grid_coarse = expand.grid(x_coarse, y_coarse)

par(mfrow=c(1,1))
plot(grid_coarse[,1], grid_coarse[,2])

# find the coarse cells in the domain
d_inter = rdist(as.matrix(domain_fine), as.matrix(grid_coarse))
d_inter[which(d_inter<1e-8)]=0

x = seq(cell_width/2,cell_width*nside,by=cell_width)
y = x
coarse_cell = expand.grid(x,y)
d_hood = rdist(matrix(c(cell_width*nside/2,cell_width*nside/2), ncol=2), coarse_cell)
d_hood[which(d_hood<1e-8)]=0

idx_coarse = rep(0, nrow(domain_fine))

for (i in 1:nrow(d_inter)){
  
  min_col = which.min(d_inter[i,])
#   print(d_inter[i,min_col])
  
  if (d_inter[i,min_col] %in% d_hood){
    idx_coarse[i] = which.min(d_inter[i,])
  }
}

domain_coarse = data.frame(grid_coarse[unique(idx_coarse), ])*rescale
colnames(domain_coarse) = c('x', 'y')

coordinates(domain_coarse) <- ~ x + y
proj4string(domain_coarse) <- CRS(proj4string(us.shp))

# get the coarse cell states
# some centroids in other states, put them back where we want them
coarse_states = as.vector(over(domain_coarse, us.shp)$STATE_NAME)
MN = c('South Dakota', 'North Dakota', 'Iowa')
coarse_states[coarse_states %in% MN] = 'Minnesota'

domain_coarse       = as.data.frame(domain_coarse)/rescale
domain_coarse$state = coarse_states

pdf('figures/coarse_domain_test.pdf')
par(mfrow=c(1,1))
plot(domain_coarse[,1]*rescale, domain_coarse[,2]*rescale)
plot(us.shp, add=T, lwd=2)
dev.off()

gridname = paste0('umw_', nside)
# write.table(domain_coarse,file=paste('data/', gridname, '_coarse.csv', sep=''),sep=",",row.names=F)
domain = domain_coarse
save(domain, N, d, d_hood, centers_pls, file=paste('data/grid/', gridname, 'by.rdata',sep=''))

############################################################################################################################
# split the domain into two pieces
############################################################################################################################

xlo_wi = min(domain_coarse[domain_coarse$state == 'Wisconsin','x'])
yhi_wi = max(domain_coarse[domain_coarse$state == 'Wisconsin','y'])

xhi_mn = max(domain_coarse[domain_coarse$state == 'Minnesota','x'])
yhi_mn = max(domain_coarse[domain_coarse$state == 'Minnesota','y'])

domW = domain_coarse[domain_coarse$state %in% c('Minnesota'), ]

# add a buffer
y_domW = unique(domW$y)
cw     = nside*cell_width

domW = rbind(domW, domain_coarse[which(domain_coarse$x < 0.421),])
domW = rbind(domW, domain_coarse[which(domain_coarse$y > 1.214),])
#domW = rbind(domW, domain_coarse[which(domain_coarse$y > 1.198),])

plot(domain_coarse$x*rescale, domain_coarse$y*rescale)
plot(us.shp, add=TRUE)
points(domW$x*rescale, domW$y*rescale, col='red', pch=19)

domW = domW[!duplicated(domW[,1:2]),]

plot(domain_coarse$x*rescale, domain_coarse$y*rescale)
plot(us.shp, add=TRUE)
points(domW$x*rescale, domW$y*rescale, col='red', pch=19)

gridname = paste0('umwW_', nside)
domain=domW
save(domain, N, d, d_hood, centers_pls, file=paste('data/grid/', gridname, 'by.rdata',sep=''))

# now for the eastern half
domE = domain_coarse[domain_coarse$state %in% c('Wisconsin', 'Michigan'), ]

# add a buffer
y_domE = unique(domE$y)
cw   = nside*cell_width

edge = cell_width*21

# domE = domain_coarse[which(domain_coarse$x > 0.421 - cw*4.5),]
domE = domain_coarse[which(domain_coarse$x > 0.421 - edge),]

plot(domain_coarse$x*rescale, domain_coarse$y*rescale)
plot(us.shp, add=TRUE)
points(domE$x*rescale, domE$y*rescale, col='red', pch=19)

domE = domE[!duplicated(domE[,1:2]),]

plot(domain_coarse$x*rescale, domain_coarse$y*rescale)
plot(us.shp, add=TRUE)
points(domE$x*rescale, domE$y*rescale, col='red', pch=19)

gridname = paste0('umwE_', nside)
domain=domE
save(domain, N, d, d_hood, centers_pls, file=paste('data/grid/', gridname, 'by.rdata',sep=''))

# 
# #######################################################################################
# 
# fit   = read_stan_csv('calibration/output/pollen_fit_7taxa_smoothed_veg_v2.csv')
# gamma = unname(summary(fit)$summary[,'mean']['gamma'])
# psi   = unname(summary(fit)$summary[,'mean']['psi'])
# 
# w <- build_weight_matrix(d, idx_cores, N, N_cores, psi)
# 
# k = 10
# 
# # indices of fine cells in a bigger coarse cell
# idx_fine = which(d[idx_cores[k],] %in% d_hood)
# 
# weights = w[k,idx_fine]
# 
# C = sum(w[k,])
#
# gamma_coarse = gamma + sum(weights)/C
# 
# save(gamma_coarse, file='prediction/r/dump/gamma_coarse_5by.rdata')

