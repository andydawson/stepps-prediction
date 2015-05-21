library(fields)
# library(rstan)
library(sp)
library(rgdal)
library(maptools)

source('r/utils/pred_helper_funs.r')
source('eda/r/utils/dataPlotFuns.r')

# us albers shape file
# us.shp <- readShapeLines('attic/r/data/map_data/us_alb.shp',
#                          proj4string=CRS('+init=epsg:3175'))
us.shp <- readOGR('data/map_data/us_alb.shp', 'us_alb')

# date of pollen_ts pull
mydate = '2014-07-22'


# reconstruction limits and bin-width
tmin = 0
tmax = 2000
int  = 100

# rescale
rescale = 1e6

# suff=''
version="v0.3"
suff    = paste(gridname, '_', version, sep='') 
# suff = '3by_v0.3_test'

states_pol = c('minnesota', 'wisconsin', 'michigan:north', 'michigan:south')
states_pls = c('minnesota', 'wisconsin', 'michigan:north', 'michigan:south')

taxa_all = toupper(c('oak', 'pine', 'maple', 'birch', 'tamarack', 'beech', 'elm', 'spruce', 'ash', 'hemlock'))
taxa_sub = toupper(c('oak', 'pine', 'maple', 'birch', 'tamarack', 'beech', 'elm', 'spruce', 'ash', 'hemlock'))

K = as.integer(length(taxa_sub) + 1)
W = K-1

# conversion tables
tree_type = read.table('data/assign_HW_CON.csv', sep=',', row.names=1, header=TRUE)
convert   = read.table('data/dict-comp2stepps.csv', sep=',', row.names=1, header=TRUE)

# PLS
pls.raw   = data.frame(read.table(file='data/western_comp_stepps_v0.3-1.csv', sep=",", row.names=NULL, header=TRUE))

# pls   = read.table(file='data/pls_UMW.csv', sep=",", row.names=NULL, header=TRUE)
# pls   = pls[-which(is.na(pls[,'oak'])),]
# pls   = pls[pls$state %in% states_pls,]
# #pls   = pls[1000:1500,]
# pls$x = pls$x
# pls$y = pls$y

colnames(pls.raw) = tolower(colnames(pls.raw))

# pull the subset of proportions
taxa.start.col = min(match(tolower(rownames(convert)), colnames(pls.raw)), na.rm=TRUE)

pls_dat  = pls.raw[,taxa.start.col:ncol(pls.raw)]
colnames(pls_dat) = as.vector(convert[match(colnames(pls_dat), tolower(rownames(convert))),1])
pls_dat_collapse  = sapply(unique(colnames(pls_dat)), 
                           function(x) rowSums( pls_dat[ , grep(x, names(pls_dat)), drop=FALSE]) )
counts = data.frame(pls_dat_collapse[,sort(colnames(pls_dat_collapse))])
meta   = pls.raw[,1:(taxa.start.col-1)]
meta   = split_mi(meta)
counts = counts[which(meta$state2 %in% states_pls),]
meta   = meta[which(meta$state2 %in% states_pls),]

centers_pls = data.frame(x=meta$x, y=meta$y)/rescale # megameters!
plot(centers_pls[,1]*rescale, centers_pls[,2]*rescale, asp=1, axes=F,  col='antiquewhite4', xlab='',ylab='', pch=19, cex=0.2)
plot(us.shp, add=T)

y_veg = convert_counts(counts, tree_type, taxa_sub)
N_pls = nrow(pls.raw)



# ##
# ## organize taxa
# ##
# 
# #taxa_sub = c('oak', 'pine', 'birch', 'maple')
# taxa_all =c("MAPLE", "ALDER", "BIRCH", "BEECH", "HICKORY", "OTHER.HARDWOOD", "CEDAR.JUNIPER", "ASH",
#             "WALNUT", "TAMARACK", "OTHER.HARDWOOD", "IRONWOOD", "SPRUCE", "PINE", "SYCAMORE", 
#             "POPLAR.TULIP.POPLAR", "OAK", "WILLOW", "BASSWOOD", "HEMLOCK", "ELM", "BLACK.GUM.SWEET.GUM", 
#             "RUMEX", "ROSE", "AMBROSIA", "CORNUS")
# # taxa_sub =c("MAPLE", "ALDER", "BIRCH", "HICKORY", "OTHER.HARDWOOD", "CEDAR.JUNIPER", "ASH",
# #             "WALNUT", "TAMARACK", "OTHER.HARDWOOD", "IRONWOOD", "SPRUCE", "PINE", "SYCAMORE", 
# #             "POPLAR.TULIP.POPLAR", "OAK", "WILLOW", "BASSWOOD", "HEMLOCK", "ELM", "BLACK.GUM.SWEET.GUM") 
# taxa_sub =c("MAPLE", "ALDER", "BIRCH", "BEECH", "HICKORY", "OTHER.HARDWOOD", "CEDAR.JUNIPER", "ASH",
#             "WALNUT", "TAMARACK", "IRONWOOD", "SPRUCE", "PINE", "SYCAMORE", 
#             "POPLAR.TULIP.POPLAR", "OAK", "WILLOW", "BASSWOOD", "HEMLOCK", "ELM", "BLACK.GUM.SWEET.GUM", 
#            "CORNUS")
# # taxa_sub = c('oak', 'pine', 'birch', 'maple', 'spruce', 'hemlock')
# taxa = c(taxa_sub, 'other')
# #taxa=taxa_all
# 
# taxa_all = c("ASH", "BEECH", "BIRCH", "ELM", "HEMLOCK", "MAPLE", "OAK", "OTHER.CONIFER", 
#              "OTHER.HARDWOOD", "PINE", "SPRUCE", "TAMARACK")
# taxa_sub =taxa_all
# taxa = taxa_all
# 
# # K = as.integer(length(taxa_sub) + 1)
# # W = K-1
# K = length(taxa)


# POLLEN

pollen_ts = read.table(paste('data/pollen_ts_', mydate, '.csv', sep=''), header=TRUE, stringsAsFactors=FALSE)
# last_date = "2015-03-20"
# pollen_ts = read.table(paste0("data/pollen/pollen_ts_", last_date, ".csv"), header=TRUE)
pollen_ts = pollen_ts[which((pollen_ts$ages <= 2500) & (pollen_ts$state %in% states_pol)),]

pollen_ts <- pollen_to_albers(pollen_ts) # reproject from lat long to Albers

# # pollen_ts = pollen_ts[which((pollen_ts[,'x'] <= xhi) & (pollen_ts[,'x'] >= xlo) & 
# #                               (pollen_ts[,'y'] <= yhi) & (pollen_ts[,'y'] >= ylo)),]
# 
# idx = which(colnames(pollen_ts) == 'oak')
# 
# ## plot the subdomain and core locations 
# par(mfrow=c(1,1))
# plot(centers_pls[,1], centers_pls[,2])
# points(pollen_ts$x*1000, pollen_ts$y*1000, col='blue', pch=19)
# us.shp <- readShapeLines('data/map_data/us_alb.shp',
#                          proj4string=CRS('+init=epsg:3175'))
# plot(us.shp, add=T, lwd=2)
# 
# plot(pollen_ts$x*1000, pollen_ts$y*1000, col='blue', pch=19)
# points(centers_pls[,1], centers_pls[,2])
# 
# us.shp <- readShapeLines('data/map_data/us_alb.shp',
#                          proj4string=CRS('+init=epsg:3175'))
# plot(us.shp, add=T, lwd=2)

##########################################################################################################################
## chunk 2: prepare pollen data; aggregate over time intervals
##########################################################################################################################

# sum counts over int length intervals
pollen_agg = build_pollen_counts(tmin=tmin, tmax=tmax, int=int, pollen_ts=pollen_ts, taxa_all=taxa_all, taxa_sub=taxa_sub)
meta_pol   = pollen_agg[[2]]
counts     = pollen_agg[[1]]

ages    = unique(sort(meta_pol$ages))
T       = length(ages) 
lag     = unname(as.matrix(dist(matrix(ages), upper=TRUE)))
N_cores = length(unique(meta_pol$id))

y       = counts[, colnames(counts) %in% taxa_sub]
y$other = rowSums(counts[, !colnames(counts) %in% taxa_sub])

# make sure columns match!
y = y[,taxa]

y       = unname(round(as.matrix(y)))


centers_polU = matrix(NA, nrow=N_cores, ncol=2)
for (i in 1:N_cores){
  id = unique(meta_pol$id)[i]
  idx = min(which(meta_pol$id == id))
  print(idx)
  centers_polU[i,] = c(meta_pol$x[idx], meta_pol$y[idx])  
}

# indices for which cells the cores fall in
# idx_cores <- build_idx_cores(centers_polU, centers_veg, N_cores)

par(mfrow=c(1,1))
plot(centers_polU[,1]*1000000, centers_polU[,2]*1000000, col='blue', pch=4, cex=1.4)
us.shp <- readShapeLines('data/map_data/us_alb.shp',
                         proj4string=CRS('+init=epsg:3175'))
plot(us.shp, add=T, lwd=2)

# save(K, N, T, N_cores, N_knots,
#      y, 
#      idx_cores, 
#      centers_pls, centers_veg, centers_polU, taxa, ages, y_veg, N_pls,
#      file=paste('r/pred/dump/pollen_pie_data.rdata',sep=""))



build_core_dat_new <- function(y, ages, centers_polU, N_cores, taxa){
  core_dat = data.frame(x=integer(0), y=integer(0), age=character(), core_sites=integer(0))
  K = length(taxa)
  
  for (i in 1:length(ages)){
    print(i)
    print((N_cores*(i-1) + 1))
    print(N_cores*i)
    
    age_rows = seq(i, length(ages)*N_cores, by=T)
    y_sub = y[age_rows,]
    
#     y_sub    = y1[(N_cores*(i-1) + 1):(N_cores*i),]
    idx_data = which(rowSums(y_sub) != 0) 
    
    print(idx_data)
    
    for (k in 1:K){
      
      core_dat = rbind(core_dat, data.frame(x     = centers_polU[idx_data,1]*1000000, 
                                            y     = centers_polU[idx_data,2]*1000000, 
                                            age   = rep(ages[i], length(idx_data)),
                                            taxa  = rep(taxa[k], length(idx_data)),
                                            counts = y_sub[idx_data,k])) 
    }
    
    
  }
  
  return(core_dat)  
}



core_dat_new<-build_core_dat_new(y, ages, centers_polU, N_cores, taxa)

# recast the data
library(reshape2)
pollen_df = dcast(core_dat_new, x+y+age~taxa, sum)
# pollen_df = pollen_df[which(pollen_df$age %in% ages_sub ),]
# make the pie plots


# veg proportions
mkprop <- function(x) {
  s = sum(x, na.rm=TRUE)
  if (s > 0) {
    y = x/s
  } else {
    y = rep(0, length(x))
  }
  y
}

pollen_df[,4:ncol(pollen_df)] = t(apply(as.matrix(pollen_df[,4:ncol(pollen_df)]), 1, mkprop))

years  <- sort(unique(pollen_df$age))
nyears <- length(years)

# 
# buffer = 0#100000
# xlo = min(centers_pls[,1]) - buffer
# xhi = max(centers_pls[,1]) + buffer
# ylo = min(centers_pls[,2]) - buffer
# yhi = max(centers_pls[,2]) + buffer

buffer = 50000#100000
xlo = min(pollen_df$x) - buffer
xhi = max(pollen_df$x) + buffer
ylo = min(pollen_df$y) - buffer
yhi = max(pollen_df$y) + buffer


#pdf('figures/pred_model/pie_maps_nov6_mn.pdf', width=6, height=14)
par(mfrow=c(1,1))
# par(mfrow=c(4,2))

postscript(paste('prediction/figures/pollen_pies/pollen_pie_plot_UMW_', 
                 paste(years[1]*100 - int/2, 'to', max(years)*100 + int/2, 'ybp_', int, 'int', sep=''), 
                 '.eps', sep=''), width=10, height=10)

add_legend = rep(FALSE, nyears)
add_legend[1] = TRUE

for (i in 1:nyears){
  
  #   time_rows = seq(i,N_cores*T, by=T)
  main_title = paste(years[i]*100 - int/2, 'to', years[i]*100 + int/2, 'ybp', sep=' ')
  
  idx_years = which(pollen_df$age == years[i])
  
  pollen_df_sub = pollen_df[idx_years,]
  
  props = pollen_df_sub[,4:ncol(pollen_df_sub)]
  centers = cbind(pollen_df_sub$x, pollen_df_sub$y)
  
#   postscript(paste('figures/pred/pollen_pies/pollen_pie_plot_UMW_', paste(years[i]*100 - int/2, 'to', years[i]*100 + int/2, 'ybp_', int, 'int', sep=''), '.eps', sep=''), width=10, height=10)
  pieMap(proportions = props, 
         centers  = centers,
         restrict = FALSE,
         inputRestricted = FALSE,
         xlim   = c(xlo, xhi),
         ylim   = c(ylo, yhi),
         radius = 22000,
         scale  = 1,
         xlab   = 'x',
         ylab   = 'y',
         add_legend = add_legend[i],
         main_title = main_title)
#   dev.off()
}
dev.off()