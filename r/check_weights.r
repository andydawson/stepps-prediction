library(ggplot2)

source('r/utils/pred_plot_funs.r')

load('r/dump/pred_data_12taxa_262cells_120knots_0to2000ypb_umw_5by_v0.3.rdata')
w_5by_umw = rowSums(w)

load('r/dump/pred_data_12taxa_459cells_78knots_0to2000ypb_umwW_3by_v0.3.rdata')
w_3by_umwW  = rowSums(w)
pc_3by_umwW = pollen_check
centers_vegW = centers_veg

load('r/dump/pred_data_12taxa_387cells_66knots_0to2000ypb_umwE_3by_v0.3.rdata')
w_3by_umwE  = rowSums(w)
pc_3by_umwE = pollen_check
centers_vegE = centers_veg

load('r/dump/pred_data_12taxa_699cells_120knots_0to2000ypb_umw_3by_v0.3.rdata')
w_3by_umw  = rowSums(w)
pc_3by_umw = pollen_check

load('r/dump/pred_data_12taxa_4199cells_79knots_0to2000ypb_umwW_1by_v0.3.rdata')
w_1by_umwW = rowSums(w)

load('r/dump/pred_data_12taxa_3608cells_69knots_0to2000ypb_umwE_1by_v0.3.rdata')
w_1by_umwE = rowSums(w)

load('r/dump/pred_data_12taxa_6341cells_120knots_0to2000ypb_umw_1by_v0.3.rdata')
w_1by_umw = rowSums(w)

summary(w_1by_umw / sum_w_pot)

summary(w_3by_umw*9 / sum_w_pot)

summary(w_5by_umw*25 / sum_w_pot)

# check halfing for 3by
limits = get_limits(centers_pls)

breaks = c(0,0.5, 0.6, 0.7, 0.8, 0.9, 1)
breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                    function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

w_all = pc_3by_umw[which((pc_3by_umw$int == TRUE) & (pc_3by_umwE$int == TRUE)), 'sum_w']
w_e   = pc_3by_umwE[which((pc_3by_umw$int == TRUE) & (pc_3by_umwE$int == TRUE)), 'sum_w']
metaE = data.frame(pc_3by_umwE[which((pc_3by_umw$int == TRUE) & (pc_3by_umwE$int == TRUE)), c('x', 'y')], prop_w=w_e/w_all)
metaE$prop_w_binned = cut(metaE$prop_w, breaks, include.lowest=TRUE, labels=FALSE)

w_all = pc_3by_umw[which((pc_3by_umw$int == TRUE) & (pc_3by_umwW$int == TRUE)), 'sum_w']
w_w   = pc_3by_umwW[which((pc_3by_umw$int == TRUE) & (pc_3by_umwW$int == TRUE)), 'sum_w']
metaW = data.frame(pc_3by_umwW[which((pc_3by_umw$int == TRUE) & (pc_3by_umwW$int == TRUE)), c('x', 'y')], prop_w=w_w/w_all)

metaALL = data.frame(pc_3by_umw[, c('x', 'y', 'sum_w')])

rescale=1e6
rast=data.frame(centers_vegE*rescale, z=rep(1,nrow(centers_vegE)))

p <- ggplot() + geom_tile(data=rast, aes(x=x, y=y, fill=factor(1))) + scale_fill_manual(values='#99CCFF')
p <- p + geom_point(data=metaE, aes(x=x*rescale, y=y*rescale, size = prop_w), colour='#003333') 
p <- add_map_albers(p, map_data=us.fort, limits=limits)
print(p)
ggsave(file='figures/weights_3by_umwE.pdf')

rast=data.frame(centers_vegW*rescale, z=rep(1,nrow(centers_vegW)))
p <- ggplot() + geom_raster(data=rast, aes(x=x, y=y, fill=factor(1))) + scale_fill_manual(values='#99CCFF')
p <- p + geom_point(data=metaW, aes(x=x*rescale, y=y*rescale, size = prop_w), colour='#003333') 
p <- add_map_albers(p, map_data=us.fort, limits=limits)
print(p)
ggsave(file='figures/weights_3by_umwW.pdf')

# rast=data.frame(centers_veg*rescale, z=rep(1,nrow(centers_veg)))
# p <- ggplot() + geom_raster(data=rast, aes(x=x, y=y, fill=factor(1))) + scale_fill_manual(values='#99CCFF')
# p <- p + geom_point(data=metaALL, aes(x=x*rescale, y=y*rescale, size = sum_w), colour='#003333') 
# p <- add_map_albers(p, map_data=us.fort, limits=limits)
# print(p)
# ggsave(file='figures/weights_3by_umwW.pdf')
