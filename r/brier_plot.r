library(ggplot2)
library(fields)

source('r/utils/pred_plot_funs.r')

# edit this file to process different runs
source('r/runs_cal.r')
load('data/comp_data_12taxa_mid_ALL_v0.3.rdata')
colnames(centers_comp) = c('x', 'y')

rescale = 1e6

# brier = list()
# for (run in runs){
#   suff_fit  = run$suff_fit
#   suff_figs = run$suff_figs
#   
#   # where to put the figures
#   subDir <- paste("figures/", suff_fit, sep='')
#   
#   load(paste0('figures/', suff_fit, '/brier_scores_', suff_figs, '.rdata'))
#   brier = rbind(brier, cbind(dat, type=rep(suff_figs, nrow(dat))))
# }
# 
# brier_sums = aggregate(resids~type, sum, data=brier)
# 
# limits = get_limits(brier[,2:3])
# 
# levels(brier$type) <- c('Base PL', 'Base G', 'Variable PL', 'Variable G')
# brier$type <- factor(brier$type, levels=c('Base G', 'Base PL', 'Variable G','Variable PL'))
# # levels(brier$type) <- factor(brier$type, levels=levels(brier$type)[c(2,1,4,3)])
# 
# p <- ggplot() + geom_tile(data=brier, aes(x=x, y=y, fill=resids)) + 
#   scale_fill_gradientn(colours=tim.colors(), name="Brier score") + coord_fixed() 
# p <- add_map_albers(p, map_data=us.fort, limits)
# p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.3)), strip.text.x = element_text(size = rel(1.0)))
# p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# p <- p + facet_grid(type~.)
# print(p)
# 
# fname = 'figures/maps_brier.pdf'
# ggsave(p, file=fname, scale=3)
# ggsave(p, file='figures/maps_brier.png', scale=3)
# sys_str = paste("pdfcrop", fname, fname, sep=' ')
# system(sys_str)

#####################################################################################################
# real Brier with pls data
#####################################################################################################

pls = read.csv('../stepps-data/data/pls_umw_v0.5.csv')

# compute the distance between comp and preds for each spatial location
# analogous to brier
compute_distance <- function(comp, preds){
  
  ntaxa  = ncol(comp)
  ncells = nrow(comp)
  
  score = rep(0, ncells)
  for (n in 1:ncells){
#     score[n] = score[n] + sum((preds[n,] - comp[n,])^2)#sqrt(sum((preds[n,] - comp[n,])^2))
    score[n] = score[n] + sqrt(sum((preds[n,] - comp[n,])^2))
  }
  
  return(score)
}

brier_fake = list()
for (run in runs){
  suff_fit  = run$suff_fit
  suff_figs = run$suff_figs
  
  load(paste0('figures/', suff_fit, '/r_mean.rdata'))
  colnames(r_mean) = taxa
  
  # check ordering
  dist_mat = rdist(as.matrix(centers_r_mean)*rescale, as.matrix(centers_comp))
  dist_mat[dist_mat < 1e-6] = 0
  idx = apply(dist_mat, 1, function(x) if (any(x==0)){ which(x == 0) } else {0})
  centers_comp = centers_comp[idx,]
  r_comp       = r_comp[idx,]
  
  brier_score = compute_distance(r_comp, r_mean)
  
  brier_fake = rbind(brier_fake, data.frame(brier_score, centers_comp, type=rep(suff_figs, nrow(centers_comp))))
}

limits = get_limits(brier_fake[,2:3])

levels(brier_fake$type) <- c('Base PL', 'Base G', 'Variable PL', 'Variable G')
brier_fake$type <- factor(brier_fake$type, levels=c('Base G', 'Base PL', 'Variable G','Variable PL'))


brier_sums = aggregate(brier_score~type, sum, data=brier_fake)
brier_sums

p <- ggplot() + geom_tile(data=brier_fake, aes(x=x, y=y, fill=brier_score)) + 
  scale_fill_gradientn(colours=tim.colors(), name="Distance") + coord_fixed() 
p <- add_map_albers(p, map_data=us.fort, limits)
p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.4)), 
                            strip.text.x = element_text(size = rel(1.4)),
                            legend.title=element_text(size = rel(1.1)),
                            legend.text=element_text(size = rel(1.0)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# p <- p + facet_grid(type~.)
p <- p + facet_wrap(~type)
print(p)

fname = 'figures/maps_l2.pdf'
ggsave(p, file=fname, width=14)#scale=1.5)
ggsave(p, file='figures/maps_l2.png', width=14)#, scale=1.5)
sys_str = paste("pdfcrop", fname, fname, sep=' ')
system(sys_str)

#####################################################################################################
# real Brier with pls data
#####################################################################################################

pls = read.csv('../stepps-data/data/pls_umw_v0.5.csv')

compute_brier <- function(counts, preds){
  
  ntaxa  = ncol(counts)
  ncells = nrow(counts)
  
  score = rep(0, ncells)
  for (n in 1:ncells){
    for (i in 1:ntaxa){
      if (counts[n,i] != 0){
        score[n] = score[n] + ( (preds[n,i] - 1)^2  + (preds[n,i] - 0)^2 * (ntaxa-1) ) * counts[n,i]
      }
    }
    
    if (sum(counts[n,]) != 0){
      score[n] = score[n]#/sum(counts[n,])
    } else {
      score[n] = NA
    }
    
  }
  
  return(score)
}

brier_pls = list()
for (run in runs){
  suff_fit  = run$suff_fit
  suff_figs = run$suff_figs

  load(paste0('figures/', suff_fit, '/r_mean.rdata'))
  colnames(r_mean) = taxa
  
  # check ordering
  dist_mat = rdist(as.matrix(centers_r_mean)*rescale, as.matrix(pls[,1:2]))
  dist_mat[dist_mat < 1e-6] = 0
  idx = apply(dist_mat, 1, function(x) if (any(x==0)){ which(x == 0) } else {0})
  pls = pls[idx,]
  

  brier_score = compute_brier(pls[,3:ncol(pls)], r_mean)
  brier_score = brier_score/sum(pls[,3:ncol(pls)])
  
  brier_pls = rbind(brier_pls, cbind(brier_score, pls[,1:2], type=rep(suff_figs, nrow(pls))))
}

limits = get_limits(brier_pls[,2:3])

levels(brier_pls$type) <- c('Base PL', 'Base G', 'Variable PL', 'Variable G')
brier_pls$type <- factor(brier_pls$type, levels=c('Base G', 'Base PL', 'Variable G','Variable PL'))

brier_sums = aggregate(brier_score~type, sum, data=brier_pls)
brier_sums

p <- ggplot() + geom_tile(data=brier_pls, aes(x=x, y=y, fill=brier_score)) + 
  scale_fill_gradientn(colours=tim.colors(), name="Brier score") + coord_fixed() 
p <- add_map_albers(p, map_data=us.fort, limits)
p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.3)), strip.text.x = element_text(size = rel(1.3)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# p <- p + facet_grid(type~.)
p <- p + facet_wrap(~type)
print(p)

fname = 'figures/maps_brier_pls.pdf'
ggsave(p, file=fname, scale=3)
ggsave(p, file='figures/maps_brier_pls.png', scale=3)
sys_str = paste("pdfcrop", fname, fname, sep=' ')
system(sys_str)

# plot the pls tree counts

p <- ggplot() + geom_tile(data=pls, aes(x=x, y=y, fill=rowSums(pls[,3:ncol(pls)]))) + 
  scale_fill_gradientn(colours=tim.colors(), name="Tree count") + coord_fixed() 
p <- add_map_albers(p, map_data=us.fort, limits)
p <- theme_clean(p) + theme(strip.text.y = element_text(size = rel(1.3)), strip.text.x = element_text(size = rel(1.0)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)

fname = 'figures/map_pls_counts.pdf'
ggsave(p, file=fname, scale=3)
