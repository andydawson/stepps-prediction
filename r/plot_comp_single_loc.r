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

# edit this file to process different runs
source('r/runs.r')

nruns    = length(runs)
r_cell   = list(length=nruns)
r_quants = array(NA, dim=c(nruns, T, K, 3))

for (i in 1:nruns){
  run=runs[[i]]  

  # suff_dat  = run$suff_dat
  # suff_fit  = run$suff_fit
  # suff_figs = run$suff_figs

  g_cell = load_cell(run, 100)

  # # where to put the figures
  # subDir <- paste("figures/", suff_fit, sep='')
  # create_figure_path(subDir)
  
  # load the data and posterior draws
  load(paste0('r/dump/', run$suff_dat, '.rdata'))

  r_cell[[i]] = build_r_cell(g_cell, N, T, K) 

  for (t in 1:T){
    r_quants[i,t,,] = t(apply(r_cell[[i]]$rc[t,,], 1, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))) 
  }
 
}


r_melt = melt(r_quants, varnames=c("run", "time", "taxon", "quant"))

p <- ggplot(r_melt) + geom_line(data=r_melt, aes(x=value, y=time, colour=run, linetype=quant))
p <- p + facet_wrap(~taxon, ncol=6)


# functions
build_r_cell <- function(g_cell, N, T, K){
  
  W       = K-1
  niter   = nrow(g_cell) 
  
  gc      = array(NA, dim=c(T, W, niter))
  rc      = array(NA, dim=c(T, K, niter))
  
  for (k in 1:W){
    print(k)
    
    g_k_cols = seq(k, T*W, by=W)
    gc[,k,]   = t(g_cell[,g_k_cols])    
  }
  
  for (i in 1:niter){
    
    if ( (i %% 200) == 0 ) { print(paste0("Iteration ", i))}
    
    sum_exp_g = rowSums(exp(gc[,,i]))
       
    rc[,,i] <- sum2one_constraint_cell(K, N, T, as.matrix(gc[,,i]), sum_exp_g) 
  }
  
  return(list(rc=rc, gc=gc))
}



# r[,,i] <- sum2one_constraint_cell(K, N, T, as.matrix(g[,,i]), sum_exp_g) 
# additive log_ratio transformation
cppFunction('
  NumericMatrix sum2one_constraint_cell(int K, int N, int T, NumericMatrix g, NumericVector sum_exp_g) {
    //std::cout << "K " << K << "; N " << N << "; T " << T << std::endl; 
    NumericMatrix r(T, K);
    for (int k = 0; k<(K-1); k++)
      for (int j = 0; j<T; j++)
        r(j,k) = exp(g(j,k))/(1+sum_exp_g(j));

    for (int j = 0; j<T; j++)
      r(j,K-1) = 1.0 / (1 + sum_exp_g[j]);
   
return r;
  }
', verbose=TRUE)

load_cell <- function(run, id) {
  fname = sprintf('output/%s.bin', run$suff_fit)
  bin      = file(fname, "rb")
  nwarmup  = readBin(bin, "integer")
  nsamples = readBin(bin, "integer")
  nparams  = readBin(bin, "integer")

  # seek to start of parameter names
  seek(bin, where=(nsamples+nwarmup)*nparams*4, origin="current")
  params    = readBin(bin, "character", n=nparams)
  par_names = readBin(bin, "character", n=nparams)

  # find first parameter
  pcol0 = 1
  for (i in 1:10) {
    if (grepl("__", params[i], fixed=TRUE) == FALSE) {
      pcol0 = i
      break
    }
  }

  params = params[pcol0:length(params)]
  par_names = par_names[pcol0:length(par_names)]

  # get T, W
  load(paste0('r/dump/', run$suff_dat, '.rdata'))
  
  # figure out which columns we want (cols)
  g_cols    = which(par_names == 'g')
  col0      = g_cols[(id-1)*T*W + 1]
  ncols     = T*W

  # save start of samples
  start = 3*4 + nwarmup*nparams*4

  # read appropriate columns
  out = array(NA, dim=c(nsamples, ncols))
  for (row in 1:nsamples) {
    seek(bin, where=start+nparams*(row-1)*4+(pcol0+col0-2)*4, origin="start")
    out[row,] = readBin(bin, "numeric", n=ncols, size=4)
  }
  colnames(out) <- params[col0:(col0+ncols-1)]
  close(bin)

  out
}