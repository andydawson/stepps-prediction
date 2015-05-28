library(Rcpp)

sourceCpp('r/utils/build_mu_g_omp.cpp')

build_mu_g <- function(post_dat, rho, eta, T, K, d, d_inter, d_knots, od, mpp, mu0) {

  N = nrow(d_inter)

  post      = post_dat$post
  par_names = post_dat$par_names

  alpha_s_start = min(which(par_names == 'alpha_s'))
  alpha_t_start = min(which(par_names == 'alpha_t'))

  if (od) {
    ones = matrix(1, nrow=N, ncol=1)
    temp = qr(ones)
    Q = qr.Q(temp)
    P = Q %*% t(Q)
  }

  mu      = post[,1,which(par_names == 'mu')]
  mu_t    = post[,1,which(par_names == 'mu_t')]
  sigma   = post[,1,which(par_names == 'sigma')]
  lambda  = post[,1,which(par_names == 'lambda')]
  alpha_s = post[,1,which(par_names == 'alpha_s')]
  alpha_t = post[,1,which(par_names == 'alpha_t')]

  build_mu_g_omp(rho, sigma, lambda, mu, mu_t, alpha_s, alpha_t, d_knots, d_inter, T, K, od, mu0, P)
}
