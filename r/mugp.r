# library(doMC)
# registerDoMC(4)

    # c++: vector<vector_d> alpha_t(W*(T-1));
    #
    # alpha_t: [W, T-1, N_knots]
    # row major: k, t, j: k*(T-1)*N_knots + (t-1)*N_knots + v;
    #
    # note that c++ code writes out as:
    #
    # for (int v = 0; v < N_knots; ++v) {
    #   for (int kt = 0; kt < (W * (T - 1)); ++kt) {
    #     vars__.push_back(alpha_t[kt][v]);
    #   }
    # }
    #
    # so in R: t, k, j: j*(T-1)*W + k*(T-1) + t-1
    #
    # this results in [T-1,W,N_knots]
    #
    # so the following is correct if t starts at 2:
    # alpha_t_cols = seq(alpha_t_start + (k-1)*(T-1) + t-1 - 1, alpha_t_start + N_knots*W*(T-1) - 1, by=W*(T-1))
    # alpha_t_pre[,,t-1] = post[1:niter,1,alpha_t_cols]
    # alpha_t = alpha_t_pre[,,t-1]
    # alpha_t_t = t(alpha_t)
    # Halpha_t[,(k-1)*(T-1) + t-1,i] = q_Qinv %*% alpha_t_t[,i]
    #
    # want to use this in, eg,
    # for k (taxa):
    #   for i (iters):
    #     for t (time):
    #       Halpha_t[,(k-1)*(T-1) + t-1,i] = q_Qinv %*% alpha_t_t[,i]
    #
    # so want to get something like:
    #   alpha_t[j,i,t,k]: post[i,1,start + j*(T-1)*W + k*(T-1) + t-1]
    #
    # let's try:
    #   alpha_t_t[,j,t,k]: post[,1,start + j*(T-1)*W + k*(T-1) + t-1]
    #




build_mu_g_k <- function(k, mu_k, mu_t_k, rho, sigma, lambda, alpha_s, alpha_t, d_knots, d_inter, P, T, W, N, N_knots, od, mu0, niter) {

  C_s <- exp(-d_knots/rho[k])
  c_s <- exp(-d_inter/rho[k])
  C_s_inv = chol2inv(chol(C_s))

  if (od) {
    cs_Csinv = c_s %*% C_s_inv - P %*% c_s %*% C_s_inv
  } else {
    cs_Csinv = c_s %*% C_s_inv
  }

  Halpha_s = array(NA, dim=c(N, niter))
  mu_g     = array(NA, dim=c(N*T, niter))

  if (mu0){
    Halpha_t = array(NA, dim=c(N, T-1, niter))
  } else {
    Halpha_t = array(NA, dim=c(N, T, niter))
  }


  for (i in 1:niter) {

    if ( (i %% 200) == 0 ) { print(paste0("Iteration ", i))}

    Halpha_s[,i] = cs_Csinv %*% alpha_s[i,]

    if (mu0) {
      mu_g[,i] = mu_k[i] + Halpha_s[,i]
    } else {
      mu_g[,i] = mu_k[i] + mu_t_k[i,1] + Halpha_s[,i]
    }

    Q <- exp(-d_knots/lambda[i])
    q <- exp(-d_inter/lambda[i])
    Q_inv = chol2inv(chol(Q))

    q_Qinv = q %*% Q_inv

    for (t in 2:T){

      if (mu0) {
        Halpha_t[,t-1,i] = q_Qinv %*% alpha_t[,i,t-1]
        mu_g[,i] = mu_k[i] + mu_t_k[i,t] + cs_Csinv %*% alpha_s[i,] + Halpha_t[,t-1,i]
      } else {
        mu_g[,i] = mu_k[i] + mu_t_k[i,t-1] + cs_Csinv %*% alpha_s[i,] + q_Qinv %*% alpha_t[,i,t-1]
      }
    }
  }

  return(list())
}

build_mu_g_parallel <- function(post_dat, rho, eta, T, K, d, d_inter, d_knots, od, mpp, mu0) {

  post      = post_dat$post
  par_names = post_dat$par_names

  N       = nrow(d_inter)
  N_knots = ncol(d_inter)
  niter   = dim(post[,1,])[1]
  W       = K-1

  ones     = matrix(1, nrow=N, ncol=1)
  mu_g     = array(NA, dim=c(N*T, W, niter))
  Halpha_s = array(NA, dim=c(N, W, niter))

  if (mu0) {
    Halpha_t = array(NA, dim=c(N, W*(T-1), niter))
  } else {
    Halpha_t = array(NA, dim=c(N, W*T, niter))
  }

  alpha_t = array(NA, dim=c(N_knots, niter, T-1))

  print("Done allocating")

  ksi    = post[,1, which(par_names == 'ksi')]

  alpha_s_start = min(which(par_names == 'alpha_s'))
  alpha_t_start = min(which(par_names == 'alpha_t'))

  if (od) {
    temp = qr(ones)
    Q = qr.Q(temp)
    P = Q %*% t(Q)
  }

  mu     = post[,1,which(par_names == 'mu')]
  mu_t   = post[,1,which(par_names == 'mu_t')]
  sigma  = post[,1,which(par_names == 'sigma')]
  lambda = post[,1,which(par_names == 'lambda')]

  res = foreach (k = 1:W) %dopar% {
  for (k in 1:W) {

    if (mu0) {
      mut_cols = seq(k, (T-1)*W, by=W)
    } else {
      mut_cols = seq(k, T*W, by=W)
    }

    alpha_s_cols = seq(alpha_s_start + k - 1, alpha_s_start + N_knots*W - 1, by=W)
    alpha_s = post[,1,alpha_s_cols]

    # XXX: AWKWARD!
    for (j in 1:N_knots) {
      for (t in 2:T) {
        alpha_t[j,,t-1] = post[,1, alpha_t_start + j*(T-1)*W + k*(T-1) + t-1 - 1]
      }
    }

    build_mu_g_k(k, mu[,k], mu_t[,mut_cols], rho, sigma[,k], lambda[,k], alpha_s, alpha_t, d_knots, d_inter, P, T, W, N, N_knots, od, mu0, niter)
  }
  # XXX: unpack res...

  return(list(mu_g=mu_g, mu=mu, mu_t=mu_t, Halpha_s=Halpha_s, Halpha_t=Halpha_t))
}

# ## lastly via OpenMP for parallel use
# build_mu_g_omp <- '
#    // assign to C++ vector
#    std::vector<double> x = Rcpp::as<std::vector< double > >(xs);
#    size_t n = x.size();
# #pragma omp parallel for shared(x, n)
#    for (size_t i=0; i<n; i++) {
#        x[i] = ::log(x[i]);
#    }
#    return Rcpp::wrap(x);
# '
# 
# settings <- getPlugin("Rcpp")
# settings$env$PKG_CXXFLAGS <- paste('-fopenmp', settings$env$PKG_CXXFLAGS)
# settings$env$PKG_LIBS <- paste('-fopenmp -lgomp', settings$env$PKG_LIBS)
# funOpenMP <- cxxfunction(signature(xs="numeric"), body=openMPCode, plugin="Rcpp", settings=settings)