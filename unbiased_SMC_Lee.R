### SMC sampler code in R
rm(list=ls())
set.seed(1)


# intial conditions
mu <- function(N) rnorm(N, mean=mu0, sd=sd0)
log_mu <- function(x) dnorm(x, mean=mu0, sd=sd0,log=T)
# unnormalised targets
log_gamma <- function(p, x) {
  -0.5*( (x-means[p])/sds[p])^2
}
log_G <- function(p, x) {
  if (p == 1) {
    res <- log_gamma(1,x) - log_mu(x)
  } else if (p <= n) {
    res <- log_gamma(p,x) - log_gamma(p-1,x)
  }
  return(res)
}
M <- function(p, x, sigma2_prop) {
  N <- length(x)
  y <- x + sqrt(sigma2_prop)*rnorm(N)
  log_ratios <- log_gamma(p, y)-log_gamma(p, x)
  accepts <- log(runif(N)) < log_ratios
  y*accepts + x*(1-accepts)
}
smc <- function(mu, M, G, n, N, mu0, sd0, sigma2_prop) {
  zetas <- matrix(0,n,N)
  w_normalised <- matrix(0,n,N)
  as <- matrix(0,n-1,N)
  log_gs <- matrix(0,n,N)
  log_Zs <- rep(0,n)
  
  zetas[1,] <- mu(N)
  log_gs[1,] <- log_G(1,zetas[1,])
  log_Zs[1] <-  log(mean(exp(log_gs[1,]-max(log_gs[1,]))))+max(log_gs[1,])
  w_normalised[1,] <- exp(log_gs[1,]-max(log_gs[1,]))/sum(exp(log_gs[1,]-max(log_gs[1,])))
  
  if(n>1){
    for (p in 2:n) {
      # simulate ancestor indices, then particles
      as[p-1,] <- sample(N,N,replace=TRUE,prob=w_normalised[p-1,])
      zetas[p,] <- M(p, zetas[p-1,as[p-1,]], sigma2_prop)
      log_gs[p,] <- log_G(p, zetas[p,])
      log_Zs[p] <- log_Zs[p-1] + log(mean(exp(log_gs[p,]-max(log_gs[p,]))))+max(log_gs[p,])
      w_normalised[p,] <- exp(log_gs[p,]-max(log_gs[p,]))/sum(exp(log_gs[p,]-max(log_gs[p,])))
    }
  }
  return(list(zetas=zetas,as=as,log_Z=log_Zs[n],w_normalised=w_normalised))
}

# biased pimh kernel
pimh_kernel <- function(smc_out){
  # obtain proposal
  smc_prop <- smc(mu, M, G, n, N, mu0, sd0, sigma2_prop)
  
  # conditionally accept proposal
  if(log(runif(1))<(smc_prop$log_Z - smc_out$log_Z)){
    return_res <- smc_prop
  }else{
    return_res <- smc_out
  }
  return(return_res)
}

# unbiased version
coupled_pimh_kernel <- function(smc_out1,smc_out2){
  # obtain proposal
  smc_prop <- smc(mu, M, G, n, N, mu0, sd0, sigma2_prop)
  u_rnd <- runif(1)
  
  # conditionally accept proposal for first chain
  if(log(u_rnd)<(smc_prop$log_Z - smc_out1$log_Z)){
    return_res1 <- smc_prop
  }else{
    return_res1 <- smc_out1
  }
  
  # conditionally accept proposal for second chain
  if(log(u_rnd)<(smc_prop$log_Z - smc_out2$log_Z)){
    return_res2 <- smc_prop
  }else{
    return_res2 <- smc_out2
  }
  
  return_res <- list(return_res1=return_res1,return_res2=return_res2)
  return(return_res)
}

# unbiased algorithm
unbiased_RB_est <- function(h, N){
  
  # set up arrays
  R <- 1e3
  logZ_arr1 <- rep(NA,R+1)
  zeta_arr1 <- matrix(NA,R+1,N)
  w_normalised_arr1 <- matrix(NA,R+1,N)
  logZ_arr2 <- rep(NA,R)
  zeta_arr2 <- matrix(NA,R,N)
  w_normalised_arr2 <- matrix(NA,R,N)
  
  # initialise
  tau <- Inf
  i <- 1
  smc_out1 <- list(log_Z=-Inf)
  smc_out2 <- list(log_Z=-Inf)
  
  # run algorithm
  smc_out1 <- pimh_kernel(smc_out1)
  logZ_arr1[1] <- smc_out1$log_Z
  zeta_arr1[1,] <- smc_out1$zetas[n,]
  w_normalised_arr1[1,] <- smc_out1$w_normalised[n,]
  
  while(is.infinite(tau)){
      # iterate for both chains
      smc_out_both <- coupled_pimh_kernel(smc_out1,smc_out2)
      smc_out1 <- smc_out_both$return_res1
      smc_out2 <- smc_out_both$return_res2
      
      # save results for two chains
      logZ_arr1[i+1] <- smc_out1$log_Z
      zeta_arr1[i+1,] <- smc_out1$zetas[n,]
      w_normalised_arr1[i+1,] <- smc_out1$w_normalised[n,]
      logZ_arr2[i] <- smc_out2$log_Z
      zeta_arr2[i,] <- smc_out2$zetas[n,]
      w_normalised_arr2[i,] <- smc_out2$w_normalised[n,]

      if(smc_out1$log_Z==smc_out2$log_Z){ tau <- i }
      i <- i+1
  }
  
  H0 <- sum(  h(zeta_arr1[1,])  * w_normalised_arr1[1,]   )
  for(i in 1:tau){
    if(i < tau){
      delta_h1 <- sum( h(zeta_arr1[i+1,])  * w_normalised_arr1[i+1,]  )
      delta_h2 <- sum( h(  zeta_arr2[i,])  * w_normalised_arr2[i,]    )
      delta_h <- delta_h1-delta_h2
    }else{
      delta_h <- 0
    }
    H0 <- H0 + delta_h
  }
  return(list(H0=H0,tau=tau))
}
