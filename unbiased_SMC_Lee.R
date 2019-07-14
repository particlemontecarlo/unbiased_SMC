### SMC sampler code in R
rm(list=ls())
set.seed(1)


# intial conditions
mu <- function(N) rnorm(N, mean=0, sd=1)
log_mu <- function(x) dnorm(x, mean=0, sd=1,log=T)
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
smc <- function(mu, M, G, n, N, sigma2_prop) {
  zetas <- matrix(0,n,N)
  as <- matrix(0,n-1,N)
  log_gs <- matrix(0,n,N)
  log_Zs <- rep(0,n)
  
  zetas[1,] <- mu(N)
  log_gs[1,] <- log_G(1,zetas[1,])
  log_Zs[1] <-  log(mean(exp(log_gs[1,]-max(log_gs[1,]))))+max(log_gs[1,])
  
  if(n>1){
    for (p in 2:n) {
      # simulate ancestor indices, then particles
      w_normalised <- exp(log_gs[p-1,]-max(log_gs[p-1,]))/sum(exp(log_gs[p-1,]-max(log_gs[p-1,])))
      as[p-1,] <- sample(N,N,replace=TRUE,prob=w_normalised)
      zetas[p,] <- M(p, zetas[p-1,as[p-1,]], sigma2_prop)
      log_gs[p,] <- log_G(p, zetas[p,])
      log_Zs[p] <- log_Zs[p-1] + log(mean(exp(log_gs[p,]-max(log_gs[p,]))))+max(log_gs[p,])
    }
  }
  return(list(zetas=zetas[n,],as=as,log_Z=log_Zs[n]))
}


## example
n <- 1
N <- 1e3
means <- 0.5#seq(4,0,length=n)
sds <- 0.1#n:1
sigma2_prop <- 2
R <- 100
log_Z_collect <- rep(NA,R)
for(r in 1:R){
  log_Z_collect[r] <- smc(mu, M, G, n, N, sigma2_prop)$log_Z
}
hist(exp(log_Z_collect))

# apply(out$zetas, 1, mean)
# apply(out$zetas, 1, sd)



# biased pimh kernel
pimh_kernel <- function(smc_out){
  # obtain proposal
  smc_prop <- smc(mu, M, G, n, N, sigma2_prop)
  
  # conditionally accept proposal
  if(log(runif(1))<(smc_prop$log_Z - smc_out$log_Z)){
    return_res <- smc_prop
  }else{
    return_res <- smc_out
  }
  return(return_res)
}
smc_out<-smc(mu, M, G, n, N, sigma2_prop)

# 
R <- 1000
logZ_arr <- rep(NA,R)
zeta_arr <- matrix(NA,R,N)
smc_out <- list(log_Z=-Inf)
for(i in 1:R){
  smc_out <- pimh_kernel(smc_out)
  logZ_arr[i] <- smc_out$log_Z
  zeta_arr[i,] <- smc_out$zetas
}
plot(logZ_arr)

# unbiased version
coupled_pimh_kernel <- function(smc_out1,smc_out2){
  # obtain proposal
  smc_prop <- smc(mu, M, G, n, N, sigma2_prop)
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


# run unbiased algorithm
R <- 1000
logZ_arr1 <- rep(NA,R)
zeta_arr1 <- matrix(NA,R,N)
logZ_arr2 <- rep(NA,R)
zeta_arr2 <- matrix(NA,R,N)
smc_out1 <- list(log_Z=-Inf)
smc_out1 <- pimh_kernel(smc_out1)
for(i in 1:R){
  # iterate for both chains
  smc_out_both <- coupled_pimh_kernel(smc_out1,smc_out2)
  smc_out1 <- smc_out_both$return_res1
  smc_out2 <- smc_out_both$return_res2
  
  # save results for two chains
  logZ_arr1[i] <- smc_out1$log_Z
  zeta_arr1[i,] <- smc_out1$zetas
  logZ_arr2[i] <- smc_out2$log_Z
  zeta_arr2[i,] <- smc_out3$zetas
}





