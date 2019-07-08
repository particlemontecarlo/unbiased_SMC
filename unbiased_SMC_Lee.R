### SMC sampler code in R
rm(list=ls())



smc <- function(mu, M, G, n, N) {
  zetas <- matrix(0,n,N)
  as <- matrix(0,n-1,N)
  gs <- matrix(0,n,N)
  log.Zs <- rep(0,n)
  
  zetas[1,] <- mu(N)
  gs[1,] <- G(1, zetas[1,])
  log.Zs[1] <- log(mean(gs[1,]))
  
  for (p in 2:n) {
    # simulate ancestor indices, then particles
    as[p-1,] <- sample(N,N,replace=TRUE,prob=gs[p-1,]/sum(gs[p-1,]))
    zetas[p,] <- M(p, zetas[p-1,as[p-1,]])
    gs[p,] <- G(p, zetas[p,])
    log.Zs[p] <- log.Zs[p-1] + log(mean(gs[p,]))
  }
  return(list(zetas=zetas,gs=gs,as=as,log.Zs=log.Zs))
}


## example
n <- 10
N <- 1e3
means = seq(4,0,length=n)
sds <- n:1

mu <- function(N) rnorm(N, mean=means[1], sd=sds[1])
gamma <- function(p, x) {
  dnorm(x, mean=means[p], sd=sds[p])
}
G <- function(p, x) {
  if (p < n) {
    res <- gamma(p+1,x)/gamma(p,x) 
  }else {
    res <- rep(1, length(x))
  }
  return(res)
}
M <- function(p, x) {
  N <- length(x)
  y <- x + rnorm(N)
  ratios <- gamma(p, y)/gamma(p, x)
  accepts <- runif(N) < ratios
  y*accepts + x*(1-accepts)
}
out <- smc(mu, M, G, n, N)
apply(out$zetas, 1, mean)
apply(out$zetas, 1, sd)
out$log.Zs



## test normaising constant estimate
n <- 5
N <- 1e3
means = seq(4,0,length=n)
sds <- n:1
gamma <- function(p, x) {
  exp(-0.5*  ((x-means[p])/sds[p])^2  )
}

R <- 100
z_ests <- rep(NA,R)
for(r in 1:R){
  out <- smc(mu, M, G, n, N)
  apply(out$zetas, 1, mean)
  apply(out$zetas, 1, sd)
  z_ests[r] <- exp(out$log.Zs)[n]
}

hist(z_ests)
abline(v=sqrt(2*pi*sds[n]^2)/sqrt(2*pi*sds[1]^2))



## ok now to have 2D
n <- 100
N <- 1e3
means = seq(1,0,length=n)
sds <- seq(from=2,to=1,length.out=n)
mu <- function(N) matrix(rnorm(D*N, mean=means[1], sd=sds[1]),nrow=2)
log_gamma <- function(p, x) {
  return(dnorm(x[1,], mean=means[p], sd=sds[p],log=T) + dnorm(x[2,], mean=means[p], sd=sds[p],log=T))
}
log_G <- function(p, x) {
  if (p < n) {
    res <- log_gamma(p+1,x)-log_gamma(p,x) 
  }else {
    res <- rep(0, length(x))
  }
  return(res)
}
M <- function(p, x) {
  N <- dim(x)[2]
  y <- x + matrix(rnorm(N*D),nrow=D)
  log_ratios <- log_gamma(p, y) - log_gamma(p, x)
  accepts <- log(runif(N)) < log_ratios
  res <- y[,accepts] + x[,1-accepts]
  return(res)
}
smc <- function(mu, M, G, D, n, N) {
  as <- matrix(0,n-1,N)
  gs <- matrix(0,n,N)
  log.Zs <- rep(0,n)
  
  zetas <- mu(N)
  log_gs_p <- log_G(1, zetas)
  log.Z.increment <- log(mean(exp(log_gs_p-max(log_gs_p)))) + max(log_gs_p)
  log.Zs[1] <- log.Z.increment
  
  for (p in 2:n) {
    # simulate ancestor indices, then particles
    zetas_pm1 <- zetas
    log_gs_pm1 <- exp(log_gs_p-max(log_gs_p))
    as[p-1,] <- sample(N,N,replace=TRUE,prob=log_gs_pm1/sum(log_gs_pm1))
    zetas <- M(p, zetas_pm1[,as[p-1,]])
    log_gs_p <- log_G(p, zetas)
    log.Z.increment <- log(mean(exp(log_gs_p-max(log_gs_p)))) + max(log_gs_p)
    log.Zs[p] <- log.Zs[p-1] + log.Z.increment
  }
  return(list(zetas=zetas,as=as,log.Zs=log.Zs))
}

R <- 100
z_ests <- rep(NA,R)
D <- 2
for(r in 1:R){
  out <- smc(mu, M, G, D,n, N)
  apply(out$zetas, 1, mean)
  apply(out$zetas, 1, sd)
  z_ests[r] <- exp(out$log.Zs)[n]
}

hist(z_ests)
abline(v=(2*pi*sds[n]^2)/(2*pi*sds[1]^2))

