### SMC sampler code in R
rm(list=ls())

### script to obtain unbiased estimates with respect to a target annealed using sequential monte carlo
# we have that the target pi(x) = exp(U(x)) where U(x)  is some target energy
# we then specify 
# i) the initial sampling distribution
# ii) the sequence of targets, in this case we will choose
#
#   pi_t(x) = pi(x)^alpha_t * pi0(x)^(1-alpha_t)
#
# so that we vary alpha between 0 and 1
# in particular we will set alpha_t = (t/(T-1))
# we will set pi0 = N(0,sigma02 * I)
# we will use D to denote the dimension of the state-space
# each step leaves pi_t invariant

D <- 1
mu0 <- 1
sigma02 <- 5
sigma2_prop <- 0.5
T_final <- 100
nparticles <- 1e3


# target
log_gamma <- function(xparticles){
  return(rowSums(-0.5*(xparticles^2)))
}

# initial distribution
r_pi_0 <- function(nparticles){
  return(matrix(rnorm(D*nparticles,mean=mu0,s=sqrt(sigma02)),ncol=D))
}
log_pi_0 <- function(xparticles){
  return(rowSums(dnorm(xparticles,mean=mu0,s=sqrt(sigma02),log=T)))
  #return(colSums(-0.5*(log(2*pi*sigma02) + (xparticles^2)/sigma02) ))
}

# sequence of targets
log_pi_t <- function(t,xparticles){
  alpha_t <- (t/(T_final-1))
  return(alpha_t*log_gamma(xparticles) + (1-alpha_t)*log_pi_0(xparticles))
}


# specify incremental weights
log_w_t <- function(t,xparticles_m1){
  return(log_pi_t(t,xparticles_m1)-log_pi_t(t-1,xparticles_m1))
}


#specify transition kernel
r_transition <- function(t,xparticles_m1){
  xparticles_prop <- xparticles_m1 + matrix(rnorm(D*nparticles,0,s=sqrt(sigma2_prop)),ncol=D)
  logd_prop <- log_pi_t(t,xparticles_prop)
  logd_init <- log_pi_t(t,xparticles_m1)
  log_acc_ratio <- pmin(0,logd_prop-logd_init)

  nparticles <- dim(xparticles_m1)[1]
  log_u_rands <- log(runif(nparticles))
  accept_mask <- log_acc_ratio>log_u_rands
  
  xparticles <- xparticles_m1
  xparticles[accept_mask,] <- xparticles_prop[accept_mask,]
  return(xparticles)
}

# run particle filter
smc_sampler <- function(){
  # initialise arrays
  logZ_arr <- rep(NA,T_final)

  # sample initial particles
  xparticles <- r_pi_0(nparticles)
  
  logW_normalised <- -rep(log(nparticles),nparticles)
  logw1 <- log_pi_0(xparticles) + logW_normalised 
  
  max_lw <- max(logw1)
  logZ_hat <- log(sum(exp(logw1-max_lw))) + max_lw
  logZ_arr[1] <- logZ_hat
  logW_normalised <- logw1 - logZ_hat 

  
  for(t in 2:T_final){

    # resample
    anc <- sample.int(nparticles, size = nparticles, replace = T, prob = exp(logW_normalised) ) 
    xparticles <- as.matrix(xparticles[anc,])
    xparticles_m1 <- xparticles
    
    logW_normalised <- -rep(log(nparticles),nparticles)
    
    # propagate particles
    xparticles <- r_transition(t-1,xparticles_m1)
    
    # reweight
    logw <- log_w_t(t,xparticles_m1)+logW_normalised
    
    # log-sum-exp to renormalise
    max_lw <- max(logw)
    logZ_hat <- log(sum(exp(logw-max_lw))) + max_lw
    logZ_arr[t] <- logZ_hat
    logW_normalised <- logw - logZ_hat
  }
  logZ <- sum(logZ_arr)
  return_res <- list(logZ=logZ,xparticles=xparticles)
  return(return_res)
}


# plot normalising constant estimate
R <- 1e3
pf_res <- list()
for(i in 1:R){
  print(i)
  pf_res[[i]] <- smc_sampler()
}
logZ_collect <- sapply(pf_res,function(x) x$logZ)
xparticle_collect <- sapply(pf_res,function(x) x$xparticles[1,])


# perform tests
# ok for unit variance target, the normalising constant should approximate 2*pi
sqrt(2*pi)
mean(exp(logZ_collect))
var(exp(logZ_collect)/R)^0.5


# run unbiased SMC
# will need to define a coupled and single PIMH kernel




