### SMC sampler code in R
rm(list=ls())
set.seed(1)

# source the 
source("unbiased_SMC.R")

# experiment settings
N <- 1e2          # number of particles
n <- 50           # number of steps in SMC
sigma2_prop <- 5  # variance of proposal 
mu0 <- 5          # mean of initial normal distribution
sd0 <- 10         # standard deviation of initial normal distribution
n_ests <- 1e3     # number of unbiased estimators


# define the log-target (default is a mixture)
mu1 <- -0.5
mu2 <- 3
log_gamma <- function(p,x){
  log((0.5*dnorm(x,mu1,1)+0.5*dnorm(x,mu2,3))^((p-1)/(n-1))) + log(dnorm(x,mu0,sd0/2)^(1-(p-1)/(n-1)))
}

# plot a few of the targets
x_plt <- seq(from=-5,to=10,length.out=100)
plot(x_plt,exp(log_gamma(n,x_plt)),type='l',ylab='density')
matplot(x_plt,exp(log_gamma(1,x_plt)),add=T,type='l')
matplot(x_plt,exp(log_gamma(n/2,x_plt)),add=T,type='l')
title('sequence of targets')

# run standard SMC
smc_out <- smc(mu, M, G, n, N, mu0, sd0, sigma2_prop)
sum(smc_out$zetas[n,]*smc_out$w_normalised[n,])

# unbiased version
h <- function(x){ x }
H0_ests <- rep(NA,n_ests)
tau_arr <- rep(NA,n_ests)
for(j in 1:n_ests){
  unbiased_res <- unbiased_RB_est(h, N)
  H0_ests[j] <- unbiased_res$H0
  tau_arr[j] <- unbiased_res$tau
}

# get confidence intervals
mean_est <- mean(H0_ests)
CI_u <- mean_est  + sqrt(2*var(H0_ests)/n_ests)
CI_l <- mean_est - sqrt(2*var(H0_ests)/n_ests)

# plot the distribution of estimators
hist(H0_ests,1e2,main='unbiased estimators')
abline(v=((mu1+mu2)/2)) # plot true value
abline(v=CI_u,col='blue') # plot true value
abline(v=CI_l,col='blue') # plot true value
abline(v=mean_est,col='darkgreen')
legend('topright',legend=c('CI','unbiased est','true'),lty=c(1,1,1),col=c('blue','darkgreen','black'))
print(sprintf('test results for SMC with confidence interval at +-3 sigma'))
print(sprintf('true value: %.3f     estimate value: %.3f    upper CI: %.3f    lower CI: %.3f',
              (mu1+mu2)/2,mean_est,CI_u,CI_l))


# plot the distribution of meeting times
hist(tau_arr,1e2,main='histogram of meeting times')




