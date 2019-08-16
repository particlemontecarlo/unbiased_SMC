### SMC sampler code in R
rm(list=ls())
set.seed(1)

# source the 
source("unbiased_SMC_Lee.R")

## test 1 - importance sampling
n <- 1
N <- 1e4
means <- 0.5
sds <- 0.5
sigma2_prop <- 2
R <- 100
mu0 <- 0
sd0 <- 1
log_Z_collect <- rep(NA,R)
for(r in 1:R){
  log_Z_collect[r] <- smc(mu, M, G, n, N, mu0, sd0, sigma2_prop)$log_Z
}
hist(exp(log_Z_collect))
abline(v=sqrt(2*pi*sds[1]^2))


## test 2 - sequential monte carlo sampling
n <- 50
N <- 1e2
means <- seq(5,1,length.out=n)
sds <- seq(5,1,length.out=n)
sigma2_prop <- 5
mu0 <- 5
sd0 <- 10
# R <- 100
smc_means <- rep(NA,R)
smc_out <- smc(mu, M, G, n, N, mu0, sd0, sigma2_prop)
rowSums(smc_out$zetas*smc_out$w_normalised)
rowSums((smc_out$zetas-means)^2*smc_out$w_normalised)^0.5


# test pimh kernel - see the bias
R <- 1e5
n_pimh <- 5
N <- 1e1
logZ_arr <- matrix(NA,R,n_pimh)
rb_est <- matrix(NA,R,n_pimh)
smc_out <- list(log_Z=-Inf)
for(i in 1:R){
  for(j in 1:n_pimh){
    smc_out <- pimh_kernel(smc_out)
    logZ_arr[i,j] <- smc_out$log_Z
    rb_est[i,j] <- sum(smc_out$zetas[n,]*smc_out$w_normalised[n,])
  }
}

plot(logZ_arr)
plt_mat <- matrix(c(rep(means[n],n_pimh),colMeans(rb_est)),nrow=2)
rb_plot <- colMeans(rb_est)
rb_avg_std <- sqrt(diag(var(rb_est))/R)
rb_u <- rb_plot+3*rb_avg_std
rb_l <- rb_plot-3*rb_avg_std
rb_ylims <- c(0.99*min(rb_plot,means[n]),1.01*max(rb_plot,means[n]))
plot(rb_plot,ylim=rb_ylims,type='l')
matplot((matrix(c(rb_u,rb_l),ncol=2)),type='l',add=T)
abline(h=means[n],col='blue')
legend('topright',legend=c('estimated','true'),lty=c(1,1),col=c('black','blue'))


# unbiased version
h <- function(x){ x }

n_ests <- 1e5#1*n_pimh*R # get estimators on a comparable budget
H0_ests <- rep(NA,n_ests)
tau_arr <- rep(NA,n_ests)
for(j in 1:n_ests){
  unbiased_res <- unbiased_RB_est(h, N)
  H0_ests[j] <- unbiased_res$H0
  tau_arr[j] <- unbiased_res$tau
}

CI_u <- mean(H0_ests)  + sqrt(2*var(H0_ests)/n_ests)
CI_l <- mean(H0_ests)  - sqrt(2*var(H0_ests)/n_ests)
abline(h=CI_u,col='darkgreen')
abline(h=CI_l,col='darkgreen')




hist(H0_ests,1e2)
hist(tau_arr,1e2)




