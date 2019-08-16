### SMC sampler code in R
rm(list=ls())
set.seed(1)

# source the 
source("unbiased_SMC.R")


# define unnormalised normal target
N <- 1e5
sigma2_prop <- 5
mu0 <- 5
sd0 <- 10
n <- 100
means <- seq(5,1,length.out=n)
sds <- seq(5,1,length.out=n)
mu1 <- -0.5
mu2 <- 3
# log_gamma <- function(p,x){
#   log(dnorm(x,mu1,1)^((p-1)/(n-1))) + log(dnorm(x,mu0,sd0/2)^(1-(p-1)/(n-1)))
# }
log_gamma <- function(p,x){
  log((0.5*dnorm(x,mu1,1)+0.5*dnorm(x,mu2,3))^((p-1)/(n-1))) + log(dnorm(x,mu0,sd0/2)^(1-(p-1)/(n-1)))
}
# plot the targets
x_plt <- seq(from=-5,to=10,length.out=100)
plot(x_plt,exp(log_gamma(n,x_plt)),type='l')
matplot(x_plt,exp(log_gamma(1,x_plt)),add=T,type='l')
matplot(x_plt,exp(log_gamma(n/2,x_plt)),add=T,type='l')

smc_out <- smc(mu, M, G, n, N, mu0, sd0, sigma2_prop)
sum(smc_out$zetas[n,]*smc_out$w_normalised[n,])




# test pimh kernel - see the bias
R <- 1e2
n_pimh <- 50

logZ_arr <- matrix(NA,R,n_pimh)
rb_est <- matrix(NA,R,n_pimh)
smc_out <- list(log_Z=-Inf)
zeta_sample <- matrix(NA,R,n_pimh)
for(i in 1:R){
  for(j in 1:n_pimh){
    smc_out <- pimh_kernel(smc_out)
    logZ_arr[i,j] <- smc_out$log_Z
    rb_est[i,j] <- sum(smc_out$zetas[n,]*smc_out$w_normalised[n,])
    #zeta_sample[i,j] <- 
  }
}

plot(logZ_arr)
plt_mat <- matrix(c(rep(means[n],n_pimh),colMeans(rb_est)),nrow=2)
rb_plot <- colMeans(rb_est)
rb_avg_std <- sqrt(diag(var(rb_est))/R)
rb_u <- rb_plot+3*rb_avg_std
rb_l <- rb_plot-3*rb_avg_std
rb_ylims <- c(0.99*min(rb_l),1.01*max(rb_u))
plot(rb_plot,ylim=rb_ylims,type='l')
matplot((matrix(c(rb_u,rb_l),ncol=2)),type='l',add=T)
abline(h=(mu1+mu2),col='blue')
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




