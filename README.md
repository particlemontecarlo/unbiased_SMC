# Unbiased Sequential Monte Carlo Samplers
The following includes code to produce unbiased Monte Carlo estimates of expectations obtained using SMC or particle methods. The code is generic enough to specify a target which can be annealed towards.
## Tests
Tests include estimating the normalising constant when the number of observations is equal to 1 and making sure that the target anneals towards the correct target. Additional tests include estimating the mean of a mixture of gaussians for which the true value is available analytically. 
## Main functions
The main functionality is presented in 'unbiased_SMC_example.R' where example code is provided both for simple SMC and unbiased SMC. The estimators for the unbiased SMC are generated sequentially. It is recommended to modify the sequence of target distributions to begin with. 




