# Unbiased Sequential Monte Carlo Samplers
The following includes code to produce unbiased Monte Carlo estimates of expectations obtained using SMC or particle methods. The code is generic enough to specify a target which can be annealed towards.
## Tests
Tests include estimating the normalising constant when the number of observations is equal to 1 and making sure that the target anneals towards the correct target. Additional tests include comparing the biased with the unbiased version to see the bias.
## Main functions
SMC




