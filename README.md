# Step-vs-Ramp
Characterize neuronal firing rate dynamic (in LIP) with Bayesian model selection 

To test whether the single-trial population spike trains are more consistent with a linear model or with a step model, global_ratio_data.m
estimates Bayes Factor or the global likelihood ratio for these two models. The linear model assumes that each trial is an instantiation
of the same non-homogeneous Poisson process with a linear rate change over time. The step model assumes that each trial is an
instantiation of non-homogeneous Poisson processes with common initial rate and common final rate and a rate jump at some random time,
which can be different across trials.

Bollimunta, A., Totten, D., and Ditterich, J. Neural dynamics of choice: Single-trial analysis of decision-related activity in the
parietal cortex. Journal of Neuroscience, 32: 12684-12701 (2012).
