Overview: This code implements a function to estimate Linear Mixed Models (LMMs).
LMMs are useful for modeling data with multiple levels of variability(In this task,
"Machine" is fixed effect ,("Worker",c("Worker","Machine") are random effect.)
The implementation includes three main functions:
1. LMMsetup: Prepares model matrices for fixed and random effects.
2. LMMprof: Evaluates the negative log-likelihood for a given set of parameters and computes the corresponding estimates.
3. lmm: The main function that estimates model parameters by optimizing the profile likelihood.
