# Teamwork：
# s2626102 Jingwen Jiang 40%, Write the main content of the code and Write comments.
# s2752993 Ziyi Ye       30%, Discuss, modify, refine the code and Write comments.
# s2646482 Chengpeng Dai 30%, Discuss, modify, refine the code and Write comments.

# Overview: This code implements a function to estimate Linear Mixed Models (LMMs).
# LMMs are useful for modeling data with multiple levels of variability(In this task,
# "Machine" is fixed effect ,("Worker",c("Worker","Machine") are random effect.)
# The implementation includes three main functions:
# 1. LMMsetup: Prepares model matrices for fixed and random effects.
# 2. LMMprof: Evaluates the negative log-likelihood for a given set of parameters and computes the corresponding estimates.
# 3. lmm: The main function that estimates model parameters by optimizing the profile likelihood.

LMMsetup <- function(form, dat, ref) {
  # Purpose: Prepares the design matrices for fixed effects (X) and random effects (Z), and calculates dimensions.
  # Inputs:
  #   - form: A formula specifying the fixed effects.
  #   - dat: The dataset containing all variables.
  #   - ref: A list specifying the random effects structure.
  # Outputs:
  #   - A list containing the fixed effect matrix (X), random effect matrix (Z), and dimensions (t2, t3).
  X <- model.matrix(form, data = dat) # Fixed effects design matrix
  Z <- NULL # Initialize the random effects design matrix
  for (block in ref) {
    # Create random effects model matrix for each block in ref
    Z_block <- model.matrix(as.formula(paste("~", paste(block, collapse = ":"), "- 1")), data = dat)
    Z <- cbind(Z, Z_block) # Combine blocks into Z
  }
  t <- calculate_counts(ref,dat)
  return(list(X = X, Z = Z, t=t)) # Return setup
}


calculate_counts <- function(ref, dat) {
  # Purpose:
  # Calculate the number of unique values for single or multiple columns defined in `ref`.
  # If `ref` contains a single column, count the unique values in that column.
  # If `ref` contains multiple columns, calculate the product of unique counts for those columns.
  
  # Input:
  # ref: A list where each element is either a single column name (as a string) 
  #      or a vector of multiple column names.
  # dat: A data frame containing the columns specified in `ref`.
  
  # Output:
  # A numeric vector where each element corresponds to the calculated count
  # for the respective element in `ref`.
  
  # Apply the calculation to each element of `ref`
  counts <- sapply(ref, function(I) {
    # If the element of `ref` (I) has only one column name
    if (length(I) == 1) {
      # Calculate the number of unique values in the corresponding column of `dat`
      return(length(unique(dat[[I]])))
    } else {
      # Calculate the product of unique values across all specified columns
      return(prod(sapply(I, function(x) length(unique(dat[[x]])))))
    }
  })
  # Return the resulting counts
  return(counts)
}


LMMprof <- function(theta, setup, y) {
  # Purpose: Computes the negative log-likelihood and estimates for given parameters.
  # Inputs:
  #   - theta: A vector of parameters (log-transformed variances).
  #   - setup: Output of LMMsetup containing design matrices and dimensions.
  #   - y: The response variable vector.
  # Outputs:
  #   - A list containing the negative log-likelihood and fixed effect estimates.
  X <- setup$X
  Z <- setup$Z
  t <- setup$t
  
  # QR decomposition of Z
  QR_decomp <- qr(Z)
  R <- qr.R(QR_decomp)
  n <- nrow(X)  # number of observations
  p <- ncol(R) #number of random effects
  
  # Variance components
  sigma <- exp(theta[1]) # Residual variance
  psi_diag <- unlist(lapply(2:length(theta), function(i) {rep(exp(theta[i])^2, t[i - 1])})) # Random effect variances
  # Construct covariance matrix
  A <- R %*% diag(psi_diag) %*% t(R) + diag(sigma^2, p)
  U <- chol(A) # Cholesky decomposition of A
  k <- forwardsolve(t(U), diag(1, nrow(A))) # Solve lower triangular system
  W_tl <- backsolve(U, k) # Inverse of A upper triangular
  W_br <- diag(sigma^(-2), n - p) # Residual weight matrix
  W_combined <- as.matrix(bdiag(W_tl, W_br)) # Combine weight matrices
  
  # Apply weights using QR decomposition
  W_combined_qy <- qr.qy(QR_decomp, W_combined)
  X_qty <- qr.qty(QR_decomp, X)
  y_qty <- qr.qty(QR_decomp, y)
  
  # Compute fixed effect estimates
  XtWX <- t(X) %*% W_combined_qy %*% X_qty
  XtWy <- t(X) %*% (W_combined_qy %*% y_qty)
  L <- chol(XtWX) # Cholesky decomposition of XtWX
  z <- forwardsolve(t(L), XtWy)
  beta_hat <- backsolve(L, z) # Compute beta estimates
  
  # Residual calculations
  residual <- y - X %*% beta_hat
  
  # Negative log-likelihood
  lw1 <- 0.5 * t(residual) %*% W_combined_qy %*% qr.qty(QR_decomp, residual) # Likelihood term 1
  lw2 <- sum(log(diag(U))) + 0.5 * (n - p) * log(sigma^2) # Likelihood term 2
  neg_log_likelihood <- lw1 + lw2
  return(list(nll = as.numeric(neg_log_likelihood), beta = beta_hat))
}

lmm <- function(form, dat, ref = list()) {
  # Purpose: Main function to estimate model parameters using optimization.
  # Inputs:
  #   - form: A formula for the fixed effects.
  #   - dat: The dataset containing all variables.
  #   - ref: A list specifying the random effects structure (default is no random effects).
  # Outputs:
  #   - A list containing the fixed effect estimates, variance estimates, and negative log-likelihood.
  setup <- LMMsetup(form, dat, ref) # Prepare model matrices
  y <- model.response(model.frame(form, data = dat)) # Extract response variable
  init_theta <- rep(1, 1+length(ref)) # Initial parameter values
  lower_bounds <- rep(-10, length(init_theta)) # Lower bounds for optimization
  upper_bounds <- rep(10, length(init_theta)) # Upper bounds for optimization
  
  # Optimize the negative log-likelihood
  opt <- optim(
    init_theta,
    fn = function(theta) LMMprof(theta, setup, y)$nll, 
    method = "L-BFGS-B",
    lower = lower_bounds,
    upper = upper_bounds
  )
  
  final_theta <- opt$par # Optimized parameters
  final_results <- LMMprof(final_theta, setup, y) # Final model evaluation
  return(list(beta = final_results$beta, theta = final_theta, neg_log_likelihood = final_results$nll))
}



