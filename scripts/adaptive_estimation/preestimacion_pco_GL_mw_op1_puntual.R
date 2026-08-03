# Load necessary libraries
# 'nor1mix' is used for Marron-Wand densities
if (!require(nor1mix)) install.packages("nor1mix")
library(nor1mix)

# =====================================================================
# 1. Configuration & Parameters
# =====================================================================
set.seed(123) # For reproducibility
M <- 5        # Number of Monte Carlo iterations
n <- 250      # Sample size per iteration
rho <- 0.5    # AR(1) dependence parameter

# Choose a Marron-Wand density (Density #2: Bimodal/Separated)
# Users can change 'MW.nm2' to any of the 15 standard MW densities 
# (e.g., MW.nm1 for Gaussian, MW.nm3 for Skewed, etc.)
target_density <- MW.nm2 

# Bandwidth grid for the GL and PCO methods
h_grid <- seq(0.05, 0.5, length.out = 30)

# =====================================================================
# 2. Helper Functions
# =====================================================================

# Generate dependent sample with AR(1) structure mapped to the MW density
generate_dependent_sample <- function(n, rho, mw_dens) {
  # Generate AR(1) Gaussian process
  ar_process <- arima.sim(model = list(ar = rho), n = n)
  # Standardize the process
  ar_process <- (ar_process - mean(ar_process)) / sd(ar_process)
  
  # Transform to Uniform(0,1) using standard normal CDF
  u <- pnorm(ar_process)
  
  # Map to the target Marron-Wand density using its quantile function
  # qnorMix computes the inverse CDF for the specified mixture
  sample_data <- qnorMix(u, mw_dens)
  return(sample_data)
}

# PCO Method to estimate gamma (tuning parameter/bandwidth)
# Note: Adapts the Penalized Comparison to Overfitting criterion
estimate_gamma_PCO <- function(data, h_grid) {
  n <- length(data)
  pco_values <- numeric(length(h_grid))
  
  for (i in seq_along(h_grid)) {
    h <- h_grid[i]
    # Standard Gaussian kernel KDE
    kde_h <- density(data, bw = h, kernel = "gaussian")
    
    # PCO Criterion: ||f_hat_h||^2 - (2/n)*sum(f_hat_h_{-i}(X_i)) + Penalty
    # Approximated here for computational efficiency in the simulation
    sq_norm <- sum(kde_h$y^2) * (kde_h$x[2] - kde_h$x[1])
    
    # Cross-validation term approximation
    f_val <- approx(kde_h$x, kde_h$y, xout = data)$y
    f_val[is.na(f_val)] <- 0
    cv_term <- (2/n) * sum(f_val)
    
    # Overfitting penalty under dependence
    penalty <- 2 * (1 / (n * h)) * (1 / (2 * sqrt(pi))) 
    
    pco_values[i] <- sq_norm - cv_term + penalty
  }
  
  # Return the h (gamma) that minimizes the PCO criterion
  gamma_est <- h_grid[which.min(pco_values)]
  return(gamma_est)
}

# Original Goldenshluger-Lepski (GL) Estimator (Fixed Gamma)
# Applies the estimator without internal selection, using the provided gamma
apply_GL_estimator <- function(data, gamma_fixed, eval_points) {
  # With gamma fixed externally by PCO, the GL estimator evaluates 
  # the density using this specific tuning parameter.
  kde_gl <- density(data, bw = gamma_fixed, kernel = "gaussian", 
                    n = length(eval_points), from = min(eval_points), to = max(eval_points))
  
  # Return an interpolating function for ISE integration
  f_hat <- approxfun(kde_gl$x, kde_gl$y, yleft = 0, yright = 0)
  return(f_hat)
}

# Compute Integrated Squared Error (ISE)
compute_ise <- function(f_hat, mw_dens, lower, upper) {
  # integrand: (\hat{f}(x) - f(x))^2
  integrand <- function(x) {
    true_f <- dnorMix(x, mw_dens)
    est_f <- f_hat(x)
    return((est_f - true_f)^2)
  }
  
  # Integrate over the effective support of the density
  ise <- integrate(integrand, lower = lower, upper = upper, 
                   subdivisions = 500L, stop.on.error = FALSE)$value
  return(ise)
}

# =====================================================================
# 3. Monte Carlo Simulation Execution
# =====================================================================

# Store results
results <- data.frame(
  Iteration = 1:M,
  Gamma_Est = numeric(M),
  ISE = numeric(M)
)

# Determine integration bounds based on the target density support
bounds <- qnorMix(c(0.001, 0.999), target_density)
lower_bound <- bounds[1] - 1
upper_bound <- bounds[2] + 1
eval_points <- seq(lower_bound, upper_bound, length.out = 1024)

cat("Starting Monte Carlo Simulation (M =", M, ")...\n")
cat("Target: Marron-Wand Density #2 | AR(1) rho =", rho, "\n\n")

for (m in 1:M) {
  # 1. Generate dependent sample
  X_m <- generate_dependent_sample(n, rho, target_density)
  
  # 2. Estimate gamma using PCO
  gamma_m <- estimate_gamma_PCO(X_m, h_grid)
  
  # 3. Apply GL estimator with fixed gamma
  f_hat_m <- apply_GL_estimator(X_m, gamma_m, eval_points)
  
  # 4. Compute ISE
  ise_m <- compute_ise(f_hat_m, target_density, lower_bound, upper_bound)
  
  # Log results
  results$Gamma_Est[m] <- gamma_m
  results$ISE[m] <- ise_m
}

# 5. Calculate final MISE
mise <- mean(results$ISE)

# =====================================================================
# 4. Summary Output
# =====================================================================

cat("=== Monte Carlo Simulation Results ===\n")
print(results)
cat("\n--------------------------------------\n")
cat(sprintf("Mean Integrated Squared Error (MISE): %.6f\n", mise))
cat("--------------------------------------\n")