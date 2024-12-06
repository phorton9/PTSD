results_simulation_3_arms <- function(parameters, mu_eval, varianceT){
  u1 <- parameters[1]  # Upper threshold for Stage 1
  l1 <- parameters[2]  # Lower threshold for Stage 1
  u2 <- parameters[3]  # Upper threshold for Stage 2
  h <- parameters[4]   # Decision boundary
  p <- parameters[5]   # Proportion of samples in Stage 1
  N <- parameters[6]   # Total sample size
  
  # Split sample size into Stage 1 and Stage 2
  m1 <- p * N
  m2 <- (1 - p) * N
  tau1 <- floor(m1 / (K + 1))  # Samples per arm in Stage 1
  
  K <- length(mu_eval)
  
  # Simulation Parameters
  iters <- 10000  # Number of iterations per configuration
  
  mu_comb <- list(
    mu_eval, #Linear
    c(0.24888, 0, 0)                         #LFC
    )
  
  # Initialize storage for results
  nsp1 <- numeric(length(mu_comb))
  gsp1 <- numeric(length(mu_comb))
  rjp1 <- numeric(length(mu_comb))
  
  # Simulation Loop
  xx <- 0  # Counter for treatment configurations
  for (mu in mu_comb) {
    # Initialize metrics for this configuration
    reject_1 <- numeric(iters)
    reject_2 <- numeric(iters)
    goodSelectiotau1 <- numeric(iters)
    goodSelectiotau2 <- numeric(iters)
    n0_samplesp <- numeric(iters)
    tau1_samplesp <- numeric(iters)
    tau2_samplesp <- numeric(iters)
    n3_samplesp <- numeric(iters)
    
    xx <- xx + 1  # Increment configuration counter
    
    for (x in 1:iters) {
      # Generate Stage 1 data
      placebo1_1 <- rnorm(tau1, sd = sqrt(varS))
      x11_1 <- rnorm(tau1, mean = mu[1], sd = sqrt(varS))
      x12_1 <- rnorm(tau1, mean = mu[2], sd = sqrt(varS))
      x13_1 <- rnorm(tau1, mean = mu[3], sd = sqrt(varS))
      
      # Calculate Stage 1 z-scores
      z11_1 <- (mean(x11_1) - mean(placebo1_1)) / sqrt(2 * varS) * sqrt(tau1)
      z12_1 <- (mean(x12_1) - mean(placebo1_1)) / sqrt(2 * varS) * sqrt(tau1)
      z13_1 <- (mean(x13_1) - mean(placebo1_1)) / sqrt(2 * varS) * sqrt(tau1)
      
      # Determine arms to keep based on z-scores
      maxZ1_1 <- max(z11_1, z12_1, z13_1)
      keep_1 <- which(abs(c(z11_1, z12_1, z13_1) - maxZ1_1) <= h)
      kept_placebo <- ifelse((maxZ1_1 > u1) | (maxZ1_1 < l1), 0, 1)
      
      # Allocate samples for Stage 2
      lenI <- length(keep_1)
      tau2 <- floor(m2 / (lenI + kept_placebo))
      placebo2_1 <- rnorm(tau2, sd = sqrt(varS))
      x21_1 <- rnorm(tau2, mean = mu[1], sd = sqrt(varS))
      x22_1 <- rnorm(tau2, mean = mu[2], sd = sqrt(varS))
      x23_1 <- rnorm(tau2, mean = mu[3], sd = sqrt(varS))
      
      # Calculate Stage 2 z-scores
      z21_1 <- (mean(c(x11_1, x21_1)) - mean(c(placebo1_1, placebo2_1))) / sqrt(2 * varS) * sqrt(tau1 + tau2)
      z22_1 <- (mean(c(x12_1, x22_1)) - mean(c(placebo1_1, placebo2_1))) / sqrt(2 * varS) * sqrt(tau1 + tau2)
      z23_1 <- (mean(c(x13_1, x23_1)) - mean(c(placebo1_1, placebo2_1))) / sqrt(2 * varS) * sqrt(tau1 + tau2)
      
      # Process outcomes for Stage 2
      kept_1 <- c(z21_1, z22_1, z23_1)[keep_1]
      if ((length(kept_1) == 1) | (maxZ1_1 < l1)) {
        goodSelectiotau1[x] <- maxZ1_1 == z11_1
      }
      if ((maxZ1_1 > l1) & (1 %in% keep_1) & (max(kept_1) == z21_1) & (length(kept_1) > 1)) {
        goodSelectiotau2[x] <- 1
      }
      
      # Update sample usage
      n0_samplesp[x] <- tau1
      tau1_samplesp[x] <- tau1
      tau2_samplesp[x] <- tau1
      n3_samplesp[x] <- tau1
      
      if (maxZ1_1 > u1) {
        reject_1[x] <- 1
        if (length(kept_1) > 1) {
          if (1 %in% keep_1) {
            tau1_samplesp[x] <- tau1 + tau2
          }
          if (2 %in% keep_1) {
            tau2_samplesp[x] <- tau1 + tau2
          }
          if (3 %in% keep_1) {
            n3_samplesp[x] <- tau1 + tau2
          }
        }
      } else if (maxZ1_1 > l1) {
        n0_samplesp[x] <- tau1 + tau2
        if (max(kept_1) > u2) {
          reject_2[x] <- 1
        }
        if (1 %in% keep_1) {
          tau1_samplesp[x] <- tau1 + tau2
        }
        if (2 %in% keep_1) {
          tau2_samplesp[x] <- tau1 + tau2
        }
        if (3 %in% keep_1) {
          n3_samplesp[x] <- tau1 + tau2
        }
      }
    }
    
    # Store results for this configuration
    nsp1[xx] <- mean(n0_samplesp + tau1_samplesp + tau2_samplesp + n3_samplesp)
    gsp1[xx] <- mean(goodSelectiotau1 + goodSelectiotau2)
    rjp1[xx] <- mean(reject_1 + reject_2)
  }
  
  nsp1[1] + 500*(1-gsp1[1]) + 500*(1-rjp1[2])
}