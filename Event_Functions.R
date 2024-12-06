library(mvtnorm)
gSelection <- function(x, mu1, var1, eps) {
  #' Calculate the probability of making a good treatment selection
  #'
  #' param x A numeric vector of six parameters: upper/lower bounds, threshold, proportion, and sample size.
  #' param mu1 A numeric vector of mean values for each group.
  #' param var1 A numeric variance for the groups.
  u1 <- x[1]  # Upper bound stage 1
  l1 <- x[2]  # Lower bound stage 1
  u2 <- x[3]  # Upper bound stage 2
  h <- x[4]   # Threshold
  p <- x[5]   # Proportion
  N <- x[6]   # Total sample size
  
  # Allocate total samples
  m1 <- p * N
  m2 <- (1 - p) * N
  
  K <- length(mu1)
  
  # Identify indices within epsilon range of max(mu1)
  L_list <- which(abs(mu1 - max(mu1)) <= eps)
  
  combined_gs <- 0
  for (L in L_list) {
    combinations_gs <- list()
    
    for (i in setdiff(1:K, L)) {
      for (x in c('greater', 'less')) {
        combinations_gs <- c(combinations_gs, list(list(i, c(), x)))
      }
    }
    
    for (i in 1:K) {
      missing <- setdiff(1:K, c(L, i))
      if (length(missing) > 1) {
        for (j in 1:length(missing)) {
          comb_j <- combn(missing, j, simplify = FALSE)
          for (comb in comb_j) {
            for (x in c('greater', 'less')) {
              combinations_gs <- c(combinations_gs, list(list(i, comb, x)))
            }
          }
        }
      } else {
        for (x in c('greater', 'less')) {
          combinations_gs <- c(combinations_gs, list(list(i, sort(missing), x)))
        }
      }
    }
    
    # Probability of making a good selection in the second stage.
    stage2 <- 0
    for (comb in combinations_gs) {
      max_arm <- comb[[1]][[1]]
      continue <- sort(unique(c(L, max_arm, comb[[2]])))
      control <- comb[[3]][[1]]
      
      comp1 <- setdiff(continue, L)
      comp2 <- setdiff(1:K, max_arm)
      
      lc0 <- length(comp1)
      
      mu <- rep(0, lc0 + K)
      cov1 <- matrix(1/2, nrow = length(mu), ncol = length(mu))
      
      tau1 <- floor(m1 / (K + 1))
      tau2 <- if (control == 'greater') {
        floor(m2 / length(continue))
      } else {
        floor(m2 / (length(continue) + 1))
      }
      
      sqX <- sqrt(tau1 / (tau1 + tau2))
      
      for (i in 1:length(mu)) {
        if (i <= lc0) {
          mu[i] <- (mu1[L] - mu1[comp1[i]]) * sqrt((tau1 + tau2) / 2)
        } else if (i <= lc0 + (K - 1)) {
          mu[i] <- (mu1[max_arm] - mu1[comp2[i - lc0]]) * sqrt(tau1 / 2)
        } else {
          mu[i] <- mu1[max_arm] * sqrt(tau1 / 2)
        }
      }
      
      # Adjust covariance matrix
      for (x in 1:lc0) {
        for (i in (lc0 + 1):(lc0 + (K - 1))) {
          if ((max_arm == L) & (comp1[x] == comp2[i - lc0])) {
            cov1[i, x] <- sqX
          } else if ((max_arm == L)) {
            cov1[i, x] <- sqX / 2
          } else if (comp1[x] == comp2[i - lc0]) {
            cov1[i, x] <- sqX / 2
          } else if ((comp1[x] == max_arm) & (comp2[i - lc0] == L)) {
            cov1[i, x] <- -sqX
          } else if (comp2[i - lc0] == L) {
            cov1[i, x] <- -sqX / 2
          } else if (comp1[x] == max_arm) {
            cov1[i, x] <- -sqX / 2
          } else {
            cov1[i, x] <- 0
          }
        }
      }
      
      for (x in 1:lc0) {
        for (i in (lc0 + K):(lc0 + K)) {
          cov1[i, x] <- if (L == max_arm) sqX / 2 else if (comp1[x] == max_arm) -sqX / 2 else 0
        }
      }
      
      for (x in (lc0 + 1):(lc0 + (K - 1))) {
        for (i in (lc0 + K):(lc0 + K)) {
          cov1[i, x] <- 1 / 2
        }
      }
      
      diag(cov1) <- 1
      cov1[upper.tri(cov1)] <- t(cov1)[upper.tri(cov1)]
      
      lower <- rep(0, length(mu))
      upper <- rep(0, length(mu))
      
      for (i in 1:length(upper)) {
        if (i <= lc0) {
          lower[i] <- 0
          upper[i] <- Inf
        } else if (i <= (lc0 + (K - 1))) {
          if (comp2[i - lc0] %in% continue) {
            lower[i] <- 0
            upper[i] <- h
          } else {
            lower[i] <- h
            upper[i] <- Inf
          }
        } else {
          lower[i] <- if (control == 'greater') u1 else l1
          upper[i] <- if (control == 'greater') Inf else u1
        }
      }
      stage2 <- stage2 + pmvnorm(lower = lower, upper = upper, mean = mu / sqrt(var1), corr = cov1)[1]
    }
    
    # Probability of making a good selection in the first stage.
    missing <- setdiff(1:K, L)
    meang <- rep(0, K-1)
    for (i in 1:(K-1)) {
      meang[i] <- (mu1[L] - mu1[missing[i]]) * sqrt(tau1 / 2)
    }
    lowerg1 <- rep(h, K-1)
    upperg1 <- rep(Inf, K-1)
    
    meang2 <- rep(0, K + (K-1))
    
    for (i in 1:K) {
      meang2[i] <- mu1[i] * sqrt(tau1 / 2)
    }
    for (i in 1:(K-1)) {
      meang2[i+K] <- (mu1[L] - mu1[missing[i]]) * sqrt(tau1 / 2)
    }
    
    lowerg2 <- c(rep(-Inf, K), rep(0, K-1))
    upperg2 <- c(rep(l1, K), rep(Inf, K-1))
    lowerg3 <- c(rep(-Inf, K), rep(h, K-1))
    upperg3 <- c(rep(l1, K), rep(Inf, K-1))
    
    covg1 <- matrix(1/2, nrow = length(meang), ncol = length(meang))
    diag(covg1) <- 1
    
    covg12 <- matrix(1/2, nrow = length(meang2), ncol = length(meang2))
    for (i in 2:K) {
      for (j in K:(K + (K-1))) {
        if ((i-1) == (j - K)) {
          covg12[j, i] <- -1/2
        } else {
          covg12[j, i] <- 0
        }
      }
    }
    diag(covg12) <- 1
    covg12[upper.tri(covg12)] <- t(covg12)[upper.tri(covg12)]
    
    gs11 <- pmvnorm(lower = lowerg1, upper = upperg1, mean = meang / sqrt(var1), corr = covg1)[1]
    gs12 <- pmvnorm(lower = lowerg2, upper = upperg2, mean = meang2 / sqrt(var1), corr = covg12)[1]
    gs13 <- pmvnorm(lower = lowerg3, upper = upperg3, mean = meang2 / sqrt(var1), corr = covg12)[1]
    
    # Probability of making a good selection in the first stage.
    stage1 <- gs11 + gs12 - gs13
  }
  combined_gs <- combined_gs + stage2+stage1
  combined_gs
}

type1 <- function(x, mu1, var1) {
  #' Calculate the maximum probability of rejecting H0 when H0 is true (type I error)
  #'
  #' param x A numeric vector of six parameters: upper/lower bounds, threshold, proportion, and sample size.
  #' param mu1 A numeric vector of mean values for each group.
  #' param var1 A numeric variance for the groups.
  u1 <- x[1]  # Upper bound stage 1
  l1 <- x[2]  # Lower bound stage 1
  u2 <- x[3]  # Upper bound stage 2
  h <- x[4]   # Threshold
  p <- x[5]   # Proportion
  N <- x[6]   # Total sample size
  eps <- 0
  
  K <- length(mu1)
  
  m1 <- p*N
  m2 <- (1-p)*N
  
  added <- c(2:K)
  new <- list()
  for (i in 1:length(added)){
    comb_i <- combn(added, i, simplify = FALSE)
    new <- c(new, comb_i)
  }
  
  new_list <- lapply(new, function(x) sort(c(x, 1)))
  
  combinations_list <- list(list(reject = 1, continue = 1))
  for (item in new_list) {
    for (i in 1:length(item)) {
      comb_i <- combn(item,i, simplify = FALSE)
      for (item_t in comb_i) {
        combinations_list <- append(combinations_list,list(list(reject = item_t, continue = item)))
      }
    }
  }
  
  ii <- 0
  
  # Probability of type I error in the second stage.
  stage2 <- 0
  for (comb in combinations_list) {
    ii <- ii+1
    cont1 <- length(comb$continue)
    lc0 <- length(comb$reject)
    tau1 <- floor(m1/(K+1))
    tau2 <- floor(m2/(cont1+1))
    
    mu <- rep(0, lc0 + (K-1) + 1)
    cov1 <- matrix(1/2, nrow = length(mu), ncol = length(mu))
    
    sqX <- sqrt(tau1/(tau1+tau2))
    
    
    for (x in 1:lc0) {
      for (i in (lc0+1):(lc0+(K-1))){
        if (comb$reject[x] == 1) {
          cov1[i,x] <- sqX/2
        }else if (comb$reject[x] == (i-(lc0-1))) {
          cov1[i,x] <- -sqX/2
        }
        else{
          cov1[i,x] <- 0
        }
      }
    }
    
    for (x in 1:lc0) {
      for (i in (lc0+K):(lc0+K)){
        if (comb$reject[x] == 1) {
          cov1[i,x] <- sqX
        }else{
          cov1[i,x] <- sqX/2
        }
      }
    }
    
    diag(cov1) <- 1
    cov1[upper.tri(cov1)] <- t(cov1)[upper.tri(cov1)]
    
    lower <- rep(0, length(mu))
    upper <- rep(0, length(mu))
    
    for (i in 1:length(upper)) {
      if (i <= lc0) {
        lower[i] <- u2
        upper[i] <- Inf
      } else if (i <= (lc0 + (K-1))) {
        if ((i-(lc0-1)) %in% comb$continue) {
          lower[i] <- 0
          upper[i] <- h
        } else{
          lower[i] <- h
          upper[i] <- Inf
        }
      } else {
        upper[i] <- u1
        lower[i] <- l1
      }
    }
    if (length(comb$reject) %% 2 == 1){
      stage2 <- stage2 + pmvnorm(lower = lower, upper = upper, mean = mu/sqrt(var1), sigma = cov1)[1]
    } else {
      stage2 <- stage2 - pmvnorm(lower = lower, upper = upper, mean = mu/sqrt(var1), sigma = cov1)[1]
    }
  }
  
  stage2 <- stage2*K
  
  # Probability of type I error in the first stage.
  stage1 <- 0
  for (i in 1:K){
    if (i == 1) {
      p1 <- 1-pnorm(u1)
    }
    else{
      lowers1 <- rep(u1,i)
      uppers1 <- rep(Inf,i)
      means1 <- rep(0,i)
      covs1 <- matrix(1/2, nrow = length(means1), ncol = length(means1))
      diag(covs1) <- 1
      p1 <- pmvnorm(lower = lowers1, upper = uppers1, mean = means1/sqrt(var1), corr = covs1)[1]
    }
    stage1 <- stage1 + ((i %% 2)*2-1)*choose(K,i)*p1
    stage2 + stage1
    
  }
  stage2 + stage1
}


power <- function(x, mu1, var1) {
  #' Calculate the probability of rejecting H0 under the Least Favorable Configuration (LFC) 
  #'
  #' param x A numeric vector of six parameters: upper/lower bounds, threshold, proportion, and sample size.
  #' param mu1 A numeric vector of mean values for each group.
  #' param var1 A numeric variance for the groups.
  u1 <- x[1]  # Upper bound stage 1
  l1 <- x[2]  # Lower bound stage 1
  u2 <- x[3]  # Upper bound stage 2
  h <- x[4]   # Threshold
  p <- x[5]   # Proportion
  N <- x[6]   # Total sample size
  eps <- 0
  
  m1 <- p * N
  m2 <- (1 - p) * N
  
  K <- length(mu1)
  
  mu_list <- c(mu1[1], rep(0, K - 1))
  
  # Probability of correctly rejecting null hypothesis in the second stage.
  stage2 <- 0
  for (L in 1:K){
    
    added <- setdiff(1:K, L)
    new <- list()
    for (i in 1:length(added)){
      comb_i <- combn(added, i, simplify = FALSE)
      new <- c(new, comb_i)
    }
    
    new_list <- lapply(new, function(x) sort(c(x, L)))
    
    combinations_list <- list(list(reject = L, continue = L))
    for (item in new_list) {
      for (i in 1:length(item)) {
        comb_i <- combn(item, i, simplify = FALSE)
        for (item_t in comb_i) {
          combinations_list <- append(combinations_list, list(list(reject = item_t, continue = item)))
        }
      }
    }
    
    ii <- 0
    for (comb in combinations_list) {
      ii <- ii + 1
      cont1 <- length(comb$continue)
      lc0 <- length(comb$reject)
      tau1 <- floor(m1 / (K + 1))
      tau2 <- floor(m2 / (cont1 + 1))
      comp1 <- setdiff(1:K, L)
      
      mu <- rep(0, lc0 + (K - 1) + 1)
      cov1 <- matrix(1 / 2, nrow = length(mu), ncol = length(mu))
      
      for (i in 1:length(mu)) {
        if (i <= lc0) {
          mu[i] <- mu_list[comb$reject[i]] * sqrt((tau1 + tau2) / 2)  # Use mu_list
        } else if (i <= lc0 + (K - 1)) {
          mu[i] <- (mu_list[L] - mu_list[comp1[i - lc0]]) * sqrt(tau1 / 2)  # Use mu_list
        } else {
          mu[i] <- mu_list[L] * sqrt(tau1 / 2)  # Use mu_list
        }
      }
      
      sqX <- sqrt(tau1 / (tau1 + tau2))
      
      for (x in 1:lc0) {
        for (i in (lc0 + 1):(lc0 + (K - 1))) {
          if (comb$reject[x] == L) {
            cov1[i, x] <- sqX / 2
          } else if (comb$reject[x] == comp1[i - lc0]) {
            cov1[i, x] <- -sqX / 2
          } else {
            cov1[i, x] <- 0
          }
        }
      }
      
      for (x in 1:lc0) {
        for (i in (lc0 + K):(lc0 + K)) {
          if (comb$reject[x] == L) {
            cov1[i, x] <- sqX
          } else {
            cov1[i, x] <- sqX / 2
          }
        }
      }
      
      diag(cov1) <- 1
      cov1[upper.tri(cov1)] <- t(cov1)[upper.tri(cov1)]
      
      lower <- rep(0, length(mu))
      upper <- rep(0, length(mu))
      
      for (i in 1:length(upper)) {
        if (i <= lc0) {
          lower[i] <- u2
          upper[i] <- Inf
        } else if (i <= (lc0 + (K - 1))) {
          if (comp1[i - lc0] %in% comb$continue) {
            lower[i] <- 0
            upper[i] <- h
          } else {
            lower[i] <- h
            upper[i] <- Inf
          }
        } else {
          upper[i] <- u1
          lower[i] <- l1
        }
      }
      if (length(comb$reject) %% 2 == 1) {
        stage2 <- stage2 + pmvnorm(lower = lower, upper = upper, mean = mu / sqrt(var1), sigma = cov1)[1]
      } else {
        stage2 <- stage2 - pmvnorm(lower = lower, upper = upper, mean = mu / sqrt(var1), sigma = cov1)[1]
      }
    }
    
  }
  stage2
  
  # Probability of correctly rejecting null hypothesis in first stage.
  stage1 <- 0
  for (i in 1:K) {
    if (i == 1) {
      p1 <- 1 - pnorm(u1, mean = mu_list[1] * sqrt(tau1 / 2) / sqrt(var1))  # Use mu_list
      p2 <- 1 - pnorm(u1)
      stage1 <- stage1 + p1 + (K - 1) * p2
    } else {
      lowers1 <- rep(u1, i)
      uppers1 <- rep(Inf, i)
      means1 <- c(mu_list[1] * sqrt(tau1 / 2), rep(0, i - 1))  # Use mu_list
      means2 <- rep(0, i)
      covs1 <- matrix(1 / 2, nrow = length(means1), ncol = length(means1))
      diag(covs1) <- 1
      p1 <- pmvnorm(lower = lowers1, upper = uppers1, mean = means1 / sqrt(var1), corr = covs1)[1]
      p2 <- pmvnorm(lower = lowers1, upper = uppers1, mean = means2 / sqrt(var1), corr = covs1)[1]
      if (i < K) {
        mult <- sum(combn(K, i)[1, ] == 1)
      } else {
        mult <- 1
      }
      stage1 <- stage1 + ((i %% 2) * 2 - 1) * (mult * p1 + (choose(K, i) - mult) * p2)
    }
  }
  stage2 + stage1
}



ESS <- function(x, mu1, var1) {
  #' Calculate the Expected Sample Size (ESS)
  #' 
  #' param x A numeric vector of six parameters: upper/lower bounds, threshold, proportion, and sample size.
  #' param mu1 A numeric vector of mean values for each group.
  #' param var1 A numeric variance for the groups.
  u1 <- x[1]  # Upper bound stage 1
  l1 <- x[2]  # Lower bound stage 1
  u2 <- x[3]  # Upper bound stage 2
  h <- x[4]   # Threshold
  p <- x[5]   # Proportion
  N <- x[6]   # Total sample size
  
  m1 <- p * N  # Sample size for stage 1
  m2 <- (1 - p) * N  # Sample size for stage 2
  K <- length(mu1) # Number of treatment arms
  tau1 <- floor(N * p / (K + 1))  # Stage 1 allocation
  
  comb_ESP <- list()
  
  # Generate all possible combinations for the expected sample size calculation
  for (L in 1:K) {
    missing <- setdiff(1:K, L)
    for (j in 1:length(missing)) {
      combinations <- combn(missing, j, simplify = FALSE)
      for (comb in combinations) {
        comb_ESP <- c(comb_ESP, list(list(L, comb)))
      }
    }
  }
  
  # Probability of continuing because of rejected null hypothesis but no treatment selection
  p_no_sel <- 0
  for (comb in comb_ESP) {
    mu <- rep(0, K)
    cov1 <- matrix(0.5, nrow = length(mu), ncol = length(mu))
    diag(cov1) <- 1
    comparisons <- comb[[2]]
    maxIndex <- comb[[1]]
    mu[1] <- mu1[maxIndex] * sqrt(tau1 / 2)
    
    missing <- setdiff(1:K, maxIndex)
    
    for (i in 1:length(missing)) {
      mu[i + 1] <- (mu1[maxIndex] - mu1[missing[i]]) * sqrt(tau1 / 2)
    }
    
    lower <- c(u1, rep(h, K - 1))
    upper <- rep(Inf, K)
    
    for (j in 1:(K - 1)) {
      if (missing[j] %in% comparisons) {
        if (missing[j] < maxIndex) {
          upper[missing[j] + 1] <- h
          lower[missing[j] + 1] <- 0
        } else {
          upper[missing[j]] <- h
          lower[missing[j]] <- 0
        }
      }
    }
    p_no_sel <- p_no_sel + pmvnorm(lower = lower, upper = upper, mean = mu / sqrt(var1), corr = cov1)[1]
  }
  
  p_no_dec <- 0  # Probability of continuing because no decision on hypothesis
  for (i in 1:K) {
    missing <- setdiff(1:K, i)
    
    mu <- rep(0, K)
    cov1 <- matrix(0.5, nrow = length(mu), ncol = length(mu))
    diag(cov1) <- 1
    
    for (x in 1:(K - 1)) {
      mu[x] <- (mu1[i] - mu1[missing[x]]) * sqrt(tau1 / 2)
    }
    mu[K] <- mu1[i] * sqrt(tau1 / 2)
    lower <- rep(0, K)
    upper <- rep(Inf, K)
    
    lower[K] <- l1
    upper[K] <- u1
    
    p_no_dec <- p_no_dec + pmvnorm(lower = lower, upper = upper, mean = mu / sqrt(var1), corr = cov1)[1]
  }
  
  tau1 * (K + 1) + (p_no_sel + p_no_dec) * m2
}


