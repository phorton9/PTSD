# Generate a New Candidate Sample for Simulated Annealing
# This function generates a new sample point within a decayed range of the current point.
# Args:
#   lower_lim: Lower bounds of the search space.
#   upper_lim: Upper bounds of the search space.
#   current_point: Current position in the search space.
#   decay_val: Decay factor controlling the range of exploration.
# Returns:
#   A new candidate point within the specified range.
next_sample_gen <- function(lower_lim, upper_lim, current_point, decay_val) {
  range <- (upper_lim - lower_lim) * decay_val
  new_u_lim <- current_point + 0.5 * range
  new_l_lim <- current_point - 0.5 * range
  
  # Ensure the new limits stay within the global bounds
  rlim_upper <- pmin(new_u_lim, upper_lim)
  rlim_lower <- pmax(new_l_lim, lower_lim)
  
  # Generate a random point within the adjusted limits
  runif(length(rlim_upper)) * (rlim_upper - rlim_lower) + rlim_lower
}

# Acceptance Probability Calculation for Simulated Annealing
# Determines whether to accept a new candidate point based on its score and the current temperature.
# Args:
#   score_candidate: Score of the new candidate point.
#   score_current_point: Score of the current point.
#   temperature: Current annealing temperature.
# Returns:
#   Acceptance probability (a value between 0 and 1).
temp_check <- function(score_candidate, score_current_point, temperature) {
  min(1, exp(-(score_candidate - score_current_point) / temperature))
}

# Simulated Annealing Optimization
# Performs the simulated annealing algorithm to optimize a given fitness function.
# Args:
#   temperatures: Vector of temperatures for the annealing schedule.
#   decays: Vector of decay values controlling exploration range at each temperature step.
#   lower_lim: Lower bounds of the search space.
#   upper_lim: Upper bounds of the search space.
#   start_point: Initial point for the search.
#   lambda1, lambda2: Penalty parameters for the fitness function.
#   mu1, var1: Problem-specific parameters for the fitness function.
#   eps: Threshold for "good" treatment selection.
#   alpha: Type I error threshold.
#   iters: Number of iterations to perform at each temperature level.
# Returns:
#   A vector containing the current score and parameters for the current point.
simulated_Annealing <- function(temperatures, decays, lower_lim, upper_lim, start_point, lambda1, lambda2, mu1, var1, eps, alpha, iters) {
  # Initialize the search
  current_point <- start_point
  cur_score <- fitness_function(current_point, lambda1, lambda2, mu1, var1, eps, alpha)
  best_score <- Inf  # Initialize with a very large value
  best_candidate <- start_point
  
  # Iterate through temperature levels
  for (i in seq_along(temperatures)) {
    for (x in seq_len(iters)) {
      # Generate a new candidate point
      new_candidate <- next_sample_gen(lower_lim, upper_lim, current_point, decays[i])
      new_score <- fitness_function(new_candidate, lambda1, lambda2, mu1, var1, eps, alpha)
      
      # Calculate acceptance probability
      accept_threshold <- temp_check(new_score, cur_score, temperatures[i])
      
      # Handle edge cases where acceptance probability is invalid
      if (is.na(accept_threshold) || is.null(accept_threshold)) {
        accept_threshold <- 0
      }
      
      # Accept or reject the new candidate
      if (runif(1) < accept_threshold) {
        current_point <- new_candidate
        cur_score <- new_score
      }
      
      # Update the best candidate if the current score is better
      if (cur_score < best_score) {
        best_score <- cur_score
        best_candidate <- current_point
      }
    }
  }
  
  # Return the optimization results
  c(current_score = cur_score, 
    u1 = current_point[1], 
    l1 = current_point[2], 
    u2 = current_point[3], 
    h = current_point[4], 
    p = current_point[5], 
    N = current_point[6])
}

source("Event_Functions.R", chdir = T)

# Penalty Function
# Applies a penalty based on the type 1 error.
penalty_function <- function(x, mu1, var1, alpha) {
  100000*max(0,type1(x,mu1,var1) - alpha)
}

# Penalized optimization function.
fitness_function <- function(x, lambda1, lambda2, mu1, var1, eps, alpha) {
  ESS(x, mu1, var1) + 
    lambda1 * (1 - gSelection(x, mu1, var1, eps)) + 
    lambda2 * (1 - power(x, mu1, var1)) + 
    penalty_function(x, mu1, var1, alpha)
}


