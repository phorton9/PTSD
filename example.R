# Simulated Annealing Setup Script
# This file contains both an example of optimizing input parameters as well as
# one simluation for three treatment arms.

# Load required functions
source("Simulated_Annealing_Functions.R", chdir = TRUE)
source("Simulation_Function.R", chdir = TRUE)

# Define upper and lower bounds for candidate generation
lower_lim <- c(2, -1, 2, 0, 0, 50)
upper_lim <- c(3, 1, 3, 1, 0.95, 500)

# Temperature decay configuration
t0 <- 1000  # Initial temperature
t_ratio <- 0.6  # Decay ratio
temp_array <- c(t0)
t <- t0

# Generate the temperature schedule
while ((t > 0.01) || length(temp_array) < 15) {
  t <- t * t_ratio
  temp_array <- c(temp_array, t)
}

# Candidate generation range decay
decay <- c(rep(1, 5), seq(1, 0.1, -0.1), rep(0.05, length(temp_array) - 15))

# Problem-specific parameters
mu <- c(0.24888, 0.24888 * 2 / 3, 0.24888 / 3)  # Mean values
ex_variance <- 0.2488  # Example variance
epsilon <- 0  # Threshold for good treatment selection
alpha <- 0.025  # Type I error threshold
lambda1 <- 500  # Penalty parameter 1
lambda2 <- 500  # Penalty parameter 2
initial_point <- c(2.69, -0.8, 2.3, 0.1, 0.5, 150)  # Starting point
nLoops <- 10  # Number of iterations per temperature level

# Run the simulated annealing algorithm
result <- simulated_Annealing(
  temperatures = temp_array,
  decays = decay,
  lower_lim = lower_lim,
  upper_lim = upper_lim,
  start_point = initial_point,
  lambda1 = lambda1,
  lambda2 = lambda2,
  mu1 = mu,
  var1 = ex_variance,
  eps = epsilon,
  alpha = alpha,
  iters = nLoops
)

# Display the result
print(result)

# Simulation
sim_results <- results_simulation_3_arms(result[2:7], mu, ex_variance)
print(sim_results)
