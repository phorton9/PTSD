# Pragmatic Two-Stage Design for Clinical Trials with Treatment Selection

This repository contains the code used for the proposed method described in the paper *Pragmatic Two-Stage Design for Clinical Trials with Treatment Selection*.

## Requirements
To run the code, you need to install the `mvtnorm` library for calculations. You can install it in R using:

```R
install.packages("mvtnorm")
```

## Files
The repository includes five files that provide functionality for optimizing design parameters and running simulations:

1. Event_Functions
  * Contains functions for each component of the optimization problem.
2. Simulated Annealing_Functions
  * Includes the candidate generator and simulated annealing functions for optimization.
3. PTSD_Simulation
  * Simulates the PTSD results used to generate Tables 2, 3, and 4 in the reference paper.
4. Simulation_Function
  * Generates an optimization score over a fixed number of iterations (default: 10,000) based on input design parameters.
5. Example
  * Provides an example for optimizing the design under a specific set of parameters and simulating results.

## Usage
- Refer to the `example` file for guidance on how to use the code for both optimization and simulation tasks.
- Customize the design parameters in the provided functions to match your specific requirements.
