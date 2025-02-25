# Load required packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(Rcpp)
sourceCpp("mcem_multi_dye.cpp")

# --- Simulation study parameters ---
num_simulations <- 50   # Number of simulation replications
K <- 2                  # Number of dyes / parameters
Tsteps <- 1000          # Number of time steps in the simulation
N_init <- c(1000, 1000)   # Initial populations
q_init <- c(0.01, 0.1) # Starting values for MCEM estimation
M <- 500                # Number of MC samples per MCEM iteration
max_iter <- 2000         # Maximum number of MCEM iterations
tol <- 1e-5             # Convergence tolerance for MCEM

# Preallocate a data frame to store the true and estimated parameters
results <- data.frame(
  sim      = 1:num_simulations,
  true_q1  = numeric(num_simulations),
  true_q2  = numeric(num_simulations),
  est_q1   = numeric(num_simulations),
  est_q2   = numeric(num_simulations)
)

# --- Simulation Study Loop ---
set.seed(123)  # for reproducibility

for (i in 1:num_simulations) {
  # 1. Sample true decay parameters uniformly from, e.g., [0.01, 0.1]
  true_q <- runif(K, min = 0.01, max = 0.1)
  
  # 2. Simulate data using your simulation function.
  #    (Assumes simulate_multi_dye() returns a list that includes N_obs)
  sim_res <- simulate_multi_dye_continuous(K, true_q, N_init, Tsteps, q_first = 0.4)
  
  # Extract the observed counts vector
  N_obs <- sim_res$N_obs
  
  # 3. Estimate the parameters using the MCEM algorithm in C++.
  
  # After simulating, truncate N_obs to remove trailing zeros.
  last_nonzero <- max(which(sim_res$N_obs > 0))
  if(last_nonzero < length(sim_res$N_obs)) {
    cat("Truncating data at time step", last_nonzero - 1, "\n")
    N_obs_trim <- sim_res$N_obs[1:last_nonzero]
    Tsteps_trim <- last_nonzero - 1  # because N_obs has length Tsteps+1
  } else {
    N_obs_trim <- sim_res$N_obs
    Tsteps_trim <- Tsteps
  }
  
  # Then, call your MCEM function with the truncated data:
  est <- mcem_multi_dye_cpp(N_obs_trim, K, Tsteps_trim, q_init, M, max_iter, tol)
#  est <- mcem_multi_dye_cpp(N_obs, K, Tsteps, q_init, M, max_iter, tol)
  est_q <- est$q_hat
  
  # 4. Record the true and estimated values.
  results$true_q1[i] <- true_q[1]
  results$true_q2[i] <- true_q[2]
  results$est_q1[i]  <- est_q[1]
  results$est_q2[i]  <- est_q[2]
  
  cat("Simulation", i, ": true_q =", round(true_q, 4), 
      "; est_q =", round(est_q, 4), "\n")
}




# Suppose your simulation study results data frame is called `results`
# and has the columns:
#   sim, true_q1, true_q2, est_q1, est_q2

# Create a new data frame that holds the sorted (min, max) pair for each simulation.
results_sorted <- data.frame(
  sim       = rep(results$sim, each = 2),
  pair      = rep(c("q1", "q2"), times = nrow(results)),
  true      = numeric(2 * nrow(results)),
  estimated = numeric(2 * nrow(results))
)

# Loop through each simulation replication and sort the true and estimated pairs.
for (i in 1:nrow(results)) {
  # True and estimated pairs (unsorted)
  true_pair <- c(results$true_q1[i], results$true_q2[i])
  est_pair  <- c(results$est_q1[i],  results$est_q2[i])
  
  # Sort each pair (lowest first, highest second)
  sorted_true <- sort(true_pair)
  sorted_est  <- sort(est_pair)
  
  # Store in the new data frame:
  idx <- (2 * i - 1):(2 * i)
  results_sorted$true[idx]      <- sorted_true
  results_sorted$estimated[idx] <- sorted_est
}

# Now create a scatter plot comparing the sorted true vs. estimated values.
p_sorted <- ggplot(results_sorted, aes(x = true, y = estimated, color = pair)) +
  geom_point(size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "True Parameter Value",
       y = "Estimated Parameter Value",
       title = "Simulation Study: True vs. Estimated Decay Parameters") +
  theme_minimal(base_size = 14)

print(p_sorted)









simulate_multi_dye_continuous <- function(K, q_true, N_init, Tsteps, dt = 1, q_first = NULL) {
  # Optionally: if you still want a special treatment for the first interval,
  # you could simulate a discrete decay at time 0 using q_first.
  # For now, we assume that after t=0, the process is fully continuous.
  
  # Convert discrete q to continuous lambda.
  # Here we assume dt=1, so lambda = -log(1 - q_true).
  lambda <- -log(1 - q_true)
  
  max_time <- Tsteps * dt
  current_time <- 0
  
  # Initialize state: each dye's count.
  state <- N_init
  
  # Record the times at which events occur and the corresponding states.
  event_times <- c(0)
  state_history <- matrix(state, nrow = 1)
  
  # Gillespie simulation loop: simulate until we exceed max_time or have no photons left.
  while (current_time < max_time && sum(state) > 0) {
    # Calculate propensities for each dye decay.
    # For dye k, the decay propensity is lambda[k] * (number of photons in dye k).
    propensities <- lambda * state
    total_propensity <- sum(propensities)
    
    # If no decays can occur, exit the loop.
    if (total_propensity <= 0) break
    
    # Draw the waiting time until the next event (exponentially distributed).
    dt_event <- rexp(1, rate = total_propensity)
    current_time <- current_time + dt_event
    
    # If the next event is beyond max_time, we stop.
    if (current_time > max_time) break
    
    # Determine which dye decays (select reaction proportional to its propensity).
    reaction_idx <- sample(1:K, size = 1, prob = propensities)
    
    # Update the state: one photon decays from the selected dye.
    state[reaction_idx] <- state[reaction_idx] - 1
    
    # Record the event time and new state.
    event_times <- c(event_times, current_time)
    state_history <- rbind(state_history, state)
  }
  
  # Now, bin the continuous process at regular time intervals (e.g., t = 0, dt, 2*dt, ..., max_time).
  # For each bin, we record the most recent state (i.e. count) prior to or at that time.
  bin_times <- seq(0, max_time, by = dt)
  N_obs <- numeric(length(bin_times))
  
  for (i in seq_along(bin_times)) {
    t_bin <- bin_times[i]
    # Find the index of the last event that occurred before or at t_bin.
    idx <- max(which(event_times <= t_bin))
    # Sum the counts over all dyes (or you can record per-dye if desired).
    N_obs[i] <- sum(state_history[idx, ])
  }
  
  # For completeness, you can also bin the latent states (per dye) if you need them.
  # For instance, record a matrix of counts for each dye at the binned times.
  Nmat <- matrix(0, nrow = length(bin_times), ncol = K)
  for (i in seq_along(bin_times)) {
    t_bin <- bin_times[i]
    idx <- max(which(event_times <= t_bin))
    Nmat[i, ] <- state_history[idx, ]
  }
  
  # Return a list similar to your original simulation function.
  list(
    K = K,
    Tsteps = length(bin_times) - 1,  # number of intervals
    q_true = q_true,
    lambda = lambda,
    N_init = N_init,
    Nmat = Nmat,      # latent states at binned times
    N_obs = N_obs,    # total observed counts at binned times
    times = bin_times # the times corresponding to each bin
  )
}
