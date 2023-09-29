library(mvtnorm)

# Function to calculate the posterior log probability p(theta|y)
calculate_posterior_log <- function(y, X, true_theta, sigma, theta_0, Sigma) {
  likelihood <- sum(log(dnorm(y, mean = X %*% true_theta, sd = sigma)))
  prior <- dmvnorm(true_theta, mean = theta_0, sigma = Sigma, log = TRUE)
  posterior <- likelihood + prior
  return(posterior)
}

# Function to calculate the inverse likelihood (1 / p(y|theta))
calculate_inverse_likelihood <- function(y, X, true_theta, sigma) {
  inverse_likelihood <- 1 / dnorm(y, mean = X %*% true_theta, sd = sigma)
  return(inverse_likelihood)
}

# Function to calculate the combined posterior log probability p(theta|y) * sum(1 / p(y|theta))
calculate_combined_posterior_log <- function(y, X, true_theta, sigma, theta_0, Sigma) {
  posterior_log <- calculate_posterior_log(y, X, true_theta, sigma, theta_0, Sigma)
  inverse_likelihood_sum <- sum(calculate_inverse_likelihood(y, X, true_theta, sigma))
  combined_posterior_log <- posterior_log + log(inverse_likelihood_sum)
  return(combined_posterior_log)
}

# Function to draw samples from the combined posterior using Metropolis-Hastings
theta_samples_from_combined_posterior <- function(y, X, num_samples, theta_0, Sigma, sigma, s.cand) {
  num_parameters <- length(theta_0)
  samples <- matrix(0, nrow = num_samples, ncol = num_parameters)
  current_theta <- theta_0
  
  for (s in 1:num_samples) {
    proposal_theta <- as.vector(rmvnorm(1, mean = current_theta, sigma = s.cand * diag(num_parameters)))
    
    # Calculate the Metropolis-Hastings acceptance ratio in log space
    combined_posterior_current_log <- calculate_combined_posterior_log(y, X, current_theta, sigma, theta_0, Sigma)
    combined_posterior_proposal_log <- calculate_combined_posterior_log(y, X, proposal_theta, sigma, theta_0, Sigma)
    
    log_ratio <- combined_posterior_proposal_log - combined_posterior_current_log
    
    if (!is.na(log_ratio) && log(runif(1)) < log_ratio) {
      current_theta <- proposal_theta
    }
    
    samples[s, ] <- current_theta
  }
  
  return(samples)
}

# Example usage:
set.seed(1)
N <- 100
X <- matrix(rnorm(N * 50), ncol = 50)
true_theta <- rep(c(0.2, -0.5), 25)
true_sigma <- 1
y <- rnorm(N, as.vector(X %*% matrix(true_theta, ncol = 1)), true_sigma)
theta_0 <- rep(0, 50)
Sigma <- diag(50)
num_samples <- 2000
initial_theta <- rep(c(0.1, -0.4), 25)
s.cand <- 0.002
theta_samples <- theta_samples_from_combined_posterior(y, X, num_samples, initial_theta, Sigma, true_sigma, s.cand)

# Function to calculate weights for each data point
calculate_weights <- function(y, X, true_theta, true_sigma) {
  inverse_likelihoods <- 1 / dnorm(y, mean = X %*% true_theta, sd = true_sigma)
  weights <- inverse_likelihoods / sum(inverse_likelihoods)
  return(weights)
}

# Function to calculate the desired quantity (predictive distribution for each data point)
calculate_mixture_estimator <- function(y, X, theta_samples, true_sigma) {
  num_samples <- nrow(theta_samples)
  num_observations <- length(y)
  
  # Initialize variables
  numerator <- numeric(num_observations)
  denominator <- numeric(num_observations)
  desired_quantity <- numeric(num_observations)
  
  for (i in 1:num_observations) {
    weighted_likelihoods <- numeric(num_samples)
    weights <- numeric(num_samples)
    
    for (s in 1:num_samples) {
      weights[s] <- calculate_weights(y[i], X[i, ], theta_samples[s, ], true_sigma)
      weighted_likelihoods[s] <- dnorm(y[i], mean = X[i, ] %*% theta_samples[s, ], sd = true_sigma) * weights[s]
    }
    
    numerator[i] <- sum(weighted_likelihoods)
    denominator[i] <- sum(weights)
    
    desired_quantity[i] <- numerator[i] / denominator[i]
  }
  
  return(desired_quantity)
}


# Calculate the desired quantity (predictive distribution) for each data point
mixture_estimator <- calculate_mixture_estimator(y, X, theta_samples, true_sigma)

mu_i_hat = mixture_estimator
mix= (log(mu_i_hat))


# Compute p(y_i|y_(-i))
compute_conditional_pdf_Y_i_given_Y_minus_i <- function(Y, X, Sigma_0, mu_0, sigma_0) {
  n <- length(Y)  # Number of observations
  pdf_values <- numeric(n)  # Initialize a vector to store PDF values
  
  for (i in 1:n) {
    # Exclude the i-th observation
    Y_minus_i <- Y[-i]
    X_minus_i <- X[-i, ]  # Corrected X for X_{-i}
    
    # Extract the i-th observation
    Y_i <- Y[i]
    X_i <- X[i, ]
    
    # Calculate Lambda_minus_i
    Lambda_minus_i <- solve(Sigma_0) + (t(X_minus_i) %*% X_minus_i) / (sigma_0^2)
    
    # Calculate eta_minus_i
    eta_minus_i <- solve(solve(Sigma_0) + (1 / (sigma_0^2)) * t(X_minus_i) %*% X_minus_i) %*%
      (solve(Sigma_0) %*% mu_0 + (1 / (sigma_0^2)) * t(X_minus_i) %*% Y_minus_i)
    
    # Calculate mean and variance of the conditional distribution
    mean_Y_i_given_Y_minus_i <- X_i %*% eta_minus_i
    variance_Y_i_given_Y_minus_i <- sigma_0^2 + t(X_i) %*% Lambda_minus_i %*% (X_i)
    
    # Compute the PDF value for Y_i
    pdf_values[i] <- dnorm(Y_i, mean = mean_Y_i_given_Y_minus_i, sd = sqrt(variance_Y_i_given_Y_minus_i))
  }
  
  return(pdf_values)
}

test <- compute_conditional_pdf_Y_i_given_Y_minus_i(y,X,Sigma,theta_0,true_sigma)
original = (log(test))

# Calculate the Mean Squared Error (MSE)
mse <- mean((original - mix)^2)
print(paste("MSE:", mse))



