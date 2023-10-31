library(mvtnorm)
library(mcmcse)
library(coda)

# Example:
set.seed(10)
N <- 50
p <- 5
X <- matrix(rnorm(N * p), ncol = p)
true_theta <- (rep(0.2, p))
true_sigma <- 2
y <- rnorm(N, as.vector(X %*% matrix(true_theta, ncol = 1)), true_sigma)
theta_0 <- as.matrix(rep(2, p))
Sigma <- diag(p)
num_samples <- 50000
initial_theta <- rep(0.1, p)
s.cand <- 0.08


# Function to calculate the posterior log probability p(theta|y)
calculate_posterior_log <- function(y, X, true_theta, sigma, theta_0, Sigma) 
{
  likelihood <- dmvnorm( y, mean = X %*% true_theta, sigma = sigma)
  prior <- dmvnorm( true_theta, mean = theta_0, sigma = Sigma)
  posterior <- likelihood * prior
  return(posterior)
}

calculate_posterior_log(y,X,true_theta,true_sigma*diag(length(y)),theta_0,Sigma)
#sigma = true_sigma*diag(length(y))

# Function to calculate the inverse likelihood (1 / p(y|theta))
calculate_inverse_likelihood <- function(y, X, true_theta, sigma) 
{
  inverse_likelihood <- numeric(length(y))
  #i = 1
  for (i in 1:length(y)) 
  {
    inverse_likelihood[i] <- (1 / dnorm( y[i], mean = X[i,] %*% true_theta, sd = sigma))
  }
  return(inverse_likelihood)
}

calculate_inverse_likelihood(y,X,true_theta, true_sigma)


# Function to calculate the combined posterior log probability p(theta|y) * sum(1 / p(y|theta))
calculate_combined_posterior_log <- function(y, X, true_theta, sigma, theta_0, Sigma) 
{
  posterior_log <- calculate_posterior_log(y, X, true_theta, sigma*diag(length(y)), theta_0, Sigma)
  
  inverse_likelihood_sum <- sum(calculate_inverse_likelihood(y, X, true_theta, sigma))
  
  combined_posterior_log <- posterior_log * (inverse_likelihood_sum)
  
  return(combined_posterior_log)
}

calculate_combined_posterior_log(y,X,true_theta,true_sigma,theta_0,Sigma)


# Function to draw samples from the combined posterior using Metropolis-Hastings
theta_samples_from_log_q_mix_theta <- function(y, X, num_samples, theta_0, Sigma, sigma, s.cand) 
{
  num_parameters <- length(theta_0)
  samples <- matrix(0, nrow = num_samples, ncol = num_parameters)
  current_theta <- theta_0
  #s = 403
  for (s in 1:num_samples) 
  {
    proposal_theta <- current_theta + rnorm(num_parameters, 0, s.cand)
    
    # Calculate the Metropolis-Hastings acceptance ratio in log space
    combined_posterior_current_log <-  calculate_combined_posterior_log(y, X, as.numeric(current_theta), sigma, theta_0, Sigma)
    combined_posterior_proposal_log <- calculate_combined_posterior_log(y, X, as.numeric(proposal_theta), sigma, theta_0, Sigma)
    
    
    log_ratio <- min(1, combined_posterior_proposal_log / combined_posterior_current_log)
    
    
    # Print the log ratio
   # cat("Step:", s, "Log Ratio:", log_ratio, "\n")
    
    if ((runif(1)) < log_ratio) 
    {
      current_theta <- proposal_theta
    }
    
    samples[s, ] <- current_theta
  }
  
  return(samples)
}


# Draw samples from the combined posterior
theta_samples <- theta_samples_from_log_q_mix_theta (y, X, num_samples, initial_theta, Sigma, true_sigma, s.cand)

dim(theta_samples)
length(unique(theta_samples[,1])) / nrow(theta_samples)
# Assuming you have theta_samples matrix with num_samples rows and num_parameters columns


# Function to calculate weights for each data point
calculate_weights <- function(y, X, true_theta, sigma) 
{
  inverse_likelihoods <- numeric(length = length(y))
  
  for(i in 1: length(y))
  {
    inverse_likelihoods[i] <-  (1 / dnorm(y[i], mean = X[i,] %*% true_theta, sd = sigma))
  }
  
  weights <- inverse_likelihoods / sum(inverse_likelihoods)
  return(weights)
}


calculate_weighted_likelihoods <- function(y, X, theta_samples, sigma) 
{
  num_samples <- nrow(theta_samples)
  num_observations <- length(y)
  #i= 100
  #vector <- numeric()
  
  weighted_likelihoods <- matrix(0,num_samples,num_observations)
  weights <- matrix(0,num_samples,num_observations)
  #s =20
  
  for (s in 1:num_samples) 
  {
    weights[s,] <- calculate_weights(y, X, theta_samples[s,], true_sigma)
    weighted_likelihoods[s,] <- (dnorm(y, mean = as.numeric(X %*% (theta_samples[s, ])), sd = sqrt(true_sigma))) * (weights[s,])
  }
  
  # Store the weighted_likelihoods for this observation in the matrix
  vector <- cbind(weighted_likelihoods,(weights))
  
  return(vector)
}

# Example usage:
# Assuming you have y, X, theta_samples, and true_sigma defined
weighted_likelihoods_matrix <- calculate_weighted_likelihoods(y, X, theta_samples, true_sigma)

vector_1 <- weighted_likelihoods_matrix
dim(vector_1)

chain_1_cov <- mcse.multi(vector_1, r = 1)$cov ## using batch means estimator
eigen(chain_1_cov)$values

num_observations <- length(y)

# vec_1_sum <- numeric()
# for (i in 1:(2*num_observations)) 
# {
#   vec_1_sum[i] <- sum(vector_1[,i])
# }
vec_1_sum <- colMeans(vector_1)
# 'delta_g_matrix' is your Δg(x) matrix
calculate_delta_g_1_matrix <- function(x) 
{
  n <- length(x) / 2  # Assuming 'x' is your input vector of length 2n
  delta_g_1_matrix <- matrix(0, nrow = n, ncol = 2 * n)  # Initialize the result matrix
  
  for (i in 1:n) 
  {
    delta_g_1_matrix[i, i] <- 1 / (x[i + n])
    delta_g_1_matrix[i, i + n] <- - x[i] / ( x[i + n])^2
  }
  
  return(delta_g_1_matrix)
}
chain_2_cov <-  calculate_delta_g_1_matrix(vec_1_sum) %*%(chain_1_cov) %*% t(calculate_delta_g_1_matrix(vec_1_sum))

calculate_vector_2 <- function(x)
{
  n <- length(x) / 2  # Assuming 'x' is your input vector of length 2n
  g_2_matrix <- matrix(0, nrow = n, ncol =  1)  # Initialize the result matrix
  
  #i = 2
  for (i in 1:n)
  {
    g_2_matrix[i, ] <- x[i] / (x[i+n])
  }
  
  return(g_2_matrix)
}
vector_2 <- calculate_vector_2(vec_1_sum)


calculate_delta_g_2_matrix <- function(x) 
{
  n <- length(x)   # Assuming 'x' is your input vector of length 2n
  delta_g_2_matrix <- matrix(0, nrow = n, ncol =  n)  # Initialize the result matrix
  
  for (i in 1:n) 
  {
    delta_g_2_matrix[i, i] <- 1 / (x[i])
  }
  
  return(delta_g_2_matrix)
}

chain_3_cov <-  calculate_delta_g_2_matrix(vector_2) %*% (chain_2_cov) %*% t(calculate_delta_g_2_matrix(vector_2))

calculate_vector_3 <- function(x)
{
  n <- length(x)  # Assuming 'x' is your input vector of length 2n
  g_3_matrix <- matrix(0, nrow = n, ncol =  1)  # Initialize the result matrix
  #i = 2
  for (i in 1:n) 
  {
    g_3_matrix[i, ] <- log(x[i])
  }
  
  return(g_3_matrix)
}
vector_3 <- calculate_vector_3(vector_2)

calculate_delta_g_3_matrix <- function(x) 
{
  n <- length(x)   # Assuming 'x' is your input vector of length 2n
  delta_g_3_matrix <- matrix(0, nrow = 1, ncol =  n)  # Initialize the result matrix
  
  for (i in 1:n)
  {
    delta_g_3_matrix[, i] <- 1
  }
  return(delta_g_3_matrix)
}

chain_4_cov <-  calculate_delta_g_3_matrix(vector_3) %*%(chain_3_cov) %*% t(calculate_delta_g_3_matrix(vector_3))

print(paste("Final estimate of Variance is: ",round(chain_4_cov,4)))

sum(vector_3)
