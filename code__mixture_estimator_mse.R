library(mvtnorm)

# Function to calculate the posterior log probability p(theta|y)
calculate_posterior_log <- function(y, X, true_theta, sigma, theta_0, Sigma) 
{
  likelihood <- dmvnorm( y, mean = as.matrix(X %*% true_theta), sigma = sigma)
  prior <- dmvnorm( true_theta, mean = theta_0, sigma = Sigma)
  posterior <- likelihood * prior
  return(posterior)
}

# Function to calculate the inverse likelihood (1 / p(y|theta))
calculate_inverse_likelihood <- function(y, X, true_theta, sigma) 
{
  inverse_likelihood <- numeric(length(y))

    for (i in 1:length(y)) 
  {
    inverse_likelihood[i] <- (1 / dnorm( y[i], mean = X[i,] %*% true_theta, sd = sqrt(sigma)))
  }
  return(inverse_likelihood)
}


# Function to calculate the combined posterior log probability p(theta|y) * sum(1 / p(y|theta))
q_mix_theta <- function(y, X, true_theta, sigma, theta_0, Sigma) 
{
  posterior <- calculate_posterior(y, X, true_theta, sigma*diag(length(y)), theta_0, Sigma)
  
  inverse_likelihood_sum <- sum(calculate_inverse_likelihood(y, X, true_theta, sigma))
  
  combined_posterior <- posterior * (inverse_likelihood_sum)
  
  return(combined_posterior)
}


# Function to draw samples from the combined posterior using Metropolis-Hastings
theta_samples_from_q_mix_theta <- function(y, X, num_samples, theta_0, Sigma, sigma, s.cand) 
{
  num_parameters <- length(theta_0)
  samples <- matrix(0, nrow = num_samples, ncol = num_parameters)
  current_theta <- theta_0

  for (s in 1:num_samples) 
  {
    proposal_theta <- current_theta + rnorm(num_parameters, 0, s.cand)
    
    # Calculate the Metropolis-Hastings acceptance ratio in log space
    q_mix_theta_current <-  q_mix_theta(y, X, as.numeric(current_theta), sigma, theta_0, Sigma)
    q_mix_theta_proposal <- q_mix_theta(y, X, as.numeric(proposal_theta), sigma, theta_0, Sigma)
    
    
    ratio <- min(1, q_mix_theta_proposal / q_mix_theta_current)
    
    
    # Print the ratio
    #cat("Step:", s, "Ratio:", ratio, "\n")
    
    if (log(runif(1)) < log(ratio)) 
    {
      current_theta <- proposal_theta
    }
    
    samples[s, ] <- current_theta
  }
  
  return(samples)
}

# Function to calculate weights for each data point
calculate_weights <- function(y, X, true_theta, sigma) 
{
  inverse_likelihoods <- numeric(0)
  
  for(i in 1: length(y))
  {
    inverse_likelihoods[i] <-  (1 / dnorm(y[i], mean = X[i,] %*% true_theta, sd = sqrt(sigma)))
  }
  
  weights <- inverse_likelihoods / sum(inverse_likelihoods)
  return(weights)
}


calculate_weighted_likelihoods <- function(y, X, theta_samples, sigma) 
{
  num_samples <- nrow(theta_samples)
  num_observations <- length(y)
  
  vector <- numeric()
  
  weighted_likelihoods <- matrix(0,num_samples,num_observations)
  weights <- matrix(0,num_samples,num_observations)
  
  for (s in 1:num_samples) 
  {
    weights[s,] <- calculate_weights(y, X, theta_samples[s,], true_sigma)
    weighted_likelihoods[s,] <- (dnorm((y), mean = as.numeric(X %*% (theta_samples[s, ])), sd = sqrt(true_sigma))) * (weights[s,])
  }
  
  # Store the weighted_likelihoods for this observation in the matrix
  vector <- cbind(weighted_likelihoods,(weights))
  
  return(vector)
}

calculate_vector_2 <- function(x)
{
  n <- length(x) / 2  # Assuming 'x' is your input vector of length 2n
  g_2_matrix <- matrix(0, nrow = n, ncol =  1)  # Initialize the result matrix
  
  for (i in 1:n) 
  {
    g_2_matrix[i, ] <- log(x[i]) - log(x[i+n])
  }
  
  return(g_2_matrix)
}


# Example:
set.seed(10)
N <- 50
p = 15
X <- matrix(rnorm(N * p), ncol = p)
true_theta <- (rep(0.2, p))
true_sigma <- 2
y <- rnorm(N, as.vector(X %*% matrix(true_theta, ncol = 1)), true_sigma)
theta_0 <- as.matrix(rep(2, p))
Sigma <- diag(p)
num_samples <- 1e5
initial_theta <- rep(0.1, p)
s.cand <- 0.04

calculate_posterior_log(y,X,true_theta,true_sigma*diag(length(y)),theta_0,Sigma)

calculate_inverse_likelihood(y,X,true_theta, true_sigma)

q_mix_theta(y,X,true_theta,true_sigma,theta_0,Sigma)

# Draw samples from the combined posterior
theta_samples <- theta_samples_from_q_mix_theta (y, X, num_samples, initial_theta, Sigma, true_sigma, s.cand)
#autocorr.plot(theta_samples)

dim(theta_samples)
length(unique(theta_samples[,1])) / nrow(theta_samples)
# Assuming you have theta_samples matrix with num_samples rows and num_parameters columns

# Example usage:
# Assuming you have y, X, theta_samples, and true_sigma defined
vector_1 <- calculate_weighted_likelihoods(y, X, theta_samples, true_sigma)

dim(vector_1)

num_observations <- length(y)
vec_1_sum <-  colMeans(vector_1)
dim(as.matrix(vec_1_sum))


vector_2 <- calculate_vector_2(vec_1_sum)
dim(vector_2)

# Calculate the desired quantity (predictive distribution) for each data point
mix = ((vector_2))

##################################################################################
##################################################################################
##################################################################################


# Compute p(y_i|y_(-i))
compute_conditional_pdf_Y_i_given_Y_minus_i <- function(Y, X, Sigma_0, mu_0, sigma_0) 
{
  n <- length(Y)  # Number of observations
  pdf_values <- numeric(n)  # Initialize a vector to store PDF values
  
  density <- numeric()  
  
  for (i in 1:n) 
  {
    # Extract the i-th observation
    Y_i <- Y[i]
    X_i <- X[i, ]
    
    # Exclude the i-th observation
    Y_minus_i <- Y[-i]
    X_minus_i <- X[-i, ]  # Corrected X for X_{-i}
    
    
    # Calculate Lambda_minus_i
    Lambda_minus_i <- solve((t(X_minus_i) %*%  X_minus_i)/(sigma_0^2) + solve( Sigma_0))
    
    # Calculate eta_minus_i
    eta_minus_i <- solve((t(X_minus_i) %*% X_minus_i)/(sigma_0^2) + solve(Sigma_0)) %*%
      (solve(Sigma_0) %*% mu_0 + (t(X_minus_i)  %*% as.matrix(Y_minus_i))/(sigma_0^2))
    
    # Calculate mean and variance of the conditional distribution
    mean_Y_i_given_Y_minus_i <- t(X_i) %*% eta_minus_i
    variance_Y_i_given_Y_minus_i <- sigma_0^2  + t(X_i) %*% Lambda_minus_i %*% (X_i)
    
    # Compute the PDF value for Y_i
    #pdf_values[i] <- dnorm(Y[i], mean = mean_Y_i_given_Y_minus_i, sd = sqrt(variance_Y_i_given_Y_minus_i),log= TRUE)
    density[i] <- dnorm(Y[i], mean = mean_Y_i_given_Y_minus_i, sd = sqrt(variance_Y_i_given_Y_minus_i))
  }
  
  return(density)
}


test <- compute_conditional_pdf_Y_i_given_Y_minus_i(y,X,Sigma,theta_0,true_sigma)
original = (log(test))

# Calculate the Mean Squared Error (MSE)
mse <- mean((mix - original)^2)
print(paste("MSE:", mse))
