###### Bronze estimator

# Function to calculate the posterior probability p(theta|y)
q_theta <- function(y, X, true_theta, sigma, theta_0, Sigma) 
{
  likelihoods <- numeric(length(y))
  n = length(y)
  
  for (i in 1:length(y)) 
  {
    likelihoods[i] <- ( dnorm( y[i], mean = X[i,] %*% true_theta, sd = sqrt(sigma)))
  }
  
  likelihood <- (prod(likelihoods))^((n-1)/n)
  prior <- dmvnorm( true_theta, mean = theta_0, sigma = Sigma)
  
  posterior <- likelihood * prior
  return(posterior)
}


# Function to draw samples from the combined posterior using Metropolis-Hastings
theta_samples_from_q_theta_bronze <- function(y, X, num_samples, theta_0, Sigma, sigma, s.cand) 
{
  num_parameters <- length(theta_0)
  samples <- matrix(0, nrow = num_samples, ncol = num_parameters)
  current_theta <- theta_0
  
  for (s in 1:num_samples) 
  {
    proposal_theta <- current_theta + rnorm(num_parameters, 0, sqrt(s.cand))
    
    # Calculate the Metropolis-Hastings acceptance ratio in log space
    q_theta_current <-  q_theta(y, X, as.numeric(current_theta), sigma, theta_0, Sigma)
    q_theta_proposal <- q_theta(y, X, as.numeric(proposal_theta), sigma, theta_0, Sigma)
    
    
    ratio <- min(1, q_theta_proposal / q_theta_current)
    
    
    # Print the log ratio
   # cat("Step:", s, "Ratio:", ratio, "\n")
    
    if ((runif(1)) < (ratio) )
    {
      current_theta <- proposal_theta
    }
    
    samples[s, ] <- current_theta
  }
  
  return(samples)
}


y_hj <- sapply(1:length(y), function(k) rnorm(1, mean = as.vector(X %*% theta_samples[k, ]), sqrt(true_sigma)))


# Function to calculate weights for each data point
calculate_wi_bronze <- function(y, X, theta_samples, sigma) {
 
  likelihoods <- numeric(length(y))
  n = length(y)
  for (s in 1:nrow(theta_samples)) 
  {
    likelihoods = (dnorm( y, mean = as.numeric(X %*% (theta_samples[s,])), sd = sqrt(true_sigma)))
      #for (i in 1:length(y)) 
      #{
      # likelihoods[i] <-  (dnorm( y[i], mean = as.numeric(X[i,] %*% theta_samples[s,]), sd = sqrt(sigma)))
      #} 
  }
  likelihood <- (prod(likelihoods))

  prior <- dmvnorm( true_theta, mean = theta_0, sigma = Sigma)
  
  posterior <- likelihoods * prior
  
  weights <- log(posterior/q_theta(y,X,true_theta,true_sigma,theta_0,Sigma))
  return(weights)
}


calculate_weighted_likelihoods_bronze <- function(y, X, theta_samples, sigma) 
{
  num_samples <- nrow(theta_samples)
  num_observations <- length(y)
  
  vector <- numeric()
  
  weighted_likelihoods <- matrix(0,num_samples,num_observations)
  weights <- matrix(0,num_samples,num_observations)
  
  for (s in 1:num_samples) 
  {
    weights[s,] <- calculate_wi_bronze(y,X,theta_samples[s,],true_sigma)
    weighted_likelihoods[s,] <- -(dnorm((y_hj), mean = as.numeric(X %*% (theta_samples[s, ])), sd = sqrt(true_sigma),log = TRUE)) %*% (weights[s,])
  }
  
  # Store the weighted_likelihoods for this observation in the matrix
  vector <- cbind(weighted_likelihoods,(weights))
  
  return(vector)
}

calculate_vector_2_bronze <- function(x)
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
set.seed(1)
N <- 50
X <- matrix(rnorm(N * p), ncol = p)
true_theta <- (rep(0.2, p))
true_sigma <- 2
y <- rnorm(N, as.vector(X %*% matrix(true_theta, ncol = 1)), true_sigma)
theta_0 <- as.matrix(rep(0, p))
Sigma <- diag(p)
num_samples <- 1e2
initial_theta <- rep(0.1, p)
s.cand <- 0.008

q_theta(y,X,true_theta,true_sigma,theta_0,Sigma)

# Draw samples from the combined posterior
theta_samples_bronze <- theta_samples_from_q_theta_bronze (y, X, num_samples, initial_theta, Sigma, true_sigma, s.cand)
dim(theta_samples_bronze)
length(unique(theta_samples[,1])) / nrow(theta_samples)

#y_hj_rep = calculate_y_hj_rep(y,X,theta_samples,true_sigma*diag(length(y)))

(calculate_wi_bronze(y,X,theta_samples_bronze,true_sigma))
vec_1 <- calculate_weighted_likelihoods_bronze(y, X, theta_samples_bronze, true_sigma)

vec_1_bronze <-  colMeans(vec_1)
dim(as.matrix(vec_1_bronze))


vec_2 <- calculate_vector_2(vec_1_bronze)
dim(vec_2)

bronze_estimator = mean((vec_2 - original)^2)




calculate_mse_for_p_bronze <- function( num_samples, p_values) {
  set.seed(10)
  mse_results <- numeric(length(p_values))
  for (i in 1:length(p_values)) {
  N <- 50
  p <- p_values[i]
  X <- matrix(rnorm(N * p), ncol = p)
  true_theta <- rep(0.2, p)
  true_sigma <- 2
  y <- rnorm(N, as.vector(X %*% matrix(true_theta, ncol = 1)), true_sigma)
  theta_0 <- as.matrix(rep(2, p))
  Sigma <- diag(p)
  initial_theta <- rep(0.1, p)
  s.cand <- 0.04

  
  
    num_samples_p <- num_samples  # You can adjust this if needed
    
    # Calculate q_theta and draw samples
    q_theta_values <- q_theta(y, X, true_theta, true_sigma, theta_0, Sigma)
    theta_samples <- theta_samples_from_q_theta_bronze(y, X, num_samples_p, initial_theta, Sigma, true_sigma, s.cand)
    
    # Calculate weighted likelihoods and vector 2
    vec_1 <- calculate_weighted_likelihoods_bronze(y, X, theta_samples, true_sigma)
    vec_1_mean <- colMeans(vec_1)
    vec_2 <- calculate_vector_2_bronze(vec_1_mean)
    
    # Calculate the desired quantity (predictive distribution) for each data point
    bronze = vec_2
    test <- compute_conditional_pdf_Y_i_given_Y_minus_i(y,X,Sigma,theta_0,true_sigma)
    original = (log(test))
    # Calculate the Mean Squared Error (MSE)
    mse_results[i] <- mean((bronze - original)^2)
  }
  
  return(data.frame(p = p_values, mse = mse_results))
}

# Example usage:
p_values <- c(5,10, 15, 20, 25, 30,35,40)  # Define the values of p for which you want to calculate MSE
mse_results_bronze <- calculate_mse_for_p_bronze(num_samples= 1e4, p_values)

# Example usage
p_values <- c(5,10,15,20,25,30)  # Add the desired values of "p" here

# Assuming df is your data frame with "p" and "mse" columns
ggplot(mse_results_mix, aes(x = p, y = mse)) +
  geom_line() +
  geom_point() +
  labs(x = "p ", y = "MSE") +
  ggtitle("MSE vs. p for bronze ")

final_df <- data.frame(p = p_values,bronze = mse_results_mix$mse,mix = mse_results_bronze$mse)

# Create a plot with both "mix" and "bronze" values
combined_plot <- ggplot(final_df, aes(x = p)) +
  geom_line(aes(y = mix, color = "Mix"), linetype = "solid") +
  geom_line(aes(y = bronze, color = "Bronze"), linetype = "dashed") +
  labs(x = "p Values", y = "MSE") +
  scale_color_manual(values = c("Mix" = "blue", "Bronze" = "red")) +
  ggtitle("MSE vs. p Values (Mix vs. Bronze)")

# Display the combined plot
print(combined_plot)

calculate_mse_for_N <- function(N_values, p) {
  set.seed(10)
  mse_results <- numeric(length(N_values))
  for (i in 1:length(N_values)) {
    N <- 50
    X <- matrix(rnorm(N * p), ncol = p)
    true_theta <- rep(0.2, p)
    true_sigma <- 2
    y <- rnorm(N, as.vector(X %*% matrix(true_theta, ncol = 1)), true_sigma)
    theta_0 <- as.matrix(rep(2, p))
    Sigma <- diag(p)
    initial_theta <- rep(0.1, p)
    
    num_samples_p <- num_samples  # You can adjust this if needed
    
    # Calculate q_theta and draw samples
    q_theta_values <- q_theta(y, X, true_theta, true_sigma, theta_0, Sigma)
    theta_samples <- theta_samples_from_q_theta_bronze(y, X, N_values[i], initial_theta, Sigma, true_sigma, s.cand)
    
    # Calculate weighted likelihoods and vector 2
    vec_1 <- calculate_weighted_likelihoods_bronze(y, X, theta_samples, true_sigma)
    vec_1_mean <- colMeans(vec_1)
    vec_2 <- calculate_vector_2_bronze(vec_1_mean)
    
    # Calculate the desired quantity (predictive distribution) for each data point
    mix = vec_2
    
    test <- compute_conditional_pdf_Y_i_given_Y_minus_i(y, X, Sigma, theta_0, true_sigma)
    original <- log(test)
    
    # Calculate the Mean Squared Error (MSE)
    mse_results[i] <- mean((mix - original)^2)
  }
  
  return(data.frame(N = N_values, mse = mse_results))
}

# Example usage:
N_values <- c(1e2,1e3,5000,1e4)  # Define the values of N for which you want to calculate MSE
p <- 25  # Fixed value of p
mse_results_df <- calculate_mse_for_N(N_values, p)
mse_results_df_bronze_diff_N = mse_results_df
# Assuming df is your data frame with "p" and "mse" columns
bronze_plot <- ggplot(mse_results_df_bronze_diff_N, aes(x = N, y = mse)) +
  geom_line() +
  geom_point() +
  labs(x = "N ", y = "MSE") +
  ggtitle("MSE vs. N for bronze for different N ")

final_df_1 <- data.frame(N = N_values,bronze = mse_results_df_bronze_diff_N$mse,mix = mse_result_df_mix_diff_N$mse)

# Create a plot with both "mix" and "bronze" values
combined_plot_1 <- ggplot(final_df_1, aes(x = N)) +
  geom_line(aes(y = mix, color = "Mix"), linetype = "solid") +
  geom_line(aes(y = bronze, color = "Bronze"), linetype = "dashed") +
  labs(x = "N Values", y = "MSE") +
  scale_y_continuous(limits = c(0,15))+
  scale_color_manual(values = c("Mix" = "blue", "Bronze" = "red")) +
  ggtitle("MSE vs. N Values (Mix vs. Bronze)")
combined_plot+combined_plot_1
