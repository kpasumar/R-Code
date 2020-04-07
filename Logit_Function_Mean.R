library(statmod) # Load this package for using gauss.quad.prob() function
library(gtools) # Load this package for using the logit function

#' Following function calculates the mean of logit-normal distribution
#' 
#' @param mu Mean of underlying normal distribution
#' @param sigma Standard deviation of underlying normal distribution
#' 
#' @return A single number representing mean of the logit-normal distribution
#' 
#' @details The calculation of w's and g(x)'s is vectorized
logitnorm_mean <- function(mu,sigma){
  v = 1/(1+ exp(-mu))           # Value passed into both shape parameters
  alpha_1 = 1/(sigma^2 * (1-v)) # Shape parameter 1
  alpha_2 = 1/(v * sigma^2) # Shape parameter 2
  # Calculate nodes and weights for Gaussian quadrature
  gqp <- gauss.quad.prob(n = 10,dist = "beta",alpha = alpha_1,beta = alpha_2)
  x <- gqp$nodes   # Extract the nodes into a vector
  w <- gqp$weights # Similarly the weights
  # Apply the function g (defined in the project description) onto the above x's
  g <- dnorm(logit(x),mean = mu,sd = sigma,log = TRUE) - log(1-x) - 
       dbeta(x,shape1 = alpha_1,shape2 = alpha_2,log = TRUE)
  # Calculate and return the mean
  answer <- sum(w*exp(g)) 
  return(answer)
}

# For testing
mu <- c(0.7,3.2,-1.1)
sigma <- c(0.8,0.1,2.3)
sapply(1:3, function(i) logitnorm_mean(mu[i],sigma[i]))


