# ========================================================================================
# Grid Approximation ------------------------------------------------------
# ========================================================================================
# Although grid approximation is inefficient at scale, it serves as an important pedagogical tool
# in understanding Bayesian updating more deeply. In the context of a simple coin toss example, 
# Grid Approximation works well. Let's build a grid approximation for the model:

# 1. Define the grid. You decide how many points to use in estimating the posterior, and then
#    you male a list of the parameter values on the grid.
# 2. Compute the value of the prior at each parameter value on the grid.
# 3. Compute the likelihood for each parameter value.
# 4. Compute the unstandardised posterior at each parameter value, by multiplying the prior by
#    the likelihood.
# 5. Lastly, standardise the posterior, by dividing each value by the sum of all values.

# Define the grid:
pGrid <- seq(from = 0, to = 1, length.out = 20)

# Define prior:
Prior <- rep(1, 20)

# Compute Likelihood at each value in grid:
Likelihood <- dbinom(x = 6, size = 9, prob = pGrid)

# Compute product of likelihood and prior:
unstdPosterior <- Likelihood * Prior

# Standardise the Posterior, making all values sum to 1:
Posterior <- unstdPosterior / sum(unstdPosterior)

# The above code makes a grid of only 20 points. Now visualise the posterior distribution:
plot(pGrid, Posterior, type = "b", 
     xlab = "Probability of Water", ylab = "Posterior Probability")
mtext("20 Points", cex = 2)

# The correct density for the the grid is determined by how accurate you want the approximation
# to be. More points means more precision.


# Now try different priors ------------------------------------------------
Prior <- ifelse(pGrid < .5, 0, 1)

# Compute Likelihood at each value in grid:
Likelihood <- dbinom(x = 6, size = 9, prob = pGrid)

# Compute product of likelihood and prior:
unstdPosterior <- Likelihood * Prior

# Standardise the Posterior, making all values sum to 1:
Posterior <- unstdPosterior / sum(unstdPosterior)

# The above code makes a grid of only 20 points. Now visualise the posterior distribution:
plot(pGrid, Posterior, type = "b", 
     xlab = "Probability of Water", ylab = "Posterior Probability")
mtext("20 Points", cex = 2)

Prior <- exp(- 5 * abs(pGrid - .5))

# Compute Likelihood at each value in grid:
Likelihood <- dbinom(x = 6, size = 9, prob = pGrid)

# Compute product of likelihood and prior:
unstdPosterior <- Likelihood * Prior

# Standardise the Posterior, making all values sum to 1:
Posterior <- unstdPosterior / sum(unstdPosterior)

# The above code makes a grid of only 20 points. Now visualise the posterior distribution:
plot(pGrid, Posterior, type = "b", 
     xlab = "Probability of Water", ylab = "Posterior Probability")
mtext("20 Points", cex = 2)

# ========================================================================================
# Quadratic Approximation -------------------------------------------------
# ========================================================================================
# A Guassian distribution is convenient because it can be completely described by only two 
# numbers: the location of its centre (mean), and its general spread (variance). A Guassian
# approximation is called "quadratic approximation" because the logoarithm of a Guassian 
# distribution forms a parabola. And a parabola is a quadratic function. So this approximation
# essentially represents any log-posterior with a parabola. 

# The procedure in R contains two steps:

# 1. Find the posterior mode. Usually accomplished by some optimisation algorithm, a procedure
#    that virtually 'climbs' the posterior distribution.
# 2. Once you find the peal of the posterior, you must estimate the curvature near the peak. 

# This curvature is sufficient to compute a quadratic equation of the entire posterior 
# distribution. In some cases, these calculations can be done analytically. 

library(rethinking)
options(mc.cores = detectCores())

# The rethinking package provides a function to compute the above, map(). MAP stands for maximum
# a posteriori, which is just a fancy name for the mode of the posterior distribution. Thus, to
# compute the quadratic approximation to the coin toss example. To use map(), you provide a 
# formula, a list of data, and a list of starting values for the parameters. The formula defines
# the likelihood and the prior:
globQuad <- map(
 alist(
  w ~ dbinom(9, p), # binomial likelihood
  p ~ dunif(0, 1) # uniform prior
 ),
 data = list(w = 6)
)
# Display summary of quadratic equation:
precis(globQuad)
# ... one can read this approximation like: assuming the posterior us Gaussian, it is maximised
# at .67, and its sigma is .16

# Since we already know the posterior, let's compare to see how good the approximation is. We
# can use the analytical approach, using dbeta function:
w <- 6
n <- 9
curve(dbeta(x, w + 1, n - w + 1), from = 0, to = 1, col = "blue")
# Quadratic approximation:
curve(dnorm(x, .67, .16), lty = "solid", add = TRUE)

# As the number of samples/observations increase, however, the quaudratic approximation gets 
# better:
globQuad <- map(
 alist(
  w ~ dbinom(36, p), # binomial likelihood
  p ~ dunif(0, 1) # uniform prior
 ),
 data = list(w = 6)
)
precis(globQuad)

w <- 6
n <- 36
curve(dbeta(x, w + 1, n - w + 1), from = 0, to = 1, col = "blue")

# ... adding quadratic approximation to plot with n = 36:
curve(dnorm(x, .17, .06), lty = "solid", add = TRUE)

# However, it is important to note that the rate of improvement in approximation as sample size
# increases varies greatly depending on the details of the model.