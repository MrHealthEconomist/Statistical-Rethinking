# ========================================================================================
# Libraries ---------------------------------------------------------------
# ========================================================================================

library(rethinking)
library(doParallel)

data("Howell1")
Data <- Howell1
str(Data)

# Envrionmental options ---------------------------------------------------
clusteredCores <- makeCluster(3)

registerDoParallel(clusteredCores)
getDoParWorkers()

# ========================================================================================
# Polynomial Regression ---------------------------------------------------
# ========================================================================================
# Let's plot the relationship between weight and height, without filtering any values:
plot(Data$height, Data$weight)
# ... Now that the data isn't filtered, unlike the previous section, the full relationship is 
# clearly curved, or *polynomial*, in nature. 

# In this context, polynomial means equations for mu[i] that add additional terms with squares, 
# cubes, and even higher powers of the predictor variable. Although there is still only one
# predictor variable in the model, so this is still a bivariate regressio. But a definition of 
# mu[i] has more parameters now. However, as always, polynomial regression is in general not a
# best practice technique. This is because polynomials are very hard to interpret. Nonetheless, 
# the following is the most common of polynomial regression equations, a parabolic model of the
# mean:

# mu[i] = alpha + beta1 * x[i] + beta2 * x[i] ^ 2

# This is a parabolic (second order) polynomial. The alpha + beta1 * x[i] part is the same 
# linear function of x in a linear regression. The additional term uses the square of x[i] to 
# construct a parabola, rather than a straight line. The new parameter beta2 measures the
# curvature of the relationship. 

# When fitting a polynomial model, it is good practice to standardise the data. Standardising 
# data leaves the mean at zero and additionally *rescales* the range of the data. This can help
# to make interpretation easier. For a standardised variable, a change of one unit is equivalent
# to a change of one sd/sigma. However, it can also be the converse, making interpretation 
# harder compared to on a natural scale. Lastly, standardising allows software to more easily 
# compute values of higher powers. 

# In order to standardise weight, all we have to do is subtract the mean and then divide by the
# standard deviation:
Data$weightStandardised <- (Data$weight - mean(Data$weight)) / sd(Data$weight)

# Let's plot height and weightStandardised to ensure no information has been lost in the 
# standardisation process:
plot(Data$height, Data$weightStandardised)
# ... the only change is on the horizontal axes...

# To fit the parabolic model, just modify the definition of mu[i]. Here is the model with very 
# weak priors:

#  h[i] ~ Normal(mu[i], sigma)
#  mu[i] = alpha + beta1 * x[i] + beta2 * x[i] ^ 2
#  alpha ~ Normal(178, 100)
#  beta1 ~ Normal(0, 10)
#  beta2 ~ Normal(0, 10)
#  sigma ~ Uniform(0, 50)

# Let's preprocess the quadratic term:
Data$weightStandardisedSqr <- Data$weightStandardised ^ 2
heightModel_5 <- map(
 alist(
  height ~ dnorm(mu, sigma),
  mu <- a + b1 * weightStandardised + b2 * weightStandardisedSqr,
  a ~ dnorm(178, 100),
  b1 ~ dnorm(0, 10),
  b2 ~ dnorm(0, 10),
  sigma ~ dunif(0, 50)
 ),
 data = Data
)
precis(heightModel_5)
# The parameter alpha is still the intercept, so it tells us the expected value of height when
# weightStandardised is zero. However, it is no longer equal to the mean height in the sample, 
# since there is no guarantee it should be a polynomial regression. And those beta1 and beta2 
# parameters are the linear and square components of the curve. However, again, this does not
# make them transparent. 

# We have to plot these model fits in order to understand what they are saying. We will 
# calculate the mean relationship and the 89% intervals of the mean and the predictions:
weightSeq <- seq(from = -2.2, 2, length.out = 30)

predData <- list(weightStandardised = weightSeq, weightStandardisedSqr = weightSeq ^ 2)
mu <- link(heightModel_5, data = predData)

muMean <- apply(mu, 2, mean)
muPI <- apply(mu, 2, PI)
simHeight <- sim(heightModel_5, data = predData)
heightPI <- apply(simHeight, 2, PI, prob = .89)

plot(Data$height ~ Data$weightStandardised, col = col.alpha("blue", alpha = .5))
lines(weightSeq, muMean)
shade(muPI, weightSeq)
shade(heightPI, weightSeq)

# Let's fit a higher order polynomial of the form:

#   h[i] ~ Normal(mu, sigma)
#   u[i] = alpha + beta1 * x[i] + beta2 * x[i] ^ 2 + beta3 * x[i] ^ 3

# Now we follow the same procedures as the previous model, just with an extra term:
Data$weightStandardisedCubed <- Data$weightStandardised ^ 3

# Fit the model:
heightModel_6 <- map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b1 * weightStandardised + b2 * weightStandardisedSqr + 
      b3 * weightStandardisedCubed,
    a ~ dnorm(178, 100),
    b1 ~ dnorm(0, 10),
    b2 ~ dnorm(0, 10),
    b3 ~ dnorm(0, 10),
    sigma ~ dunif(0, 50)
  ),
  data = Data
)

# Plot the model:
plot(Data$height, Data$weightStandardisedCubed)

weightSeq <- seq(from = -10.5, to = 6, length.out = 30)

predData <- list(weightStandardised = weightSeq, weightStandardisedSqr = weightSeq ^ 2, 
                 weightStandardisedCubed = weightSeq ^ 3)

mu <- link(heightModel_6, data = predData)

muMean <- apply(mu, 2, mean)
muPI <- apply(mu, 2, PI, prob = .89)
simHeight <- sim(heightModel_6, data = predData)
heightPI <- apply(simHeight, 2, PI, prob = .89)

plot(Data$height ~ Data$weightStandardised, 
     col = col.alpha("blue", alpha = .5),
     xlab = "Weight (Standardised)",
     ylab = "Height",
     main = "Polynomial fit of Height vs. z-Weight")
lines(weightSeq, muMean)
shade(muPI, weightSeq, col = col.alpha("red", alpha = .4))
shade(heightPI, weightSeq, col = col.alpha("grey", alpha = .4))
# ... this cubic curve is even more flexible than the parabola, so it fits the data even better.
# However, it is not clear whether these models make a lot of sense. Although they might be good
# geocentric fits, there are other models that could be better, particularly for interpretation.

# What happens if you want to convert the standardised z-score data back to the original scale 
# of the data. All that is really needed is to first turn off the horizontal axis when you plot
# the raw data:
plot(height ~ weightStandardised, data = Data, xaxt = "n", 
     col = col.alpha("skyblue", alpha = .6))
# ... the xaxt argument turns off the horizontal axis. Then you explicitly construct the axis, 
# using the axis function:
at <- c(-2, -1, 0, 1, 2)
labels <- at * sd(Data$weight) + mean(Data$weight)
axis(side = 1, at = at, labels = round(labels, 1))
# The first line above defines the location of the labels, in standardised units. The second 
# line then takes those units and converts them back to the original scale. The third line 
# draws the axis.

# End file ----------------------------------------------------------------
stopCluster(clusteredCores)