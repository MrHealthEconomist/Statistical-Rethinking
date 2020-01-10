# ==========================================================================================
# Libraries ---------------------------------------------------------------
# ==========================================================================================
library(rethinking)
data("Howell1")
head(Howell1)
library(doParallel)


# Environment settings ----------------------------------------------
clusteredCores <- makeCluster(3)
registerDoParallel(clusteredCores)
getDoParWorkers()


# ==========================================================================================
# Exercise 4H2 ------------------------------------------------------------
# ==========================================================================================
Data <- Howell1[Howell1$age < 18, ]
dim(Data)

plot(height ~ weight, data = Data)

sd(Data$height)
summary(Data$height)

dens(Data$height)
# We'd like to say how knowing the values in x can help us describe or predict the values in h. 
# Thus, to get weight into the model in this way, we define the mean mu as a function of the 
# values in x.
set.seed(77)
LinearBayesModel <- map(
 alist(
  height ~ dnorm(mu, sigma), # likelihood
  
  mu <-  a + b * weight, # linear model
  
  # Median child height to be 111 with a sigma 
  # of 50 for variance:
  a ~ dnorm(120, 100), # alpha prior
  
  # Why have a prior with a mean zero (Beta)? This 
  # prior places just as much probability below zero 
  # as it does above zero:
  b ~ dnorm(0, 10), # beta prior
  
  # Sigma to be uniform for equal Prob Mass:
  sigma ~ dunif(0, 50) # sigma prior
 ),
 data = Data
)
precis(LinearBayesModel)

# Sample Posterior:
postSamples <- extract.samples(object = LinearBayesModel, n = 20)

# ... display raw data:
plot(Data$weight, Data$height,
     xlim = range(Data$weight), ylim = range(Data$height),
     col = col.alpha("navyblue", .6),
     pch = 16,
     xlab = "Weight", ylab = "Height")
# ... plot lines of uncertainty:
for (i in 1:20) {
 abline(a = postSamples$a[i], b = postSamples$b[i],
        col = col.alpha("grey", .4), lwd = 2)
}

# Plot of regression contours and intervals for linear model -----------------
# ... remember the interval incorporates the uncertainty in both the slope Beta and intercept 
# alpha.

# Define sequence of weights to predict for height:
weightSeq <- seq(from = 0, to = 50)

# Extract the intervals for the mu (mean) distribution:
mu <- link(LinearBayesModel, data = data.frame(weight = weightSeq))

muMean <- apply(mu, 2, mean)
muHDPI <- apply(mu, 2, HPDI, prob = .89)

plot(Data$weight, Data$height,
     xlim = range(Data$weight), ylim = range(Data$height),
     col = col.alpha("navyblue", .6),
     pch = 16,
     xlab = "Weight", ylab = "Height")

# MAP line, i.e. the average mu (mean) for each weight value or *distribution of the average 
# value*:
lines(weightSeq, muMean, lty = "dashed", lwd = 2)
# HPDI region of uncertainty for the average, mu:
shade(muHDPI, weightSeq, 
      col = col.alpha("darkgrey", alpha = .5))

# Prediction intervals ----------------------------------------------------
simHeight <- sim(LinearBayesModel, data = list(weight = weightSeq))
# check dimensions...
str(simHeight)

# Summarise the simulated height distributions:
heightPI <- apply(simHeight, 2, PI, prob = .89)

# Wrap up: the MAP line, the shaded region of the 89% posterior region of plausible mu, and 
# the boundaries of the simulated heights that the model expects (according to the assumptions
# inherent in the modelling process):
plot(Data$weight, Data$height,
     col = col.alpha("skyblue", alpha = .6), lwd = 3, pch = 16,
     xlab = "Weight", ylab = "Height", 
     main = "Bayesian Regression Model: Height vs Weight")
lines(weightSeq, muMean, lty = "solid", lwd = 3, col = "black")
shade(muHDPI, weightSeq, 
      col = col.alpha("navyblue", alpha = .5))
# PI region for our simulated, *predicted* heights:
shade(heightPI, weightSeq, col = col.alpha("darkgrey", alpha = .2))

# It seems a polynomial model of height vs weight will fit the data better, as height and weight 
# *visually* seem to share an exponential relationship...

# ==========================================================================================
# Exercise 4H3 ------------------------------------------------------------
# ==========================================================================================
# The task is to model the relationship between height and the natural logarithm of weight. This
# exercise uses the entire Howell1 data set. The fit of the model is the quadratic 
# approximation indicated by:

#   h[i] ~ N(mu[i], sigma),
#   mu[i] = alpha + betalog(w[i])
#   alpha ~ N(178, 100),
#   beta ~ N(0, 100),
#   sigma ~ U(0, 50)

# ... where h[i] is the height of the individual i and w[i] is the weight of individual i.

# Load & Transform Data ---------------------------------------------------
data("Howell1")
Data <- Howell1
dim(Data)
plot(Data$weight, Data$height, 
     col = col.alpha("skyblue", alpha = .6), pch = 16)

# Standardise weight:
Data$weightStnd <- (Data$weight - mean(Data$weight)) / sd(Data$weight)
# Inspect information loss:
plot(Data$weightStnd, Data$height, 
     col = col.alpha("skyblue", alpha = .6), pch = 16) # None...

# Quadratic form of weight:
Data$weightSqr <- (Data$weightStnd) ^ 2

# Model:
heightLogModel <- map(
  alist(
    height ~ dnorm(mu, sigma), # likelihood
    mu <- a + b * weightSqr, # linear log model
    a ~ dnorm(178, 100), # alpha prior
    b ~ dnorm(0, 100), # beta prior
    sigma ~ dunif(0, 50)
  ),
  data = Data
)
precis(heightLogModel, corr = TRUE)



# ... display raw data:
plot(Data$weightSqr, Data$height,
     xlim = range(Data$weightSqr), ylim = range(Data$height),
     col = col.alpha("skyblue", alpha = .7), pch = 16,
     xlab = "z-score Weight", ylab = "Height")
# ... plot simulated line fits:
for (i in 1:50) {
  abline(a = postSamples$a[i], b = postSamples$b[i],
         col = col.alpha("black", .1), lwd = 4)
}

# Generate predictions and visualise -----------------

# Define sequence of weight values:
weightSeq <- seq(from = -2.5, to = 2.5, length.out = 50)

# Generate data to predict relationship with height:
predData <- list(weightSqr = (weightSeq) ^ 2)
# Extract mu distribution intervals:
mu <- link(heightLogModel, data = predData)

muMean <- apply(mu, 2, mean)
muPI <- apply(mu, 2, PI, prob = .89)

simHeight <- sim(heightLogModel, data = predData)
heightPI <- apply(simHeight, 2, PI, prob = .89)

# Plot raw data:
plot(Data$weightStnd, Data$height, 
     col = col.alpha("navyblue", .6), pch = 16,
     xlab = "z-score Weight", ylab = "height",
     main = "Standardised Weight vs Height: a Linear Bayesian model")
lines(weightSeq, muMean, 
      col = col.alpha("black", alpha = .9), lwd = 2)
shade(muPI, weightSeq, 
      col = col.alpha("red", alpha = .1))
shade(heightPI, weightSeq,
      col = col.alpha("darkgrey", alpha = .4))

# So, it *visually* indicates that the model would be a better fit with an additional, higher 
# order (cubic) term...

# End file ----------------------------------------------------------------
stopCluster(cl = clusteredCores)