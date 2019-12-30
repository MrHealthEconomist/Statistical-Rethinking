# ========================================================================================
# Libraries ---------------------------------------------------------------
# ========================================================================================
library(rethinking)
library(tidyverse)

options(mc.cores = detectCores())

# ========================================================================================
# Sampling from a Grid-Approximate Posterior ------------------------------
# ========================================================================================
# Let's regenerate the coin tossing example.

# Probability grid:
pGrid <- seq(from = 0, to = 1, length.out = 1000)

# Uniform prior:
Prior <- rep(1, 1000)

# Likelihood:
Likelihood <- dbinom(x = 6, size = 9, prob = pGrid)

# Unstandardised posterior:
UnstdPosterior <- Likelihood * Prior

# Normalised posterior:
Posterior <- UnstdPosterior / sum(UnstdPosterior)

# Now we wish to draw 10000 samples from the posterior. Imagine the posterior is a bucket full of
# parameter values, with numbers such as .1, .7, .5, 1, etcetera. Within the bucket, each value
# exists in proportion to its posterior probability, such that values near the peak are much more
# common than those in the tails. Thus, the sampling can be thought of as scooping out 10000 
# parameter values from the bucket. Provided that the 'bucket' is well mixed, the resulting 
# samples will have the same proportions as the exact posterior density. Therefore, the 
# individual values of *p* will appear in our samples in proportion to the posterior 
# plausibility/probability of each value. Remember from Hastie et al. that this is how bagging 
# is similar with Bayesian sampling methods; hence replacement.  

# Here's how we can do it:
samples <- sample(pGrid, prob = Posterior, size = 1e4, replace = TRUE)
# Note that 1e4 is 10000

# We can plot the results:
plot(samples, 
     col =  alpha(colour = "skyblue", alpha = .5), 
     pch = 16)

plot(density(samples))
# The author also provides a function for plotting the density estiamtes of the sampled 
# distribution in the rethinking package:
dens(samples, 
     col = "blue", lwd = 3)

# What this has done so far is crudely replicate the posterior density. Below will be to use 
# samples to describe and understand the posterior.

# ========================================================================================
# Sampling to Summarise ---------------------------------------------------
# ========================================================================================
# Summarising Posterior sampling outputs can be thought of as comprising three main themes, that
# encompass questions about:

# 1. intervals of defined boundaries;
# 2. intervals of defined probability mass; and
# 3. point estimates

# Intervals of defined boundaries -----------------------------------------
# For example, suppose you are asked for the posterior probability that the proportion of heads
# is less than .5. Using the grid-approximate posterior, you can just add up all the 
# probabilities, where the corresponding parameter value is less than .5. 

# For example, add up posterior probability where p < .5:
sum(Posterior[pGrid < .5])
# ... so about 17% of the posterior probability is below .5. However, since grid approximation 
# is not practical in general, it is not always so easy. Once there is more than one parameter
# in the posterior ditribution, this some becomes to analytically complex and cumbersome. 
# Therefore, let's see how to perform the same calculation, using the samples from the posterior.
# This approach, unlike the previous one above, does generalise to complex models with theta > 1.
# All you have to do is similarly add up all of the samples below .5, but then also divide the
# resulting count by the total number of samples (n). Thus, this finds the 'frequency' (remember:
# 10000) of parameter values below the .5 threshold:
sum(samples < .5) / 1e4
# ... and that's nearly the same answer as the grid approximation approach; although the answer
# will not be exactly the same, because the exact samples drawn from the posterior will be 
# different. 

# Using the same approach, you can also ask how much posterior probability lies between .5 and
# .75:
sum(samples > .5 & samples < .75) / 1e4
# ... so about 61% percent of the posterior probability lies between .5 & .75. 

# Intervals of defined mass -----------------------------------------------
# More commonly refered to as 'credibility intervals'.

# Suppose for example, you want to know the boundaries of the lower 80% posterior probability. 
# You know the interval starts at p = 0. To find out where it stops, think of the samples as
# data and ask where the 80th percentile lies:
quantile(samples, probs = .8)
# Similarly, the middle 80% lies between the 10th and 90th percentiles:
quantile(samples, probs = c(.1, .9))

# These intervals do a good job of communicating the shape of a distribution, as long as the 
# distribution isn't too assymetrical. However, in terms of supporting inferences about which 
# parameters are consistent with the data. Consider a posterior consistent with observing three
# heads in three tosses and a uniform/flat prior. It is highly skewed, having its maximum value
# at the boundary, p = 1:
pGrid <- seq(from = 0, to = 1, length.out = 1000)
Prior <- rep(1, 1000)
Likelihood <- dbinom(x = 3, size = 3, prob = pGrid)
UnstdPosterior <- Likelihood * Prior
Posterior <- UnstdPosterior / sum(UnstdPosterior)
samples <- sample(pGrid, size = 1e4, replace = TRUE, prob = Posterior)
# The rethinking package provides a convenient function to compute CI's:
PI(samples, prob = .5)
# ... note that this interval assigns 25% of the probability mass above and below the interval.
# So it provides the central 50% probability. But in this example, it ends up excluding the most
# probable parameter values, near p = 1. Therefore, in terms of describing the shape of the 
# posterior distribution the percentile interval can be misleading.

# In contrast, the Highest Posterior Density Interval (HDPI) is the narrowest interval 
# containing the specified probability mass. Simply, it can be thought of in this way: if you 
# want an interval that best represents the parameter values most consistent with the data, then
# you want the densest of these intervals; that's what the HDPI is. Again, the rethinking 
# package provides a function for this:
HPDI(samples, prob = .5)
# ... the HDPI captures the parameters with the highest posterior probability, as well as being
# noticeably narrower: .16 in width compared to .23 for the percentile interval.

# Note however, that the two types of interval summaries are only different because the 
# distribution is skewed. The HDPI also has some disadvantages. It is more computationally 
# intensive than standard PI/CI and suffers from greater *simulation variance*, which means 
# that it is sensitive to how many samples you draw from the posterior. CI's are also easier to
# interpret for non-Bayesian practitioners. Thus, intervals are often best used to communicate
# the shape of a distribution; many non-bayesian commit the fallacy of viewing frequentist CI's
# as probabilties - they are not! 

# Overall, if the choice of interval type makes a big difference, then you shouldn't be using
# intervals to summarise the posterior. If choice of interval leads to different inferences,
# then you'd be better off plotting the entire posterior distribution. 

# Point estimates ---------------------------------------------------------
# The final common summary task for the posterior is to produce point estimates. Given the 
# entire posterior distribution, what value should we report. The Bayesian parameter estiamte 
# is precisely the entire posterior distribution, which is a function that maps each unique 
# parameter value onto a plausibility value. Thus, it's hardly ever necessary to choose a point
# estimate. However, if you must produce a point estimate, you'll have to ask and answer 
# questions regarding the objectives of the analysis and the characteristics of the posterior.

#  For example, consider the coing tossing example again, in which we observe three heads out of
# three tosses. Let's consider the alternative point estimates. First, its very common for 
# scientists to report the parameter value with the highest probability, called a *maximum a 
# posteriori* (MAP) estimate. One can easily compute the MAP in this example:
pGrid[which.max(Posterior)]
# ... or if you instead have samples from the posterior, you can still approximate the same 
# point:
chainmode(samples, adj = .01)
# But why is this point, the mode, interesting? Why not report the posterior mean or median?
mean(samples)
median(samples)

# Once principled way to go beyond using the entire posterior as the estimate is to choose a loss
# function. A loss function is a rule that tells you the cost assoicated with using any 
# particular point estimate. It is important to note that different loss functions imply 
# different point estimates. What a loss function imposes a penalisation proportional to the
# distance from the true value from the estimated value:
sum(Posterior * abs(.5 - pGrid))
# All the above code does is compute the weighted loss average, where each loss is weighted by 
# its corresponding posterior probability. There is a trick for repeating this calculation for
# every possible decision, using the sapply function:
lossFunction <- sapply(pGrid, function(d) sum(Posterior * abs(d - pGrid)))
# ... Now lossFunction contains a list of loss values, for each possible decision, correspinding 
# to the values in pGrid.

# From here, it's easy to find the parameter value that minimises the loss: 
pGrid[which.min(lossFunction)]
# ... and this is actually the posterior median, the parameter value that splits the posterior
# density such that half of the moss is above it and half below it. 

# The quadratic loss function (d - p)^2, is the mean loss function.

# ========================================================================================
# Sampling to simulate prediction -----------------------------------------
# ========================================================================================
# Generating implied observations from a model is useful for at least four distinct reasons:

# 1. Model checking - after a model is fit to real data, it is worth simulating implied 
#    observations, to check both whether the fit worked correclty and to investigate model
#    behaviour. 
# 2. Software validation - in order to be sure that our modeling fitting software is working,
#    it helps to simulate observations under a known model and then attempt o recover the 
#    values of the parameters the data were simulated.
# 3. Research design - if you can simulate observations from your hypothesis, then you can 
#    evaluate whether the research design can be effective. In a narrow sense, this means doing
#    *power analysis*, but the possibilities are much broader.
# 4. Forecasting - estimates can be used to simulate new predictions, for new cases and future
#    observations. These forecasts can be useful as applied prediction, but also for model 
#    criticism and revision. 

# Another important aspect regarding Bayesian modelling is that given a realised observation, 
# the likelihood function says how plausible the observation is. And given only the parameters, 
# the likelihood defines a distribution of possible observations that we can sample from, to 
# simulate observation. In this way, Bayesian modelling is generative. We can call such 
# simulated data 'Dummy Data'. With the coin tossing example, the dummy data arises from a 
# binomial likelihood function. For example, suppose we observe two tosses; there are only 
# three possible observations, which are 0 H, 1 H, and 2 H. Let's use p = .5:
dbinom(0:2, size = 2, prob = .7)

# Now we're going to simulate observations, using their likelihoods. This is done by sampling 
# from the distribution described above. You could use the sample function to do this, but R 
# provides convenient sampling functions for all ordinary probability distributions. So a single
# dummy observation of Heads can be sampled with:
rbinom(1, size = 2, prob = .7)
# ... the 'r' in rbinom stands for 'random'. It can also generate more than one simulation at a
# time. For example, a set of 10 simulations can be made by:
rbinom(10, size = 2, prob = .7)

# Let's generate 100000 dummy observations, just to verify that each value (0, 1, or 2) appears
# in proportion to its likelihood:
dummyHeads <- rbinom(1e5, size = 2, .7)
table(dummyHeads) / 1e5
# ... and those values are very close to the analytically calculated binomial likelihoods above. 
# The difference in values is due to simulation variance.

# Only two tosses of the coin is not much of a sample, however. So now, let's simulate the same
# sample size as before; 9 tosses:
dummyHeads <- rbinom(1e5, size = 9, prob = .7)
simplehist(dummyHeads, xlab = "Dummy Heads Count")
# Notice that most of the time, the expected observation does not contain heads in its true 
# proportion, i.e. .7. That is the nature of observation: there is a one-to-many relationship 
# between data and data-generating process.

# ========================================================================================
# Model checking ----------------------------------------------------------
# ========================================================================================
# Model checking means (1) ensuring the model fitting worked correctly and (2) evaluating the 
# adequacy of a model for some purpose. Since Bayesian models are *generative*, i.e. they are
# able to simulate observations as well as estimate parameters from observations, once you 
# condition a model on data, you can simulate to examine the model's empirical expectations. 
# One way to assess the implied distribution of outcomes for each value of p, we can use
# a posterior predictive distribution.

# To simulate predicted observations for a single value of p, e.g. p = .6, one can use the 
# rbinom to generate random binomial samples for the coin toss example: 
Heads <- rbinom(1e4, size = 9, prob = .6)
# ... the above generates 10000 simulated predictions of 9 coin tosses, assuming that p = .6. 
# The predictions are stored as counts of heads, so the theoretical minimum is zero, and the 
# theoretical maximum is 9. 
simplehist(Heads)
# All you need to propogate parameter uncertainty into these predictions is by replacing the 
# p = .6 value with samples from the posterior:
Heads <- rbinom(1e4, size = 9, prob = samples)
simplehist(Heads)
# ... since the sampled values appear in proportion to their posterior probabilities, the 
# resulting simulated observations are averaged over the posterior. 

# Note: it is important to remember that working with samples transforms a problem of integral 
# calculus into a problem of data summary. Additionally, to sum, posterior predictive checks 
# combine uncertainty about parameters with uncertainty about outcomes - as described by the 
# assumed likelihood function. 

# ========================================================================================
# Practice ----------------------------------------------------------------
# ========================================================================================
pGrid <- seq(from = 0, to = 1, length.out = 1000)
Prior <- rep(1, 1000)
Likelihood <- dbinom(6, size = 9, prob = pGrid)
UnstdPosterior <- Likelihood * Prior
Posterior <- UnstdPosterior / sum(UnstdPosterior)

set.seed(100)
samples <- sample(pGrid, prob = Posterior, size = 1e4, replace = TRUE)

# How much posterior probability lies below p = .2?
sum(samples < .2) / 1e4

# How much posterior probability lies above p = .8?
sum(samples > .8) / 1e4

# How much posterior probability lies between p = .2 & p = .8?
sum(samples > .2 & samples < .8) / 1e4

# 20% of the posterior lies below which value of p?
quantile(samples, probs = .2)

PI(samples, prob = .66)

HPDI(samples, prob = .66)


# Another exercise with 8 Heads out of 15 Tosses:
Likelihood <- dbinom(8, size = 15, prob = pGrid)

UnstdPosterior <- Likelihood * Prior
Posterior <- UnstdPosterior / sum(UnstdPosterior)

set.seed(100)
samples <- sample(pGrid, prob = Posterior, size = 1e5, replace = TRUE)

HPDI(samples, prob = .9)

Heads <- rbinom(1e5, size = 15, prob = samples)
simplehist(Heads)


Heads <- rbinom(1e4, size = 9, prob = .6)
# ... the above generates 10000 simulated predictions of 9 coin tosses, assuming that p = .6. 
# The predictions are stored as counts of heads, so the theoretical minimum is zero, and the 
# theoretical maximum is 9. 
simplehist(Heads)
# All you need to propogate parameter uncertainty into these predictions is by replacing the 
# p = .6 value with samples from the posterior:
Heads <- rbinom(1e4, size = 9, prob = samples)
simplehist(Heads)
table(Heads) / 1e5

# End file ----------------------------------------------------------------