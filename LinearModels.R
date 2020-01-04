# ========================================================================================
# Libraries ---------------------------------------------------------------
# ========================================================================================
library(rethinking)
library(parallel)

options(mc.cores = detectCores())

# ========================================================================================
# Linear Models -----------------------------------------------------------
# ========================================================================================
# Under a probability interpretation, necessary for a Bayesian procedure, linear regression uses 
# a Guassian distribution to describe a model's uncertainty about some measurement of interest. 
pOutcome <- replicate(1000, sum(runif(16, -1, 1)))
hist(pOutcome)
plot(density(pOutcome))

# Normal by addition ------------------------------------------------------
# The above simulates random walks on a distribution according to a binomial outcome, which 
# converges to a normal distribution. The more 'steps' or draws sampled from the posterior, the 
# closer the distribution converges to the dnorm. 

dens(pOutcome, 
     col = "lightblue", lwd = 4)

# Normal by multiplication ------------------------------------------------
# Here's another way to get a normal distribution. Suppose the growth rate of an organism is 
# influenced by a dozen loci, each with several alleles that code for more growth. Suppose that 
# all of these loci interact with one another, such that each increase the growth of the 
# organism by a percentage. This means that their effects are multiplicative, rather than 
# additive. E.g. we can sample a random growth rate for this example:
prod(1 + runif(12, 0, .1))
# ... the above code samples 12 random numbers between 1.0 and 1.1, each representing a 
# proportional increase in growth. Thus 1.0 means no additional growth, and 1.1 means a 10% 
# increase in growth. The product of all 12 sampled numbers is then computed and returned as an
# output. What distribution will these random products produce? Let's generate 10000 of them and 
# see:
growth <- replicate(10000, prod(1 + runif(12, 0, .1)))
dens(growth, norm.comp = TRUE, 
     col = "skyblue", lwd = 3)

# Normal by log-multiplication --------------------------------------------
# Large deviates that are multiplied together do not produce Guassians, but they do tend to 
# produce Guassian distributions on the log-scale. For example:
logBig <- replicate(10000, log(prod(1 + runif(12, 0, .5))))
dens(logBig)

# ========================================================================================
# A Guassian model of height ----------------------------------------------
# ========================================================================================
# There will be two parameters describing the distributions shape, the mean and the sd. Bayesian
# updating will allow us to consider every possible combination of values for mu and sigma, and 
# to score each combination by its relative plausibility, in light of the data. These relative 
# plausibilities are the posterior probabilities for each combination of values (mu, sigma).

data("Howell1")
# The data contained in Howell1 are partial census data for the Dobe area !Kung San, 

Data <- Howell1
str(Data)
Data$height

# All we want for now is the heights of adults from the Data of 18 years or older:
Data_2 <- Data[Data$age >= 18, ]
str(Data_2)

# The model ---------------------------------------------------------------
# Our goal is to model these values using a Guassian distribution. First, let's plot the 
# distribution:
dens(Data_2$height)

# It's reasonable for the moment to adopt a stance that the model's likelihood should be a 
# Guassian. However, not that only gawking (viewing) the data to try and decide how to model them
# is usually not a good idea. 

# To complete the model, we are going to need some priors. The parameters to be estiamted are 
# both mu and sigma, so we need a prior Pr(mu, sigma), the joint probability for all the 
# parameters. In most cases, priors are specified independently for each parameter, which amounts
# suming Pr(mu, sigma) = Pr(mu)Pr(sigma). Then we can write:

#  h[i] ~ N(mu, sigma) # likelihood
#  mu ~ N(178, 20) # mu prior
#  sigma ~ Uniform(0, 50) # sigma prior

# The prior for mu is a broad Guassian prior, centered on 178cm, with a 95% probability between
# 178 +- 40. Why 178cm? The author is 178cm tall. And the range from 138cm to 218cm encompasses
# a large range of plausible mean heights for human populations. However, in many regression 
# problems, using prior information is typically more subtle, because parameters don't always 
# have such clear physical meaning. 

# Whatever the prior, it's a very good idea to plot your priors, so you can have a sense of the
# assumption they build into the model:
curve(dnorm(x, 178, 20), from = 100, to = 250)

# The sigma prior is a truly flat prior, a uniform one, that functions just to constrain sigma to
# have positive probability between zero and 50cm. 
curve(dnorm(x, 0, 50), from = -10, to = 60)
# ... s standard deviation ike sigma must be positive, so bounding it at zero makes sense. 

# You can simulate heights by sampling from the prior. Remember that every posterior is also 
# potentially a prior for a subsequent analysis, so you can process priors just like posteriors.

# Sample mu heights:
sampleMu <- rnorm(1e4, 178, 20)
# Sample height sd:
sampleSigma <- runif(1e4, 0, 50)
# Sample prior heights:
priorHeights <- rnorm(1e4, sampleMu, sampleSigma)

# Plot prior heights:
dens(priorHeights)
# ... the density plot shows a vaguely bell-shaped density with thick tails. It is the expected
# distribution of heights, averaged over the prior.

# Grid approximation of the posterior -------------------------------------
# The use of brute force calculations for the posterior is not, if at all, usually possible.
# Thus, this technique rather offers insight into what the 'target' it is Bayesian computation.
muList <- seq(from = 140, to = 160, length.out = 200)
sigmaList <- seq(from = 4, to = 9, length.out = 200)

Post <- expand.grid(mu = muList, sigma = sigmaList)

Post$LL <- sapply(1:nrow(Post), function(i) sum(dnorm(
 Data_2$height,
 mean = Post$mu[i],
 sd = Post$sigma[i], 
 log = TRUE
)))

Post$prod <- Post$LL + dnorm(Post$mu, 178, 20, log = TRUE) + 
 dunif(Post$sigma, 0, 50, log = TRUE)

Post$prob <- exp(Post$prod - max(Post$prod))

# You can inspect this posterior distribution, located in Post$prob, using a variety of plotting
# commands:
contour_xyz(Post$mu, Post$sigma, Post$prob)

image_xyz(Post$mu, Post$sigma, Post$prob)

# Sampling from the Posterior ---------------------------------------------
# Since there are two parameters, and we want to sample combinations of them, we first randomly
# sample row numbers in Post in proportion to the values in Post$prob. Then we pull out the 
# parameter values on those randomly sampled rows:
sampleRows <- sample(1:nrow(Post), size = 1e4, replace = TRUE, prob = Post$prob)

sampleMu <- Post$mu[sampleRows]
sampleSigma <- Post$sigma[sampleRows]

# You end up with 10000 samples, with replacement, from the posterior for the height data:
plot(sampleMu, sampleSigma, cex = .5, pch = 16, col = col.alpha(acol = "blue", .1))

# Now that you have these samples, you can describe the distribution of confidence in each 
# combination of mu and sigma by summarising samples. For example, to describe the marginal 
# posterior densities of mu and sigma, all we need is:
dens(sampleMu)
dens(sampleSigma)
# ... the jargon 'marginal' here means averaging over the other parameters. 

# To summarise the widths of these densities with highest posterior intervals:
HPDI(sampleMu)
HPDI(sampleSigma)

# Other summary statistics:
mean(sampleMu)
mean(sampleSigma)
median(sampleMu)
median(sampleSigma)

# Before moving on to using quadratic approximation, it is worth repeating the analysis of 
# height data above with only a fraction of the original data. This will demonstrate that, in
# principle, the posterior is not always so Guassian in shape. There's no trouble with the mean
# For a Guassian likelihood and a Guassian prior on mu, the posterior distribution is always 
# Guassian as well, regardless of sample size. It is the standard deviation (sigma) that causes
# problems.

# The reasons for the posterior of sigma tending to have a long right-tail are complex. But a 
# useful way to concieve of the problem is that variance must be positive. Thus, there must be 
# more uncertainty about how large the variance or sigma is than about how small it is. 

# Let's quickly analyse only 20 of the heights from the data to reveal this issue. To sample 20
# random heights from the original data:
Data_3 <- sample(Data_2$height, size = 20)

# Now we repeat the code from the previous section, modified to focus on 20 heights:
muList <- seq(from = 150, to = 170, length.out = 200)
sigmaList <- seq(from = 4, to = 20, length.out = 200)
post2 <- expand.grid(mu = muList, sigma = sigmaList)

post2$LL <- sapply(1:nrow(post2), function(i) 
 sum(dnorm(Data_3, mean = post2$mu[i], sd = post2$sigma[i], log = TRUE)))

post2$prod <- post2$LL + dnorm(post2$mu, 178, 20, log = TRUE) + 
 dunif(post2$sigma, 0, 50, log = TRUE)

post2$prob <- exp(post2$prod - max(post2$prod))

sampleRows_2 <- sample(1:nrow(post2), size = 1e4, prob = post2$prob, replace = TRUE)

sampleMu_2 <- post2$mu[sampleRows_2]
sampleSigma_2 <- post2$sigma[sampleRows_2]

plot(sampleMu_2, sampleSigma_2, cex = .5, 
     col = col.alpha(rangi2, .1),
     xlab = "mu", ylab = "sigma", pch = 16)
dens(sampleSigma_2, norm.comp = TRUE, 
     col = "skyblue", lwd = 4)

# Fitting the model with map ----------------------------------------------
# Our interest in quadratic approximation is that it is a handy way to quickly make inferences
# about the shape of the posterior. The posterior's peak will lie at the maximum a posteriori 
# estimate (MAP), and we can get a useful image of the posterior's shape by using the quadratic
# approximation of the posterior distribution at this peak. Think MLE.

# To find the values of mu and sigma that maximise the posterior probability, we can use map, a 
# command in the rethinking package. The way that map works is by using the model definition. 
# Each line in the definition has a corresponding definition in the form of R cod. The 'engine' 
# inside map then uses these definitions to define the posterior probability at each combination
# of parameter values. Then it can climb the posterior distribution and find the peak, i.e. the 
# MAP. It is thus declarative. 

#  h[i] ~ N(mu, sigma), in R = height ~ dnorm(mu, sigma)
#  mu ~ N(178, 20), in R = mu ~ dnorm(178, 20)
#  sigma ~ Uniform(0, 50), in R = sigma ~ dunif(0, 50)

# We have to place the code equivalents into a list:
modelString <- alist(
 height ~ dnorm(mu, sigma),
 mu ~ dnorm(178, 20),
 sigma ~ dunif(0, 50)
)
# ... note that alist() does not evaluate the embedded code you put in it, unlike the normal 
# list() function. Thus, when you define a list of formulas in R, you should use alist(), so
# the code isn't executed. However, when you define a list of start values for parameters, you
# should use list(), so that code like mean(Data_2$height) will be evaluated to a numeric value.

# Now we fit the model to the dataframe:
heightModel <- map(modelString, 
                   data = Data_2)

# Now we can look at the fit *maximum a posteriori* model:
precis(heightModel)

# These numbers provide Guassian approximations for each parameters *marginal* distribution. This
# means the plausibility of each value of mu, after averaging over the plausibilities of each 
# value of sigma, is given by  Guassian distribution with a mu of 154.61 and a sigma of .41.

# The priors we used before are very weak, both because they are nearly flat and because there 
# is so much data. Let's, therefore, see what a more informative prior does. We will provide a 
# more informative prior for mu. All we are going to do is change the standard deviation of the
# prior to .1, so that it is a very narrow prior:
heightModel_2 <- map(
 alist(
  height ~ dnorm(mu, sigma), # likelihood
  mu ~ dnorm(178, .1), # mu prior
  sigma ~ dunif(0, 50) # sigma prior
 ),
 data = Data_2
)
precis(heightModel_2)
# Notice that the estimate for mu has hardly moved off the prior. The prior was very 
# concentrated around 178. So this is not suprising. But also notice that the estimate for sigma
# has changed quite a lot, even though we didn't change its prior at all. Once the sampler is
# certain that the mean is near 178 - as the prior insists - then the sampler has to estimate 
# sigma conditional on that fact. This results in a different posterior for sigma, even though
# all we changed is prior information about the other parameter.  

# Sampling from a map fit -------------------------------------------------
# Just like a mu ans sigma are sufficient to describe a one-dimensional Guassian distribution,
# a list of means and a matrix of variances and covariances are sufficient to describe a multi-
# dimensional Guassian distribution. To see the matrix of variances and covariances for model 1:
vcov(heightModel)
# ... the above is thus a variance-covariance matrix. It is the multi-dimensional glue of a 
# quadratic approximation, because it tells us how each parameter relates to every other 
# parameter in the posterior. A variance-covariance matrix can be factored into two elements:

# 1) a vector of variances for the parameters; and
# 2) a correlation matrix that tells us how changes in any parameter lead to correlated chnages
#    in the others.

diag(vcov(heightModel))
cov2cor(vcov(heightModel))
# The two-element vector in the output is the list of variances. The two-by-two matrix in the
# output is the correlation matrix. Each entry shows the correlation, bounded between -1 & +1, 
# for each pair of parameters. 

# However, how do we get samples from this multi-dimensional posterior? We sample vectors of 
# values from a multi-dimensional Guassian distribution:
posteriorSample <- extract.samples(heightModel, n = 1e4)
head(posteriorSample)
# Each value is a sample from th posterior, so the mu and sigma for each column will be very 
# close to the MAP values from before. You can comfirm this by summarising:
precis(posteriorSample)
plot(posteriorSample, 
     col = col.alpha("skyblue", alpha = .3), pch = 16)

# Let's peak under the hood of the extract.samples function. The work is done by the multi-
# dimensional version of rnorm, mvrnorm, mvrnorm simulates random vectors of multivariate 
# Guassian values; unlike rnorm which simulates single, random values. Here is how to do what
# extract.samples() does in one step:
library(MASS)

posteriorSample <- mvrnorm(n = 1e4, mu = coef(heightModel), Sigma = vcov(heightModel))
head(posteriorSample)

# Food for thought: getting sigma right -----------------------------------
# The quadratic assumption for sigma can be problematic. A conventional way to improve the 
# situation is the estimate log(sigma) instead. Why does this help? The posterior distribution
# of sigma will often not be Guassian, and the distribution of its logarithm can be much closer
# to Guassian. Thus, if we impose the quadratic approximation on the logarithm, rather than the
# sigma itself, we can get a better approximation of the uncertainty. Here's how you can do it
# using map:
heightLogModel <- map(
  alist(
    height ~ dnorm(mu, exp(log_sigma)), # likelihood
    mu ~ dnorm(178, 20), # mu prior
    log_sigma ~ dnorm(2, 10) # sigma prior
  ), data = Data_2
)
# ... notice the prior for log_sigma. Since log_sigma is continuous now, it can have a Guassian 
# prior. 

# To get the distribution of sigma, you just need to use the same exp() function as in the model
# definition to get back to the natural scale:
posteriorSamples <- extract.samples(heightLogModel)
sigma <- exp(posteriorSamples$log_sigma)
sigma

# ========================================================================================
# Adding a Predictor ------------------------------------------------------
# ========================================================================================
# Let's now take a look at how height in this data covaries with weight. Let's plot height 
# against weight to get an idea of the strength of covariation between the two variables:
plot(Data_2$height ~ Data_2$weight,
     col = "navyblue", lwd = 2,
     xlab = "Weight",
     ylab = "Height",
     main = "Covariation between Height & Weight")

# The linear model strategy -----------------------------------------------
# The strategy is to make the parameter for the mean of a Guassian distribution, mu, into a 
# linear function of the predictor variable and other, new parameters that we invent. This 
# strategy is often simply called the linear model. The linear model strategy assumes that the
# predictor variable has a perfectly constant and additive relationship to the mean of the 
# outcome. In Bayesian analysis, we/the computational software considers every possible 
# combination of the parameter values. For each combination of values, the machine computes the
# posterior probability, which is a measure of relative plausibility, given the model and the 
# data. Thus, the posterior distribution ranks the infinite possible combinations of parameter
# values by their logical plausibility. 

# So. let x be the notation for the column of weight measurements. Now we have a predictor 
# variable x, which is a list of measures of the same length as h. We'd like to say how knowing
# the values in x can help us describe or predict the values in h. Thus, to get weight into the
# model in this way, we define the mean mu as a function of the values in x. Thus, we can 
# describe the model with the following notation:

# likelihood:
#   h[i] ~ N(mu, sigma) 

# linear model:
#   mu[i] = alpha + beta * x[i] 

# alpha prior:
#   alpha ~ N(178, 100)

# Beta prior:
#   Beta ~ Normal(0, 10) 

# sigma prior:
#   sigma ~ Uniform(0, 50)

# Note: remember that the indexing of values means that it is dependent upon unique a predictor
# value, i.e. each row value e.g. for mu.

# The mean mu is thus no longer a parameter to be estimated. Rather, mu[i] is constructed from
# other parameters, alpha and beta, and the predictor variable x. This line is not a stochastic 
# relationship, because the definition of mu[i] is deterministic, not probabilistic! Thus, this
# means that mu is defined by the fact that once we know alpha, beta, and x[i] we also know 
# mu[i]. 

# The parameters alpha and beta are more mysterious. Where did they come from? We made them up.
# The parameters mu and sigma are necessary and sufficient to describe a Guassian distribution.
# But alpha and beta are instead devices we invent in order to manipulate mu, allowing it to 
# vary systematically across each value of the predictor variable in the data.

# One way to think about these made-up parameters is to think of them as targets learning. Each
# parameter is something that must be described in the posterior density. For this example,

# mu[i] = alpha + beta * x[i]

# ... tells us that we are asking to questions about the mean of the outcome:

# 1. What is the expected height, when x[i] = 0? The parameter alpha answers this question. 
#    Alpha is the intercept.

# 2. What is the change in expected height, when x[i] changes by 1 unit? This is the beta 
#    coefficient.

# The remaining lines in model define priors for the parameters to be estimated: alpha, beta, 
# and sigma. All of these are weak priors, leading to inferences that will echo non-Bayesian
# methods of model fitting, such as maximum likelihood. Why have a prior with a mean zero 
# (Beta)? This prior places just as much probability below zero as it does above zero, and when 
# beta = 0, weight has no relationship to height. This is thus seen as a conservative assumption.
# Such a prior will pull probability mass towards zero, leading to more conservative estimates
# than a perfectly flat prior. However, note that a Guassian prior witha sigma of 10 is still 
# very weak, so the amount of conservatism it induces will be very. As you make the sigma in 
# this prior smaller, the amount of shrinkage towards zero increases and your model produces
# more and more conservative estimates about the relationship between h and w. 

# However, do we think that there's just as much chance that the relationship between h and w is
# negative as that it is positive? Of course not. In this context, such a silly prior is 
# harmless, because there is a lot of data which overrides the prior. But in other contexts, 
# our model might need a nudge in the 'right' direction. 

# Fitting the model -------------------------------------------------------
# To fit the linear model, all we need to do is modify the code from the previous section by  
# incorporating the model for the mean and be sure to add our start parameters to the start 
# list. Therefore, the code for the model definition is:

#  height ~ dnorm(mu, sigma)
#  mu <- a + b * weight
#  a ~ dnorm(156, 100)
#  b ~ dnorm(0, 10)
#  sigma ~ dunif(0, 50)

# Notice that the lienar model in the R code uses the assignment operator.

# Model fit:
heightModel_3 <- map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <-  a + b * weight,
    a ~ dnorm(178, 100),
    b ~ dnorm(0, 10),
    sigma ~ dunif(0, 50)
  ),
  data = Data_2
)
# Note that the parametee mu is no longer really a parameter here, because it has been replaced 
# by the linear model, a + b * weight. In the start list, mu is replaced by a, which starts at 
# the overall mean from the last observation, just like mu used to. We also have to add the new
# parameter b to the list, and as is usually a conservative first guess, we start the slope at
# 0, and constrain the maximum to 50. However, note that starting b at a value of zero is 
# different to having b's prior with a mean of zero. The values in the start list don't alter
# the posterior probabilities, while priors definitely do. 

# Interpreting the model fit ----------------------------------------------
precis(heightModel_3, corr = TRUE)
# As can be shown by the correlation matrix, alpha and beta are almost perfectly negatively 
# correlated. Right now, this is harmless. However, in more complex models, strong 
# correlations like this can make it difficult to fit the model to the data. Therefore, we would
# want to utilise spme engineering tricks to avoid it, when possible. 

# The first trick is called centering. This is the process of subtracting the mean of the 
# variable from each value. For example, to create a centered version of the weight variable:
Data_2$weightCentered <- Data_2$weight - mean(Data_2$weight)

# Now, let's refit the model and see what this gains us:
heightModel_4 <- map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b * weightCentered,
    a ~ dnorm(178, 100),
    b ~ dnorm(0, 10),
    sigma ~ dunif(0, 50)
  ),
  data = Data_2
)
precis(heightModel_4, corr = TRUE)

# Thus, the centering has enable for an easier interpretation of the intercept, as now, when 
# mean value of the predictor is zero, it means that the intercept also means: the expected value
# of the outcome is y when the predictor is at its average value x. i.e. roughly 155.

# Plotting posterior inference against data -------------------------------
# Let's start by superimposing just the MAP values over the height and weight data. Then we'll
# slowly add more and more information to the prediction plots, until we've used the entire 
# posterior distribution. 

# To superimpose the MAP values for mean height over the actual data:
plot(height ~ weightCentered, data = Data_2, 
     col = "lightblue", lwd = 2)
abline(a = coef(heightModel_4)["a"], b = coef(heightModel_4)["b"], 
       col = "grey", lwd = 2, lty = "dashed")

# Adding uncertainty around the mean --------------------------------------
# The MAP is just the posterior mean. Plots of the MAP, like above, are useful for getting the
# impression of the magnitude of the estimated influence of a variable, like weight, on an 
# outcome, like height. However, they do a poor job of communicating uncertainty. Therefore, how
# can we plot the uncertainty of the posterior onto the plot? Together, alpha and beta define a 
# line. And so we could sample a bunch of lines from the posterior distribution. Let's, 
# therefore, extract from samples from th posterior:
postSamples <- extract.samples(heightModel_4)
head(postSamples)
# or
postSamples[1:5, ]
# Each row is a correlated random sample from the joint posterior of all three parameters. The
# paired values of alpha and beta on each row define a line, i.e. an intercept and coefficient.
# The scatter of these lines is meaningful, as it alters our confidence in the relationship
# between the predictor and the outcome. 

# Let's begin with simulating the posterior by using the first 10 cases of the Data. The
# following code extracts the first 10 cases and re-estimates the model:
N <- 10
Data_N <- Data_2[1:N, ]

heightModel_N <- map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b * weight,
    a ~ dnorm(178, 100),
    b ~ dnorm(0, 10),
    sigma ~ dunif(0, 50)
  ),
  data = Data_N
)

# Now, let's plot 20 of these lines, to see what the uncertainty looks like.

# First, extract 20 samples from the posterior:
postSamples <- extract.samples(heightModel_N, n = 20)

# Then display the raw data and sample size:
plot(Data_N$weight, Data_N$height,
     xlim = range(Data_N$weight), ylim = range(Data_N$height),
     col = rangi2, xlab = "Weight", ylab = "Height")
mtext(concat("N" = "", N))

# Lastly, plot the lines with transparency:
for (i in 1:20) {
  abline(a = postSamples$a[i], b = postSamples$b[i],
         col = col.alpha("black", .30))
}

# Plotting regression intervals and contours ------------------------------
# It is, however, much more common to see the uncertainty displayed by plotting an interval or
# contour around the MAP regression line. Here's how to plot an interval around the regression 
# line.

# The interval incorporates the uncertainty in both the slope Beta and intercept alpha, at the
# same time. To understand how this works, focus for a moment on a single weight value, say 
# 50 kilograms. You can quickly make a list of 10000 value for mu for an individual who weights
# 50 kg's, by using the samples from the posterior:
muWeight50 <- postSamples$a + postSamples$b * 50
muWeight50

# The above is code for the equation: mu[i] = alpha + beta * x[i]. In this case, x[i] is defined
# as 50, rather than for each individual row index. Thus, it results in a vector of predicted 
# means, one for each random sample from the posterior. Since joint alpha and beta went into 
# computing each, the variation across those means incorporates the uncertainty in and 
# correlation between both parameters. Let's plot the density for this vector of means:
dens(muWeight50, col = rangi2, lwd = 2, xlab = "mu | weight = 50")
# ... since the components of mu have distributions, so too does mu. And since the distributions
# of alpha and beta are Guassian, so to is the distribution of mu (adding Guassian distributions
# together always produces a Gaussian distribution). 

# Since the posterior for mu is a distribution, you can find intervals for it:
HPDI(muWeight50)

# That's good so far, but we need to repeat the above calculation for every weight value on the
# horizontal axes, not just when it is 50kg. Thus, we want to draw 89% HPDIs around the MAP 
# slope. 

# The author provides a simple function link() in order to do this. What link() does is take the
# map model, sample from the posterior distribution, and then compute mu for each case in the 
# data and sample from the posterior distribution:
mu <- link(heightModel_4)
str(mu)
# ... you end up with a matrix of values for mu. Each row is a sample from the posterior. The
# default is 1000 samples. Each column is a case (row) in the data. There are 352 rows in Data_2
# correpsonding to 352 individuals. So there are 352 columns in the mu matrix. 

# Now what can we do with this matrix? Lots of things. The link function provides a posterior 
# of mu for each case we feed it. So above, we have a distribution of mu for each individual
# in the original data. We actually want something slightly different: a distribution for mu
# for each unique weight value on the horizontal axis:

# ... first, we define a sequence of weights to compute predictions for...
weightSeq <- seq(from = 25, to = 70, by = 1)

# ... then we use link to compute mu for each sample from the posterior and for each weight
# defined in weightSeq...
mu <- link(heightModel_3, data = data.frame(weight = weightSeq))
str(mu)
# ... and now there are only 46 columns in mu because we have defined 46 different values for 
# weight. 

# To visualise what we have got here, let's plot the distribution of mu values at each height:
plot(Data_2$height ~ Data_2$weight, type = "n", 
     xlab = "Weight", ylab = "Height")
# ... note here that the "n" type option hides the raw data

# ... then loop over the samples and plot each mu value:
for (i in 1:100) {
  points(weightSeq, mu[i, ], 
         pch = 16, col = col.alpha("navyblue", .1))
}

# The final step is to summarise the distribution for each weight value. We'll use apply, which
# applies a function of your choice to a matrix:
muMean <- apply(mu, 2, mean)
muHPDI <- apply(mu, 2, HPDI, prob = .89)
# ... you can plot these summaries on top of the raw data too:
plot(Data_2$height ~ Data_2$weight, 
     col = col.alpha("navyblue", .6),
     xlab = "Weight", ylab = "Height", main = "Regression with Uncertainty")
# ... plot the MAP line, i.e. the mean mu for each weight value:
lines(weightSeq, muMean, lty = "dashed", lwd = 2)
# ... lastly, plot a shaded region representing the 89% HPDI:
shade(muHPDI, weightSeq, col = col.alpha("darkgrey", alpha = .4))

# It is worth noting that depsite the tight intervals, this is *conditional* on th eassumptions 
# of the model used. Thus, it is helpful to think of the output of our model saying: 
# *conditional on the assumption that height and weight are related by a straight line fit, then
# this is the most plausible line, and these are its plausible bounds.

# How the link function works ---------------------------------------------
# The link function uses the formula you provided to fit the model in order to compute the value
# of the linear model, which it does for each sample from the posterior, for each case of data. 
# To accomplish this manually, we have to code the following:
post <- extract.samples(heightModel_3)

muLink <- function(weight) {
  muLink <- post$a + post$b * weight
  return(muLink)
}

mu <- sapply(weightSeq, muLink)

muMean_2 <- apply(mu, 2, mean)
muHPDI_2 <- apply(mu, 2, HPDI)

# As you can see, they are pretty much the same:
plot(muMean, muMean_2)
plot(muHPDI, muHPDI_2)

# Prediction intervals ----------------------------------------------------
# This section deals with generating an 89% prediction interval for actual heights, not just 
# the average height, mu. This entails incorporating the standard deviation, sigma, and its
# uncertainty too. Omitting priors for brevity, remember that the statistical model is:

# h[i] ~ Normal(mu, sigma)
# mu[i] = alpha + beta * x[i]

# What we have done so far is using samples from the posterior to visualise uncertainty in mu[i],
# the linear model of the mean. However, *actual* predictions of heights depend also upon the
# stochastic/probabilistic definition in the first line. The Guassian distribution for h[i] 
# tells us that the model expects observed heights to be distributed around mu. Additionally,
# the spread around mu is governed by sigma. Thus, this suggests that we need to incorporate
# sigma intot he predictions. 

# Imagine simulating heights. For any unique weight value, you sample from a Guassian 
# distribution with the correct mu for that weight, using the correct value of sigma sampled
# from the posterior. If you do this for every sample from the posterior, for every weight value 
# of interest, you end up with a collection of simulated heights that embody the uncertainty in
# the posterior as well as the uncertainty in Guassian likelihood:
simHeight <- sim(heightModel_3, data = list(weight = weightSeq))
str(simHeight)
# ... this matrix is much like the earlier one, mu, but it contains simulated heights, not 
# distributions of plausible average height, mu. 

# We can summarise these simulated heights in the same way that we summarised the distributions 
# of mu, by using apply:
heightPI <- apply(simHeight, 2, PI, prob = .89)
# ... heightPI contains the 89% posterior prediction interval of observable (according to the
# model and the assumptions embedded) heights, across the values of weightSeq.

# let's now plot everything we've done so far: 1) the MAP line, 2) the shaded region of the 
# 89% posterior region of plausible mu, and 3) the boundaries of the simulated heights that the
# model expects.
# First, we plot the raw data:
plot(Data_2$weight, Data_2$height, 
     col = col.alpha("skyblue", alpha = .5), lwd = 3, pch = 16,
     xlab = "Weight", ylab = "Height", 
     main = "Height vs Weight with Uncertainty")

# Then draw the MAP line:
lines(weightSeq, muMean, lwd = 2, lty = "dashed")

# Now draw the HPDI region of uncertainty for the average, mu:
shade(muHPDI, weightSeq, 
      col = col.alpha("navyblue", alpha = .3))

# Lastly, we draw the PI region for our simulated, *predicted* heights:
shade(heightPI, weightSeq, )

# The wide, grey shaded region represents the area within which the model expects to find 89% 
# of actual heights in the population, given each value of weight. Notice that this wide 
# interval is jagged, representing the simulation variance for the predicted heights. The 
# jaggedness will lessen with increased sample size. 

# We have thus encountered both uncertainty in the parameter values and uncertainty in the 
# sampling process. The posterior distribution is a ranking of the relative plausibilities of
# every possible combination of parameter values. The distribution of simulated outcomes, like
# height, is instead a distribution that includes sampling variation from some process that 
# generates Guassian random variables. This sampling variation is still a model assumption.
# Note that it is possible to view the Guassian likelihood as a purely epistemological 
# assumption (a device for estimating the mean and variance of a variable), rather than an 
# ontological assumption about what future data will look like. In that case, it may not make
# complete sense to simulate outcomes (predict).

# How the sim function operates -------------------------------------------
# Like the link function, it is useful to understand how the sim function operates. For every 
# distribution like dnorm, there is a companion simulation function. For the Guassian 
# distribution, the companion s rnorm, and it simulates sampling from a Guassian distribution.
# What we want R to do is simulate a height for each set of samples, and to do this for each 
# value of weight:
post <- extract.samples(heightModel_3)

weightSeq <- 25:70

simHeight <- sapply(weightSeq, function(weight)
  rnorm(
    n = nrow(post),
    mean = post$a + post$b * weight,
    sd = post$sigma))
heightPI <- apply(simHeight, 2, PI, prob = .89)

shade(heightPI, weightSeq, 
      col = col.alpha("yellow", alpha = .5))

# End file ----------------------------------------------------------------