# This is an R script that performs Gibbs sampling.
# We are going to build an empirical distribution that
# approximates a normal distribution descriged below.
# This example is taken from Gelman's red book, page 288.

# Load the MASS library so that we have access to mvrnorm.
# library('MASS');

# Set up the mean and the covariance matrix.

# The mean:
X1.mean = 0;
X2.mean = 0;

X.mean = matrix( c(X1.mean, X2.mean), 2, 1);

# The variance:
X1.var = 1;
X2.var = 1;
X1X2.cov = 0.5;
X1X2.corr = X1X2.cov * (X1.var * X2.var)^(-0.5);

X.var = matrix( c(
  X1.var, X1X2.cov,
  X1X2.cov, X2.var
  ), 2, 2);

# METHOD 1: Direct Sampling #

# We can generate an empirical distribution by drawing
# directly from a multivariate normal random variable.

# data = mvrnorm(1000, X.mean, X.var);

# METHOD 2: Gibbs Sampling #

# We can generate an empirical distribution by Gibbs sampling.
# We need to know the form of the conditional distributions.
# The basic idea: one has a vector of random variables (Z_i).
# If I know the conditional distribution for each
# Z_i | (Z_j)_{j \neq i} then I can approximate the joint
# density of (Z_i) by sampling from each Z_i conditioned upon
# the most recently sampled data Z_j where j \neq i.  If I construct
# an empirical distribution from this data I converge to the joint
# density of (Z_i).

# Number of iterations.
N = 1000;

# We record the data in data.
data = matrix(0, N, 2);

X1.draw = 0;
X2.draw = 0;

for (i in 1:N){
  # Fix X2 and sample from X1.
  X1gX2.mean = X1.mean + X1X2.corr * (X2.draw - X2.mean);
  X1gX2.var  = 1 - X1X2.corr^2;
  X1.draw = rnorm(1, X1gX2.mean, X1gX2.var^0.5);
  # Now Fix X1 and smaple from X2.
  X2gX1.mean = X2.mean + X1X2.corr * (X1.draw - X1.mean);
  X2gX1.var  = 1 - X1X2.corr^2;
  X2.draw = rnorm(1, X2gX1.mean, X1gX2.var^0.5);
  # Store the data.
  data[i, 1] = X1.draw;
  data[i, 2] = X2.draw;
}

# Check the sample mean and sample variance to see if our
# approximation works.
# mean(data);
# var(data);

# To visualize the data:
# plot(data);
