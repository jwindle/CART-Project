# This R script simulates the posterior distribution of the model
#   Y ~ N(0, V/phi) 
#   phi ~ Gamma.
# The gamma distribution is conjugate to the normal distribution
# when mixing the precision parameter.  Thus a closed form
# solution of the posterior distribution exists.  The calculations
# are outlined in the documentation.

# Load the mvtnorm library for dmvnorm.
print("Make sure to load 'mvtnorm' with library('mvtnorm').");
# library('mvtnorm');

# A gamma distribution has a shape parameter n and a scale
# parameter d.  We will always be working with parameters
# of the form n/2 and d/2 to match our analysis.
n = 3.0;  # shape
d = 0.5;  # rate
# Make sure to look up help('dgamma') to see what R is doing.
# I spent hours trying to figure out why this script didn't work
# only to find that I was specifying the scale parameter instead
# of the the rate parameter.

# We need to construct the covariance matrix V in our model.
# Let's just arbitrarily define that for now.
dim = 2;
V = matrix( c(
  1.0, 0.5,
  0.5, 1.0
  ), dim, dim);

# The observed data...
# Set the actual precision.
tau = 10.0;
# Produce synthetic data.
y = rmvnorm(1, matrix(0, dim, 1), V / tau);

# The number of iterations for our simulation.
N = 10000;

# Pull from our prior N times.
phi = rgamma(N, n/2, rate=d/2);

# We now want to approximate the density of the prior.
# We will use the histogram data structure.
# Check out help('hist') for more information.
# If we let H = hist(...) then H$mids gives us the midpoints
# of our bins and H$density gives us the histogram evaluated
# at H$mids.

# Minimum number of bins.
K = 41;

# The histogram.  We do want to normalize.  We do not want to plot now.
H = hist(phi, breaks=K, freq=FALSE, plot=FALSE);

# Number of midpoints;
L = length(H$mids);

# Now calculate our proportional posterior.
# The derivation is given in the documentation.
emp_posterior = 1:L;
for ( i in 1:L ) {
  emp_posterior[i] = dmvnorm(y, matrix(0, dim, 1), V/H$mids[i]) * H$density[i];
}

# We can plot the prior density with
# plot(H$mids, H$density);

# We can plot the posterior with
# plot(H$mids, emp_posterior);

# We can also calculate the proportional posterior ``explicitly.''
# I did this initially as a check upon my calculation above.
gamma_dens = dgamma(H$mids, n/2, rate=d/2)
posterior = 1:L;
for ( i in 1:L ) {
  posterior[i] = dmvnorm(y, matrix(0, dim, 1), V/H$mids[i]) * gamma_dens[i];
}

# To check everything...

# The posterior distribution should be
# Gamma(n_post/2,d_post/2) where
n_post = n + dim;
d_post = d + y %*% solve(V) %*% t(y);

# To get an idea let's plot the prior and the posterior.
# mygrid = seq(0, 5, 0.1);
# plot(mygrid, dgamma(mygrid, n/2, rate=d/2) );
# plot(mygrid, dgamma(mygrid, n_post/2, rate=d_post/2) );

##############
## PLOTTING ##
##############

# To plot everything at once:
# Check out help('par').
par(mfrow=c(2,2));

# Analytic Prior
plot(H$mids, gamma_dens, type="l");
title("Analytic Prior");

# Empirical Prior
plot(H$mids, H$density, type="l");
title("Empirical Prior");

# Analytic Posterior
plot(H$mids, dgamma(H$mids, n_post/2, rate=d_post/2) , type="l");
title("Analytic Posterior");

# Empirical Posterior
plot(H$mids, emp_posterior, type="l");
title("Empirical Posterior");

# To write a plot to file.
# postscript('file.ps');
# ... plotting commands ...
# I copied and pasted from above.
# par(mfrow=c(2,2));
# Analytic Prior
# plot(H$mids, gamma_dens, type="l");
# title("Analytic Prior");
# Empirical Prior
# plot(H$mids, H$density, type="l");
# title("Empirical Prior");
# Analytic Posterior
# plot(H$mids, dgamma(H$mids, n_post/2, rate=d_post/2) , type="l");
# title("Analytic Posterior");
# Empirical Posterior
# plot(H$mids, emp_posterior, type="l");
# title("Empirical Posterior");
# dev.off();
