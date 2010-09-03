# We have generated the data.  Now we want to calculate
# the posterior distribution.  We really do not need
# to use Monte Carlo here.

# We need to have defined the covariance matrix K beforehand.

# Generate prior distribution.  This can be done analytically
# or with Monte Carlo.

# The prior parameters.
n = 3.0;
d = 2.0;

# The grid.
# Make sure you start at a point bigger than zero to avoid singularities.
mygrid = seq(0.05, 3, 0.05);

# The number of points in our grid.
mysize = length(mygrid);

# The prior density.
gamma.prior = dgamma(mygrid, n/2, rate=d/2);

# Now find the posterior distribution given all of our observations.

# First set aside some space for the posterior.
gamma.post = gamma.prior;

# The dimension of our observation.
N = nrow(data);

# The number of ``time steps,'' that is the number of observations.
T = ncol(data);

# We assume that Q has already been defined so that
# for the normal random variable Y | sig ~ N(0, sig^2 Q).
# Calculate the posterior.
gamma.post = log(gamma.prior);
for (i in 1:N){
  for (k in 1:mysize){
    gamma.post[k] = gamma.post[k] + dmvnorm(data[i,], rep(0,T), Q / gamma.prior[k], log=TRUE);
  }
  # gamma.post = gamma.post / (sum(gamma.post) * 0.1);
}

gamma.post2 = log(gamma.prior);
# gamma.post2 = rep(0, mysize);
for (i in 1:N){
  for (k in 1:mysize){
    gamma.post2[k] = gamma.post2[k] + 0.5 * T * log(gamma.prior[k]) - 0.5 * gamma.prior[k] * t(data[i,]) %*% solve(Q) %*% data[i,];
  }
}

# Exponentiate and normalize.
gamma.post = exp(gamma.post);
gamma.post2 = exp(gamma.post2);
# gamma.post = gamma.post / sum(gamma.post);

# Now calculate our posterior using some analysis.  These derivations
# can be found in the notes by J. Scott.
n_post = n + N*T;

# I made the mistake of inverting the wrong matrix here.
# We need to use the variance of Y marginally.  That is we
# must use Q from Y ~ N(0, sig Q).
d_post = d;
for(i in 1:N){
  d_post = d_post + t(data[i,]) %*% solve(Q) %*% data[i,];
}

# Plot the data.
par(mfrow=c(2,3));
plot(mygrid, gamma.prior);
plot(mygrid, gamma.post);
#plot(mygrid, gamma.post2);
plot(mygrid, dgamma(mygrid, n_post/2, rate=d_post/2) );
