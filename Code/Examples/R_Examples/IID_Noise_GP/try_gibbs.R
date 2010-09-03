# We want to implement the Gibbs sampler for our toy problem.
# Make sure the prior parameters and the covariance matrix
# match any other script you might use for comparison.

# ASSUMPTIONS!!!
# Currently, we are assuming that there is only ONE data point.

# Function for covariance matrix.
# We have chosen to use the squared exponential covariance function.
#CovFunc <- function(t1, t2, k1, k2, k3){
#  k1 * exp( -1*(t1-t2)*(t1-t2)/(2*k2) ) + k3 * (t1==t2);
#}

# Size of ``time series.''
#T = 2;

# Generate covariance matrices.

# F | sig ~ sig N(0,K).
#K = matrix(0, T, T);
#for (i in 1:T){
#  for (j in 1:T){
#    K[i,j] = CovFunc(i,j, 0.5, 1, 0.5);
#  }
#}

# Y | F, sig ~ sig N(0, I).
Id = diag( rep(1,T) );

# The joint covariance matrix is given by (Y F) ~ N(0,A)
# where A = [Id + K   K 
#              K^T    K ];
A = matrix(0, 2*T, 2*T);
# It would be nice to find a cleaner way to express this in code.
A[1:T,1:T] = Id + K;
A[1:T,(T+1):(2*T)] = K;
A[(T+1):(2*T), 1:T] = t(K);
A[(T+1):(2*T), (T+1):(2*T)] = K;

# We need to use the inverse of A.
InvA = solve(A);

# We want to look at the precision phi = 1 / sig^2
# instead of the variance sig^2.  This transformation
# ensures that phi is from a gamma distribution, which
# is conjugate in this case.  Our prior parameters for
# phi ~ Gamma(n/2, d/2).
n.prior = 3.0;
d.prior = 2.0;

# Given one data point our joint posterior is
# p(sig, f | y).  We can simulate this by sampling
# p(f | y, sig) and then p(sig | f, y) repeatedly.
# y = matrix(data[1,]);

# Now we want to do our Gibbs sampling.
samples = 5000;

phi.data = matrix(0, samples + 1, 1);

# I have no reason for doing this, but it seems reasonable.
phi.data[1] = 0.5; # rgamma(1, n.prior/2, d.prior/2);

# As mentioned above we need two conditional distributions,
# p(f | y, sig) and p(sig | f, y).  We constructed the matrix
# A above because the joint distribution of (y,f) is given by
# p(y, f| sig) ~ N(0, sig^2 A).  From J. Scott's notes we know
# that p(\sig | y, f) is then given by Gamma(n.post/d, d.post/2)
# where n.post and d.post are defined as below.  We must be
# careful with the conditional distribution p(f | y, sig).  Again,
# according to J. Scott's notes we know that f ~ N(m, \sig^2 C)
# where m and C are defined below.  This is something that requires
# a bit of elbow grease to show, though the calculation is
# straightforward.  I introduced new variables to match
# J. Scott's notes.  J. Scott's F = Id in this case.

# To match J. Scott's notes...
R = K;
V = Id;
Q = R + V;
InvQ = solve(Q);
C = R - R %*% InvQ %*% R;

# Now we actually iterate.
for (i in 1:samples){
  # f ~ p(f | y, sig) = N(m, C / phi) where
  m = R %*% InvQ %*% y;
  # C is defined above.
  # Notice that we could speed up our calculation by
  # storing some of the above computataions.
  # Sample f...
  f = rmvnorm(1, m, C / phi.data[i]);
  # phi ~ p(phi | f, y) ~ Gamma(n.post/2, d.post/2);
  # Let x = (y, f).
  x = matrix( c(y,f) );
  # Then the posterior parameters are given below.
  # Note that there are 2T components to x = (y,f).
  n.post = n.prior + 2 * T;
  d.post= d.prior + t(x) %*% InvA %*% x;
  # print( c(1.0/sqrt(phi.data[i]), f %*% y, d.post-d.prior) );
  phi.data[i+1] = rgamma(1, n.post/2, rate=d.post/2);
}

# We can plot our data with
# hist(phi.data, breaks=40, prob=TRUE);

# Let's compare our Gibbs sampler to the solution
# we have calculated by hand, which can be found
# in J. Scott's notes.

# Now calculate our posterior using some analysis.  These derivations
# can be found in the notes by J. Scott.
n.post = n.prior + T;

# I made the mistake of inverting the wrong matrix here.
# We need to use the variance of Y marginally.  That is we
# must use Q from Y ~ N(0, sig Q).
d.post = d.prior;
d.post = d.post + t(y) %*% solve(Q) %*% y;

# Now set up a grid and calculate the posterior on that grid.
mygrid = seq(0.05, 4, 0.05);
gamma.post = dgamma(mygrid, n.post/2, rate=d.post/2);

# Now plot.
par(mfrow=c(1,3));
plot(mygrid, gamma.post);
hist(phi.data, breaks=40, prob=TRUE);
hist(check.data[,1], breaks=40, prob=TRUE);

# To compary some summary statistics look at the mean and variance.
# The mean...
print( c(mean(phi.data), n.post/d.post, mean(check.data[,1])) );
# The variance...
print( c(var(phi.data), 2*n.post/d.post^2, var(check.data[,1])) );
