# Generate data for the GP with IID Noise.

# Make sure to load the necessary library.
# library(mvtnorm)

# Function for covariance matrix.
# We have chosen to use the squared exponential covariance function.
CovFunc <- function(t1, t2, k1, k2, k3){
  k1 * exp( -1*(t1-t2)*(t1-t2)/(2*k2) ) + k3 * (t1==t2);
}

# Size of ``time series.''
T = 40;

# Generate covariance matrix for p(f) ~ N(0, sig^2 K).
K = matrix(0, T, T);
for (i in 1:T){
  for (j in 1:T){
    K[i,j] = CovFunc(i,j, 0.5, 1, 0.5);
  }
}

# Generate covariance matrix for p(y) ~ N(0, sig^2 Q)
# where Q = I + K.  We perturb each element by
# our time series independently.
Q = diag(1,T) + K;

# The var scale we want.
v = 2.0;

# Inverting this to get get our analog of phi.
phi = 1/v;

# And generating our scaled matrix.
vQ = v*Q;

# Now generate our independent draws...

# Number of data points.
N = 1;

# Draw N data points from Y ~ N(0, v Q).
y = rmvnorm(N, rep(0, T), vQ);

y = matrix(y, T, 1);

# To write the data...
# write(data, "y.data", ncolumns=T);
