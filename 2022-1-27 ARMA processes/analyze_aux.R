

# Goal: Understand intuitively which kinds of time series are stationary (at least for their expectations)
# e.g., I don't understand why this holds for a moving-average process?

# Note about simts pkg: 
#  The error "Need to supply initial values within the ts.model object"
#  means you specified arguments that are out of possible range


# https://cran.r-project.org/web/packages/simts/vignettes/vignettes.html
library(simts)


n = 10000


# White noise
x = gen_gts( n, WN(sigma2 = 1) )
mean(x)
plot(x)

# AR1, positive autocorrelation
# specification in package docs:
# X_t = φ X_{t-1} + \varepsilon_t
x = gen_gts( n, AR1(phi = 0.95, sigma2 = 1) )
mean(x)
plot(x)

# ARMA(1): special case with constant moving avg (should be just like AR(1))
# X_t = ∑_{j = 1}^p φ_j X_{t-j} + ∑_{j = 1}^q θ_j \varepsilon_{t-j} + \varepsilon_t
# "ar" argument is the phi's (provide only one to have an ARMA(1) process)
# "ma" is the theta's
x = gen_gts( n, ARMA(ar = .95, ma = 0, sigma2 = 1) )
mean(x)
plot(x)

# ARMA(1): special case with (almost) no autocorrelation
x = gen_gts( n, ARMA(ar = 0.01, ma = 0.9, sigma2 = 1) )
mean(x)
plot(x)
# @I don't understand why this doesn't increase its average over time? 


# ARMA(1): autocorrelated and moving average is positive
x = gen_gts( n, ARMA(ar = 0.95, ma = 0.9, sigma2 = 1) )
mean(x)
plot(x)




