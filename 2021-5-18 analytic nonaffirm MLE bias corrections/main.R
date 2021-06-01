
# Goal: Try to correct small-sample bias using the Firth 1993 
# approach of doing Bayesian estimation with Jeffreys prior
# (e.g., Zhou paper)
# and other analytic bias corrections as in Godwin paper

# Remember: You'll need to use canonical parameterization for Jeffreys! 

# PRELIMINARIES -----------------------------

library(here)
setwd(here())
source("helper_SAPH.R")


library(metafor)
library(weightr)
library(ggplot2)
library(dplyr)
# data-wrangling packages
library(dplyr)
library(ggplot2)
library(data.table)
library(tidyverse)
library(fastDummies)
# meta-analysis packages
library(metafor)
library(robumeta)
# other
library(xtable)
library(testthat)
# for this project
library(truncdist)
#library(ExtDist)
library(gmm)  # https://stackoverflow.com/questions/63511986/error-package-or-namespace-load-failed-for-gmm-in-dyn-loadfile-dllpath-dl
library(tmvtnorm)
library(TruncatedNormal)
library(mle.tools)
library(fitdistrplus)
library(Deriv)

# **note that mle.tools package will also calculate Fisher info in case we 
#   want to try the Jeffreys correction
# see expected.varcov 

# 2021-5-31: TRY CORDEIRA CORRECTION AS IN GODWIN -----------------------------

# ~ Check my theory -----------------------------

# check the derivatives and partial derivatives I got on iPad

# simple package example
func2 <- function(vec){
  vec[1]^2 + vec[2]^2
}
# first derivatives wrt vec[1] and wrt x[2], evaluated at (2,3)
jacobian(func2, c(2,3))  
# second derivatives wrt vec[1] and wrt x[2] and the mixed ones, evaluated at (2,3)
hessian(func2, c(2,3))  # second derivative evaluated at x=4


# params: (mu, sigma)
# NOTE parameterization in terms of sigma rather than sigma^2
myLPDF = function(params, .x, .crit) {
  mu = params[1]
  sigma = params[2]
  -log( sqrt(2*pi) * sigma) - (.x - mu)^2 / (2 * sigma^2 ) -
    pnorm( q = .crit, mean = mu, sd = sigma, log=TRUE)
}
myLPDF( c(0, .5), 0, 1.96)

# same, but separating parameters
myLPDF2 = function(mu, sigma, .x, .crit) {
  -log( sqrt(2*pi) * sigma) - (.x - mu)^2 / (2 * sigma^2 ) -
    pnorm( q = .crit, mean = mu, sd = sigma, log=TRUE)
  
}

# second derivatives wrt mu and sigma, evaluated at (2,3)
hessian(myLPDF, c(2,3), .x = 0, .crit = 1.96 )

# ~~ Check my first derivatives -----------------------------

# Mills ratio for right truncation
mills = function(params, .crit) {
  mu = params[1]
  sigma = params[2]
  uStar = (.crit - mu)/sigma
  dnorm(uStar) / pnorm(uStar) 
}

# mine
myJacobian = function(params, .x, .crit) {
  
  mu = params[1]
  sigma = params[2]
  
  # entry 1: derivative wrt mu
  J1 = ( (.x - mu) / sigma^2 ) +  # this part matches Deriv!
    mills(params = params, .crit = .crit) * (1/sigma)
  
  # entry 2: derivative wrt sigma
  J2 = (-1/sigma) + ( (.x - mu)^2 / sigma^3 ) +
    mills(params = params, .crit = .crit) * (.crit - mu)/sigma^2
  
  return( matrix( c(J1, J2), nrow = 1 ) )
}

# set up some values to plug in
params = c(-0.1, .2)
x = 0
crit = 3
mu = params[1]
sigma = params[2]
myLPDF(params, .x = x, .crit = crit)

# dl/dmu matches :)
( myJ = myJacobian( params = params,
                    .x = x, 
                    .crit = crit) )

library(Deriv)
J1 = Deriv(myLPDF2, "mu")
expect_equal( J1( mu = params[1], sigma = params[2], .x = x, .crit = crit ),
              myJ[1] )


# dl/dsigma matches :)
J2 = Deriv(myLPDF2, "sigma")
expect_equal( J2( mu = params[1], sigma = params[2], .x = x, .crit = crit ),
              myJ[2] )


# ~~ Check my second derivatives -----------------------------


myHessian = function(params, .x, .crit) {
  
  mu = params[1]
  sigma = params[2]
  mills = mills(params = params, .crit = .crit)
  uStar = (.crit - mu)/sigma
  
  # entry 1,1: dl/dmu^2
  H11 = (-1/sigma^2) + (1/sigma^2) * ( mills^2 + uStar * mills )
  
  # entry 2,2: dl/dsigma^2
  H22 = (1/sigma^2) - ( 3 * (.x - mu)^2 / sigma^4 ) +
    ( mills^2 + uStar * mills ) * ( ( .crit - mu ) / sigma^2 )^2 -
    mills * ( 2 * (.crit - mu) / sigma^3 )
  
  # entry 1,2 and 2,1: dl / dmu dsigma
  H12 = ( -2 * (.x - mu) / sigma^3 ) +
    ( mills^2 + uStar * mills ) * ( ( .crit - mu ) / sigma^3 ) -
    mills / sigma^2
  
  return( matrix( c(H11, H12, H12, H22), nrow = 2 ) )
}


# dl/dmu^2 matches :)
( myH = myHessian( params = params,
                   .x = x, 
                   .crit = crit) )

H11 = Deriv(J1, "mu")
expect_equal( H11( mu = params[1], sigma = params[2], .x = x, .crit = crit ),
              myH[1,1] )


# dl/dsigma^2 matches :)
H22 = Deriv(J2, "sigma")
expect_equal( H22( mu = params[1], sigma = params[2], .x = x, .crit = crit ),
              myH[2,2] )

# dl/dmu dsigma matches :)
H12 = Deriv(J1, "sigma")
expect_equal( H12( mu = params[1], sigma = params[2], .x = x, .crit = crit ),
              myH[1,2] )



# ~~ Check my third derivatives -----------------------------
# these ones are needed for Cordeiro correction, but not Jeffreys prior

# entry: a number like "121" to say which derivative we want
#  where parameter 1 is mu and parameter 2 is sigma
myThirdDerivs = function(params, .x, .crit, .entry) {
  
  mu = params[1]
  sigma = params[2]
  
  # terms that will show up a lot
  mills = mills(params = params, .crit = .crit)
  uStar = (.crit - mu)/sigma
  termA = mills^2 + uStar * mills
  termB = 2*mills + uStar
  termC = termA * termB - mills
  
  ### third erivatives wrt mu
  # entry 111: dl/dmu^3
  if ( .entry == 111 ) {
    return( (1/sigma^3) * termC )
  }
  
  # entry 221=212=122: dl/(dsigma^2 dmu)
  if ( .entry %in% c(221, 212, 122) ) {
    return( (1/sigma^4) * ( 6*(.x - mu) - 2*(.crit - mu)*termA +
                              (.crit - mu)^2/sigma * termC +
                              2*mills*sigma ) )
  }
  
  # entry 121=211=112: dl/(dmu dsigma dmu)
  if ( .entry %in% c(121, 211, 112) ) {
    return( (1/sigma^3) * ( 2 - 2*termA + uStar*termC ) )
  }
  
  # entry 222: dl/(dsigma^3)
  if ( .entry == 222 ) {
    return( (1/sigma^4) * ( (-2*sigma) + (12/sigma)*(.x-mu)^2 -
              ( 4*(.crit - mu)^2/sigma )*termA +
              ( (.crit - mu)^3/sigma^2 )*termC +
              6*(.crit - mu)*mills + 2*( (.crit - mu)^2/sigma )*termA ) )
    }
  
  # H11 = (-1/sigma^2) + (1/sigma^2) * ( mills^2 + uStar * mills )
  # 
  # # entry 2,2: dl/dsigma^2
  # H22 = (1/sigma^2) - ( 3 * (.x - mu)^2 / sigma^4 ) +
  #   ( mills^2 + uStar * mills ) * ( ( .crit - mu ) / sigma^2 )^2 -
  #   mills * ( 2 * (.crit - mu) / sigma^3 )
  # 
  # # entry 1,2 and 2,1: dl / dmu dsigma
  # H12 = ( -2 * (.x - mu) / sigma^3 ) +
  #   ( mills^2 + uStar * mills ) * ( ( .crit - mu ) / sigma^3 ) -
  #   mills / sigma^2
  
}

# check entry 111
( mine = myThirdDerivs( params = params,
                        .x = x, 
                        .crit = crit,
                        .entry = 111 ) )

func = Deriv(H11, "mu")
expect_equal( func( mu = params[1],
                    sigma = params[2],
                    .x = x,
                    .crit = crit ),
              mine )

# check entry 221 and equivalent ones
( mine = myThirdDerivs( params = params,
                        .x = x, 
                        .crit = crit,
                        .entry = 221 ) )

func = Deriv(H22, "mu")
expect_equal( func( mu = params[1],
                    sigma = params[2],
                    .x = x,
                    .crit = crit ),
              mine )

# 122
func = Deriv(H12, "sigma")
expect_equal( func( mu = params[1],
                    sigma = params[2],
                    .x = x,
                    .crit = crit ),
              mine )


# check entry 121
( mine = myThirdDerivs( params = params,
                        .x = x, 
                        .crit = crit,
                        .entry = 121 ) )

func = Deriv(H12, "mu")
expect_equal( func( mu = params[1],
                    sigma = params[2],
                    .x = x,
                    .crit = crit ),
              mine )

# check entry 222
( mine = myThirdDerivs( params = params,
                        .x = x, 
                        .crit = crit,
                        .entry = 222 ) )

func = Deriv(H22, "sigma")
expect_equal( func( mu = params[1],
                    sigma = params[2],
                    .x = x,
                    .crit = crit ),
              mine )

# all derivatives match; woohoo! 
#@SHOULD PROBABLY SAVE THIS AS SEPARATE FILE AND CLEAN UP THE 
# THREE FUNCTIONS INTO A SINGLE FUNCTION THAT GIVES JUST THE DERIVATIVE
#(OF WHATEVER ORDER) THAT YOU ASKED FOR



# 2021-5-18: TRY COX-SNELL BIAS CORRECTION -----------------------------

# ~ Toy example from package -----------------------------

# simple example from mle.tools::coxsnell.bc
# not the Jeffreys correction, but rather a different correction that also uses 
#  something related to the Fisher info
set.seed(1)
## Normal distribution
pdf <- quote(1 / (sqrt(2 * pi) * sigma) * exp(-0.5 / sigma ^ 2 * (x - mu) ^ 2))
lpdf <- quote(- log(sigma) - 0.5 / sigma ^ 2 * (x - mu) ^ 2)
x <- rnorm(n = 100, mean = 0.0, sd = 1.0)

( mu.hat <- mean(x) )
( sigma.hat = sqrt((length(x) - 1) * var(x) / length(x)) )

coxsnell.bc(density = pdf, logdensity = lpdf, n = length(x), parms = c("mu", "sigma"),
            mle = c(mu.hat, sigma.hat), lower = '-Inf', upper = 'Inf')

# ~ Random generates straight from trunc dist -----------------

# try one example using parameters that yielded very biased estimates in 
#   2021-5-3 folder, Expt 3.3

p = data.frame( Mu = 0.1,
                T2 = 0.25,
                m = 500,
                t2w = 0.25,
                se = .5,
                
                Nmax = 1,
                hack = "affirm",
                
                k = 20, 
                k.hacked = 0 )

crit = qnorm(.975)


# withot bias correction:
# meanHat = 0.4442739
# varHat = 3.349863

sim.reps = 1

meanHat = c()
varHat = c()
tstats = c()

# set up PDF and score for Cox-Snell thing
# example for plain-vanilla Normal:
# pdf <- quote(1 / (sqrt(2 * pi) * sigma) * exp(-0.5 / sigma ^ 2 * (x - mu) ^ 2))
# lpdf <- quote(- log(sigma) - 0.5 / sigma ^ 2 * (x - mu) ^ 2)


mean = .Mu / .se,
#@doesn't yet use the delta-method thing
sd = sqrt( (1/.se^2) * (.T2 + .t2w + .se^2) )


myPDF = "1 / (sqrt(2 * pi) * sigma) * exp(-0.5 / sigma ^ 2 * (x - mu) ^ 2) * pnorm( q = crit, mean = mu, sd = sigma)"

myLPDF = "-log(sigma) - 0.5/sigma^2 * (x - mu)^2 + pnorm( q = crit, mean = mu, sd = sigma, log=TRUE)"

# ?D:
# The internal code knows about the arithmetic operators +, -, *, / and ^, and the single-variable functions exp, log, sin, cos, tan, sinh, cosh, sqrt, pnorm, dnorm, asin, acos, atan, gamma, lgamma, digamma and trigamma, as well as psigamma for one or two arguments (***but derivative only with respect to the first***). (Note that only the standard normal distribution is considered.)

for ( i in 1:sim.reps ) {
  
  if ( i %% 25 == 0 ) cat( paste("\n\nStarting rep", i ) )
  
  x2 = rtrunc(n = 20,
              spec = "norm",
              mean = p$Mu / p$se,
              sd = sqrt( (1/p$se^2) * (p$T2 + p$t2w + p$se^2) ),
              b = crit)
  
  # get MLEs
  mle.fit = mle.tmvnorm( as.matrix(x2, ncol = 1), lower=-Inf, upper=crit)
  mles = coef(mle.fit)
  
  # Cox-Snell bias correction
  #@BE CAREFUL ABOUT SD VS. VAR THING!
  coxsnell.bc(density = myPDF,
              logdensity = myLPDF,
              n = length(x2),
              parms = c("mu", "sigma"),
              mle = c(mles[1], mles[2]),
              lower = '-Inf',
              upper = 'Inf')
  
  
  
  meanHat = c(meanHat, mles[1])
  varHat = c(varHat, mles[2])
  tstats = c(tstats, x2)
}


mean(meanHat); p$Mu/p$se
mean(varHat); (1/p$se^2) * (p$T2 + p$t2w + p$se^2)

# should be normal
hist(meanHat)

# trunc normal
hist(tstats)



# 2021-5-18: TRY JEFFREYS BIAS CORRECTION -----------------------------

# ~ Is Zhou's variance wrong? -----------------------------

# In Appendix A, Zhou calculated the Fisher info using a variance that seems to be wrong. 
#  Simulate to confirm that I'm correct about the variance mistake

mu = 1.1
sig = 0.7
u = .8  # cut point for LEFT-truncated normal

( uStar = (u - mu)/sig )

# Mills ratio for LEFT truncation (i.e., alpha-function in Zhou)
mills = function(.uStar) {
  dnorm(.uStar) / ( 1 - pnorm(.uStar ) )
}


mills(uStar)

# Zhou's is clearly WRONG
( VarZhou = sig^2 * ( 1 + uStar*mills(uStar) ) )

# mine is from Wikipedia:
# https://en.wikipedia.org/wiki/Truncated_normal_distribution#One_sided_truncation_(of_lower_tail)[4]
( VarMine =  sig^2 * ( 1 + uStar*mills(uStar) - mills(uStar)^2 ) )
# THIS IS RIGHT


vartrunc( spec = "norm",
          mean = mu,
          sd = sig,
          a = u,
          b = Inf )



# simulate
x2 = rtrunc(n = 5000,
            spec = "norm",
            mean = mu,
            sd = sig,
            a = u)

var(x2)
# exactly matches vartrunc









