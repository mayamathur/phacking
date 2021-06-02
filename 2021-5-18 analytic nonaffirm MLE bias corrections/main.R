
# This script goes with the theory in iPad PDF "SAPH - trunc Normal Jeffreys prior"

# Goal:
#  - Try the analytic bias corrections in Godwin paper (based on Cordeiro)
#  - Could also try the Firth 1993 approach of doing Bayesian estimation with
#   Jeffreys prior (e.g., Zhou paper)

# Remember: You'll need to use canonical parameterization for Jeffreys! 

# PRELIMINARIES -----------------------------

library(here)

# general helpers
setwd(here())
source("helper_SAPH.R")

# helpers for the theory herein
setwd(here("2021-5-18 analytic nonaffirm MLE bias corrections"))
source("2021-5-18 helper.R")

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

options(scipen=99)

# **note that mle.tools package will also calculate Fisher info
# in case we 
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

# second derivatives wrt mu and sigma, evaluated at (2,3)
hessian(myLPDF, c(2,3), .x = 0, .crit = 1.96 )

# ~~ Check my first derivatives -----------------------------

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

J1 = Deriv(myLPDF2, "mu")
expect_equal( J1( mu = params[1], sigma = params[2], .x = x, .crit = crit ),
              myJ[1] )


# dl/dsigma matches :)
J2 = Deriv(myLPDF2, "sigma")
expect_equal( J2( mu = params[1], sigma = params[2], .x = x, .crit = crit ),
              myJ[2] )


# ~~ Check my second derivatives -----------------------------


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



# ~~ Check my expected Fisher info -----------------------------

# because I simplified these expressions manually

#bm
( myF = expectFisher( params = params,
                      .crit = crit) )

# check entry 11
# nothing stochastic, so:
expect_equal( -myH[1,1], myF[1,1 ])

# check entry 22
# do this by directly plugging in the trunc normal moment into 
# unsimplified entry from Hessian function
mills = mills(params = params, .crit = crit)
uStar = (crit - mu)/sigma
truncNormalVar = sigma^2*(1 - uStar*mills - mills^2) 

F22 = -( (1/sigma^2) - ( 3 * truncNormalVar / sigma^4 ) +
  ( mills^2 + uStar * mills ) * ( ( crit - mu ) / sigma^2 )^2 -
  mills * ( 2 * (crit - mu) / sigma^3 ) )

expect_equal(F22, myF[2,2])

# entry 12
truncNormalMean = mills*sigma
F12 = -( ( -2 * truncNormalMean / sigma^3 ) +
  ( mills^2 + uStar * mills ) * ( ( crit - mu ) / sigma^3 ) -
  mills / sigma^2 )

expect_equal(F12, myF[1,2], tol = 0.0000000001)
expect_equal(F12, myF[2,1], tol = 0.0000000001)


# interesting that off-diagonal entries are so close to 0
round(myF, 3)


# ~~ Check my third derivatives -----------------------------


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









