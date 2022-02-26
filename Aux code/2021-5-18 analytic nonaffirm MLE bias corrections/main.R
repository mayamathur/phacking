

#BM: SEARCH "@". CONFUSED ABOUT WHY I'M GETTING NEGATIVES IN THE INVERSE FISHER INFO EVEN FOR N=1. AND YET THE FISHER INFO ENTRIES CHECK OUT AGAINST DERIV().


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
library(Hmisc)
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
# but it doesn't work with complicated functions


# 2021-6-4: TRY JEFFREYS BIAS CORRECTION -----------------------------

# set up some values to plug in
params = c(0, 1) # mu, sigma
#crit = qnorm(.975)
crit = 999  # no truncation
mu = params[1]
sigma = params[2]


( myF = expectFisher( params = params,
                      .crit = crit) )

expect_equal( det(myF), 
              myF[1,1]*myF[2,2] - myF[1,2]*myF[2,1] )

# this is the Jeffreys prior for these particular parameter values
Jeffreys(params = params, .crit = crit)  # for standard untruncated normal, equal to 1/sqrt(2)

# compare to Jeffreys prior for untruncated normal with unknown mean and variance: 1/sigma^2
# https://stats.stackexchange.com/questions/156199/jeffreys-prior-for-normal-distribution-with-unknown-mean-and-variance
1/(sigma^2)

# ~ Get intuition for the prior --------------------------------

### Hold constant sigma at its true value and vary mu, and plot prior as fn of mu
# plot just the prior at different values
dp = data.frame( mu = seq(-2, 2, .05),
                 sigma = sigma )

dp = dp %>% rowwise() %>%
  mutate( prior = Jeffreys(params = c(mu, sigma), .crit = crit) )


ggplot( data = dp, 
        aes(x = mu, y = prior) ) +
  geom_line()

### Now vary both mu and sigma
dp = expand.grid( mu = seq(-2, 2, .1),
                  sigma = c(0.001, 0.1, 0.5, 1, 2, 3, 5) )

dp$sigma.pretty = paste( "sigma = ", dp$sigma, sep = "" )

dp = dp %>% rowwise() %>%
  mutate( prior = Jeffreys(params = c(mu, sigma), .crit = crit) )

# **in project log
ggplot( data = dp, 
        aes(x = mu, y = prior) ) +
  geom_line() +
  facet_wrap( sigma.pretty ~.,
              scales = "free_y" )

### Now also plot the unnormalized posterior for a given observation from distribution

x = 1.9

dp %>% rowwise() %>%
  mutate(joint.post = joint_post( params = c(mu, sigma),
                                  .xVec = x,
                                  .crit = crit ) )


# ~ Code up MCMC from scratch --------------------------------

# try one example using parameters that yielded very biased estimates in 
#   2021-5-3 folder, Expt 3.3 (k=20)
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




# test drive
set.seed(451)
x = rtrunc(n = 20,
           spec = "norm",
           mean = p$Mu / p$se,
           sd = sqrt( (1/p$se^2) * (p$T2 + p$t2w + p$se^2) ),
           b = crit)

# MLEs for comparison
mle.fit = mle.tmvnorm( as.matrix(x, ncol = 1), lower=-Inf, upper=crit)
( mles = coef(mle.fit) )
# targets:
p$Mu / p$se
sqrt( (1/p$se^2) * (p$T2 + p$t2w + p$se^2) )

# true parameters
joint_post( params = c(p$Mu / p$se, sqrt( (1/p$se^2) * (p$T2 + p$t2w + p$se^2) ) ),
            .xVec = x,
            .crit = crit )

#bm: wanted to look into off-the-shelf options to do the Bayesian estimation
# like brms
# or alternatively think about getting mean of posterior analytically


# optim( par = c(mles[1], mles[2]),
#        # negative because optim gets the minimum
#        fn = function(.par) -joint_post(params = .par,
#                                         .xVec = x,
#                                         .crit = crit) )




# # from earlier
# for ( i in 1:sim.reps ) {
#   
#   if ( i %% 25 == 0 ) cat( paste("\n\nStarting rep", i ) )
#   
#   x = rtrunc(n = 20,
#               spec = "norm",
#               mean = p$Mu / p$se,
#               sd = sqrt( (1/p$se^2) * (p$T2 + p$t2w + p$se^2) ),
#               b = crit)
#   
#   # get MLEs
#   mle.fit = mle.tmvnorm( as.matrix(x2, ncol = 1), lower=-Inf, upper=crit)
#   mles = coef(mle.fit)
#   
#   # Bayesian
#   
#   
#   
#   meanHat = c(meanHat, mles[1])
#   varHat = c(varHat, mles[2])
#   tstats = c(tstats, x2)
# }


mean(meanHat); p$Mu/p$se
mean(varHat); (1/p$se^2) * (p$T2 + p$t2w + p$se^2)



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

# # set up some values to plug in
# params = c(-0.1, .2)
# x = 0
# crit = 3

# another set of values
params = c(0, 1)
x = 0.5
crit = 0.5
mu = params[1]
sigma = params[2]


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
( myF = expectFisher( params = params,
                      .crit = crit) )
solve(myF)

# check entry 11
# nothing stochastic, so:
expect_equal( -myH[1,1], myF[1,1 ])

# check entry 22
# do this by directly plugging in the trunc normal moment into 
# unsimplified entry from Hessian function
Mills = mills(params = params, .crit = crit)
uStar = (crit - mu)/sigma
truncNormalMoment = sigma^2*(1 - uStar*Mills) 

F22 = -( (1/sigma^2) - ( 3 * truncNormalMoment / sigma^4 ) +
           ( Mills^2 + uStar * Mills ) * ( ( crit - mu ) / sigma^2 )^2 -
           Mills * ( 2 * (crit - mu) / sigma^3 ) )
expect_equal(F22, myF[2,2])

# another way to check: take the innards of H22 function from Deriv and replace the single (.x - mu)^2 term with its expectation
H22_expect = function (mu, sigma, .crit) 
{
  
  Mills = mills(params = params, .crit = crit)
  uStar = (crit - mu)/sigma
  
  # Zhou's second moment
  secondMoment = sigma^2*(1 - uStar*Mills) 
  .e1 <- sigma^2
  .e2 <- (2 * .e1)^2
  .e3 <- dnorm(.crit, mu, sigma)
  .e4 <- mu - .crit
  .e5 <- pnorm(.crit, mu, sigma, TRUE)
  
  # straight from H22():
  # 4 * ((.x - mu)^2 * (1 - 16 * (sigma^4/.e2))/.e2) -
  #   ((((.crit - mu)/sigma)^2 - (2 + .e3 * .e4/.e5)) * .e3 * .e4/.e5 - 
  #                                                       1)/.e1
  
  4 * ( secondMoment * (1 - 16 * (sigma^4/.e2))/.e2) -
    ((((.crit - mu)/sigma)^2 - (2 + .e3 * .e4/.e5)) * .e3 * .e4/.e5 - 
       1)/.e1
}

expect_equal( -H22_expect(mu, sigma, crit), myF[2,2] )




# entry 12
F12 = -( ( -2 * 0 / sigma^3 ) +
           ( Mills^2 + uStar * Mills ) * ( ( crit - mu ) / sigma^3 ) -
           Mills / sigma^2 )

expect_equal(F12, myF[1,2], tol = 0.0000000001)
expect_equal(F12, myF[2,1], tol = 0.0000000001)


# another way to check: take the innards of H22 function from Deriv and replace the single (.x - mu)^2 term with its expectation
H12_expect = function (mu, sigma, .crit) 
{
  .e1 <- dnorm(.crit, mu, sigma)
  .e2 <- pnorm(.crit, mu, sigma, TRUE)
  
  # straight from H22
  # ((((.crit - mu)/sigma)^2 - (1 + .e1 * (mu - .crit)/.e2)) * 
  #     .e1/.e2 - 2 * ((.x - mu)/sigma^2))/sigma
  
  ((((.crit - mu)/sigma)^2 - (1 + .e1 * (mu - .crit)/.e2)) * 
      .e1/.e2 - 2 * ((0)/sigma^2))/sigma
}


expect_equal( -H12_expect(mu, sigma, crit), myF[1,2] )



# ~~ Compare my expected Fisher info to Zhou's -----------------------------

# IMPORTANT: Theirs is for a LEFT-truncated normal, so will only
# coincide with mine when crit = mu (by symmetry)
# and I know they have mistake with variance


( ZhouF = expectFisherZhou(params, crit) )
# matches when crit = mu :)


# ~~ Simulation checks -----------------------------


### Get MLEs and see if their variances match the Fisher info
sim.reps = 1000
n = 20

muHat = c()
varHat = c()

for ( i in 1:sim.reps ) {
  
  if ( i %% 25 == 0 ) cat( paste("\n\nStarting rep", i ) )
  
  # right-truncation case
  x = rtrunc(n = n,
             spec = "norm",
             mean = mu,
             sd = sigma,
             b = crit)
  # get MLEs
  mle.fit = mle.tmvnorm( as.matrix(x, ncol = 1), lower=-Inf, upper=crit)
  
  
  # # left-truncation case
  # x = rtrunc(n = n,
  #            spec = "norm",
  #            mean = mu,
  #            sd = sigma,
  #            a = crit)
  # # get MLEs
  # mle.fit = mle.tmvnorm( as.matrix(x, ncol = 1), lower=crit, upper=Inf)
  
  
  mles = coef(mle.fit)
  
  muHat = c(muHat, mles[1])
  varHat = c(varHat, mles[2])
}

sigmaHat = sqrt(varHat)
var(muHat)
var(sigmaHat)

# expected variances of MLEs
solve(ZhouF)/n  # doesn't seem to match at all?
solve(myF)/n



### Simulate and get observed Fisher each time, then take mean to approximate
#  the expected Fisher

# right truncation case
x = rtrunc(n = 1000,
           spec = "norm",
           mean = mu,
           sd = sigma,
           b = crit)


d = data.frame(x=x)

d = d %>% rowwise() %>%
  # this expression flattens the Hessian's entries into columns
  mutate( as.data.frame( matrix( as.numeric( myHessian(params = params,
                                                       .x = x,
                                                       .crit = crit) ), nrow = 1 ) )
  )


-colMeans(d[2:5])  # empirical expected Fisher info
myF  # should match IF you did right truncation

# should match myF entries
-H22_expect(mu, sigma, crit)
-H12_expect(mu, sigma, crit)
#*MYHESSIAN SEEMS ALMOST RIGHT, BUT H12_EXPECT = MYFISHER SEEM WRONG
myF
#bm

# target: 2 for H22


# symmetrical left truncation case: Zhou
# reverse mu and sigma
x = rtrunc(n = 8000,
           spec = "norm",
           mean = mu,
           sd = sigma,
           a = crit)

d = data.frame(x=x)
d = d %>% rowwise() %>%
  # this expression flattens the Hessian's entries into columns
  mutate( J11 = zhou_J11(mu = mu, 
                         sigma = sigma,
                         .x = x,
                         .crit = crit),
          J12 = zhou_J12(mu = mu, 
                         sigma = sigma,
                         .x = x,
                         .crit = crit),
          J22 = zhou_J22(mu = mu, 
                         sigma = sigma,
                         .x = x,
                         .crit = crit) )

-colMeans(d[2:4])  # empirical expected Fisher info
( ZhouF = expectFisherZhou(params, crit) )
# these match! :)



# it really seems like Fisher info should be the same for left-trunc vs. right-trunc
#bm

# check symmetry in general
extrunc(spec = "norm", mean = 1, sd = 0.5, a = 1.5, b = Inf)
extrunc(spec = "norm", mean = -1, sd = 0.5, a = -Inf, b = -1.5)

vartrunc(spec = "norm", mean = 1, sd = 0.5, a = 1.5, b = Inf)
vartrunc(spec = "norm", mean = -1, sd = 0.5, a = -Inf, b = -1.5)



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

# ~~ Check my derivatives OF expected Fisher -----------------------------

# these are for Cordeiro correction

# change parameters so the entries don't evaluate to 0
# set up some values to plug in
params = c(5, 3)
x = 1
crit = 3
mu = params[1]
sigma = params[2]
uStar = (crit-mu)/sigma
uStar2 = (crit-mu)/sigma^2
Mills = mills(params = params, .crit = crit)

# ~~~ Check intermediate pieces --------------------

# check derivative of termA wrt sigma
termA_func = function(.mu, .sigma, .crit) {
  mills = mills(params = c(.mu, .sigma), .crit = .crit)
  uStar = (.crit - .mu)/.sigma
  mills^2 + uStar * mills
}

func = Deriv(termA_func, ".sigma")


termA = termA_func(.mu = mu,
                   .sigma = sigma, 
                   .crit = crit)

# intermediate
# handy derivatives
d1 = uStar2*termA  # derivative of termA wrt sigma
d2 = -uStar2 # derivative of uStar wrt sigma

# intermediates that all match
(2*Mills)*d1 + uStar*d1 + d2*Mills 
(2*Mills + uStar)*d1 + Mills*d2

mine = ( 2*uStar2*Mills + (uStar^2/sigma) )*termA - uStar2*Mills


expect_equal( func( .mu = params[1],
                    .sigma = params[2],
                    .crit = crit ),
              mine )


# ~~~ Check entry 111 ----------------
( mine = expectFisherDerivs( params = params,
                             .x = x, 
                             .crit = crit,
                             .entry = 111 ) )


F11_func = function(.mu, .sigma, .crit) {
  mills = mills(params = c(.mu, .sigma), .crit = .crit)
  uStar = (.crit - .mu)/.sigma
  termA = mills^2 + uStar * mills
  (1/sigma^2) * ( 1 - termA )
}

func = Deriv(F11_func, ".mu")

expect_equal( func( .mu = params[1],
                    .sigma = params[2],
                    .crit = crit),
              mine )


# ~~ Check entry 121 ----------------
( mine = expectFisherDerivs( params = params,
                             .x = x, 
                             .crit = crit,
                             .entry = 121 ) )

F12_func = function(.mu, .sigma, .crit) {
  mills = mills(params = c(.mu, .sigma), .crit = .crit)
  uStar = (.crit - .mu)/.sigma
  termA = .mills^2 + .uStar * .mills
  (mills/.sigma^2) * (3 - mills*uStar - uStar^2)
}

func = Deriv(F12_func, ".mu")

# correct :)
expect_equal( func( .mu = params[1],
                    .sigma = params[2],
                    .crit = crit ),
              mine )

# ~~ Check entry 221 ----------------

( mine = expectFisherDerivs( params = params,
                             .x = x, 
                             .crit = crit,
                             .entry = 221 ) )

F22_func = function(.mu, .sigma, .crit) {
  mills = mills(params = c(.mu, .sigma), .crit = .crit)
  uStar = (.crit - .mu)/.sigma
  termA = mills^2 + uStar * mills
  (1/.sigma^2) * (2 - uStar*mills - 3*mills^2 - mills^2*uStar^2 - mills*uStar^3)
}

func = Deriv(F22_func, ".mu")

expect_equal( func( .mu = params[1],
                    .sigma = params[2],
                    .crit = crit ),
              mine )


# ~~ Check entry 112 ----------------

( mine = expectFisherDerivs( params = params,
                             .x = x, 
                             .crit = crit,
                             .entry = 112 ) )

func = Deriv(F11_func, ".sigma")

expect_equal( func( .mu = params[1],
                    .sigma = params[2],
                    .crit = crit ),
              mine )



# innards of expected Fisher info function for reference:
# mu = params[1]
# sigma = params[2]
# 
# # terms that will show up a lot
# mills = mills(params = params, .crit = .crit)
# uStar = (.crit - mu)/sigma
# termA = mills^2 + uStar * mills
# # termB = 2*mills + uStar
# # termC = termA * termB - mills
# 
# # entry 11
# F11 = (1/sigma^2) * ( 1 - termA )
# 
# # entry 22
# F22 = (1/sigma^2) * (2 - uStar*mills - 3*mills^2 - mills^2*uStar^2 - mills*uStar^3)
# 
# # entry 12=21
# F12 = (mills/sigma^2) * (3 - mills*uStar - uStar^2)
# 
# return( matrix( c(F11, F12, F12, F22), nrow = 2 ) )


# ~ 2021-6-14: Try Godwin/Cordeiro for a single choice of parameters ----------------------


# ~~ Matrix version: -----------------------------------------

# scenario 1: truncate at the mean
# empirical bias (1000 sim reps): 0.980, 0.09
# theory bias: 0.009, -0.19
params = c(0, 1)
mu = params[1]
sigma = params[2]
# BE CAREFUL TO CHOOSE CORRECT SIGN OF CRIT BASED ON TRUNCATION DIRECTION
crit = -0.5  # shouldn't be biased if truncation point is extreme in the appropriate direction
n = 20

# # scenario 2: extremely high trunc point
# #  sanity check because this should be unbiased
# # empirical bias (1000 sim reps): -0.004, -0.04
# # theory bias: 0.000, 0.0125
# params = c(0, 1)
# mu = params[1]
# sigma = params[2]
# crit = 99  # shouldn't be biased if truncation point is super high
# n = 20

# std normal mu bias:
# for crit = -5: 0.02
# crit = 0: 0.17
# crit = 5: 0.002

# empirical bias:
reps = 5000
muHat = c()
sigmaHat = c()
# as sanity check, I confirmed that letting n -> infinity eliminates the bias
for ( i in 1:reps ) {
  
  if ( i %% 25 == 0 ) cat( paste("\n\nStarting rep", i ) )
  
  
  # RIGHT-TRUNC
  # x = rtrunc(n = n,
  #             spec = "norm",
  #             mean = params[1],
  #             sd = params[2],
  #             b = crit)
  # 
  # # get MLEs
  # mle.fit = mle.tmvnorm( as.matrix(x, ncol = 1),
  #                        lower = -Inf,
  #                        upper = crit)
  
  # LEFT-TRUNC
  x = rtrunc(n = n,
             spec = "norm",
             mean = params[1],
             sd = params[2],
             a = crit)
  
  # get MLEs
  mle.fit = mle.tmvnorm( as.matrix(x, ncol = 1),
                         lower = crit,
                         upper = Inf)
  
  mles = coef(mle.fit)
  
  muHat = c(muHat, mles[1])
  sigmaHat = c(sigmaHat, sqrt(mles[2]) )
}

# emprical bias:
mean(muHat)
( empBias = c( mean(muHat) - params[1],
               mean(sigmaHat) - params[2] ) )


# Godwin bias correction
# **but note that in real life, we'd have to plug in the MLEs rather than the true parameters
res = godwinBiasMatrix(params = params,
                       n = n,
                       crit = crit)
( theoryBias = res$bias )

# sanity check
# compare empirical variances of MLEs to Fisher info
# should probably match if untruncated, but not sure whether bias
#  will mess this up in truncated case
res$Kninv[1,1]; var(muHat)
res$Kninv[1,2]; cov(muHat, sigmaHat)
res$Kninv[2,2]; var(sigmaHat)
# comparisons to empirical variances:
#  - truncated at mean, n = 1000 - MATCHES
#  - truncated at mean, n = 20 - DOESN'T MATCH AT ALL
#  - untruncated, n = 20 - MATCHES

#bm: 
# look into:
#  - truncated at mean, n = 1000: Fisher info behaves well but bias correction is wrong



# ~~ Check the bias matrix by implementing each entry manually ------------

# Godwin Eq. (6)

### Bias of mean estimate
p = nrow(res$Kn)  # number of parameters


biasedParamIndex = 1

Kn = n * expectFisherZhou2(mu = mu,
                      sigma = sigma,
                      .crit = crit)

Kninv = solve(Kn)

expect_equal( as.vector(res$Kninv), as.vector(Kninv) )


outerSum = 0
for ( i in 1:p ) {
  
  # k^{si}
  term1 = Kninv[biasedParamIndex,i]
  
  innerSum = 0
  
  for (j in 1:p) {
    for (l in 1:p) {
      
      # k_{ij}^{l}
      term2 = -1 * n * zhou_expect_fisher_derivs( mu = params[1],
                                                  sigma = params[2],
                                                  .crit = crit,
                                                  .entry = as.numeric( paste( c(i,j,l), collapse = "" ) ) ) 
      
      
      
      # k_{ijl} 
      term3 = 0.5 * n * zhou_third_derivs( mu = params[1],
                                           sigma = params[2],
                                           .crit = crit,
                                           .entry = as.numeric( paste( c(i,j,l), collapse = "" ) ) )
      
      # k^{jl}
      term4 = Kninv[j,l]
      
      innerSum = innerSum + (term2 - term3)*term4
      
    }
  }
  
  outerSum = outerSum + (term1*innerSum)
}

outerSum


expect_equal( as.numeric(res$bias[1]), outerSum)

### Bias of sigma estimate


biasedParamIndex = 2

outerSum = 0
for ( i in 1:p ) {
  
  # k^{si}
  term1 = Kninv[biasedParamIndex,i]
  
  innerSum = 0
  
  for (j in 1:p) {
    for (l in 1:p) {
      
      # k_{ij}^{l}
      term2 = -1 * n * zhou_expect_fisher_derivs( mu = params[1],
                                                  sigma = params[2],
                                                  .crit = crit,
                                                  .entry = as.numeric( paste( c(i,j,l), collapse = "" ) ) ) 
      
      
      
      # k_{ijl} 
      term3 = 0.5 * n * zhou_third_derivs( mu = params[1],
                                           sigma = params[2],
                                           .crit = crit,
                                           .entry = as.numeric( paste( c(i,j,l), collapse = "" ) ) )
      
      # k^{jl}
      term4 = Kninv[j,l]
      
      innerSum = innerSum + (term2 - term3)*term4
      
    }
  }
  
  outerSum = outerSum + (term1*innerSum)
}

outerSum


expect_equal( as.numeric(res$bias[2]), outerSum)


# ~~ Estimate only mu: -----------------------------------------

# try comparing to bias correction calculated for only mu
#  when sigma is treated as known

#bm
( a_111 = get_A_entry(params = params,
                      n = n,
                      .crit = crit,
                      .entry = 111) )



( F11 = zhou_F11(mu = mu,
                 sigma = sigma,
                 .crit = crit) )

# bias
( theoryBias = a_111 / (n*F11)^2 )


# empirical bias when treating sigma as fixed:
reps = 3000
muHat = c()
for ( i in 1:reps ) {
  
  if ( i %% 25 == 0 ) cat( paste("\n\nStarting rep", i ) )
  
  
  # RIGHT-TRUNC
  # x = rtrunc(n = n,
  #             spec = "norm",
  #             mean = params[1],
  #             sd = params[2],
  #             b = crit)
  # 
  # # get MLEs
  # mle.fit = mle.tmvnorm( as.matrix(x, ncol = 1),
  #                        lower = -Inf,
  #                        upper = crit)
  
  # LEFT-TRUNC
  x = rtrunc(n = n,
             spec = "norm",
             mean = params[1],
             sd = params[2],
             a = crit)
  
  # get MLEs
  mle.fit = mle.tmvnorm( as.matrix(x, ncol = 1),
                         lower = crit,
                         fixed = list(sigma_1.1 = sigma),
                         upper = Inf)
  
  mles = coef(mle.fit)
  
  muHat = c(muHat, mles[1])
}

# empirical bias:
#mean(muHat)
( empBias = c( mean(muHat) - params[1] ) )
theoryBias

# n=20: matches very closely!




# ~~ Estimate only sigma: -----------------------------------------

# try comparing to bias correction calculated for only sigma
#  when sigma is treated as known


( a_222 = get_A_entry(params = params,
                      n = n,
                      .crit = crit,
                      .entry = 222) )



( F22 = zhou_F22(mu = mu,
                 sigma = sigma,
                 .crit = crit) )

# bias
( theoryBias = a_222 / (n*F22)^2 )


# empirical bias when treating sigma as fixed:
reps = 3000
sigmaHat = c()
for ( i in 1:reps ) {
  
  if ( i %% 25 == 0 ) cat( paste("\n\nStarting rep", i ) )
  
  
  # RIGHT-TRUNC
  # x = rtrunc(n = n,
  #             spec = "norm",
  #             mean = params[1],
  #             sd = params[2],
  #             b = crit)
  # 
  # # get MLEs
  # mle.fit = mle.tmvnorm( as.matrix(x, ncol = 1),
  #                        lower = -Inf,
  #                        upper = crit)
  
  # LEFT-TRUNC
  x = rtrunc(n = n,
             spec = "norm",
             mean = params[1],
             sd = params[2],
             a = crit)
  
  # get MLEs
  mle.fit = mle.tmvnorm( as.matrix(x, ncol = 1),
                         lower = crit,
                         fixed = list(mu_1 = mu),
                         upper = Inf)
  
  mles = coef(mle.fit)
  
  sigmaHat = c(sigmaHat, mles[2])
}

# empirical bias:
#mean(muHat)
( empBias = c( mean(sigmaHat) - params[2] ) )
theoryBias




# # 2021-5-18: TRY COX-SNELL BIAS CORRECTION -----------------------------
# 
# # ~ Toy example from package -----------------------------
# 
# # simple example from mle.tools::coxsnell.bc
# # not the Jeffreys correction, but rather a different correction that also uses 
# #  something related to the Fisher info
# set.seed(1)
# ## Normal distribution
# pdf <- quote(1 / (sqrt(2 * pi) * sigma) * exp(-0.5 / sigma ^ 2 * (x - mu) ^ 2))
# lpdf <- quote(- log(sigma) - 0.5 / sigma ^ 2 * (x - mu) ^ 2)
# x <- rnorm(n = 100, mean = 0.0, sd = 1.0)
# 
# ( mu.hat <- mean(x) )
# ( sigma.hat = sqrt((length(x) - 1) * var(x) / length(x)) )
# 
# coxsnell.bc(density = pdf, logdensity = lpdf, n = length(x), parms = c("mu", "sigma"),
#             mle = c(mu.hat, sigma.hat), lower = '-Inf', upper = 'Inf')
# 
# # ~ Random generates straight from trunc dist -----------------
# 
# # try one example using parameters that yielded very biased estimates in 
# #   2021-5-3 folder, Expt 3.3
# 
# p = data.frame( Mu = 0.1,
#                 T2 = 0.25,
#                 m = 500,
#                 t2w = 0.25,
#                 se = .5,
#                 
#                 Nmax = 1,
#                 hack = "affirm",
#                 
#                 k = 20, 
#                 k.hacked = 0 )
# 
# crit = qnorm(.975)
# 
# 
# # withot bias correction:
# # meanHat = 0.4442739
# # varHat = 3.349863
# 
# sim.reps = 1
# 
# meanHat = c()
# varHat = c()
# tstats = c()
# 
# # set up PDF and score for Cox-Snell thing
# # example for plain-vanilla Normal:
# # pdf <- quote(1 / (sqrt(2 * pi) * sigma) * exp(-0.5 / sigma ^ 2 * (x - mu) ^ 2))
# # lpdf <- quote(- log(sigma) - 0.5 / sigma ^ 2 * (x - mu) ^ 2)
# 
# 
# mean = .Mu / .se,
# #@doesn't yet use the delta-method thing
# sd = sqrt( (1/.se^2) * (.T2 + .t2w + .se^2) )
# 
# 
# myPDF = "1 / (sqrt(2 * pi) * sigma) * exp(-0.5 / sigma ^ 2 * (x - mu) ^ 2) * pnorm( q = crit, mean = mu, sd = sigma)"
# 
# myLPDF = "-log(sigma) - 0.5/sigma^2 * (x - mu)^2 + pnorm( q = crit, mean = mu, sd = sigma, log=TRUE)"
# 
# # ?D:
# # The internal code knows about the arithmetic operators +, -, *, / and ^, and the single-variable functions exp, log, sin, cos, tan, sinh, cosh, sqrt, pnorm, dnorm, asin, acos, atan, gamma, lgamma, digamma and trigamma, as well as psigamma for one or two arguments (***but derivative only with respect to the first***). (Note that only the standard normal distribution is considered.)
# 
# for ( i in 1:sim.reps ) {
#   
#   if ( i %% 25 == 0 ) cat( paste("\n\nStarting rep", i ) )
#   
#   x2 = rtrunc(n = 20,
#               spec = "norm",
#               mean = p$Mu / p$se,
#               sd = sqrt( (1/p$se^2) * (p$T2 + p$t2w + p$se^2) ),
#               b = crit)
#   
#   # get MLEs
#   mle.fit = mle.tmvnorm( as.matrix(x2, ncol = 1), lower=-Inf, upper=crit)
#   mles = coef(mle.fit)
#   
#   # Cox-Snell bias correction
#   #@BE CAREFUL ABOUT SD VS. VAR THING!
#   coxsnell.bc(density = myPDF,
#               logdensity = myLPDF,
#               n = length(x2),
#               parms = c("mu", "sigma"),
#               mle = c(mles[1], mles[2]),
#               lower = '-Inf',
#               upper = 'Inf')
#   
#   
#   
#   meanHat = c(meanHat, mles[1])
#   varHat = c(varHat, mles[2])
#   tstats = c(tstats, x2)
# }
# 
# 
# mean(meanHat); p$Mu/p$se
# mean(varHat); (1/p$se^2) * (p$T2 + p$t2w + p$se^2)
# 
# # should be normal
# hist(meanHat)
# 
# # trunc normal
# hist(tstats)
# 
# 
