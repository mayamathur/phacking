
# Try to set up an MCMC algorithm to sample from posterior with p-hacking :)

library(here)
library(crch)
library(truncdist)
library(testthat)
library(dplyr)
library(ggplot2)
library(doParallel)
library(foreach)

setwd(here("2021-4-3 Bayesian experiments"))
source("helper.R")



# simulate meta-analysis of hacked studies
# set parameters globally so that the Metropolis sampler can use the same parameters
Z = 0
hack = "affirm"
N.max = 20

d = sim_meta(k = 20,
             N.max = N.max,
             Z = Z,
             T2 = 0.5,  # same as tau^2
             t2 = 0,
             hack = hack,
             NSpreference = "last")
hist(d$Zhati, breaks = 20)


# quick test
# plug in true values
log_joint_post( .Z = Z,
                .Zi = d$Zi,
                .Ni = d$Ni,
                .Zhat = d$Zhati,
                .tau = sqrt(0.5),
                .tauw = 0,
                
                hack = "affirm",
                N.max = 20 )


# LOOK AT JOINT POST FOR VARIOUS CHOICES OF PARAMS -------------------------------------------------------------------

dp = data.frame( .Z = seq( -5, 5, 0.1 ) )

dp = dp %>% rowwise() %>%
  # for now, give it true values of the params other than the one we're manipulating
  mutate( post = log_joint_post( .Z = .Z,
                                 .Zi = d$Zi,
                                 .Ni = d$Ni,
                                 .Zhat = d$Zhati,
                                 .tau = sqrt(0.5),
                                 .tauw = 0,
                                 
                                 hack = "affirm",
                                 N.max = 20 ) )


ggplot( data = dp,
        aes(x = .Z,
            y = post) ) + 
  # true Z (should be the max)
  geom_vline(xintercept = Z, color = "red") +
  geom_line() +
  theme_classic()


# looks good :)

#bm: Ready to try writing the sampler using this joint posterior and the code in HW5 :)


# METROPOLIS -------------------------------------------------------------------

# number of sampling iterates
B = 1000

scalarPars = c("Z", "tau", "tauw")

vectorPars = c("Zi", "Ni")

# for tracking parameters that we've tried
( triedPars = list( scalars = data.frame( matrix(NA,
                                                 nrow = B,
                                                 ncol = length(scalarPars) ) ),
                    
                    # vector params stored as df with a column for each study in meta-analysis
                    Zi = data.frame( matrix(NA,
                                            nrow = B,
                                            ncol = nrow(d) ) ),
                    
                    Ni = data.frame( matrix(NA,
                                            nrow = B,
                                            ncol = nrow(d) ) ) ) )

names(triedPars$scalars) = scalarPars

# set starting values
triedPars$scalars[1,] = c(0.5, 0.5, 0.5)
triedPars$Zi[1,] = 0.5
triedPars$Ni[1,] = 10

# SD of jumping dist for parameters
jump.SD = 0.5

# number studies 
k = nrow(d)



# at each step, compute log of joint posterior, jumping dist, and ratio 
# jumping dist must be symmetric
for (i in 1:B) {
  

    ### Draw "Primary" Parameters
    # i.e., from marginals
    # last iteration's draws for scalars
    lastScalarPars = triedPars$scalars[i,]
  
  # generate new candidates
  # jumping dist is independent normals centered around last iteration's draw
  
  candZ = rnorm( n = 1,
                 mean = lastScalarPars$Z,
                 sd = jump.SD )
  
  # draw variance components from truncated normal
  # asymmetric proposal dist is fine if you use Metropolis-Hastings rather than Metropolis
  # see Gelman pg 279
  # the difference is the way the acceptance ratio is calculated later
  candtau = rtrunc( n = 1,
                    spec = "norm",
                    mean = lastScalarPars$tau,
                    sd = jump.SD,
                    a = 0)
  
  candtauw = rtrunc( n = 1,
                     spec = "norm",
                     mean = lastScalarPars$tauw,
                     sd = jump.SD,
                     a = 0 )
  
  # for the acceptance ratio
  # log-prob of making the jump we actually made
  termA.J = log( dnorm( candZ,
                        mean = lastScalarPars$Z,
                        sd = jump.SD ) ) +
    
    log( dtrunc( candtau,
                 spec = "norm",
                 mean = lastScalarPars$tau,
                 sd = jump.SD,
                 a = 0) ) + 
    
    log( dtrunc( candtauw,
                 spec = "norm",
                 mean = lastScalarPars$tauw,
                 sd = jump.SD,
                 a = 0) )
  
  #...vs log-prob of making the opposite jump
  termB.J = log( dnorm( lastScalarPars$Z,
                        mean = candZ,
                        sd = jump.SD ) ) +
    
    log( dtrunc( lastScalarPars$tau,
                 spec = "norm",
                 mean = candtau,
                 sd = jump.SD,
                 a = 0) ) + 
    
    log( dtrunc( lastScalarPars$tauw,
                 spec = "norm",
                 mean = candtauw,
                 sd = jump.SD,
                 a = 0) )
  
  
  ### Draw "Intermediate" Parameters (Zi, Ni) Conditional on Primaries
  # random intercept for each study
  ( candZi = rnorm( n = k,
                    mean = candZ,
                    sd = candtau ) )
  
  # Ni
  SDi = sqrt(candtauw^2 + 1)  # marginal SD of the Z-statistics
  crit = qnorm(.975)
  Psig.pos = 1 - pnorm( (crit - candZi) / SDi )
  Psig.neg = pnorm( (-crit - candZi) / SDi )
  
  if ( hack == "signif" ) power_i = Psig.pos + Psig.neg
  if ( hack == "affirm" ) power_i = Psig.pos
  
  # is the final draw affirmative or not?
  if ( hack == "signif" ) favored = abs(candZi) > qnorm(0.975)
  if ( hack == "affirm" ) favored = candZi > qnorm(0.975)
  
  
  # truncated geometric (max = maxN)
  # dtrunc isn't vectorized
  candNi = vapply( X = seq( 1:k ),
                   FUN = function(i) rtrunc( n = 1,
                                             spec = "geom",
                                             prob = power_i[i],
                                             b = N.max ),
                   FUN.VALUE = -99)
  #@ to avoid impossible choices of Ni, anytime the final draw was not favored, 
  # Ni should be equal to N.max
  candNi[ favored == FALSE ] = N.max
  
  
  # calculate new posterior using the candidates
  # use log densities to avoid craziness
  newPost = log_joint_post( .Z = candZ,
                            .tau = candtau,
                            .tauw = candtauw,
                            .Zi = candZi,
                            .Ni = candNi,
                            
                            # **here's where the observed data come in
                            .Zhat = d$Zhati,
                            
                            
                            
                            hack = "affirm",
                            N.max = 20 )
  
  #@ is this right? force acceptance of first set of params
  if ( i == 1 ) oldPost = -9999999
  
  # ratio of new vs. old parameter values' log-posteriors
  # with Metropolis-Hastings adjustment for asymmetric proposal
  # Gelamn pg 279
  termA = newPost - termA.J
  termB = oldPost - termB.J
  
  # r > 1 is good
  ( r = exp( termA - termB ) )
  

  # probability of accepting this step
  ( p = ifelse(r<1, r, 1) )
  accept = rbinom(n=1, size=1, p=p)
  
  if (accept) {  # automatically accept improvement steps
    
    triedPars$scalars[i+1, "Z"] = candZ
    triedPars$scalars[i+1, "tau"] = candtau
    triedPars$scalars[i+1, "tauw"] = candtauw 
    
    triedPars$Zi[i+1,] = candZi
    triedPars$Ni[i+1,] = candNi

  } else {
    triedPars$scalars[i+1,] = triedPars$scalars[i,]
    triedPars$Zi[i+1,] =  triedPars$Zi[i,]
    triedPars$Ni[i+1,] = triedPars$Ni[i,]
  }
}





triedPars$scalars


plot( triedPars$scalars$Z, type = "line" )




