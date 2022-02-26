
# Try to set up an MCMC algorithm to sample from posterior with p-hacking :)

library(here)
library(crch)
library(truncdist)
library(testthat)
library(dplyr)
library(ggplot2)
library(doParallel)
library(foreach)
library(rjags)

setwd(here("2021-4-3 Bayesian experiments"))
source("helper.R")



# simulate meta-analysis of hacked studies
# set parameters globally so that the Metropolis sampler can use the same parameters
Z = 0
hack = "affirm"
N.max = 20
tau = sqrt(.5)
tauw = 0
k = 500

d = sim_meta(k = k,
             N.max = N.max,
             Z = Z,
             T2 = tau^2, 
             t2 = tauw^2,
             hack = hack,
             NSpreference = "last")
hist(d$Zhati, breaks = 20)


# quick test
# plug in true values
goodPost = log_joint_post( .Z = Z,
                #.Zi = d$Zi,
                .Zi = candZi,
                
                .Ni = d$Ni,
                .Zhat = d$Zhati,
                .tau = tau,
                .tauw = 0,

                hack = "affirm",
                N.max = 20 )

# vs. some incorrect parameter values
badPost = log_joint_post( .Z = 0.5,
                #.Zi = d$Zi,
                .Zi = as.numeric( triedPars$Zi[i,] ),
                #.Zi = rep(10, k),
                
                .Ni = d$Ni,
                .Zhat = d$Zhati,
                .tau = tau,
                .tauw = 0,
                
                hack = "affirm",
                N.max = 20 )

# want this > 0
goodPost - badPost

# when I use the CANDIDATE Zi values, the bad posterior (bottom) is better 
# but when I use the TRUE Zi values, the new posterior (top) is better, as it should be
# think issue might be that the initial Zi values are all exactly at their mean, which makes the posterior super good
#@do I need to draw these intermediates for the first iterate as well?

# LOOK AT JOINT POST FOR VARIOUS CHOICES OF PARAMS -------------------------------------------------------------------


# ~ Try different values of Z -------------------------------------------------------------------
dp = data.frame( .Z = seq( -5, 5, 0.1 ) )

dp = dp %>% rowwise() %>%
  # for now, give it true values of the params other than the one we're manipulating
  mutate( post = log_joint_post( .Z = .Z,
                                 .Zi = d$Zi,
                                 .Ni = d$Ni,
                                 .Zhat = d$Zhati,
                                 .tau = .tau,
                                 .tauw = 0,
                                 
                                 hack = "affirm",
                                 N.max = 20 ) )


ggplot( data = dp,
        aes(x = .tau,
            y = post) ) + 
  # true parameter (should be the max)
  geom_vline(xintercept = tau, color = "red") +
  geom_line() +
  theme_classic()

# should be close
tau; dp$.tau[ which.max(dp$post) ]



# ~ Try different values of tau -------------------------------------------------------------------
dp = data.frame( .tau = seq( 0, 2, 0.1 ) )

dp = dp %>% rowwise() %>%
  # for now, give it true values of the params other than the one we're manipulating
  mutate( post = log_joint_post( .Z = Z,
                                 .Zi = d$Zi,
                                 .Ni = d$Ni,
                                 .Zhat = d$Zhati,
                                 .tau = .tau,
                                 .tauw = 0,
                                 
                                 hack = "affirm",
                                 N.max = 20 ) )


ggplot( data = dp,
        aes(x = .tau,
            y = post) ) + 
  # true Z (should be the max)
  geom_vline(xintercept = , color = "red") +
  geom_line() +
  theme_classic()


# METROPOLIS -------------------------------------------------------------------

# number of sampling iterates
B = 500

init_sampler()
triedPars$scalars[1,]
triedPars$Zi[1,]
triedPars$Ni[1,]

res = run_sampler()

# look at last row of scalar parameters
res$scalars[B,]
plot( res$scalars$Z, type = "line" )
#@even when I give it tau, tauw, chain still gets stuck

# look at acceptance probabilities and ratios
mean(res$accept)

summary(res$r)
table(res$r < 0.0000000001)

# not working - see #bm in helper.R
# another way to debug: assume true values of the taus rather than estimating them

# TRY IN JAGS -------------------------------------------------------------------



model_string <- "model{

    # Likelihood
    for(i in 1:k){
      
      #@how to determine which things go in lkl vs in priors?
      # temp only
      # no truncation or anything
      Zhat[i] ~ dnorm( Z[i], 0^2 + 1 )
      
      # Zhat[i] ~ Zhat_pdf( .Zhat , 
      #                           .Zi, 
      #                           .Ni,
      #                           .t2w,
      #                           
      #                           N.max,
      #                           hack = 'affirm')
    
    }
  
  # intermediates
    for(i in 1:k){
      Z[i] ~ dnorm( Z, 0.5^2 ) # normal random effects
      N[i] ~ dgeom( prob = power[i] ) # not truncated
    }
  
    # Prior for Z
    Z ~ dnorm(0,0.0001)
  
    # # Prior for tau. Half-Cauchy is a truncated standard t w/ 1 degree of freedom
    # tau ~ dt(0,1,1) T(0,)
    # tauw ~ dt(0,1,1) T(0,)
  

    # Update power_i
    for(i in 1:k){
         SD[i] = sqrt(0^2 + 1)  # marginal SD of the Z-statistics
          crit = qnorm(.975)
          Psig.pos[i] = 1 - pnorm( (crit - Z[i]) / SD[i] )
          power[i] = Psig.pos

    }
  
    }"

# Compile the model in JAGS
model <- rjags::jags.model(textConnection(model_string), 
                           data = list(Zhat = d$Zhati, N.max = N.max, k = k),
                           inits = list(Z = 0.5),
                           n.chains=1,
                           quiet=TRUE)
#@hitting syntax errors in the dgeom line and not sure why

# # Draw samples
# # First use function 'update' to draw 10,000 warm-up (burnin) samples
# stats::update(model, n.iter=burn, progress.bar="none") # Burnin for 10,000 samples
# 
# # Next use the function 'coda.samples' to produce the next 10,000
# # samples which will ultimately used to approximate the posterior.
# full.samples <- rjags::coda.samples(model, 
#                                     variable.names=c("mu","theta","tau","rho","gamma0","gamma1"), 
#                                     n.iter=nmc, progress.bar="none")
# 
# 
# 
# 





