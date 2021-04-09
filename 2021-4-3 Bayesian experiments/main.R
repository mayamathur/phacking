
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
tau = sqrt(.5)
k = 500

d = sim_meta(k = k,
             N.max = N.max,
             Z = Z,
             T2 = tau^2, 
             t2 = 0,
             hack = hack,
             NSpreference = "last")
hist(d$Zhati, breaks = 20)


# # quick test
# # plug in true values
# log_joint_post( .Z = Z,
#                 .Zi = d$Zi,
#                 .Ni = d$Ni,
#                 .Zhat = d$Zhati,
#                 .tau = tau,
#                 .tauw = 0,
#                 
#                 hack = "affirm",
#                 N.max = 20 )




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
res = run_sampler()

# look at last row of scalar parameters
res$scalars[B,]
plot( res$scalars$Z, type = "line" )

# look at acceptance probabilities and ratios
mean(res$accept)

summary(res$r)
table(res$r < 0.0000000001)



# uh-oh...now it stays stuck for long periods without changing



