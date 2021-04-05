
# Try to set up an MCMC algorithm to sample from posterior with p-hacking :)

library(here)
library(crch)
library(truncdist)
library(testthat)
library(dplyr)
library(ggplot2)

setwd(here("2021-4-3 Bayesian experiments"))
source("helper.R")



# simulate meta-analysis of hacked studies
Z = 0.5
d = sim_meta(k = 20,
             N.max = 20,
             Z = Z,
             T2 = 0.5,
             t2 = 0,
             hack = "affirm",
             NSpreference = "last")
hist(d$Zhati, breaks = 20)


# quick test
# plug in true values
log_joint_post( .Z = Z,
                .Zi = d$Zi,
                .Ni = d$Ni,
                .Zhat = d$Zhati,
                .T2 = 0.5,
                .t2w = 0,
                
                hack = "affirm",
                N.max = 20 )


# LOOK AT JOINT POST FOR VARIOUS CHOICES OF PARAMS -------------------------------------------------------------------

# for now, 
dp = data.frame( .Z = seq( -5, 5, 0.1 ) )

dp = dp %>% rowwise() %>%
  # for now, give it true values of the params other than the one we're manipulating
  mutate( post = log_joint_post( .Z = .Z,
                                 .Zi = d$Zi,
                                 .Ni = d$Ni,
                                 .Zhat = d$Zhati,
                                 .T2 = 0.5,
                                 .t2w = 0,
                                 
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

# for tracking parameters that we've tried
B = 10

scalarPars = c("Z", "T2", "t2w")

vectorPars = c("Zi", "Ni")


# initialized with starting values
( triedPars = list( scalars = data.frame( matrix(NA,
                                                 nrow = B,
                                                 ncol = length(scalarPars) ) ),
                    
                    # vector params: a column for each study in meta-analysis
                    Zi = data.frame( matrix(NA,
                                            nrow = B,
                                            ncol = length(Zhats) ) ) ) )

names(triedPars$scalars) = scalarPars








