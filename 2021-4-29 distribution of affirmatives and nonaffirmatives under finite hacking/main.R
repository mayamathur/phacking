
# Goal: I think that the distribution of the affirmative and nonaffirmative results
#  is still truncated-t even when the hacking is finite rather than infinite. 
#  That would really help with doing 2-step bias corrections without assuming 
#  an N_max. Let's check if this is true. :)

# Conclusion: YES. This does hold! 


# PRELIMINARIES ------------------------------

# borrow helper fns from previous sims
library(here)
setwd( here("2021-4-20 simple sims with partial try-to-Nmax hacking") )
source("helper.R")

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

# SIMULATE DATA ------------------------------

# parameters also needed later
p = data.frame( Mu = 1,
                T2 = 0.1,
                m = 50,
                t2w = .1,
                se = 1,
                
                Nmax = 1,
                hack = "affirm",
                
                k = 2000,
                k.hacked = 2000 )


# simulate a huge dataset, all finitely hacked
d = sim_meta(Nmax = p$Nmax,
             Mu = p$Mu,
             T2 = p$T2,
             m = p$m,
             t2w = p$t2w,
             se = p$se,
             hack = p$hack,
             return.only.published = FALSE,
             
             k = p$k,
             k.hacked = p$k.hacked )


# dataset of only published results
dp = d %>% filter(Di == 1)
dim(dp)


# dataset of only published, hacked results
# (here, all results are hacked)
dph = dp %>% filter(hack == "affirm")

# # save for later
# setwd( here("2021-4-29 distribution of affirmatives and nonaffirmatives under finite hacking") )
# # add in the parameters
# fwrite( cbind( dph, p ), "sim_meta_Nmax1_hugek.csv")


# read back in
setwd( here("2021-4-29 distribution of affirmatives and nonaffirmatives under finite hacking") )
dph = fread("sim_meta_Nmax1_hugek.csv")

# since Nmax is really large in this dataset, we always get an affirmative result to report


# DATA SIMULATION SANITY CHECKS ------------------------------

summary(d$N)

table(d$hack)

nrow(d)
nrow(dp)

length(unique(d$study))


# all results, sorted by hacking status and publication status
t = d %>%
  group_by(hack, Di) %>%
  summarise( n(),
             k = length(unique(study)),
             mean(affirm),
             mean(mui),
             var(mui),
             mean(yi),
             var(yi))
t



# CHECK AGAINST TRUNCATED DISTRIBUTION ------------------------------


# # distribution I THINK the results should follow
# #  (i.e., using the real T2, Mu, and se rather than sample estimates)
# dph = dph %>%
#   rowwise() %>%
#   mutate( hackedExpTrue =  extrunc( spec = "t",
#                                     ncp = p$Mu / sqrt( p$T2 + p$t2w + viTrue ),
#                                     df = m-1,
#                                     a = tcrit ),
#           hackedMargVarTrue = vartrunc( spec = "t",
#                                       ncp = p$Mu / sqrt( p$T2 + p$t2w + viTrue ),
#                                       df = m-1,
#                                       a = tcrit ) ) 


# calculate real moments of truncated distribution
# because studies all have same viTrue and m, all will have same trunc dist
hackedExpTrue = extrunc( spec = "t",
                         ncp = p$Mu / sqrt( p$T2 + p$t2w + unique( d$viTrue ) ),
                         df = p$m-1,
                         a = unique( dph$tcrit ) )

hackedMargVarTrue = vartrunc( spec = "t",
                              ncp = p$Mu / sqrt( p$T2 + p$t2w + unique( d$viTrue ) ),
                              df = p$m-1,
                              a = unique( dph$tcrit ) )


# should be very close
# when Nmax = 1 (k=2000 but 345 observed), it's 2.561771 vs. 2.542194
# when Nmax = 2 (k=2000 but 670 observed), it's 2.561771 vs. 2.548238
# when Nmax = 200 (k=2000), it's 2.561771 vs. 2.527792
hackedExpTrue; t$`mean(yi)`[t$Di == TRUE]



# should be very close
# when Nmax = 1 (k=2000 but 345 observed), it's 0.2349009 vs. 0.266968
# when Nmax = 2 (k=2000), it's 0.2349009 vs. 0.2866817
# when Nmax = 200 (k=2000), it's 0.2349009 vs. 0.2623937
hackedMargVarTrue; t$`var(yi)`[t$Di == TRUE]
# seems pretty close

# so the Nmax doesn't seem to matter, but either way the expectation is a little off

# SANITY CHECK: SIMULATE DIRECTLY FROM TRUNCATED DISTRIBUTION ------------------------------

# this one looks just like the trunc-t, as expected

# this is a sanity check for the above
x = rtrunc( n=2000,
            spec = "t",
            ncp = p$Mu / sqrt( p$T2 + p$t2w + p$se^2 ),
            df = p$m-1,
            a = unique( dph$tcrit ) )

# yes, these match almost exactly
# so it doesn't seem to be extrunc/vartrunc's fault
hackedExpTrue; mean(x)
hackedMargVarTrue; var(x)


# SANITY CHECK: SIMULATE UNHACKED STUDIES AND LOOK AT DIST OF AFFIRMATIVES ------------------------------

# this one looks just like the trunc-t, as expected

# this one takes ~10 min
# # simulate a huge dataset, all unhacked
# dUH = sim_meta(Nmax = p$Nmax,
#              Mu = p$Mu,
#              T2 = p$T2,
#              m = p$m,
#              t2w = p$t2w,
#              se = p$se,
#              hack = p$hack,
#              return.only.published = FALSE,
#              
#              k = p$k,
#              k.hacked = 0 )  # none are hacked


# # save for later
# setwd( here("2021-4-29 distribution of affirmatives and nonaffirmatives under finite hacking") )
# # add in the parameters
# temp = cbind( dUH, p )
# temp$k.hacked = 0
# fwrite( temp, "sim_meta_no_hack_hugek.csv")

# read back in
setwd( here("2021-4-29 distribution of affirmatives and nonaffirmatives under finite hacking") )
dUH = fread("sim_meta_no_hack_hugek.csv")


# matches what I expected
mean(dUH$yi[ dUH$affirm == TRUE])

# also pretty much matches
# at least it's not too small as we were seeing before
# in this case, published ones are the just the LAST
#  draws of each study set
mean(dUH$yi[ dUH$affirm == TRUE & dUH$Di == 1])



# TRY TO DEBUG ------------------------------

# only 1 draw per study set
# so basically this just simulates 8000 times from N(mu, T2 + t2w + se^2)
# and then labels every affirmative as Di =  1

for ( i in 1:8000 ) {
  newRow = sim_one_study_set(Nmax = 1,
                             Mu = p$Mu,
                             T2 = p$T2,
                             m = p$m,
                             t2w = p$t2w,
                             se = p$se,
                             hack = "affirm",
                             return.only.published = FALSE)
  
  if ( i == 1 ) fake = newRow else fake = rbind(fake, newRow)
}


dim(fake)
expect_equal(  )

#bm

# THIS REPRODUCES THE PROBLEM
hackedExpTrue; mean(fake$yi[ fake$Di == 1])
hackedMargVarTrue; var(fake$yi[ fake$Di == 1])

# vs. all affirmatives
# which is exactly the same thing because Nmax = 1
hackedExpTrue; mean(fake$yi[ fake$affirm == 1])
hackedMargVarTrue; var(fake$yi[ fake$affirm == 1])


# STILL HAS PROBLEM, EVEN WITH HACK = "NO": ???
# this basically simulates 8000 times from N(mu, T2 + t2w + se^2)
#  and then marks every study as Di=1
for ( i in 1:8000 ) {
  newRow = sim_one_study_set(Nmax = 1,
                             Mu = p$Mu,
                             T2 = p$T2,
                             m = p$m,
                             t2w = p$t2w,
                             se = p$se,
                             hack = "no",  # **only change from the above
                             return.only.published = FALSE)
  
  if ( i == 1 ) fake2 = newRow else fake2 = rbind(fake2, newRow)
}

dim(fake2) 
table(fake2$Di)  # all should be published

# STILL HAS THE PROBLEM!! 
hackedExpTrue; mean(fake2$yi[ fake2$affirm == 1])
hackedMargVarTrue; var(fake2$yi[ fake2$affirm == 1])



# EARLIER VERSION:
# simulate a huge dataset, all unhacked
fake3 = sim_meta(Nmax = 1,
                 Mu = p$Mu,
                 T2 = p$T2,
                 m = p$m,
                 t2w = p$t2w,
                 se = p$se,
                 hack = "affirm",
                 return.only.published = FALSE,
                 
                 k = 8000,
                 k.hacked = 0 )  # none are hacked

dim(fake3)
table(fake3$Di)  # all should be published

extrunc( spec = "t",
         ncp = p$Mu / sqrt( p$T2 + p$t2w + unique( fake3$viTrue ) ),
         df = p$m-1,
         a = unique( fake3$tcrit ) )

vartrunc( spec = "t",
          ncp = p$Mu / sqrt( p$T2 + p$t2w + unique( fake3$viTrue ) ),
          df = p$m-1,
          a = unique( fake3$tcrit ) )


# STILL HAS THE PROBLEM!! ??????
mean(fake3$yi[ fake3$affirm == 1])
var(fake3$yi[ fake3$affirm == 1])

#bm


# simplify even further
# just call make_one_draw directly
for ( i in 1:8000 ) {
  
  mui = p$Mu + rnorm(mean = 0,
                     sd = sqrt(p$T2),
                     n = 1)
  
  newRow = make_one_draw( mui = mui, 
                          sd.y = p$se * sqrt(p$m),
                          Mu = p$Mu,
                          T2 = p$T2,
                          m = p$m,
                          t2w = p$t2w,
                          se = p$se )
  
  if ( i == 1 ) fake4 = newRow else fake4 = rbind(fake4, newRow)
}

dim(fake4) 

# STILL HAS THE PROBLEM!! WTF??
fake4$affirm = fake4$yi > 0 & fake4$pval < 0.05
mean(fake4$yi[ fake4$affirm == 1])
var(fake4$yi[ fake4$affirm == 1])

extrunc( spec = "t",
         ncp = p$Mu / sqrt( p$T2 + p$t2w + unique(p$se)^2 ),
         df = p$m-1,
         a = unique( fake3$tcrit ) )

vartrunc( spec = "t",
         ncp = p$Mu / sqrt( p$T2 + p$t2w + unique(p$se)^2 ),
         df = p$m-1,
         a = unique( fake3$tcrit ) )


# EVEN SIMPLER AGAIN: SIMULATE FROM SCRATCH
for ( i in 1:8000 ) {
  
  mui = p$Mu + rnorm(mean = 0,
                     sd = sqrt(p$T2),
                     n = 1)
  
  sd.y = p$se * sqrt(p$m)
  
  
  # true mean for draw n (based on within-study heterogeneity)
  muin = rnorm(mean = mui,
               sd = sqrt(p$t2w),
               n = 1)
  
  # draw subject-level data from this study's population effect
  y = rnorm( mean = muin,
             sd = sd.y,
             n = p$m)
  
  
  test = t.test(y,
                alternative = "two.sided")
  
  pval = test$p.value
  tstat = test$statistic
  vi = test$stderr^2  # ESTIMATED variance
  
  newRow = data.frame(pval = pval,
                     tstat = tstat,
                     tcrit = qt(0.975, df = p$m-1),
                     mui = mui,
                     muin = muin,
                     yi = mean(y),
                     vi = vi,
                     viTrue = sd.y^2 / p$m,  # true variance; will equal p$se^2
                     m = p$m ) 
  
  if ( i == 1 ) fake5 = newRow else fake5 = rbind(fake5, newRow)
  
}

# if (hack == "signif") success = (pval < 0.05)
# if (hack == "affirm") success = (pval < 0.05 & mean(y) > 0)

# still wrong, especially the variance...
fake5$affirm = fake5$yi > 0 & fake5$pval < 0.05
expect_equal( fake5$affirm, fake5$tstat > qt(0.975, df = p$m-1))
mean(fake5$yi[ fake5$affirm == 1])
var(fake5$yi[ fake5$affirm == 1])

# marginal moments are fine
mean(fake5$yi)
var(fake5$yi); p$T2 + p$t2w + p$se^2

# could I have one of the arguments to extrunc wrong?
#**I don't think it's the package's fault since it was correct when I simulated directly from trunc dist

vartrunc( spec = "t",
          ncp = p$Mu / sqrt( p$T2 + p$t2w + unique(p$se)^2 ),
          df = p$m-1,
          a = unique( fake3$tcrit ) )

#bm

# this is exactly right
mean(fake3$vi)

# BM: WHAT A MYSTERY. THIS LAST EXAMPLE "EVEN SIMPLER AGAIN: SIMULATE FROM SCRATCH" IS 
# ESPECIALLY WEIRD BECAUSE I'M JUST MAKING T-DISTRIBUTED SAMPLES AND LOOKING AT THE MEAN AND 
#  VARIANCE OF THE AFFIRMATIVE ONES. THE MARGINAL MOMENTS (NOT CONDIITONING ON AFFIRM STATUS)
#  ARE WHAT I EXPECTED TO SEE, BUT THE CONDITIONAL ONES ARE WRONG.




