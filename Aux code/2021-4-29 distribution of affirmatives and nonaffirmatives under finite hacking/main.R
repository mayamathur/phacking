

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
                T2 = 0.25,
                m = 500,
                t2w = .25,
                se = .5,
                
                Nmax = 200,
                hack = "affirm",
                
                k = 10000,
                k.hacked = 10000 )


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

summary(d$viTrue)
summary(d$vi)

# add in the parameters that aren't already in dataset
shortParams = p[ , !names(p) %in% names(d) ]
d = cbind( d, shortParams )

# dataset of only published results
dp = d %>% filter(Di == 1)
dim(dp)


# dataset of only published, hacked results
# (here, all results are hacked)
dph = dp %>% filter(hack == "affirm")

# # save for later
setwd( here("2021-4-29 distribution of affirmatives and nonaffirmatives under finite hacking") )
fwrite( d, "sim_meta_Nmax1_hugek_hugem_allhacked_all_studies.csv")


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



# ~ EARLIER VERSION: -------------------------------
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


# EVEN SIMPLER AGAIN: SIMULATE FROM SCRATCH  -------------------------------
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
          a = unique( fake5$tcrit ) )


# this is exactly right
mean(fake3$vi)

# BM: WHAT A MYSTERY. THIS LAST EXAMPLE "EVEN SIMPLER AGAIN: SIMULATE FROM SCRATCH" IS 
# ESPECIALLY WEIRD BECAUSE I'M JUST MAKING T-DISTRIBUTED SAMPLES AND LOOKING AT THE MEAN AND 
#  VARIANCE OF THE AFFIRMATIVE ONES. THE MARGINAL MOMENTS (NOT CONDIITONING ON AFFIRM STATUS)
#  ARE WHAT I EXPECTED TO SEE, BUT THE CONDITIONAL ONES ARE WRONG.



# # recalculate affirmative indicator, accounting for T2 + t2w
# fake5$tstat2 = fake5$yi / sqrt( p$T2 + p$t2w + unique(p$se)^2 )
# fake5$pval2 = 2 * ( 1 - pt( abs( fake5$tstat2 ), df = p$m - 1 ) )
# 
# fake5$yi / sqrt( p$T2 + p$t2w + unique(p$se)^2
# 
# 
# mean(fake5$yi[ fake5$pval2 < 0.05 & fake5$yi > 0 ])
# var(fake5$yi[ fake5$affirm == 1])


# TRUNCATE AT DIFFERENT VALUE, NOT DEPENDING ON P-VALUE -------------------------------

# truncate at a different value, not depending on p-value
mean(fake5$yi[ fake5$pval2 < 0.05 & fake5$yi > 0 ])

vartrunc( spec = "t",
          ncp = p$Mu / sqrt( p$T2 + p$t2w + unique(p$se)^2 ),
          df = p$m-1,
          a = 1 )

var(fake5$yi[ fake5$yi > 1 ])




# still wrong!! 
# **so it's nothing to do with the t-test itself


# NO HETEROGENEITY -------------------------------

for ( i in 1:8000 ) {
  
  sd.y = p$se * sqrt(p$m)

  # draw subject-level data from this study's population effect
  y = rnorm( mean = p$Mu,
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


# **needed to calculate t-stat with vi instead of true p$se^2
fake5$tstat2 = fake5$yi / sqrt(fake5$vi)

fake5$pval2 = 2 * ( 1 - pt( abs( fake5$tstat2 ), df = p$m - 1 ) )



extrunc( spec = "t",
          ncp = p$Mu / unique(p$se),
          df = p$m-1,
          a = 1 )
mean(fake5$tstat2[ fake5$tstat2 > 1 ])

vartrunc( spec = "t",
          ncp = p$Mu / unique(p$se),
          df = p$m-1,
          a = 1 )

var(fake5$tstat2[ fake5$tstat2 > 1 ])

# WORKS
# **needed to calculate t-stat with vi instead of true p$se^2
# and needed to compare the t-stats, not the yis, to truncated t


# now try truncating at affirm status
# seems close? WORKS?
crit = qt( p = 0.975,
           df = p$m - 1 )

extrunc( spec = "t",
         ncp = p$Mu / unique(p$se),
         df = p$m-1,
         a = crit )
mean(fake5$tstat2[ fake5$tstat2 > crit ])

vartrunc( spec = "t",
          ncp = p$Mu / unique(p$se),
          df = p$m-1,
          a = crit )

var(fake5$tstat2[ fake5$tstat2 > crit ])


# PUT BACK HETEROGENEITY  -------------------------------

# increase heterogeneity here to exacerbate problem, if there is one
# and reduce se a little
# and increase m to make normal approx really good
p$T2 = 0.25
p$t2w = 0.25
p$se = 0.2
p$m = 200
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


# **this is the "real" t-stat with heterogeneity
fake5$tstat2 = fake5$yi / sqrt(p$T2 + p$t2w + fake5$vi)

# confirm how the regular t-stat is calculated
expect_equal( fake5$yi / sqrt(fake5$vi), fake5$tstat )


# 1.98
extrunc( spec = "t",
         ncp = p$Mu / sqrt(p$T2 + p$t2w + unique(p$se)^2),
         df = p$m-1,
         a = 1 )

# extrunc( spec = "t",
#          ncp = p$Mu / unique(p$se),
#          df = p$m-1,
#          a = 1 )

# too low!
mean(fake5$tstat2[ fake5$tstat2 > 1 ])
# now I'm consistently getting 1

vartrunc( spec = "t",
          ncp = p$Mu / sqrt(p$T2 + p$t2w + unique(p$se)^2),
          df = p$m-1,
          a = 1 )

var(fake5$tstat2[ fake5$tstat2 > 1 ])


# now look at the "fake" t-stats that don't incorporate heterogeneity

vartrunc( spec = "t",
          ncp = p$Mu / sqrt(p$T2 + p$t2w + unique(p$se)^2),
          df = p$m-1,
          a = 1 )

var(fake5$tstat[ fake5$tstat > 1 ])





vartrunc( spec = "norm",
          mean = p$Mu / sqrt(p$T2 + p$t2w + unique(p$se)^2,
          sd = sqrt(p$T2 + p$t2w + unique(p$se)^2)
          df = p$m-1,
          a = 1 )

vartrunc( spec = "norm",
          mean = p$Mu,
          sd = sqrt(p$T2 + p$t2w + unique(p$se)^2)
          df = p$m-1,
          a = 1 )


# ~ Normal approx -----------------------------------
# should be close
# ~~ true mean of t-stats calculated WITH heterogeneity --------
p$Mu / sqrt(p$T2 + p$t2w + unique(p$se)^2); mean(fake5$tstat2)
1; var(fake5$tstat2)

# normal approx
vartrunc( spec = "norm",
          mean = p$Mu / sqrt(p$T2 + p$t2w + unique(p$se)^2),
          sd = 1,
          a = 1 )

var(fake5$tstat2[ fake5$tstat2 > 1 ])
# WORKS - THIS IS PRETTY CLOSE


# ~~ now try the CALCULATED t-stats --------------
se = unique(p$se)
p$Mu / se; mean(fake5$tstat)
# **because the CALCULATED t-stats scale by the incomplete (marginal) variance
# much closer, but true variance still a bit larger
# probably because of sampling error in se estimate
(1/se^2) * (p$T2 + p$t2w + se^2); var(fake5$tstat)

# check if that's the reason by calculating "t-stats" using true SE instead of estimated
# **YES
# so we could either ignore this source of variability or could try to delta-method it?
var(fake5$yi / se)

# x1: yi
# x2: vi
library(msm)
correctedSE = deltamethod( g = ~ x1/sqrt(x2),
             mean = c(p$Mu, se^2),
             cov = matrix( c( p$T2 + p$t2w + se^2, 0, 0, var(fake5$vi) ),
                           nrow = 2 ) )
# **a bit closer, but still not quite
correctedSE^2; var(fake5$tstat)


# normal approx
vartrunc( spec = "norm",
          mean = p$Mu / se,
          sd = sqrt( (1/se^2) * (p$T2 + p$t2w + se^2) ),
          a = 1 )

# with delta method variance
vartrunc( spec = "norm",
          mean = p$Mu / se,
          sd = correctedSE,
          a = 1 )

# still a bit too high
var(fake5$tstat[ fake5$tstat > 1 ])

# # not sure how to use t-dist here...
# vartrunc( spec = "t",
#           ncp = (p$Mu / correctedSE),
#           df = p$m-1,
#           a = 1 )

#bm




# TRY FULL SIM AGAIN USING ABOVE REALIZATIONS ------------------------------

# read back in
# note: uses different parameters from sims above
setwd( here("2021-4-29 distribution of affirmatives and nonaffirmatives under finite hacking") )
# # keep only 1 draw per study
# d = fread("sim_meta_no_hack_hugek.csv")
# d = d[ d$Di == 1, ]

# try instead with huge m to eliminate t vs. normal issues
d = fread("sim_meta_Nmax1_hugek_hugem_all_studies.csv")


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



# ~ CHECK AGAINST TRUNCATED DISTRIBUTION ------------------------------

# **use normal based on what I learned above

# first look at t-stats calculated using real vi to avoid vi estimation error
d$tstat2 = d$yi / sqrt(d$viTrue)



Mu = unique(d$Mu)
T2 = unique(d$T2)
t2w = unique(d$t2w)
m = unique(d$m)
se = unique(d$se)


Mu / se; mean(d$tstat)
# **because the CALCULATED t-stats scale by the incomplete (marginal) variance
# much closer, but true variance still a bit larger
# probably because of sampling error in se estimate
(1/se^2) * (T2 + t2w + se^2); var(d$tstat); var(d$tstat2)


# x1: yi
# x2: vi
library(msm)
correctedSE = deltamethod( g = ~ x1/sqrt(x2),
                           mean = c(Mu, se^2),
                           cov = matrix( c( T2 + t2w + se^2, 0, 0, var(d$vi) ),
                                         nrow = 2 ) )
# **a bit closer, but still not quite
correctedSE^2; var(d$tstat); var(d$tstat2)


# empirical moments with cutoff
var(d$tstat[ d$tstat > cutoff])
var(d$tstat2[ d$tstat2 > cutoff])

mean(d$tstat[ d$tstat > cutoff])
mean(d$tstat2[ d$tstat2 > cutoff])

# normal approx without delta method
cutoff = 2
vartrunc( spec = "norm",
          mean = p$Mu / se,
          sd = sqrt( (1/se^2) * (p$T2 + p$t2w + se^2) ),
          a = cutoff )

extrunc( spec = "norm",
         mean = p$Mu / se,
         sd = sqrt( (1/se^2) * (p$T2 + p$t2w + se^2) ),
         a = cutoff )

# with delta method variance
# **this one is much better than above one in this case
# probably because m isn't very large
vartrunc( spec = "norm",
          mean = Mu / se,
          sd = correctedSE,
          a = cutoff )

extrunc( spec = "norm",
         mean = Mu / se,
         sd = correctedSE,
         a = cutoff )



# t-distribution instead (as in the first sims in this script)
extrunc( spec = "t",
         ncp = Mu / sqrt( T2 + t2w + se^2 ),
         df = m-1,
         a = cutoff )

vartrunc( spec = "t",
          ncp = Mu / sqrt( T2 + t2w + se^2 ),
          df = m-1,
          a = cutoff )

# t-dist with delta method variance
# makes very little difference to t-distribution
#   even though it made a pretty big difference to normal
extrunc( spec = "t",
         ncp = Mu / correctedSE,
         df = m-1,
         a = cutoff )

vartrunc( spec = "t",
          ncp = Mu / correctedSE,
          df = m-1,
          a = cutoff )


# NORMAL VS. T ASYMPTOTIC EQUIVALENCE ---------------------------


# shouldn't these be the same?
extrunc( spec = "norm",
         mean = Mu/se,
         sd = correctedSE,
         a = -Inf,
         b = Inf )

# this is definitely like a bug with package:
# answer should obviously be 1 still
extrunc( spec = "norm",
         mean = Mu,
         sd = correctedSE,
         a = -999,
         b = 999 )

# x = rtrunc( n=2000,
#             spec = "norm",
#             mean = Mu / se,
#             sd = correctedSE,
#             a = -999,
#             b = 999 )
# mean(x)

extrunc( spec = "norm",
         mean = Mu / se,
         sd = correctedSE,
         a = ,
         b = Inf )




# df -> infty; no truncation
# **Wikipedia: Limit of noncentral t as df -> infty is N(ncp, 1)

# this definitely DOES need to divide by se in the mean
vartrunc( spec = "norm",
         mean = Mu / se,
         sd = correctedSE,
         a = -Inf,
         b = Inf )
# yes; with large m, this is exactly right
var(d$tstat)
var(d$tstat2)


# this will always have variance of 1 
# so how are we supposed to use it?
vartrunc( spec = "t",
          ncp = Mu / correctedSE,
          df = 1000,
          a = -Inf,
          b = Inf )

# I think we'd need to actually calculate the t-stats differently:
var(d$tstat / correctedSE)
# this works

# # calculated the ESTIMATED corrected SE for the rescaling
# d %>% rowwise() %>%
#   mutate( correctedSE.est = deltamethod( g = ~ x1/sqrt(x2),
#                                          mean = c(Mu, d$vi),
#                                          cov = matrix( c( p$T2 + p$t2w + sqrt(d$vi)^2, 0, 0, var(d$vi) ),
#                                                        nrow = 2 ) ) )
# 

d$tstatRescaled = d$tstat / correctedSE


# now look at truncation
## Method 1: Normal
# yes, looks good
cutoff = 2
extrunc( spec = "norm",
         mean = Mu / se,
         sd = correctedSE,
         a = cutoff )
mean( d$tstat[d$tstat > cutoff] )

vartrunc( spec = "norm",
         mean = Mu / se,
         sd = correctedSE,
         a = cutoff )
var( d$tstat[d$tstat > cutoff] )


# lower part of truncation
extrunc( spec = "norm",
         mean = Mu / se,
         sd = correctedSE,
         b = cutoff )
mean( d$tstat[d$tstat < cutoff] )

vartrunc( spec = "norm",
          mean = Mu / se,
          sd = correctedSE,
          b = cutoff )
var( d$tstat[d$tstat < cutoff] )

## Method 2: truncated t
# DOES NOT WORK
extrunc( spec = "t",
          ncp = Mu / correctedSE,
          df = 1000,
          a = cutoff )
#*note use of RESCALED t-stat here
mean( d$tstatRescaled[d$tstatRescaled > cutoff] )

vartrunc( spec = "t",
          ncp = Mu / correctedSE,
          df = 1000,
          a = cutoff )
var( d$tstatRescaled[d$tstatRescaled > cutoff] )



# ~ 1. NO HACKING, NMAX = 1 ---------------------------------

# read back in
# note: uses different parameters from sims above
setwd( here("2021-4-29 distribution of affirmatives and nonaffirmatives under finite hacking") )

# try instead with huge m to eliminate t vs. normal issues
d = fread("sim_meta_Nmax1_hugek_hugem_all_studies.csv")


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


d$tstat2 = d$yi / sqrt(d$viTrue)

dph = d[ d$Di == TRUE, ]
mean(dph$affirm)

Mu = unique(d$Mu)
T2 = unique(d$T2)
t2w = unique(d$t2w)
m = unique(d$m)
se = unique(d$se)


library(msm)
correctedSE = deltamethod( g = ~ x1/sqrt(x2),
                           mean = c(Mu, se^2),
                           cov = matrix( c( T2 + t2w + se^2, 0, 0, var(d$vi) ),
                                         nrow = 2 ) )

crit = unique(d$tcrit)


## Method 1: Normal
# WORKS
extrunc( spec = "norm",
         mean = Mu / se,
         sd = correctedSE,
         a = crit )
mean(d$tstat[ d$tstat > crit] )
mean(dph$tstat )

vartrunc( spec = "norm",
          mean = Mu / se,
          sd = correctedSE,
          a = crit )
var(d$tstat[ d$tstat > crit] )


# lower part of truncation
extrunc( spec = "norm",
         mean = Mu / se,
         sd = correctedSE,
         b = crit )
mean( d$tstat[d$tstat < crit] )

vartrunc( spec = "norm",
          mean = Mu / se,
          sd = correctedSE,
          b = crit )
var( d$tstat[d$tstat < crit] )



# ~ 2. ALL HACKED, NMAX = 10 ---------------------------------

# read back in
# note: uses different parameters from sims above
setwd( here("2021-4-29 distribution of affirmatives and nonaffirmatives under finite hacking") )

# try instead with huge m to eliminate t vs. normal issues
d = fread("sim_meta_Nmax10_hugek_hugem_allhacked_all_studies.csv")


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


d$tstat2 = d$yi / sqrt(d$viTrue)

dph = d[ d$Di == TRUE, ]
table(dph$affirm)

Mu = unique(d$Mu)
T2 = unique(d$T2)
t2w = unique(d$t2w)
m = unique(d$m)
se = unique(d$se)


library(msm)
correctedSE = deltamethod( g = ~ x1/sqrt(x2),
                           mean = c(Mu, se^2),
                           cov = matrix( c( T2 + t2w + se^2, 0, 0, var(d$vi) ),
                                         nrow = 2 ) )

crit = unique(d$tcrit)


## Method 1: Normal
# NOW DOES NOT WORK
extrunc( spec = "norm",
         mean = Mu / se,
         sd = correctedSE,
         a = crit )
mean( dph$tstat )
mean(d$tstat[ d$tstat > crit] )

vartrunc( spec = "norm",
          mean = Mu / se,
          sd = correctedSE,
          a = crit )
var(d$tstat[ d$tstat > crit] )


# lower part of truncation
extrunc( spec = "norm",
         mean = Mu / se,
         sd = correctedSE,
         b = crit )
mean( d$tstat[d$tstat < crit] )

vartrunc( spec = "norm",
          mean = Mu / se,
          sd = correctedSE,
          b = crit )
var( d$tstat[d$tstat < crit] )


# NOW BOTH MOMENTS ARE WRONG...
# BM: TRY TO FIGURE THIS OUT!




# ~ 3. ALL HACKED, BUT NMAX = 1 ---------------------------------

# read back in
# note: uses different parameters from sims above
setwd( here("2021-4-29 distribution of affirmatives and nonaffirmatives under finite hacking") )

# try instead with huge m to eliminate t vs. normal issues
d = fread("sim_meta_Nmax1_hugek_hugem_allhacked_all_studies.csv")


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


d$tstat2 = d$yi / sqrt(d$viTrue)

dph = d[ d$Di == TRUE, ]
table(dph$affirm)

Mu = unique(d$Mu)
T2 = unique(d$T2)
t2w = unique(d$t2w)
m = unique(d$m)
se = unique(d$se)


library(msm)
correctedSE = deltamethod( g = ~ x1/sqrt(x2),
                           mean = c(Mu, se^2),
                           cov = matrix( c( T2 + t2w + se^2, 0, 0, var(d$vi) ),
                                         nrow = 2 ) )

crit = unique(d$tcrit)


## Method 1: Normal
# NOW DOES NOT WORK
extrunc( spec = "norm",
         mean = Mu / se,
         sd = correctedSE,
         a = crit )
mean( dph$tstat )
mean(d$tstat[ d$tstat > crit] )

vartrunc( spec = "norm",
          mean = Mu / se,
          sd = correctedSE,
          a = crit )
var( dph$tstat )


# lower part of truncation
extrunc( spec = "norm",
         mean = Mu / se,
         sd = correctedSE,
         b = crit )
mean( d$tstat[d$tstat < crit] )

vartrunc( spec = "norm",
          mean = Mu / se,
          sd = correctedSE,
          b = crit )
var( d$tstat[d$tstat < crit] )


# ~ 4. ALL HACKED, NMAX = 200 ---------------------------------

# so every study set should get an affirmative

# read back in
# note: uses different parameters from sims above
setwd( here("2021-4-29 distribution of affirmatives and nonaffirmatives under finite hacking") )

# try instead with huge m to eliminate t vs. normal issues
d = fread("sim_meta_Nmax1_hugek_hugem_allhacked_all_studies.csv")


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


d$tstat2 = d$yi / sqrt(d$viTrue)

dph = d[ d$Di == TRUE, ]
table(dph$affirm)

Mu = unique(d$Mu)
T2 = unique(d$T2)
t2w = unique(d$t2w)
m = unique(d$m)
se = unique(d$se)


library(msm)
correctedSE = deltamethod( g = ~ x1/sqrt(x2),
                           mean = c(Mu, se^2),
                           cov = matrix( c( T2 + t2w + se^2, 0, 0, var(d$vi) ),
                                         nrow = 2 ) )

crit = unique(d$tcrit)


## Method 1: Normal
# NOW DOES NOT WORK
extrunc( spec = "norm",
         mean = Mu / se,
         sd = correctedSE,
         a = crit )
mean( dph$tstat )
mean(d$tstat[ d$tstat > crit] )

vartrunc( spec = "norm",
          mean = Mu / se,
          sd = correctedSE,
          a = crit )
var( dph$tstat )


# lower part of truncation
extrunc( spec = "norm",
         mean = Mu / se,
         sd = correctedSE,
         b = crit )
mean( d$tstat[d$tstat < crit] )

vartrunc( spec = "norm",
          mean = Mu / se,
          sd = correctedSE,
          b = crit )
var( d$tstat[d$tstat < crit] )






# spin idea: 
# 0. simulate Nmax = 1 but hack = "affirm", which should be exactly the same as what I already have
# 1. simulate Nmax = huge so that every single set has an affirmative (eliminate seleciton issues due to sets that never get an affirmative)



# NEXT UP: THINK ABOUT WHERE WE ARE WRT THE PREVIOUS SIM RESULTS USING T DISTRIBUTION. 
# SHOULD I ACTUALLY CHANGE ANYTHING ABOUT THAT?
# OVERALL GOAL WAS TO LOOK AT DIST OF STATS UNDER FINITE VS. INFINITE HACKING. 
# SO SHOULD PROBABLY SEE IF BIAS IS SAME AT DIFFERENT LEVELS OF NMAX. 
# ALSO LOOK AT EXPECTATION OF NONAFFIRMATIVES.

# not sure how we're supposed to use truncated t because it always has variance of 1 
# asymptotically...I guess we'd need to rescale 



# bm


# going to try same scenario, but taking away within-study heterogeneity



# NUMBERED SIMULATION EXPERIMENTS -------------------------
# as in table in project log


results.dir = here( "2021-4-29 distribution of affirmatives and nonaffirmatives under finite hacking/Numbered simulation experiments" )

# ~ EXPT 2: ALL HACKED, NMAX = 10 -------------------------

x = quick_sim( .p = data.frame( Mu = 1,
                       T2 = 0.25,
                       m = 500,
                       t2w = .25,
                       se = .5,
                       
                       Nmax = 10,
                       hack = "affirm",
                       
                       k = 10000,
                       k.hacked = 10000,
                       
                       .results.dir = results.dir,
                       sim.name = "expt_2_dataset" ) )

x$res


# ~ EXPT 4: ALL HACKED, NMAX = 200 -------------------------

x = quick_sim( .p = data.frame( Mu = 1,
                                T2 = 0.25,
                                m = 500,
                                t2w = .25,
                                se = .5,
                                
                                Nmax = 200,
                                hack = "affirm",
                                rho = 0,
                                
                                k = 10000,
                                k.hacked = 10000,
                                
                                .results.dir = results.dir,
                                sim.name = "expt_4_dataset" ) )

x$res



# ~ EXPT 5: ALL HACKED, NMAX = 10, T2 = t2w = 0 -------------------------

# try without heterogeneity, still Nmax = 10
# THIS ONE IS PERFECT
x = quick_sim( .p = data.frame( Mu = 1,
                                T2 = 0,
                                m = 500,
                                t2w = 0,
                                se = .5,
                                
                                Nmax = 10,
                                hack = "affirm",
                                
                                k = 10000,
                                k.hacked = 10000,
                                
                                .results.dir = results.dir,
                                sim.name = "expt_5_dataset" ) )

x$res


# ~ EXPT 6: ALL HACKED, NMAX = 10, T2 > 0, t2w = 0 -------------------------

# this messes it up!!
x = quick_sim( .p = data.frame( Mu = 1,
                                T2 = 0.25,
                                m = 500,
                                t2w = 0,
                                se = .5,
                                
                                Nmax = 10,
                                hack = "affirm",
                                
                                k = 10000,
                                k.hacked = 10000,
                                
                                
                                sim.name = "expt_6_res" ),
               
               .results.dir = results.dir )

setwd(results.dir)
load("expt_6_res")
x = returnList
d = x$d

Mu = unique(d$Mu)
T2 = unique(d$T2)
t2w = unique(d$t2w)
m = unique(d$m)
se = unique(d$se)



x$correctedSE^2

# first draw of each study set
( t1 = d %>% filter( !duplicated(study) ) %>%
  summarise( n(),
             var(mui),
             var(muin),
             var(yi),
             mean(tstat),
             var(tstat) ) )

# all draws from each study set
( t2 = d %>% 
  summarise( n(),
             var(mui),
             var(muin),
             var(yi),
             mean(tstat),
             var(tstat) ) )


# **IMPORTANT:
# marginal mean(tstat) and var(tstat) match Mu/se and correctedSE^2
# ONLY when I look at JUST the first draw per study set 
# if I use all draws per study set, mean(tstat) << Mu/se and var(tstat) slightly < correctedSE^2

# try calculating truncated expectation by cheating
# i.e., using the empirical marginal mean of ALL tstats
extrunc( spec = "norm",
         mean = t2$`mean(tstat)`,
         sd = sqrt(2),
         a = crit )

vartrunc( spec = "norm",
         mean = t2$`mean(tstat)`,
         sd = sqrt(2),
         a = crit )

x$res
# this doesn't really work


t3 = d %>% filter(!duplicated(study)) %>%
  group_by(N < 10) %>%
  summarise( mean(tstat),
             var(tstat) )
t3

extrunc( spec = "norm",
         mean = t3$`mean(tstat)`[t3$`N < 10` == 1],
         #sd = sqrt(2),
         sd = sqrt( t3$`var(tstat)`[t3$`N < 10` == 1] ),
         a = crit )

vartrunc( spec = "norm",
          mean = t3$`mean(tstat)`[t3$`N < 10` == 1],
          #sd = sqrt(2),
          sd = sqrt( t3$`var(tstat)`[t3$`N < 10` == 1] ),
          a = crit )

x$res

# this also doesn't work...not sure why...


# look at published ones only
# **this nicely illustrates the issue of selecting on mui effectively
( t4 = d %>% 
    filter(Di ==1) %>%
    summarise( n(),
               mean(mui),
               var(mui),
               mean(muin),
               var(muin),
               var(yi),
               mean(tstat),
               var(tstat) ) )

# **note that the distribution of muin after p-hacking selection is NOT normal
#  but rather skewed, which means any normal approximation below is not going to 
#  work perfectly
hist( d$muin[ d$Di == 1 ] )

correctedSE = deltamethod( g = ~ x1/sqrt(x2),
                           mean = c(Mu, se^2),
                           cov = matrix( c( t4$`var(muin)` + se^2, 0, 0, var(d$vi) ),
                                         nrow = 2 ) )

extrunc( spec = "norm",
         mean = t4$`mean(muin)` / se,
         #sd = sqrt(2),
         sd = correctedSE,  # need to scale for t-stat
         a = crit )

vartrunc( spec = "norm",
          mean = t3$`mean(tstat)`[t3$`N < 10` == 1],
          #sd = sqrt(2),
          sd = sqrt( t3$`var(tstat)`[t3$`N < 10` == 1] ),
          a = crit )


# calculate expectation of each one individually
dp = d %>% filter(Di == 1) %>%
  rowwise() %>%
  mutate( theoryExp_in = extrunc( spec = "norm",
                                   mean = muin / se,
                                   #sd = sqrt(2),
                                   sd = se,  # no T2 in here
                                   a = crit ),
          
          # **GOOD ONE - save it! 
          theoryExp_i = extrunc( spec = "norm",
                                  mean = mui / se,
                                  #sd = sqrt(2),
                                  sd = sqrt( (1/se^2) * (t2w + se^2) ),  
                                  a = crit ),
          
          varExp_i = vartrunc( spec = "norm",
                                 mean = mui / se,
                                 #sd = sqrt(2),
                                 sd = sqrt( (1/se^2) * (t2w + se^2) ),  
                                 a = crit ) )

mean(dp$theoryExp_in)
mean(dp$theoryExp_i)  # this one is HELLA close!!!! **save this one
mean(dp$tstat)

mean(dp$varExp_i)  # this is still off, probably because it's not really normal
var(dp$tstat)




# look at published nonaffirmatives
dn = d %>% filter( Di == 1 & affirm == FALSE )



# ~ EXPT 7: ALL HACKED, NMAX = 10, T2 = 0, t2w > 0 -------------------------

# # this one is CORRECT!!
# x = quick_sim( .p = data.frame( Mu = 1,
#                                 T2 = 0,
#                                 m = 500,
#                                 t2w = 0.25,
#                                 se = .5,
# 
#                                 Nmax = 10,
#                                 hack = "affirm",
# 
#                                 k = 10000,
#                                 k.hacked = 10000,
# 
# 
#                                 sim.name = "expt_7_res" ),
# 
#                .results.dir = results.dir )

setwd(results.dir)
load("expt_7_res")
x = returnList
d = x$d

# first draw of each study set
( t1 = d %>% filter( !duplicated(study) ) %>%
    summarise( n(),
               var(mui),
               var(muin),
               var(yi),
               mean(tstat),
               var(tstat) ) )

# all draws from each study set
( t2 = d %>% 
    summarise( n(),
               var(mui),
               var(muin),
               var(yi),
               mean(tstat),
               var(tstat) ) )




