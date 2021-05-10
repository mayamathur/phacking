
# Much of this is taken from "reanalyze_hagger_only.R"

# Naive Hagger MA estimate: 0.68
# Replication estimate: 0.05

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


prepped.data.dir = "~/Dropbox/Personal computer/Miscellaneous/Journal article library/Reproducibility:replicability/Kvarven comparing meta-analysis and multisite replications/Reanalyze Kvarven/Prepped Hagger data (meta and replication)"
res.dir = "~/Dropbox/Personal computer/Miscellaneous/Journal article library/Reproducibility:replicability/Kvarven comparing meta-analysis and multisite replications/Reanalyze Kvarven/Hagger comparison results"

setwd(prepped.data.dir)
dm = read.csv("prepped_hagger_meta_data.csv")
dr = read.csv("prepped_hagger_rep_data.csv")
dd = read.csv("prepped_dang_meta_data.csv")

# dd with just the crossing-out tasks
dd.cl = dd[dd$IV == "CL",]


# P-HACKING ADJUSTMENT IN HAGGER META-ANALYSIS -----------------------------


# very similar to correct_meta_phack1, but allows for different SEs across studies



# published affirmatives only
dmn = dm[ dm$affirm == FALSE, ]
nrow(dmn)

dmn$tstat = dmn$yi / sqrt(dmn$vi)
summary(dmn$tstat)

crit = qnorm(.975)


### MLEs from trunc normal ###
# these are the MLEs of the *t-stats*
#@ IMPORTANT: for convenience, this is using the normal distribution, 
#  so won't work well for small m
mle.fit = mle.tmvnorm( X = as.matrix(dmn$tstat, ncol = 1),
                       lower = -Inf,
                       upper = crit)
#summary(mle.fit)
mles = coef(mle.fit)

#bm: it's not able to get CIs and the MLEs are kind of nuts
# maybe try the replications for comparison, since they shouldn't have any p-hacking?

# get Wald CI a different way
tstat.mu.SE = attr( summary(mle.fit), "coef" )[ "mu_1", "Std. Error" ]
tstat.mu.CI = c( mles[1] - tstat.mu.SE * qnorm(0.975),
                 mles[1] + tstat.mu.SE * qnorm(0.975) )

# # pretty good :)
# mles[1]; p$Mu/p$se
# mles[2]; (1/p$se^2) * (p$T2 + p$t2w + p$se^2)

# rescale MLEs to represent effect sizes rather than tstats
Mhat = mles[1] * .p$se
# **use Vhat to represent MARGINAL heterogeneity (i.e., T2 + t2w)
Vhat = ( mles[2] * .p$se^2 ) - .p$se^2

# and rescale CI limits
MhatCI = tstat.mu.CI * .p$se

### Sanity checks: Moments of published nonaffirms vs. theory ###
# check that moments are what we expect
# without delta method:

theoryExpTstat = extrunc(spec = "norm",
                         mean =.p$Mu /.p$se,
                         #@doesn't use the delta-method thing
                         sd = sqrt( (1/.p$se^2) * (.p$T2 +.p$t2w +.p$se^2) ),
                         b = crit )

theoryVarTstat = vartrunc(spec = "norm",
                          mean =.p$Mu /.p$se,
                          #@doesn't use the delta-method thing
                          sd = sqrt( (1/.p$se^2) * (.p$T2 +.p$t2w +.p$se^2) ),
                          b = crit )

# delta-method version (not checked and seems not to work):
# library(msm)
# # https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
# sd.y = .p$se * sqrt(.p$m)
# viSE = sqrt( 2 * sd.y^4 / (.p$m-1) )
# correctedSE = deltamethod( g = ~ x1/sqrt(x2),  # the t-stat
#                            mean = c(.p$Mu, .p$se^2),
#                            cov = matrix( c( .p$T2 + .p$t2w + .p$se^2, 0, 0, viSE^2 ),
#                                          nrow = 2 ) )
# 
# extrunc(spec = "norm",
#         mean =.p$Mu /.p$se,
#         sd = correctedSE,
#         b = crit )


# FROM EARLIER COMPARISON -----------------------------

# ~ REANALYZE HAGGER'S META-ANALYSIS -----------------------------


# plain meta-analysis
meta.m = rma.uni( yi = yi,
                  vi = vi,
                  data = dm,
                  method = "REML",
                  knha = TRUE )
summary(meta.m)

##### Funnel Plots ######
setwd(res.dir)
pdf( file = paste( "funnel_hagger_meta", sep = "" ) )
funnel.rma(meta.m,
           level = c(0.95),
           legend = TRUE,
           main = "Hagger meta-analysis",
           refline = 0)
dev.off()


# ** reported corrected muhat: 0.45 (vs. uncorrected 0.68)
#  in original analysis, 3PSM corrected to 0.50, so our result isn't that different
( meta.wtr = weightfunct( effect = dm$yi,
                          v = dm$vi ) )



##### Normality #####
dm$ens = calib_ests(yi = dm$yi,
                    sei = sqrt(dm$vi))
# looks fairly non-normal
setwd(res.dir)
pdf( file = paste( "calib_ests_meta", sep = "" ) )
plot( stats::density(dm$ens),
      main = "Hagger meta-analysis" )
dev.off()


##### Worst-Case Meta-Analysis ######
temp = dm[ dm$affirm == 0, ]
# ** worst-case: 0.30
library(robumeta)
robu( yi ~ 1,
      data = temp,
      studynum = 1:nrow(temp),
      var.eff.size = vi,
      modelweights = "HIER",
      small = TRUE )




# ~ REANALYZE HAGGER'S REPLICATIONS -----------------------------

# plain meta-analysis
meta.r = rma.uni( yi = yi,
                  vi = vi,
                  data = dr,
                  method = "REML",
                  knha = TRUE )
summary(meta.r)

##### Funnel Plots ######
setwd(res.dir)
pdf( file = paste( "funnel_hagger_rep", sep = "" ) )
funnel.rma(meta.r,
           level = c(0.95),
           legend = TRUE,
           main = "Hagger replications",
           refline = 0)
dev.off()


##### Normality #####
dr$ens = calib_ests(yi = dr$yi,
                    sei = sqrt(dr$vi))
# looks fairly non-normal
setwd(res.dir)
pdf( file = paste( "calib_ests_rep", sep = "" ) )
plot( stats::density(dr$ens),
      main = "Hagger replications" )
dev.off()


# ~ REANALYZE DANG META-ANALYSIS -----------------------------

# keep only crossing-out-letter tasks

# plain meta-analysis
# correct number of studies per Dang's Table 1, but our estimate is a bit higher than reported 0.58
#  even when using DL
meta.d = rma.uni( yi = yi,
                  vi = vi,
                  data = dd.cl,
                  method = "REML",
                  knha = TRUE )
summary(meta.d)

##### Funnel Plots ######
setwd(res.dir)
pdf( file = paste( "funnel_dang_meta", sep = "" ) )
funnel.rma(meta.d,
           level = c(0.95),
           legend = TRUE,
           main = "Dang meta-analysis",
           refline = 0)
dev.off()


# ** barely changed: 0.65 (tau=0.42)
( meta.wtr = weightfunct( effect = dd.cl$yi,
                          v = dd.cl$vi ) )



##### Normality #####
dd.cl$ens = calib_ests(yi = dd.cl$yi,
                       sei = sqrt(dd.cl$vi))

setwd(res.dir)
pdf( file = paste( "calib_ests_meta", sep = "" ) )
plot( stats::density(dd.cl$ens),
      main = "Dang meta-analysis" )
dev.off()


##### Worst-Case Meta-Analysis ######
temp = dd.cl[ dd.cl$affirm == 0, ]

# # robumeta has SVD problem here
# library(robumeta)
# robu( yi ~ 1,
#       data = temp,
#       studynum = 1:nrow(temp),
#       var.eff.size = vi,
#       modelweights = "HIER",
#       small = TRUE )

# point estimate still 0.19! 
rma.uni( yi = yi,
         vi = vi,
         data = temp,
         method = "REML",
         knha = TRUE )

