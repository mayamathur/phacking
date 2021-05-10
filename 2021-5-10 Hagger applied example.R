
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

.obj = correct_meta_phack2(yi = dm$yi,
                    vi = dm$vi)

# t-stat MLE is HUGE
.obj$sanityChecks$tstatMeanMLE

# **plot shows that the MLE is so big because the nonaffirmative t-stats are actually left-skewed
# very interesting
plot_trunc_densities(.obj)

hist(dm$yi/sqrt(dm$vi), breaks = 20)
hist(.obj$data$tstat, breaks = 20)

.obj$sanityChecks$tstatMeanMLE  # 56!!!
# quite close to the mean in all the replications
mean( dm$yi/sqrt(dm$vi) )  # 2.45

.obj$metaCorr  # corrected: 16.8
.obj$metaNaive  # naive: 0.68


# P-HACKING ADJUSTMENT IN HAGGER REPLICATIONS -----------------------------


.obj = correct_meta_phack2(yi = dr$yi,
                           vi = dr$vi)

# this one is much more reasonable
.obj$sanityChecks$tstatMeanMLE
# quite close to the mean in all the replications
mean( dr$yi/sqrt(dr$vi) )

plot_trunc_densities(.obj)

.obj$sanityChecks$tstatMeanMLE  # 0.25
# quite close to the mean in all the replications
mean( dr$yi/sqrt(dr$vi) )  # 0.26

.obj$metaCorr  # corrected: 0.05
.obj$metaNaive  # naive: also 0.05


# P-HACKING ADJUSTMENT IN AWR META-ANALYSIS -----------------------------

setwd("~/Dropbox/Personal computer/Independent studies/2020/Meta-regression metrics (MRM)/Applied example/Prepped data")

d3 = fread("mathur_data_prepped.csv")

.obj = correct_meta_phack2(yi = d3$yi,
                           vi = d3$vi)

# this one is much more reasonable
.obj$sanityChecks$tstatMeanMLE  # 1.09
# quite close to the mean in all the replications
mean( d3$yi/sqrt(d3$vi) )  # 1.44

plot_trunc_densities(.obj)


.obj$metaCorr  # 0.29
.obj$metaNaive  # 0.21

# Conclusion: The more left-skewed the nonaffirmative t-stats are, the
#  larger the corrected estimate will be 


# P-HACKING ADJUSTMENT IN SIMULATED META-ANALYSIS -----------------------------

d = sim_meta(Nmax = 1,  
             Mu = 5,  
             T2 = 2,  
             
             # study parameters, assumed same for all studies:
             m = 500,  
             t2w = 0.25,  
             se = 0.5,  
             
             rho = 0,  
             
             hack = "affirm",  
             
             k = 1000,  
             k.hacked = 0 )

dp = d[ d$Di == TRUE, ]
table(dp$affirm)


.obj = correct_meta_phack2(yi = dp$yi,
                           vi = dp$vi)


.obj$sanityChecks$tstatMeanMLE  
mean( d$yi/sqrt(d$vi) )  

plot_trunc_densities(.obj)


.obj$metaCorr  
.obj$metaNaive 










# # FROM EARLIER COMPARISON -----------------------------
# 
# # ~ REANALYZE HAGGER'S META-ANALYSIS -----------------------------
# 
# 
# # plain meta-analysis
# meta.m = rma.uni( yi = yi,
#                   vi = vi,
#                   data = dm,
#                   method = "REML",
#                   knha = TRUE )
# summary(meta.m)
# 
# ##### Funnel Plots ######
# setwd(res.dir)
# pdf( file = paste( "funnel_hagger_meta", sep = "" ) )
# funnel.rma(meta.m,
#            level = c(0.95),
#            legend = TRUE,
#            main = "Hagger meta-analysis",
#            refline = 0)
# dev.off()
# 
# 
# # ** reported corrected muhat: 0.45 (vs. uncorrected 0.68)
# #  in original analysis, 3PSM corrected to 0.50, so our result isn't that different
# ( meta.wtr = weightfunct( effect = dm$yi,
#                           v = dm$vi ) )
# 
# 
# 
# ##### Normality #####
# dm$ens = calib_ests(yi = dm$yi,
#                     sei = sqrt(dm$vi))
# # looks fairly non-normal
# setwd(res.dir)
# pdf( file = paste( "calib_ests_meta", sep = "" ) )
# plot( stats::density(dm$ens),
#       main = "Hagger meta-analysis" )
# dev.off()
# 
# 
# ##### Worst-Case Meta-Analysis ######
# temp = dm[ dm$affirm == 0, ]
# # ** worst-case: 0.30
# library(robumeta)
# robu( yi ~ 1,
#       data = temp,
#       studynum = 1:nrow(temp),
#       var.eff.size = vi,
#       modelweights = "HIER",
#       small = TRUE )
# 
# 
# 
# 
# # ~ REANALYZE HAGGER'S REPLICATIONS -----------------------------
# 
# # plain meta-analysis
# meta.r = rma.uni( yi = yi,
#                   vi = vi,
#                   data = dr,
#                   method = "REML",
#                   knha = TRUE )
# summary(meta.r)
# 
# ##### Funnel Plots ######
# setwd(res.dir)
# pdf( file = paste( "funnel_hagger_rep", sep = "" ) )
# funnel.rma(meta.r,
#            level = c(0.95),
#            legend = TRUE,
#            main = "Hagger replications",
#            refline = 0)
# dev.off()
# 
# 
# ##### Normality #####
# dr$ens = calib_ests(yi = dr$yi,
#                     sei = sqrt(dr$vi))
# # looks fairly non-normal
# setwd(res.dir)
# pdf( file = paste( "calib_ests_rep", sep = "" ) )
# plot( stats::density(dr$ens),
#       main = "Hagger replications" )
# dev.off()
# 
# 
# # ~ REANALYZE DANG META-ANALYSIS -----------------------------
# 
# # keep only crossing-out-letter tasks
# 
# # plain meta-analysis
# # correct number of studies per Dang's Table 1, but our estimate is a bit higher than reported 0.58
# #  even when using DL
# meta.d = rma.uni( yi = yi,
#                   vi = vi,
#                   data = dd.cl,
#                   method = "REML",
#                   knha = TRUE )
# summary(meta.d)
# 
# ##### Funnel Plots ######
# setwd(res.dir)
# pdf( file = paste( "funnel_dang_meta", sep = "" ) )
# funnel.rma(meta.d,
#            level = c(0.95),
#            legend = TRUE,
#            main = "Dang meta-analysis",
#            refline = 0)
# dev.off()
# 
# 
# # ** barely changed: 0.65 (tau=0.42)
# ( meta.wtr = weightfunct( effect = dd.cl$yi,
#                           v = dd.cl$vi ) )
# 
# 
# 
# ##### Normality #####
# dd.cl$ens = calib_ests(yi = dd.cl$yi,
#                        sei = sqrt(dd.cl$vi))
# 
# setwd(res.dir)
# pdf( file = paste( "calib_ests_meta", sep = "" ) )
# plot( stats::density(dd.cl$ens),
#       main = "Dang meta-analysis" )
# dev.off()
# 
# 
# ##### Worst-Case Meta-Analysis ######
# temp = dd.cl[ dd.cl$affirm == 0, ]
# 
# # # robumeta has SVD problem here
# # library(robumeta)
# # robu( yi ~ 1,
# #       data = temp,
# #       studynum = 1:nrow(temp),
# #       var.eff.size = vi,
# #       modelweights = "HIER",
# #       small = TRUE )
# 
# # point estimate still 0.19! 
# rma.uni( yi = yi,
#          vi = vi,
#          data = temp,
#          method = "REML",
#          knha = TRUE )
# 
