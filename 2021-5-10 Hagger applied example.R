
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


correct_meta_phack2(yi = dm$yi,
                    vi = dm$vi)


#bm: it's not able to get CIs and the MLEs are kind of nuts
# maybe try the replications for comparison, since they shouldn't have any p-hacking?








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

