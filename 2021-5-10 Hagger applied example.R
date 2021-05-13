
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
plot_trunc_densities(.obj, showAffirms = TRUE)


hist(dm$yi/sqrt(dm$vi), breaks = 20)
hist(.obj$data$tstat, breaks = 20)

.obj$sanityChecks$tstatMeanMLE  # 56!!!
# quite close to the mean in all the replications
mean( dm$yi/sqrt(dm$vi) )  # 2.45

.obj$metaCorr  # corrected: 16.8
.obj$metaNaive  # naive: 0.68

# ~ Sanity check: No truncation ---------------------

# use a really high cutoff so that it's just the regular MLE

.obj = correct_meta_phack2(yi = dm$yi,
                           vi = dm$vi,
                           crit = 100)

.obj$sanityChecks$tstatMeanMLE  # **now matches the mean nicely

.obj$metaCorr  # corrected: 16.8
.obj$metaNaive  # naive: 0.68

# **This is very informative: I think the issue with the truncated MLE is that the distribution is
#  really light-tailed
plot_trunc_densities(.obj, showAffirms = FALSE) +
   # superimpose the actual critical value
   geom_vline(xintercept = 1.96, color = "gray")



# try an intermediate cutoff
.obj = correct_meta_phack2(yi = dm$yi,
                           vi = dm$vi,
                           crit = 3)

.obj$sanityChecks$tstatMeanMLE  # now matches the mean nicely
plot_trunc_densities(.obj, showAffirms = FALSE)



# ~ Explore doing MLE for each nonaffirm individually ---------------------

d = dm
d$tstat = d$yi / sqrt(d$vi)
d$affirm = d$tstat > crit


# published affirmatives only
dpn = d[ d$affirm == FALSE, ]

# MLE for EACH t-stat individually
dpn %>% rowwise() %>%
   mutate( )


#@THINK ABOUT: Currently this allows the within-study variance (t2w + SE^2)
#  to be estimated. We could alternatively set t2w = 0, for example, using the "fixed"
#  argument of mle.tmvnorm
one_nonaffirm_mle = function(.tstat,
                             .se = NA, # of yi, not tstat
                             .t2w = NA,
                             .crit
) {
   
   if (length(.tstat) > 1) stop("fn not intended for more than 1 observation")
   
   # if we're given the components of the WITHIN-study variance,
   #  treat as fixed rather than estimating it
   if ( !is.na(.se) & !is.na(.t2w) ) Vw = (1/.se^2) * (.t2w + .se^2)
   
   if ( exists("Vw") ) {
      mle.fit = mle.tmvnorm( X = as.matrix(.tstat, ncol = 1),
                             lower = -Inf,
                             upper = .crit,
                             # fixed parameters must be named exactly as
                             #  coef(mle.fit) returns them
                             fixed = list(sigma_1.1 = Vw) )
   } else{
      mle.fit = mle.tmvnorm( X = as.matrix(.tstat, ncol = 1),
                             lower = -Inf,
                             upper = .crit)
   }
   
   
   mles = coef(mle.fit)
   
   return( data.frame( tstat.mu.hat = as.numeric(mles[1]),
                     muiHat = as.numeric(mles[1]) * .se,
                       VwHat = as.numeric(mles[2]) ) )
   
}

one_nonaffirm_mle( .tstat = dpn$tstat[1],
                   .se = sqrt(dpn$vi[1]),
                   .t2w = 0,
                   .crit = qnorm(.975) )

one_nonaffirm_mle( .tstat = dpn$tstat[1],
                   # .se = sqrt(dpn$vi[1]),
                   # .t2w = 0,
                   .crit = qnorm(.975) )

one_nonaffirm_mle( .tstat = dpn$tstat[1],
                   .se = sqrt(dpn$vi[1]),
                   .t2w = 0,
                   .crit = qnorm(.975) )


dpn[1,]


# mle.fit = mle.tmvnorm( X = as.matrix(dpn$tstat[1], ncol = 1),
#                        lower = -Inf,
#                        fixed = list(sigma_1.1 = 1),
#                        upper = crit)
# 
# mle.fit = mle.tmvnorm( X = as.matrix(dpn$tstat[1], ncol = 1),
#                        lower = -Inf,
#                        fixed = list(mu_1 = 0),
#                        upper = crit)
# 
# 
# mles = coef(mle.fit)
# dpn$tstat[1]
# 
# 
# dpn$vi[1]

# P-HACKING ADJUSTMENT IN HAGGER REPLICATIONS -----------------------------


.obj = correct_meta_phack2(yi = dr$yi,
                           vi = dr$vi)

# this one is much more reasonable
.obj$sanityChecks$tstatMeanMLE
# quite close to the mean in all the replications
mean( dr$yi/sqrt(dr$vi) )

plot_trunc_densities(.obj)
plot_trunc_densities(.obj, showAffirms = TRUE)

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
plot_trunc_densities(.obj, showAffirms = TRUE)

.obj$metaCorr  # 0.29
.obj$metaNaive  # 0.21

# **Conclusion: The more left-skewed the nonaffirmative t-stats are, the
#  larger the corrected estimate will be 



# ~ Try doing it with only preregistered studies ----------------------
d = d3 %>% filter(qual.prereg2 == 1)


.obj = correct_meta_phack2(yi = d$yi,
                           vi = d$vi)

.obj$sanityChecks$tstatMeanMLE
mean(d$yi/sqrt(d$vi))  # from all studies

.obj$metaCorr  
.obj$metaNaive 

#@ strange because the tstat MLE is less than the overall one, 
# yet metaCorr > metaNaive
# I think that means that the rescaling approximation doesn't work very well?
cor(d$yi, d$vi)



# look at all of them by setting crit = 100
.obj = correct_meta_phack2(yi = d$yi,
                           vi = d$vi,
                           crit = 100)


# these match exactly :)
.obj$metaCorr  
.obj$metaNaive  

.obj$sanityChecks$tstatMeanMLE; mean(d$yi/sqrt(d$vi))
plot_trunc_densities(.obj, showAffirms = TRUE)



# P-HACKING ADJUSTMENT IN SIMULATED META-ANALYSIS -----------------------------

d = sim_meta(Nmax = 1,  
             Mu = 0,  
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

plot_trunc_densities(.obj, showAffirms = TRUE)


.obj$metaCorr  
.obj$metaNaive 

#hmm...seems pretty noisy?
#bm

# SEPARATE MLES FOR EACH NONAFFIRMATIVE -----------------------------
























