
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

# other
library(devtools)
#devtools::install_github('jtleek/tidypvals')


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


# try the full RTMA via weightr (c.f. "2021-7-4 connection to selection models")
( m1 = weightfunct( effect = dm$yi[ dm$affirm == FALSE ],
                    v = dm$vi[ dm$affirm == FALSE ],  
                    steps = c(0.025, 1),
                    weights = c(0,1), # weight such that ONLY nonaffirmatives are published
                    table = TRUE ) )

# adjusted point estimate
m1[[2]]$par[2]




# just for comparison:
# allow weightr to find the weights and use all studies
( m1 = weightfunct( effect = dm$yi,
                    v = dm$vi,  
                    steps = c(0.025, 1),
                    #weights = c(1,0), # weight such that ONLY nonaffirmatives are published
                    table = TRUE ) )

# eta
1/m1[[2]]$par[3]


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


# try the full RTMA via weightr (c.f. "2021-7-4 connection to selection models")
d3$affirm = ( d3$yi / sqrt(d3$vi) ) > qnorm(.975)
( m1 = weightfunct( effect = d3$yi[ d3$affirm == FALSE ],
                    v = d3$vi[ d3$affirm == FALSE ],  
                    steps = c(0.025, 1),
                    weights = c(0,1), # weight such that ONLY nonaffirmatives are published
                    table = TRUE ) )

# adjusted point estimate
m1[[2]]$par[2]


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





# just for comparison:
# allow weightr to find the weights and use all studies
( m1 = weightfunct( effect = d3$yi,
                    v = d3$vi,  
                    steps = c(0.025, 1),
                    #weights = c(1,0), # weight such that ONLY nonaffirmatives are published
                    table = TRUE ) )

# eta
1/m1[[2]]$par[3]


# P-HACKING ADJUSTMENT IN SIMULATED META-ANALYSIS -----------------------------

# no hacking
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


# weightr RTMA
( m1 = weightfunct( effect = dp$yi[ dp$affirm == FALSE ],
                    v = dp$vi[ dp$affirm == FALSE ],
                    steps = c(0.025, 1),
                    weights = c(0,1), # weight such that ONLY nonaffirmatives are published
                    table = TRUE ) )

# adjusted point estimate: -0.20
m1[[2]]$par[2]


# just for comparison:
# allow weightr to find the weights and use all studies
( m1 = weightfunct( effect = dp$yi,
                    v = dp$vi,  
                    steps = c(0.025, 1),
                    #weights = c(1,0), # weight such that ONLY nonaffirmatives are published
                    table = TRUE ) )

# eta
1/m1[[2]]$par[3]




# ~ Try Efron nonparametric method  ---------------------

#bm
# lower and upper trunc limits for each yi
dp$yi.cutoff = dp$tcrit * sqrt(dp$vi)


library(double.truncation)
NPMLE(u.trunc = dp$yi.cutoff,
      y.trunc = dp$yi,
      v.trunc = rep(99, nrow(dp) ) )

# package example from Efron paper
y.trunc=c(0.75, 1.25, 1.50, 1.05, 2.40, 2.50, 2.25)
u.trunc=c(0.4, 0.8, 0.0, 0.3, 1.1, 2.3, 1.3)
v.trunc=c(2.0, 1.8, 2.3, 1.4, 3.0, 3.4, 2.6)
NPMLE(u.trunc,y.trunc,v.trunc)

NPMLE( u.trunc, y.trunc, rep(99, length(y.trunc) ) )


# SEPARATE MLES FOR EACH NONAFFIRMATIVE -----------------------------


# ~ Explore ---------------------


# for the troublesome Hagger example:
d = dm
d$tstat = d$yi / sqrt(d$vi)
d$affirm = d$tstat > crit


# published affirmatives only
dpn = d[ d$affirm == FALSE, ]

# MLE for the first one assuming no within-study heterogeneity
one_nonaffirm_mle( .tstat = dpn$tstat[1],
                   .se = sqrt(dpn$vi[1]),
                   .t2w = 0,
                   .crit = qnorm(.975) )

# MLE for the first one allowing arbitrary within-study heterogeneity
one_nonaffirm_mle( .tstat = dpn$tstat[1],
                   # .se = sqrt(dpn$vi[1]),
                   # .t2w = 0,
                   .crit = qnorm(.975) )

# c.f. observed t-stat
dpn[1,]



# ~ Plot single-study MLEs: t2w = 0 ---------------------

dpn = dpn %>% rowwise() %>%
   mutate( one_nonaffirm_mle(.tstat = tstat,
                             .se = sqrt(vi),
                             .t2w = 0,
                             .crit = qnorm(.975) ) )

#@SAVE THIS IN PROJECT LOG! VERY INTERESTING AND INFORMATIVE.
# t-stats that are close to the truncation point receive by far the most
#  upward correction
#  small ones get little upward correction
ggplot( data = dpn, 
        aes(x = tstat,
            y = tstat.mu.hat) ) +
   geom_abline(slope = 1, 
               intercept = 0) +
   geom_point(shape = 21) +
   theme_bw()

# also look at the effect sizes
# with inverse-variance sizing
ggplot( data = dpn, 
        aes(x = yi,
            y = muiHat,
            size = 1/vi),
        ) +
   
   geom_abline(slope = 1, 
               intercept = 0) +
   geom_point(shape = 21) +
   theme_bw()


# try meta-analyzing the corrected ones
rma.uni( yi = dpn$muiHat, 
         vi = dpn$vi,  #@ultimately will need to include SE of the MLE itself, I think
         method = "REML",
         knha = TRUE)

# 2.82: still very large
#bm
# @meta-analysis MLE was 16.8
# vs. naive: 0.68


# ~ Plot single-study MLEs: t2w unrestricted ---------------------

#bm: seems like these are the same as the above ones?
# maybe also look at how the single-study corrections work in another
#  applied example

dpn = dpn %>% rowwise() %>%
   mutate( one_nonaffirm_mle(.tstat = tstat,
                             # .se = sqrt(vi),
                             # .t2w = 0,
                             .crit = qnorm(.975) ) )

#@SAVE THIS IN PROJECT LOG! VERY INTERESTING AND INFORMATIVE.
# t-stats that are close to the truncation point receive by far the most
#  upward correction
#  small ones get little upward correction
ggplot( data = dpn, 
        aes(x = tstat,
            y = tstat.mu.hat) ) +
   geom_abline(slope = 1, 
               intercept = 0) +
   geom_point(shape = 21) +
   theme_bw()

# also look at the effect sizes
# with inverse-variance sizing
ggplot( data = dpn, 
        aes(x = yi,
            y = muiHat,
            size = 1/vi),
) +
   
   geom_abline(slope = 1, 
               intercept = 0) +
   geom_point(shape = 21) +
   theme_bw()


# try meta-analyzing the corrected ones
rma.uni( yi = dpn$muiHat, 
         vi = dpn$vi,  #@ultimately will need to include SE of the MLE itself, I think
         method = "REML",
         knha = TRUE)

# 2.82: still very large
# @meta-analysis MLE was 16.8
# vs. naive: 0.68


# HUGE CORPUSES INSTEAD OF META-ANALYSES -----------------------------

# Goal: try to avoid small-k issues

# ~ Ioannidis corpus ----------------------------

# IMPORTANT: with this corpus, not clear if people want POSITIVE signs because they're just 
#  scraped, so probably best to consider significance for this 

setwd("~/Dropbox/Personal computer/Independent studies/BOD (believability of disciplines)/bod_git")
s = read.csv("Data from Ioannidis/full_data.csv")
dim(s)  # expected: 23509

# they've coded them such that all signs are positive
table(s$t > 0)
table(s$D > 0)

# Dr, tr are the signed versions

s$affirm = ( s$Dr / s$SE ) > qnorm(.975)
table(s$affirm)

s$signif = abs( s$Dr / s$SE ) > qnorm(.975)
table(s$signif)

# ~~ ECDF of p-values ----------------------------
plot1 = ggplot( data = s,
                aes( x = pval ) ) +
   geom_vline(xintercept = 0.05, color = "red", lwd = 1.5) +
   #geom_vline(xintercept = 0.975, color = "red", lwd = 1.5) +
   geom_histogram( binwidth = 0.025 ) +
   xlab("Two-sided p-value") + 
   theme_classic() +
   theme( panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20, face = "bold"),
          plot.title=element_text(size=24, face = "bold")) +
   ggtitle("Ioannidis corpus")

plot1


# ~~ Get MLEs meta-analysis-style -----------------------

# corpus is too large for memory limits, so let's try downsampling 
# 10^3 is okay locally
set.seed(451)
s2 = s[ sample( nrow(s),
                 size = 10^3,
                 replace = FALSE ), ] 

.obj = correct_meta_phack2(yi = s2$Dr,
                           vi = s2$SE^2)



.obj$sanityChecks$tstatMeanMLE  # 1.09


plot_trunc_densities(.obj)
plot_trunc_densities(.obj, showAffirms = TRUE)

.obj$metaCorr  
.obj$metaNaive # 0.02

# number of nonaffirmatives
.obj$sanityChecks$kNonaffirmPub

# try the full RTMA via weightr (c.f. "2021-7-4 connection to selection models")

( m1 = weightfunct( effect = s2$Dr[ s2$affirm == FALSE ],
                    v = s2$SE[ s2$affirm == FALSE ]^2,  
                    steps = c(0.025, 1),
                    weights = c(0,1), # weight such that ONLY nonaffirmatives are published
                    table = TRUE ) )

# adjusted point estimate: -0.20
m1[[2]]$par[2]


# ~~ Hacking for significance rather than affirmative status -----------------------

( m1 = weightfunct( effect = s2$Dr[ s2$signif == FALSE ],
                    v = s2$SE[ s2$signif == FALSE ]^2,  
                    steps = c(0.025, 0.975, 1),
                    weights = c(0,1,0), # weight such that ONLY nonaffirmatives are published
                    table = TRUE ) )

# adjusted point estimate: -0.31
m1[[2]]$par[2]



# ~~ C.f. weightr when it can estimate the weights itself -----------------

# just for comparison:
# allow weightr to find the weights and use all studies
( m1 = weightfunct( effect = s2$Dr,
                    v = s2$SE^2,  
                    steps = c(0.025, 0.975, 1),
                    #weights = c(1,0), # weight such that ONLY nonaffirmatives are published
                    table = TRUE ) )

# eta
1/m1[[2]]$par[3]









