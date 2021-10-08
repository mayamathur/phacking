
# ~ PRELIMINARIES ----------------------------------------------------

rm(list=ls())


# data-wrangling packages
library(here)
library(dplyr)
library(tibble)
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
library(Deriv)
library(mosaic)
library(hpa)
library(pracma)
library(truncnorm)
library(tmvtnorm)
library(Hmisc)
library(truncdist)


prepped.data.dir = "~/Dropbox/Personal computer/Miscellaneous/Journal article library/Reproducibility:replicability/Kvarven comparing meta-analysis and multisite replications/Reanalyze Kvarven/Prepped Hagger data (meta and replication)"
res.dir = "~/Dropbox/Personal computer/Miscellaneous/Journal article library/Reproducibility:replicability/Kvarven comparing meta-analysis and multisite replications/Reanalyze Kvarven/Hagger comparison results"


#code.dir = here("2021-10-7 check RTMA Jeffreys theory")

code.dir = here()
setwd(code.dir)
source("helper_SAPH.R")



# ~ SET PARAMETERS AND SIMULATE META-ANALYSIS ----------------------------------------------------

# Mu = 1
# T2t = 2
# m = 50
# se = 
# 
# Nmax = 1
# 
# d = sim_meta(Nmax = .p$Nmax,
#              Mu = .p$Mu,
#              T2 = .p$T2,
#              m = .p$m,
#              t2w = .p$t2w,
#              se = .p$se,
#              hack = .p$hack,
#              return.only.published = FALSE,
#              rho = .p$rho,
#              
#              k = .p$k,
#              k.hacked = .p$k.hacked )

# use Hagger replications, for which MLE was horrible

setwd(prepped.data.dir)
dm = read.csv("prepped_hagger_meta_data.csv")

# .obj = correct_meta_phack2(yi = dm$yi,
#                            vi = dm$vi)
# 
# # t-stat MLE is HUGE
# .obj$sanityChecks$tstatMeanMLE
# 
# # **plot shows that the MLE is so big because the nonaffirmative t-stats are actually left-skewed
# # very interesting
# plot_trunc_densities(.obj)

yi = dm$yi
sei = sqrt(dm$vi)


# ~ CHECK VS. NUMERICAL DIFFERENTIATION ----------------------------------------------------


### entry 1 - MATCHES :)
get_D1_num = Deriv(joint_ll_2, ".Mu")

# compare to theoretical one at various (.Mu, .T2t)
res = expand_grid( Mu = c(-1, 0.5, 1),
                   T2t = c(0.1, 1, 2) )

res = res %>% rowwise() %>%
  mutate( num = get_D1_num( .yi = yi, .sei = sei, .Mu = Mu, .T2t = T2t ),
          theory = get_D1( .yi = yi, .sei = sei, .Mu = Mu, .T2t = T2t ) )

expect_equal( res$num, res$theory )



### entry 11 

# first check an intermediate quantity, d Zi.tilde / d tau - MATCHES
temp1 = Deriv( get_Zi_tilde, ".Tt")( .yi, .sei, .Mu = 1, .Tt = 1 ) 
.Mu = 1
.Tt = 1
temp2 = -0.5*(.Tt^2 + .sei^2)^(-3/2) * (.yi - .Mu) * 2*.Tt
expect_equal( temp1, temp2)


get_D12_num = Deriv( get_D1_num, ".Tt")



# compare to theoretical one at various (.Mu, .T2t)
res = expand_grid( Mu = c(-1, 0.5, 1),
                   Tt = c(0.1, 1, 2) )

res = res %>% rowwise() %>%
  mutate( num = get_D12_num( .yi = yi, .sei = sei, .Mu = Mu, .Tt = Tt ),
          theory = get_D12( .yi = yi, .sei = sei, .Mu = Mu, .Tt = Tt ) )

expect_equal( res$num, res$theory )






