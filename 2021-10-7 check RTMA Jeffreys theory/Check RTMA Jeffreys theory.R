
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


# ~~ Derivative 1 --------------------
# MATCHES :)

get_D1_num = Deriv(joint_ll_2, ".Mu")

# compare to theoretical one at various (.Mu, .Tt)
res = expand_grid( Mu = c(-1, 0.5, 1),
                   Tt = c(0.1, 1, 2) )

res = res %>% rowwise() %>%
  mutate( num = get_D1_num( .yi = yi, .sei = sei, .Mu = Mu, .Tt = Tt ),
          theory = get_D1( .yi = yi, .sei = sei, .Mu = Mu, .Tt = Tt ) )

expect_equal( res$num, res$theory )


# ~~ Derivative 12 --------------------

### Check all of D12
# DOESN'T MATCH :(

get_D12_num = Deriv( get_D1_num, ".Tt")

# compare to theoretical one at various (.Mu, .T2t)
res = expand_grid( Mu = c(-1, 0.5, 1),
                   Tt = c(0.1, 1, 2) )

res = res %>% rowwise() %>%
  mutate( num = get_D12_num( .yi = yi, .sei = sei, .Mu = Mu, .Tt = Tt ),
          theory = get_D12( .yi = yi, .sei = sei, .Mu = Mu, .Tt = Tt ) )

expect_equal( res$num, res$theory )

# debug
.yi = yi; .sei = sei; .Mu = 1; .Tt = 2; .T2t = .Tt^2
get_D12( .yi = yi, .sei = sei, .Mu = Mu, .Tt = Tt )
get_D12_num( .yi = yi, .sei = sei, .Mu = Mu, .Tt = Tt )

gamma.i = get_gamma_i( .yi, .sei, .Mu, .Tt, .crit )
Zi.tilde = (.yi - .Mu) / sqrt(.T2t + .sei^2)

# innards of get_D12_num
.e2 <- .sei^2 + .Tt^2
.e3 <- .yi - .Mu
.e4 <- sqrt(.e2)
.e5 <- .e3/.e4  
.e6 <- pnorm(.e5)
.e7 <- .e2^2
.e8 <- dnorm(.e5, 0, 1)
-(0.5 * sum(4 * (.Tt * .e3/.e7)) + sum(.Tt * ((.e6/.e4 - 
                                                 .e3 * .e8/.e2)/(.e6 * .e4)^2 - .e3^2/(.e7 * .e6 * .e4)) * 
                                         .e8))


.e5 == Zi.tilde
.e8/.e6 == gamma.i

# check partial term - MATCHES :)
-0.5 * sum(4 * (.Tt * .e3/.e7)) == -2*sum( .Tt*(.yi - .Mu)/(.T2t + .sei^2)^(2) )


### Intermediate quantity: d Zi.tilde / d tau
# MATCHES :)
# .Mu = 1
# .Tt = 2
temp1 = Deriv( get_Zi_tilde, ".Tt")( .yi, .sei, .Mu = .Mu, .Tt = .Tt ) 
temp2 = -0.5*(.Tt^2 + .sei^2)^(-3/2) * (.yi - .Mu) * 2*.Tt
expect_equal( temp1, temp2)

### Intermediate quantity: d gamma_i / d tau
# .Mu = 1
# .Tt = 2
# .T2t = .Tt^2
temp1 = Deriv( get_gamma_i, ".Tt")( .yi, .sei, .Mu = .Mu, .Tt = .Tt ) 
temp2 = get_D_gammai_wrt_tau( .yi, .sei, .Mu = .Mu, .Tt = .Tt ) 
expect_equal( temp1, temp2)


### Intermediate quantity: d term2 / d tau
# term2 = gamma.i * (.T2t + .sei^2)^(-1/2)
# BREAKS!

get_D12_term2_wrt_tau_num = Deriv( term2, ".Tt")

temp1 = get_D12_term2_wrt_tau( .yi, .sei, .Mu = .Mu, .Tt = .Tt ); sum(temp1)
temp2 = get_D12_term2_wrt_tau_num( .yi, .sei, .Mu = .Mu, .Tt = .Tt ); sum(temp2)
expect_equal( temp1, temp2)



# check innards of get_D12_term2_wrt_tau_num
gamma.i = get_gamma_i( .yi, .sei, .Mu, .Tt, .crit )
Zi.tilde = (.yi - .Mu) / sqrt(.T2t + .sei^2)

.e2 <- .sei^2 + .Tt^2
.e3 <- .yi - .Mu
.e4 <- sqrt(.e2)
.e5 <- .e3/.e4  # Zi.tilde
.e6 <- pnorm(.e5)
# ADDED SUM JUST FOR CHECKING PURPOSES
sum( .Tt * (2 * (.e3 * (.e5 + dnorm(.e5)/.e6) * dnorm(.e5, 0, 1)/.e6) - get_gamma_i(.yi, .sei, .Mu, .Tt, .crit)/.e4)/.e2 )

dnorm(.e5)/.e6 == gamma.i
.e5 == Zi.tilde

# simplify the above (from get_D12_term2_wrt_tau_num) to see what differs
sum( .Tt * ( 2 * ( (.yi - .Mu) * ( Zi.tilde*gamma.i + gamma.i^2 ) ) - gamma.i*(.sei^2 + .Tt^2)^(-1/2) ) / (.sei^2 + .Tt^2) )

# simplify further
# first part of this has an extra "2" coefficient and also ^(-1) instead of ^(-2)
# ****BM: WTF??????
sum( ( 2*(.yi - .Mu)*(.sei^2 + .Tt^2)^(-1)*( Zi.tilde*gamma.i + gamma.i^2 )*.Tt ) -
       gamma.i*(.sei^2 + .Tt^2)^(-3/2) * .Tt )



# mine again
yi = .yi
sei = .sei
Mu = .Mu
Tt = .Tt
T2t = .Tt^2

#rm( list = c(".yi", ".sei", ".Mu", ".Tt", ".T2t") ) 



gamma.i = get_gamma_i( .yi = yi, .sei = sei, .Mu = Mu, .Tt = Tt )

# "true" Z-score (including .T2t)
Zi.tilde = (yi - Mu) / sqrt(T2t + sei^2)


term2.1 = get_D_gammai_wrt_tau( .yi = yi, .sei = sei, .Mu = Mu, .Tt = Tt )
#***mine
line1 = term2.1*(T2t + sei^2)^(-1/2) - 0.5*(T2t + sei^2)^(-3/2) * 2*Tt * gamma.i

sum(line1)

