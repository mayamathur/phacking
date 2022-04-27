

# This script checks partial derivatives wrt tau (not the total variance) for Jeffreys prior.
# It checks numerical derivatives from R (used in init_stan_model) to those from Mathematica 
#  (see 2022-4-25 Check Mathematica using R).
# It also compares theory I did manually ("iPad" below) to Mathematica/R.


# IMPORTANT: THROUGHOUT, CRIT IS THE 1.96 VALUE, NOT THE CUTOFF FOR EACH STUDY (CRIT*SEI)

# HELPERS ------------------

library(truncnorm)

# data-wrangling packages
library(dplyr)
library(tibble)
library(ggplot2)
library(data.table)
library(stringr)
library(tidyverse)
library(fastDummies)
# meta-analysis packages
library(metafor)
library(robumeta)
# other
library(here)
library(xtable)
library(testthat)
library(Deriv)

setwd("~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Linked to OSF (SAPH)/Code (git)/Sherlock code")
source("helper_SAPH.R")


# ~ New likelihood, keeping .tau and .sei separate ------------------
simple_lli = function(.yi,
                      .sei,
                      .mu,
                      .tau,  
                      .crit ) {
  
  # # # TEST
  # .mu=0.1
  # .tau=1
  # .sei=0.5
  # .crit=1.96
  # .yi=0.5
  
  Si = sqrt(.tau^2 + .sei^2)
  alphaU = (.crit*.sei - .mu)/Si
  
  # sanity check:
  # mine = -log( Si* sqrt(2*pi) ) - 0.5 * Si^(-2) * (.yi-.mu)^2 - log( pnorm(alphaU) )
  # termA = dnorm( .yi,
  #                mean = .mu,
  #                sd = Si,
  #                log = TRUE)
  # 
  # termB = pnorm( q = alphaU,
  #                mean = 0,
  #                sd = 1,
  #                log.p = TRUE )
  # 
  # termA - termB
  # 
  # expect_equal(termA - termB, mine) 
  # 
  # # also check vs. dtruncnorm:
  # expect_equal( log( dtruncnorm(x = .yi,
  #            a = -Inf,
  #            b = .crit*.sei,
  #            mean = .mu,
  #            sd = Si) ), mine )
  
  # same as "mine" above
  -log( Si* sqrt(2*pi) ) - 0.5 * Si^(-2) * (.yi-.mu)^2 - log( pnorm(alphaU) )
}


# ~ New prior using the body of derivative functions defined below ------------------
# closely translated from stan
prior = function(mu, tau, k, sei, tcrit) {
  
  # this will be the TOTALS for all observations
  fishinfototal = matrix( 0, nrow = 2, ncol = 2 )
  
  # build a Fisher info matrix for EACH observation
  for (i in 1:k) {
    
    # for this observation
    fishinfo = matrix( NA, nrow = 2, ncol = 2 )
    
    # from body of R's get_D11_num:
    e2 = sei[i]^2 + tau^2
    e3 = sqrt(e2)
    e5 = sei[i] * tcrit[i] - mu
    e6 = e5/e3
    e7 = dnorm(e6, 0, 1)
    # Stan version:
    # e7 = exp( normal_lpdf(e6 | 0, 1) )
    e8 = pnorm(e6)
    #e8 = exp( normal_lcdf(e6 | 0, 1 ) )
    kmm = -(1/e2 - (e5/(e2 * e8) + e7 * e3/(e8 * e3)^2) * e7/e3)
    
    # from body of R's get_D12_num:
    e2 = sei[i]^2 + tau^2
    e3 = sqrt(e2)
    e5 = sei[i] * tcrit[i] - mu
    # e6 is scaled critical value:
    e6 = e5/e3
    e7 = pnorm(e6)
    # e7 = exp( normal_lcdf(e6 | 0, 1 ) )
    e8 = e2^2
    e9 = dnorm(e6, 0, 1)
    #e9 = exp( normal_lpdf(e6 | 0, 1) )
    
    # my own expectation of .yi - .mu:
    expectation1 = -sqrt(sei[i]^2 + tau^2) * e9/e7
    kms = -(tau * (((e7/e3 - e5 * e9/e2)/(e7 * e3)^2 - e5^2/(e8 *
                                                               e7 * e3)) * e9 + 2 * ((expectation1)/e8)))
    
    
    # from body of R's get_D22_num:
    e1 = tau^2
    e3 = sei[i]^2 + e1
    e5 = sei[i] * tcrit[i] - mu
    e6 = sqrt(e3)
    # e7 is scaled crit value:
    e7 = e5/e6
    e8 = pnorm(e7)
    # e8 = exp( normal_lcdf(e7 | 0, 1 ) )
    e9 = dnorm(e7, 0, 1)
    # e9 = exp( normal_lpdf(e7 | 0, 1 ) )
    e10 = e5 * e9
    e11 = e8 * e6
    e13 = e10/e11
    # *replace this one with its expectation:
    # e15 = (.yi - .mu)^2/e3
    # expectation of (.yi - .mu)^2:
    expectation2 = (sei[i]^2 + tau^2)*(1 - e7 * e9/e8)
    e15 = (expectation2)/e3
    
    kss = (e13 + e15 - (e1 * (e5 * ((e8/e6 - e10/e3)/e11^2 -
                                      e5^2/(e3^2 * e8 * e6)) * e9 + 2 * ((e13 + 2 * e15 -
                                                                            1)/e3)) + 1))/e3
    
    
    fishinfo[1,1] = -kmm
    fishinfo[1,2] = -kms
    fishinfo[2,1] = -kms
    fishinfo[2,2] = -kss
 
    # add the new fisher info to the total one
    fishinfototal = fishinfototal + fishinfo
  }
  
  return( list(Efish = fishinfototal,
               det = det(fishinfototal),
               prior = sqrt( det(fishinfototal) ) ) )
}


# MATHEMATICA VS. R: CHECK MY FN "lkl" IN MATHEMATICA ------------------

mu=0.1
tau=1
sei=0.5
crit=1.96
y=0.5

check = log( dtruncnorm( x = y,
                         mean = mu,
                         sd = sqrt(tau^2 + sei^2),
                         a = -Inf,
                         b = s*crit ) )

# Mathematica: -0.851658
expect_equal( check, -0.851658, tol = 0.001)

# also check my simple R fn
expect_equal( check,
              simple_lli(.yi = y,
                         .sei = sei,
                         .mu = mu,
                         .tau = tau,
                         .crit = crit)
              ,
              tol = 0.001)


# MATHEMATICA VS. R: PARTIAL DERIVATIVES -------------------------------

# all get_D_XX fns below match Mathematica! :)
# c.f. "2022-4-25 Fisher info check Mathematica"

# ~ D11 ------------------------------------------
setwd("~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Linked to OSF (SAPH)/Code (git)/Sherlock code")
source("helper_SAPH.R")


mu=0.1
tau=1
sei=0.5
crit=1.96
y=0.5


get_D1_num = Deriv(simple_lli, ".mu")

get_D11_num = Deriv(get_D1_num, ".mu")

expect_equal( get_D11_num(.yi = y,
                          .sei = sei,
                          .mu = mu,
                          .tau = tau,
                          .crit = crit),
              -0.453676,
              tol = 0.001)


# ~ D12 ------------------------------------------

get_D12_num = Deriv(get_D1_num, ".tau")

expect_equal( get_D12_num(.yi = y,
                          .sei = sei,
                          .mu = mu,
                          .tau = tau,
                          .crit = crit),
              -0.535173,
              tol = 0.001)

# ~ D22 ------------------------------------------

get_D2_num = Deriv(simple_lli, ".tau")

get_D22_num = Deriv(get_D2_num, ".tau")

# vs. Mathematica:
expect_equal( get_D22_num(.yi = y,
                          .sei = sei,
                          .mu = mu,
                          .tau = tau,
                          .crit = crit),
              0.0974377,
              tol = 0.001)



# PLOT PRIOR CONTOURS ------------------------------------

# ~ At various (mu, tau) --------------------------------------
sei = c(0.534207749800569, 0.454756736242383, 0.725584888897609, 1.03488267282425, 
        0.295282394138296, 0.449484477194935, 1.05236633781316, 0.300110117200045, 
        0.61423644667798, 0.435493032872841)
mu = 0.5
tau = 0.28
k = length(sei)
tcrit = rep(2, k)

dp = expand_grid(.mu = seq(-2, 2, .05),
                 .tau = seq(0, 2, .05))

# make plotting dataframe by calculating log-prior for a grid of values (mu, sigma)
dp = dp %>%
  rowwise() %>%
  mutate( prior.val = prior(mu = .mu,
                            tau = .tau,
                            k = k,
                            sei = sei,
                            tcrit = tcrit )$prior )


# check again for NA values occurring when EFisher is NA and remove them
table(is.na(dp$prior.val))

# set up colors for contours
get_colors = colorRampPalette( c("lemonchiffon1", "chocolate4") )
myColors = get_colors(n=15)  # chose 11 based on errors from ggplot if it was fewer

# make plot
p = ggplot( data = dp, 
            aes(x = .mu,
                y = .tau,
                z = prior.val) ) +
  
  geom_contour_filled() +
  
  # close, but not enough colors
  scale_fill_manual(values = myColors) +
  
  geom_contour(color = "white") +
  
  xlab( bquote(mu) ) +
  ylab( bquote(tau) ) +
  
  geom_vline( xintercept = 0, lty = 2 ) +
  
  scale_y_continuous(breaks = seq( min(dp$.tau), max(dp$.tau), 0.25),
                     limits = c( min(dp$.tau), max(dp$.tau) ) ) +
  
  theme_bw(base_size = 16) +
  theme(text = element_text(face = "bold"),
        axis.title = element_text(size=20),
        legend.position = "none")

p



# ~ Cross-section at one value of mu --------------------------------------
ggplot( data = dp %>% filter(.mu == 0.5),
        aes(x = .tau,
            y = prior.val) ) +
  geom_line() +
  
  xlab( bquote(tau) ) +
  ylab( "Prior" ) +
  
  ggtitle( paste("Mu = 0.5" ) ) +
  
  geom_vline( xintercept = 0, lty = 2 ) +
  
  scale_x_continuous(breaks = seq( min(dp$.tau), max(dp$.tau), 0.25),
                     limits = c( min(dp$.tau), max(dp$.tau) ) ) +
  
  theme_bw(base_size = 16) +
  theme(text = element_text(face = "bold"),
        axis.title = element_text(size=20),
        legend.position = "none")





# CONTOUR PLOT FOR DIFFERENT VALUES OF SEI ------------------------------


# the SE distribution
x = 0.02 + rexp(n = 200, rate = 3)

dp = expand_grid(.mu = seq(-2, 2, .05),
                 .tau = seq(0.05, 2, .05),
                 .sei = c(0.01, 0.2, 0.4, 2) )
tcrit = 1.96

# make plotting dataframe by calculating log-prior for a grid of values (mu, sigma)
dp = dp %>%
  rowwise() %>%
  mutate( prior.val = prior(mu = .mu,
                            tau = .tau,
                            k = 1,
                            sei = .sei,
                            tcrit = tcrit )$prior,
          
          prior.val.old = prior_Si(mu = .mu,
                                   Si = sqrt(.tau^2 + .sei^2),
                                   k = 1,
                                   sei = .sei,
                                   tcrit = tcrit )$prior )

# check again for NA values occurring when EFisher is NA and remove them
table(is.na(dp$prior.val)); table(is.na(dp$prior.val.old))

View( dp %>% filter(is.na(prior.val) ) )

#**if both sei and tau are tiny, Fisher info can be all negative infinities
#*# because of numerical issues with the dnorm and pnorm terms being nearly 0
prior(mu = 2,
      tau = 0.05,
      k = 1,
      sei = 0.02,
      tcrit = tcrit )

# set up colors for contours
get_colors = colorRampPalette( c("lemonchiffon1", "chocolate4") )
myColors = get_colors(n=15)  # chose 11 based on errors from ggplot if it was fewer

# make plot
dp.temp = dp %>% filter(.sei == 2)

p = ggplot( data = dp.temp, 
            aes(x = .mu,
                y = .tau,
                z = prior.val) ) +
  
  geom_contour_filled() +
  
  # close, but not enough colors
  scale_fill_manual(values = myColors) +
  
  geom_contour(color = "white") +
  
  xlab( bquote(mu) ) +
  ylab( bquote(tau) ) +
  
  geom_vline( xintercept = 0, lty = 2 ) +
  
  # scale_y_continuous(breaks = seq( min(dp.temp$.tau), max(dp.temp$.tau), 0.25),
  #                    limits = c( min(dp.temp$.tau), max(dp.temp$.tau) ) ) +
  
  theme_bw(base_size = 16) +
  theme(text = element_text(face = "bold"),
        axis.title = element_text(size=20),
        legend.position = "none")  +
  
  facet_wrap( as.factor(.sei) ~ .,
              scales = "free" )

p


