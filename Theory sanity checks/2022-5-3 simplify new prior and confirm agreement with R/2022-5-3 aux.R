

# This script checks partial derivatives wrt tau (not the total variance) for Jeffreys prior.
# It checks numerical derivatives from R (used in init_stan_model) to those from Mathematica 
#  (see 2022-4-20 Check Mathematica using R).
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


setwd("~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (not on OSF)/Aux code/**2022-5-3 simplify new prior and confirm agreement with R")
source("2022-5-3 helper.R")


# GET READY TO MANUALLY SIMPLIFY PRIOR ---------------------------

# set up some test values
#tau = 0
tau = 0.2
mu = 1.65
sei = 0.53
tcrit = 2
i = 1


# the following simplifications that should always hold regardless of tau, etc.
Si = sqrt( tau^2 + sei[i]^2 )
cz = (sei[i] * tcrit[i] - mu) / Si
dnor = dnorm(cz) # can't use log on this in case it's negative
pnor = pnorm(cz)
r = dnor/pnor


# SIMPLIFY KSS TERM ---------------------------

# what follows is from the body of prior():
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



expect_equal( e15, (1 - cz*r) )

expect_equal( e13, cz*r)

expect_equal( e10, cz*Si*dnorm(cz) )


# start simplifying subparts of kss

### Term A: (e8/e6 - e10/e3)/e11^2
(e8/e6 - e10/e3)/e11^2

( (pnor/Si) - cz*Si*dnor*Si^(-2) ) / ( pnor^2 * Si^(2) )

(1/Si) * (pnor - cz*dnor) / (Si^2 * pnor^2)
# ****
termA = Si^(-3) * pnor^(-1) * (1 - cz*r)  # MATCHES! 


# intermediate check
my.kss = (1 - (e1 * (e5 * (termA -
                             e5^2/(e3^2 * e8 * e6)) * e9 + 2 * ((e13 + 2 * e15 -
                                                                   1)/e3)) + 1))/e3
expect_equal( kss, my.kss)

### Term B: e5^2/(e3^2 * e8 * e6)
e5^2/(e3^2 * e8 * e6)

cz^2*Si^2 / (Si^4 * pnor * Si)
# ***
termB = cz^2 * pnor^(-1) * Si^(-3)


### Term C: ((e13 + 2 * e15 - 1)/e3)
(e13 + 2 * e15 - 1)/e3
termC = (1 - cz*r)/Si^2

### Term D: The huge, messy chunk

(e1 * (e5 * ((e8/e6 - e10/e3)/e11^2 -
               e5^2/(e3^2 * e8 * e6)) * e9 + 2 * ((e13 + 2 * e15 -
                                                     1)/e3)) + 1)

e1 * (e5 * (termA - termB) * e9 + 2 * termC) + 1

(tau^2)*( cz*Si*dnor*(termA - termB) + 2*termC ) + 1

(tau^2)*cz*Si*dnor*(termA - termB) + 2*(tau^2)*termC + 1

(tau^2)*cz*Si^(-2)*r * ( (1-cz*r) - cz^2 ) + 2*(tau^2)*termC + 1

( tau^2 * Si^(-2) ) * ( cz*r*( (1-cz*r) - cz^2 ) + 2*(1-cz*r) ) + 1

termD = -( tau^2 * Si^(-2) ) * ( cz^3*r + cz^2*r^2 + cz*r - 2 ) + 1


# intermediate check
my.kss = (1 - (termD))/e3
expect_equal( kss, my.kss)


### Whole kss

expect_equal( 1-termD,
              ( tau^2 * Si^(-2) ) * ( cz^3*r + cz^2*r^2 + cz*r - 2 ) )

# FINAL SIMPLIFICATION :) 
# WOOHOO
my.kss = ( tau^2 * Si^(-4) ) * ( cz^3*r + cz^2*r^2 + cz*r - 2 )

expect_equal(my.kss, kss)



# SIMPLIFY KMS TERM ---------------------------


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


### Term A

(e7/e3 - e5 * e9/e2)/(e7 * e3)^2

( pnor*Si^(-1) - cz*Si*dnor*(Si^-2) ) / ( pnor^2 * Si^(2) )

pnor^(-1)*Si^(-3) - cz*Si^(-3)*r*pnor^(-1)

termA = pnor^(-1)*Si^(-3) * (1 - cz*r)


# intermediate check
my.kms = -(tau * (( termA - e5^2/(e8 * e7 * e3) ) * e9 + 2 * ((expectation1)/e8)))
expect_equal( kms, my.kms)


### Term B

e5^2/(e8 * e7 * e3)

termB = Si^(-3) * cz^2 * pnor^(-1)

# intermediate check
my.kms = -( tau * ( ( termA - termB ) * e9 + 2 * ( (expectation1)/e8 ) ) )
expect_equal( kms, my.kms)


### Term C

( termA - termB ) * e9

termC = Si^(-3)*r*(1 - cz*r - cz^2)


# intermediate check
my.kms = -( tau * ( termC + 2 * ( (expectation1)/e8 ) ) )
expect_equal( kms, my.kms)


### Whole thing

-( tau * ( termC + 2 * ( (expectation1)/e8 ) ) )

-tau * ( termC + 2*expectation1/e8 )

-tau*Si^(-3)*r*( -1 - cz*r - cz^2 )

my.kms = tau*Si^(-3)*r*( cz^2 + cz*r + 1 )
expect_equal( kms, my.kms)


# SIMPLIFY KMM TERM ---------------------------

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



### Term A

(e5/(e2 * e8) + e7 * e3/(e8 * e3)^2)


termA = Si^(-1)*pnor^(-1)*(cz + r)


# intermediate check
my.kmm = -(1/e2 - termA * e7/e3)
expect_equal( kmm, my.kmm)

### Whole thing

Si^(-2) + Si^(-2)*r*(cz+r)

my.kmm = Si^(-2)*(cz*r + r^2 - 1)

expect_equal( kmm, my.kmm)


# CONFIRM AGREEMENT AT DIFFERENT VALUES ---------------------------


sei = c(0.534207749800569, 0.454756736242383, 0.725584888897609, 1.03488267282425, 
        0.295282394138296, 0.449484477194935, 1.05236633781316, 0.300110117200045, 
        0.61423644667798, 0.435493032872841)
# mu = 0.5
# tau = 0.28
k = length(sei)
tcrit = rep(2, k)

dp = expand_grid(.mu = seq(-2, 2, .05),
                 .tau = seq(0, 2, .05))


# check behavior for first row
prior(mu = dp$.mu[1],
      tau = dp$.tau[1],
      k = k,
      sei = sei,
      tcrit = tcrit )$prior

prior_simp(mu = dp$.mu[1],
           tau = dp$.tau[1],
           k = k,
           sei = sei,
           tcrit = tcrit )$prior

# make plotting dataframe by calculating log-prior for a grid of values (mu, sigma)
dp = dp %>%
  rowwise() %>%
  mutate( prior1 = prior(mu = .mu,
                         tau = .tau,
                         k = k,
                         sei = sei,
                         tcrit = tcrit )$prior,
          
          prior2 = prior_simp(mu = .mu,
                              tau = .tau,
                              k = k,
                              sei = sei,
                              tcrit = tcrit )$prior, )


# prior is still sometimes NA at tau=0 values 
table(is.na(dp$prior1))
View( dp %>% filter(is.na(prior1)))
which( is.na(dp$prior1) )

# look at agreement
expect_equal( dp$prior1, dp$prior2, tol = 0.0001 )


# CONTOUR PLOT ---------------------------

# set up colors for contours
get_colors = colorRampPalette( c("lemonchiffon1", "chocolate4") )
myColors = get_colors(n=15)  # chose 11 based on errors from ggplot if it was fewer


p = ggplot( data = dp, 
            aes(x = .mu,
                y = .tau,
                # simplified prior
                z = prior2) ) +
  
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


# ~ At what value of tau does it peak for various choices of sei? ---------------------


sei = c(0.534207749800569, 0.454756736242383, 0.725584888897609, 1.03488267282425, 
        0.295282394138296, 0.449484477194935, 1.05236633781316, 0.300110117200045, 
        0.61423644667798, 0.435493032872841)

prior_contour_plot(sei)

# larger seis => favors larger taus
prior_contour_plot(sei = 0.05)  # favors tau = 0.05
prior_contour_plot(sei = 1)  # favors tau = 0.7

# favors tau = 0.05
prior_contour_plot(sei = runif(n=50, min = 0.02, max = 1.5) )

# exponential seis, as in sims
# favors tau = 0.10
p = prior_contour_plot(sei = 0.1 + rexp(n = 100, rate = 1.5) )
p + ggtitle("sei ~ 0.1 + Exp(1.5)")

# favors tau = 0.05
p = prior_contour_plot(sei = rbeta(n = 100, 2, 5) )
p + ggtitle("sei ~ Beta(2,5)")


