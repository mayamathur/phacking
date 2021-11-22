
# Goal: Try to get Jeffreys prior numerically because the expectations will be hard


# ~ PRELIMINARIES ----------------------------------------------------

#rm(list=ls())


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


# ~ GET EXPECTED FISHER NUMERICALLY ----------------------------------------------------


.Mu = 0.5
.Tt = 0.2
.yi = yi
.sei = sei

E_fisher_RTMA_OLD( .sei = sei, .Mu = 0.5, .Tt = 0.2 )


# ~ SANITY CHECK: PRIOR SHD BE SPECIAL CASE OF TNE PRIOR ----------------------------------------------------

# if sei -> 0, priors should agree



Mu = 1
Tt = 2
k = 1000

# by setting se super small, we're effectively just truncating at yi > 0 
#  since all studies with yi > 0 are significant 
dp = sim_meta(Nmax = 1,
             Mu = Mu,
             T2 = Tt^2,
             m = 100,
             t2w = 0,
             se = 0.0001,
             hack = "affirm",
             return.only.published = TRUE,
             rho = 0,

             k = k,
             k.hacked = k )

hist(dp$yi)
kpub = nrow(dp)

( EFish.TNE = E_fisher_TNE( .mu = Mu,
                          .sigma = Tt, 
                          .n = kpub,
                          .a = -99,
                          .b = 0 ) )

# this one isn't even def'd
( Efish.RTMA = E_fisher_RTMA_OLD( .sei = sqrt(dp$vi),
                              .Mu = Mu,
                              .Tt = Tt ) )

# hmmm...these don't agree :(

# # look at its guts
# integrand_Dij(i = 1, 
#               j = 1, 
#               .yi = dp$yi[1],
#               .sei = sqrt(dp$vi[1]),
#               .Mu = Mu,
#               .Tt = Tt)
# 
# 
# integrand_Dij(i = 1, 
#               j = 1, 
#               .yi = -99,  # here's the problem!
#               .sei = sqrt(dp$vi[1]),
#               .Mu = Mu,
#               .Tt = Tt)
# 
# integrand_Dij(i = 1, 
#               j = 1, 
#               .yi = 0,  # here's the problem!
#               .sei = sqrt(dp$vi[1]),
#               .Mu = Mu,
#               .Tt = Tt)
# 
# 
# integrate( function(x) integrand_Dij(i = 1, 
#                                      j = 1, 
#                                      .yi = x,
#                                      .sei = sqrt(dp$vi[1]),
#                                      .Mu = Mu,
#                                      .Tt = Tt),
#            lower = -99,
#            upper = 99 )$value
# 
# # test case
# integrand_Dij(i = 1, 
#               j = 1, 
#               .yi = 3,
#               .sei = 0.0001,
#               .Mu = 1,
#               .Tt = 2)



# ~~ Alternative: Use TNE version for each one ----------------------------------------------------

#bm
# idea of using just one of these!!

# Fisher info for each observation separately, based on its unique SE
E_fishers = lapply( X = as.list( dp$vi ),
                    FUN = function(v) {
                      E_fisher_TNE( .mu = Mu,
                                    .sigma = sqrt( Tt^2 + v ), 
                                    .n = 1,
                                    .a = -99,
                                    .b = 0 )
                    })


# add all the matrices entrywise
# https://stackoverflow.com/questions/11641701/sum-a-list-of-matrices
( Efish.TNE.2 = Reduce('+', E_fishers) )

expect_equal(Efish.TNE.2, Efish.TNE)

# E_fisher_TNE( .mu = Mu,
#               .sigma = Tt, 
#               .n = 1,
#               .a = -99,
#               .b = 0 )

# now maybe try a case where SEs are different across studies


# ~ SANITY CHECK: AGREEMENT WITH MLE ----------------------------------------------------

# check that not using the prior agrees with mle

# you already have the usePrior argument :)


# ~ SANITY CHECK: PLOT PRIOR ----------------------------------------------------



# # STRAIGHT FROM TNE
# if ( redo.contour.plot == TRUE ) {
#   # figure out how low we can go with sigma before getting NAs in lkl
#   # because geom_contour will behave badly if there are too many NAs (a few is okay)
#   lprior(.sei = 1,
#          .Mu = 0,
#          .Tt = 3)
#   
#   # set parameters for all plot facets
#   prop.retained = .5
#   n = 20
#   # underlying mu and sigma, only used for calculating cut points on raw scale
#   mu = 0
#   sigma = 1
#   
#   dp = expand.grid( .n = c(n),
#                     .mu = seq(-2, 8, .1),
#                     .sigma = seq(0.5, 2, .1),
#                     .trunc.type = c("Single", "Symmetric double", "Asymmetric double", "No truncation") )
#   
#   # not facetting on n anymore
#   # dp$n.pretty = paste( "n = ", dp$.n, sep = "" )
#   # dp$n.pretty = factor( dp$n.pretty, levels = c("n = 10",
#   #                                               "n = 50",
#   #                                               "n = 1000") )
#   
#   # calculate cutpoints that all retain same proportion ASSUMING that truth is mu=0, sigma=1
#   # only need to calculate once for each truncation type
#   ( single.cuts = calculate_cutpoints(.trunc.type = "Single",
#                                       .prop.retained = prop.retained) )
#   ( double.symm.cuts = calculate_cutpoints(.trunc.type = "Symmetric double",
#                                            .prop.retained = prop.retained) )
#   ( double.asymm.cuts = calculate_cutpoints(.trunc.type = "Asymmetric double",
#                                             .prop.retained = prop.retained) )
#   
#   
#   # set cutpoints on raw scale equal to the ones on the Z-scale
#   # because, again, we assumed mu=0, sigma=1 above to hold constant prop.retained
#   dp$a = dp$b = NA
#   dp$a[ dp$.trunc.type == "Single" ] = single.cuts$Za
#   dp$a[ dp$.trunc.type == "Symmetric double" ] = double.symm.cuts$Za
#   dp$a[ dp$.trunc.type == "Asymmetric double" ] = double.asymm.cuts$Za
#   #dp$a[ dp$.trunc.type == "Asymmetric double" ] = -4 #@TEST ONLY
#   
#   dp$b[ dp$.trunc.type == "Single" ] = single.cuts$Zb
#   dp$b[ dp$.trunc.type == "Symmetric double" ] = double.symm.cuts$Zb
#   dp$b[ dp$.trunc.type == "Asymmetric double" ] = double.asymm.cuts$Zb
#   #dp$b[ dp$.trunc.type == "Asymmetric double" ] = -0.25 #@TEST ONLY
#   
#   dp$a[ dp$.trunc.type == "No truncation" ] = -Inf
#   dp$b[ dp$.trunc.type == "No truncation" ] = Inf
#   
#   table( dp$.trunc.type, round(dp$a, 2) )
#   table( dp$.trunc.type, round(dp$b, 2) )
#   
#   # prettify trunc type variable
#   dp$trunc.type.pretty = paste( dp$.trunc.type,
#                                 " [",
#                                 round(dp$a, 2),
#                                 ", ",
#                                 round(dp$b, 2),
#                                 "]",
#                                 sep = "" )
#   current.levels = unique(dp$trunc.type.pretty)
#   dp$trunc.type.pretty = factor( dp$trunc.type.pretty,
#                                  levels = c( current.levels[ grepl(pattern = "No truncation", x = current.levels) ],
#                                              current.levels[ grepl(pattern = "Single", x = current.levels) ],
#                                              current.levels[ grepl(pattern = "Symmetric", x = current.levels) ],
#                                              current.levels[ grepl(pattern = "Asymmetric", x = current.levels) ] ) 
#                                  levels( factor(dp$trunc.type.pretty) )
#                                  
#                                  # make plotting dataframe by calculating log-prior for a grid of values (mu, sigma)
#                                  dp = dp %>%
#                                    rowwise() %>%
#                                    mutate( lprior = lprior_Jeffreys( .pars = c(.mu, .sigma),
#                                                                      par2is = "sd",
#                                                                      .n = .n, 
#                                                                      
#                                                                      .a = a,
#                                                                      .b = b)  )
#                                  
#                                  # check again for NA values occurring when EFisher is NA and remove them
#                                  table(is.na(dp$lprior))
#                                  dp = dp %>% filter( !is.na(lprior) )
#                                  
#                                  # set up colors for contours
#                                  get_colors = colorRampPalette( c("lemonchiffon1", "chocolate4") )
#                                  myColors = get_colors(n=11)  # chose 11 based on errors from ggplot if it was fewer
#                                  
#                                  # make plot
#                                  p = ggplot( data = dp, 
#                                              aes(x = .mu,
#                                                  y = .sigma,
#                                                  z = lprior) ) +
#                                    
#                                    geom_contour_filled() +
#                                    
#                                    # close, but not enough colors
#                                    scale_fill_manual(values = myColors) +
#                                    
#                                    geom_contour(color = "white") +
#                                    
#                                    xlab( bquote(mu) ) +
#                                    ylab( bquote(sigma) ) +
#                                    
#                                    geom_vline( xintercept = 0, lty = 2 ) +
#                                    
#                                    facet_wrap( trunc.type.pretty ~.,
#                                                scales = "fixed" ) +
#                                    
#                                    scale_y_continuous(breaks = seq( min(dp$.sigma), max(dp$.sigma), 0.5),
#                                                       limits = c( min(dp$.sigma), max(dp$.sigma) ) ) +
#                                    
#                                    theme_bw(base_size = 16) +
#                                    theme(text = element_text(face = "bold"),
#                                          axis.title = element_text(size=20),
#                                          legend.position = "none")
#                                  
#                                  
#                                  my_ggsave( name = "jeffreys_prior_contours.pdf",
#                                             .width = 8,
#                                             .height = 8,
#                                             .results.dir = results.dir,
#                                             .overleaf.dir.general = overleaf.dir.figs )
# }  # end if(redo.plots == TRUE)
# 




# ~ GET NLPOSTERIOR FOR JEFFREYS ----------------------------------------------------




