

# PRELIMINARIES ------------------------------

# borrow helper fns from previous sims
library(here)
setwd( here() )
source("helper_SAPH.R")

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


# NUMBERED SIMULATION EXPERIMENTS -------------------------
# as in table in project log


results.dir = here( "2021-5-10 clean sims about affirm moments/Numbered simulation datasets" )



# ~ EXPT 1: NONE HACKED, NMAX = 1 -------------------------

# rerun this
x = quick_sim( .p = data.frame( Mu = 1,
                                T2 = 0.25,
                                m = 500,
                                t2w = .25,
                                se = .5,
                                
                                Nmax = 1,
                                hack = "affirm",
                                rho = 0,
                                
                                k = 10000,
                                k.hacked = 0,
                                
                                
                                sim.name = "expt_1_res" ),
               .results.dir = results.dir )


x$res


# alternatively, re-load saved object, including dataset
setwd(results.dir)
load("expt_1_res")
d = x$d


# ~ EXPT 2: ALL HACKED, NMAX = 10 -------------------------

x = quick_sim( .p = data.frame( Mu = 1,
                                T2 = 0.25,
                                m = 500,
                                t2w = .25,
                                se = .5,
                                
                                Nmax = 10,
                                hack = "affirm",
                                rho = 0,
                                
                                k = 10000,
                                k.hacked = 10000,
                                
                                
                                sim.name = "expt_2_res" ),
               .results.dir = results.dir )


x$res


# alternatively, re-load saved object, including dataset
setwd(results.dir)
load("expt_2_res")
d = x$d


# ~ EXPT 3: ALL HACKED, BUT NMAX = 1 -------------------------

x = quick_sim( .p = data.frame( Mu = 1,
                                T2 = 0.25,
                                m = 500,
                                t2w = .25,
                                se = .5,
                                
                                Nmax = 1,
                                hack = "affirm",
                                rho = 0,
                                
                                k = 10000,
                                k.hacked = 10000,
                                
                                
                                sim.name = "expt_3_res" ),
               .results.dir = results.dir )


x$res


# alternatively, re-load saved object, including dataset
setwd(results.dir)
load("expt_3_res")
d = x$d


# ~ EXPT 4: ALL HACKED, NMAX = 200 -------------------------

x = quick_sim( .p = data.frame( Mu = 1,
                                T2 = 0.25,
                                m = 500,
                                t2w = .25,
                                se = .5,

                                Nmax = 200,
                                hack = "affirm",
                                rho = 0,

                                k = 10000,
                                k.hacked = 10000,

                                sim.name = "expt_4_res" ),
               
               .results.dir = results.dir )

x$res

# alternatively, re-load saved object, including dataset
setwd(results.dir)
load("expt_4_res")
d = x$d



# ~ EXPT 5: ALL HACKED, NMAX = 10, T2 = t2w = 0 -------------------------

# try with T2 = 0, still Nmax = 10
# THIS ONE IS PERFECT
x = quick_sim( .p = data.frame( Mu = 1,
                                T2 = 0,
                                m = 500,
                                t2w = 0,
                                se = .5,
                                
                                Nmax = 10,
                                hack = "affirm",
                                rho = 0,
                                
                                k = 10000,
                                k.hacked = 10000,
                                
                                sim.name = "expt_5_res" ),
               
               .results.dir = results.dir )

x$res


# alternatively, re-load saved object, including dataset
setwd(results.dir)
load("expt_5_res")
d = x$d


# ~ EXPT 6: ALL HACKED, NMAX = 10, T2 > 0, t2w = 0 -------------------------

# having T2 > 0 messes it up!!
x = quick_sim( .p = data.frame( Mu = 1,
                                T2 = 0.25,
                                m = 500,
                                t2w = 0,
                                se = .5,
                                
                                Nmax = 10,
                                hack = "affirm",
                                rho = 0,
                                
                                k = 10000,
                                k.hacked = 10000,
                                
                                
                                sim.name = "expt_6_res" ),
               
               .results.dir = results.dir )

# alternatively, re-load saved object, including dataset
setwd(results.dir)
load("expt_6_res")
d = x$d

Mu = unique(d$Mu)
T2 = unique(d$T2)
t2w = unique(d$t2w)
m = unique(d$m)
se = unique(d$se)



x$correctedSE^2

# first draw of each study set
( t1 = d %>% filter( !duplicated(study) ) %>%
    summarise( n(),
               var(mui),
               var(muin),
               var(yi),
               mean(tstat),
               var(tstat) ) )

# all draws from each study set
( t2 = d %>% 
    summarise( n(),
               var(mui),
               var(muin),
               var(yi),
               mean(tstat),
               var(tstat) ) )


# **IMPORTANT:
# marginal mean(tstat) and var(tstat) match Mu/se and correctedSE^2
# ONLY when I look at JUST the first draw per study set 
# if I use all draws per study set, mean(tstat) << Mu/se and var(tstat) slightly < correctedSE^2

# try calculating truncated expectation by cheating
# i.e., using the empirical marginal mean of ALL tstats
extrunc( spec = "norm",
         mean = t2$`mean(tstat)`,
         sd = sqrt(2),
         a = crit )

vartrunc( spec = "norm",
          mean = t2$`mean(tstat)`,
          sd = sqrt(2),
          a = crit )

x$res
# this doesn't really work


t3 = d %>% filter(!duplicated(study)) %>%
  group_by(N < 10) %>%
  summarise( mean(tstat),
             var(tstat) )
t3

extrunc( spec = "norm",
         mean = t3$`mean(tstat)`[t3$`N < 10` == 1],
         #sd = sqrt(2),
         sd = sqrt( t3$`var(tstat)`[t3$`N < 10` == 1] ),
         a = crit )

vartrunc( spec = "norm",
          mean = t3$`mean(tstat)`[t3$`N < 10` == 1],
          #sd = sqrt(2),
          sd = sqrt( t3$`var(tstat)`[t3$`N < 10` == 1] ),
          a = crit )

x$res

# this also doesn't work...not sure why...


# look at published ones only
# **this nicely illustrates the issue of selecting on mui effectively
( t4 = d %>% 
    filter(Di ==1) %>%
    summarise( n(),
               mean(mui),
               var(mui),
               mean(muin),
               var(muin),
               var(yi),
               mean(tstat),
               var(tstat) ) )

# **note that the distribution of muin after p-hacking selection is NOT normal
#  but rather skewed, which means any normal approximation below is not going to 
#  work perfectly
hist( d$muin[ d$Di == 1 ] )

correctedSE = deltamethod( g = ~ x1/sqrt(x2),
                           mean = c(Mu, se^2),
                           cov = matrix( c( t4$`var(muin)` + se^2, 0, 0, var(d$vi) ),
                                         nrow = 2 ) )

extrunc( spec = "norm",
         mean = t4$`mean(muin)` / se,
         #sd = sqrt(2),
         sd = correctedSE,  # need to scale for t-stat
         a = crit )

vartrunc( spec = "norm",
          mean = t3$`mean(tstat)`[t3$`N < 10` == 1],
          #sd = sqrt(2),
          sd = sqrt( t3$`var(tstat)`[t3$`N < 10` == 1] ),
          a = crit )


# calculate expectation of each one individually
dp = d %>% filter(Di == 1) %>%
  rowwise() %>%
  mutate( theoryExp_in = extrunc( spec = "norm",
                                  mean = muin / se,
                                  #sd = sqrt(2),
                                  sd = se,  # no T2 in here
                                  a = crit ),
          
          # **GOOD ONE - save it! 
          theoryExp_i = extrunc( spec = "norm",
                                 mean = mui / se,
                                 #sd = sqrt(2),
                                 sd = sqrt( (1/se^2) * (t2w + se^2) ),  
                                 a = crit ),
          
          varExp_i = vartrunc( spec = "norm",
                               mean = mui / se,
                               #sd = sqrt(2),
                               sd = sqrt( (1/se^2) * (t2w + se^2) ),  
                               a = crit ) )

mean(dp$theoryExp_in)
mean(dp$theoryExp_i)  # this one is HELLA close!!!! **save this one
mean(dp$tstat)

mean(dp$varExp_i)  # this is still off, probably because it's not really normal
var(dp$tstat)




# look at published nonaffirmatives
dn = d %>% filter( Di == 1 & affirm == FALSE )



# ~ EXPT 7: ALL HACKED, NMAX = 10, T2 = 0, t2w > 0 -------------------------

# this one is CORRECT!!
x = quick_sim( .p = data.frame( Mu = 1,
                                T2 = 0,
                                m = 500,
                                t2w = 0.25,
                                se = .5,

                                Nmax = 10,
                                hack = "affirm",
                                rho = 0,

                                k = 10000,
                                k.hacked = 10000,

                                sim.name = "expt_7_res" ),

               .results.dir = results.dir )

# alternatively, re-load saved object, including dataset
setwd(results.dir)
load("expt_7_res")
d = x$d

# first draw of each study set
( t1 = d %>% filter( !duplicated(study) ) %>%
    summarise( n(),
               var(mui),
               var(muin),
               var(yi),
               mean(tstat),
               var(tstat) ) )

# all draws from each study set
( t2 = d %>% 
    summarise( n(),
               var(mui),
               var(muin),
               var(yi),
               mean(tstat),
               var(tstat) ) )


# ~ EXPT 8: ALL HACKED, NMAX = 200, T2 = 0 -------------------------

x = quick_sim( .p = data.frame( Mu = 1,
                                T2 = 0,
                                m = 500,
                                t2w = .25,
                                se = .5,
                                
                                Nmax = 200,
                                hack = "affirm",
                                rho = 0,
                                
                                k = 10000,
                                k.hacked = 10000,
                                
                                sim.name = "expt_8_res" ),
               
               .results.dir = results.dir )

x$res

# alternatively, re-load saved object, including dataset
setwd(results.dir)
load("expt_8_res")
d = x$d


# first draw of each study set
( t1 = d %>% filter( !duplicated(study) ) %>%
    summarise( n(),
               var(mui),
               var(muin),
               var(yi),
               mean(tstat),
               var(tstat) ) )

# all draws from each study set
( t2 = d %>% 
    summarise( n(),
               var(mui),
               var(muin),
               var(yi),
               mean(tstat),
               var(tstat) ) )


