

library(testthat)
library(ggplot2)
library(dplyr)

code.dir = "~/Dropbox/Personal computer/Independent studies/*Inchoate/Sensitivity analysis for p-hacking/Code (git)"
results.dir = "~/Dropbox/Personal computer/Independent studies/*Inchoate/Sensitivity analysis for p-hacking/Results from R"
overleaf.dir = "~/Dropbox/Apps/Overleaf/P-hacking/figures"
  
setwd(code.dir)
source("helper.R")


################################ 2020-6-3: TEST DRIVE ON A META-ANALYSIS ################################ 

# get the Anderson data from SAPB-E
setwd("~/Dropbox/Personal computer/Independent studies/Sensitivity analysis for publication bias (SAPB)/Private data component/Anderson data")

da = read.csv("anderson_prepped.csv")

# proportion significant
table(da$pval<0.05)

# remember to convert MLE from z-scale back to "real" scale by multiplying by sei
# then we have Fisher's z again in this case



# corrected

# interesting to think about this case,
#  where the significant studies just get the truncated likelihood:
correct_meta_phack( yi = da$yi,
                    sei = da$sei,
                    # sensitivity parameters assumed constant across studies for now
                    n = 1,
                    t2w = 0 )


# using only the significant ones
correct_meta_phack( yi = da$yi[da$pval<0.05],
                      sei = da$sei[da$pval<0.05],
                      # sensitivity parameters assumed constant across studies for now
                      n = 1,
                      t2w = 0 )


################################ 2020-6-3: LOOK AT HOPEFULLY-CONSERVATIVE MLE FROM IPAD ################################ 


##### Individual Examples of Solving for N to Attenutate By Various Amounts #####

zn = 3
t2w = 1

# heuristics to understand choices above:
# p-value 
2 * ( 1 - pnorm( abs(zn) ) )

# 95% of true effects for this study are in (approximately):
zn - sqrt(t2w); zn + sqrt(t2w)



mle4( zn = zn,  # the observed one
      n = 100,
      t2w = 0 )

mle4( zn = zn,  # the observed one
      n = 40,
      t2w = t2w )


mle4( zn = zn,  # the observed one
      n = 50,
      t2w = 0 )



##### Plot For Various Values #####

# this version treats as our observed "data" the actual zn AND the fact that all n-1 previous draws were N.S.
# keeps n as a sensitivity parameter, but assumes rho = 0 because I think that is conservative

# vary the p-value and calculate the yi for that p-value
dp = expand.grid( #zn = c(1.96, 2.1, 2.5),
                  zn = c(1.96, 2.2, 2.8, 4),
                  n = seq(1,200,1),
                  t2w = c(0, .5, 1, 2) )

dp = dp %>% rowwise() %>%
  mutate( mle = mle4(zn,  # the observed one
                     n,
                     t2w = t2w,
                     include.signif.term = TRUE,
                     select.tails = 2) )

dp$zn.pretty = paste( "Observed zn = ", dp$zn, sep = "" )
dp$t2w.pretty = paste( "t2w = ", dp$t2w, sep = "" )

##### Make the Plot #####
colors = c("red", "orange", "black", "purple")

ggplot( data = dp,
        aes( x = n, 
             y = mle,
             color = t2w.pretty ) ) +
  
  # # mark when z1 becomes significant
  # geom_vline( aes(xintercept = 1.96),
  #             lty = 2,
  #             color = "red") +
  # 
  # # mark when z1 becomes > z2 itself
  # geom_vline( aes(xintercept = z2),
  #             lty = 2,
  #             color = "black") +
  
  # threshold for meaningful effect size?
  # if SD = 1
  geom_hline( aes(yintercept = .1),
              lty = 2,
              color = "red") +

  
geom_line(lwd = 1.1) +
  
  theme_bw(base_size = 18) +
  
  xlab("n (number of draws)") +
  ylab("MLE") +
  
  facet_grid(~ zn.pretty ) +
  
  scale_color_manual(values = colors,
                     name = "n (fixed number of draws)" ) +
  
  scale_y_continuous( limits = c( 0, 0.25 ) )


# ggsave( "2020-6-3_mle_given_previous_zstat_signif.pdf",
#         width = 15,
#         height = 10)

################################ 2020-5-26: LOOK AT LKL GIVEN AUTOCORRELATION (NORMAL DIST) ################################ 


# vary the p-value and calculate the yi for that p-value
dp = expand.grid( z2 = c(1.96, 2.1, 2.5, 4),
                  z1 = seq(-4, 4, .05),  # to use marginal MLE function, needs to be p<0.05
                  rho = c(0, 0.1, .5, .95) )

dp = dp %>% rowwise() %>%
  mutate( mle = mle3(z2 = z2,
                     z1 = z1, 
                     rho = rho,
                     # CAN TOGGLE THIS 
                     include.signif.term = TRUE,
                     select.tails = 2) )

dp$z2.pretty = paste( "Current Z = ", dp$z2, sep = "" )

##### Make the Plot #####
colors = c("red", "orange", "black", "purple")

ggplot( data = dp,
        aes( x = z1, 
             y = mle,
             color = as.factor(rho) ) ) +
  
  # mark when z1 becomes significant
  geom_vline( aes(xintercept = 1.96),
              lty = 2,
              color = "red") +
  
  # mark when z1 becomes > z2 itself
  geom_vline( aes(xintercept = z2),
              lty = 2,
              color = "black") +
  
  # marginal MLE is just z2 itself
  geom_hline( aes(yintercept = z2),
              lty = 2,
              color = "black") +
  
  # # cutoff for tstat1 to be less than 
  # geom_vline( aes( xintercept = tcrit1),
  #             lty = 1,
  #             color = "black") +
  
  # # value of the current t-stat
  # geom_vline( aes( xintercept = tstat2),
  #             lty = 1,
  #             color = "red") +
  
  geom_line(lwd = 1.1) +
  
  theme_bw(base_size = 18) +
  
  xlab("Previous Z-stat") +
  ylab("MLE") +
  
  facet_grid(~ z2.pretty ) +
  
  scale_color_manual(values = colors,
                     name = bquote(rho) ) +
  
  scale_y_continuous( limits = c(-10,20) )


ggsave( "2020-5-29_mle_given_previous_zstat_signif.pdf",
        width = 15,
        height = 10)

# it seems like as long as z2 > 2.2 approximately, assuming rho=0 is conservative
#  as long as the previous z-stat was nonsignificant

# what if we used the conditional distribution given only the knowledge that previous stat was 
#  n.s.? Would that help because we'd be integrating over the distribution of previous ones,
# so would include some smaller values?


##### Try to Understand #####

# counter-example to conservatism! 
# here, conditioning on a smaller, highly correlated previous z-stat
#  makes the MLE become negative 
# this doesn't seem to make sense?
z2 = 2
z1 = 1.9
include.signif.term = TRUE
select.tails = 2

mle3( z2 =  z2,
      z1 = z1,
      rho = 0,
      include.signif.term = include.signif.term,
      select.tails = select.tails)

mle3( z2 =  z2,
      z1 = z1,
      rho = 0.5,
      include.signif.term = include.signif.term,
      select.tails = select.tails)

# but this even seems to happen if previous z-stat is LARGER than current one?
mle3( z2 =  2,
      z1 = 2.5,
      rho = 0,
      include.signif.term = TRUE )

# **what is it about the significance term that causes this behavior?
# try going inside the MLE function and seeing how the significance term responds to different
#  choices
# * if rho = 0, then negative Z is bad for the lkl because the dnorm term is really small
#   (since we observed a positive Zn), AND bad because then P(signif) might be pretty large
#  because it could be significant in the other direction
# * but if rho is big and Z_{n-1} is also pretty big, then the conditional mean might be
#  close to 0 rather than negative (i.e., equal to Z itself), and that is GOOD for lkl because
#  then P(signif) is small AND dnorm isn't super small either


ll3( th = 0,
     z2 = 2,  # the current one
     z1 = 2.5,  # the previous one
     t2w = 0, 
     rho = 0,
     include.signif.term = FALSE )



################################ 2020-5-26: LOOK AT LKL GIVEN AUTOCORRELATION (T DIST) ################################ 


##### Lkl and MLE of Current T-stat Given Previous #####

log_lkl2( tstat2 = 2.5,  
          se2 = 1, 
          tstat1 = -1, 
          se1 = 1,
          df = 100,
          t2w = 0, 
          cv = 0,
          th = 1 )
# -2.017 if tstat1 = -1
# -2.0399 if tstat1 = 1

# why the crazy peak at tstat2 = 2, tstat1 = -2?
mle2( tstat2 = 2,  
      se2 = 1, 
      tstat1 = -2, 
      se1 = 1,
      df = 50,
      t2w = 0, 
      cv = 0.8 )

log_lkl2( tstat2 = 2,  
          se2 = 1, 
          tstat1 = -2, 
          se1 = 1,
          df = 50,
          t2w = 0, 
          cv = 0.8,
          th = 1 )


##### Plot It #####

# * Preliminarily, seems like the MLE under autocorrelation is always larger than the 
#  MLE under independence (good for conservatism!) as long as the previous t-stat is less
#  than the current one, which makes hella sense and will hold in our application because
#  the previous one was N.S. and the current one is significant
# * Obviously the t-stats have to be POSITIVELY autocorrelated for this to make sense

# vary the p-value and calculate the yi for that p-value
dp = expand.grid( tstat1 = seq(-3, 5, .5),
                  tstat2 = c(2.5, 4),  # to use marginal MLE function, needs to be p<0.05
                  cv = c(0, .8),
                  df = 50,
                  se2 = 1)

dp$pval2 = 2 * (1 - pt( abs(dp$tstat2), df = dp$df ) )


# critical value for tstat1
dp$tcrit1 = qt( p = 0.975, df = dp$df )

dp$tstat2.pretty = paste( "tstat2 = ", dp$tstat2, sep = "" )


# wide version of data
dpw = dp %>% rowwise() %>%
  mutate( mle.cond = mle2( tstat2 = tstat2,  
                           se2 = 1, 
                           tstat1 = tstat1, 
                           se1 = 1,
                           df = df,
                           t2w = 0, 
                           cv = cv ),
          mle.marg = mle(m = df + 1, # sample size
                         pval = pval2,  # or could use the yi, alternatively
                         SDi = se2 * sqrt(df+1),  # SD, *not* SE
                         regular.MLE = FALSE),
          mle.naive = mle(m = df + 1, # sample size
                         pval = pval2,  # or could use the yi, alternatively
                         SDi = se2 * sqrt(df+1),  # SD, *not* SE
                         regular.MLE = TRUE) )


# reshape into long version for plotting joy
library(tidyr)
detach( "package:plyr", unload = TRUE )
group.vars = names(dp)[ !names(dp) %in% c("mle.cond", "mle.marg", "mle.naive") ]
dp = dp %>% gather( key = "mle.type", value = "mle.value", -c(group.vars) )

# since se2 = 1, the naive mle should be damn close to the tstat itself
dpw %>% group_by(tstat2, mle.type) %>%
  filter(cv == 0) %>%
  summarise( pval2 = pval2[1],
             mean(mle.value) )

# * the conditional MLE under positive autocorrelation seems to always be larger than the 
# marginal one, whenever previous tstat was smaller than the current one
temp = dpw %>% filter(tstat1 < tstat2 & cv > 0 )
mean( temp$mle.cond < temp$mle.marg )
mean( temp$mle.cond < temp$mle.marg )

ind = which( temp$mle.cond < temp$mle.marg )

##### Make the Plot #####
colors = c("red", "orange", "black", "purple")

ggplot( data = dp,
        aes( x = tstat1, 
             y = mle.value,
             color = as.factor(cv),
             lty = as.factor(mle.type) ) ) +
  
  # the second tstat itself
  # expect this to be the MLE if SE2 = 1 and cv=0
  # since then tstat2 = the point estimate
  geom_hline( aes(yintercept = tstat2),
             lty = 1,
             color = "black") +

  # cutoff for tstat1 to be significant
  geom_vline( aes( xintercept = tcrit1),
             lty = 1,
             color = "black") +
  
  # value of the current t-stat
  geom_vline( aes( xintercept = tstat2),
              lty = 1,
              color = "red") +
  # 
  # # the alternative cutoff
  # geom_vline(xintercept = .005, 
  #            lty = 2,
  #            color = "red") +
  
  geom_line(lwd = 1.1) +
  
  theme_bw(base_size = 18) +
  
  xlab("Previous t-stat") +
  ylab("MLE") +
  
    facet_grid(~ tstat2.pretty ) +

    scale_color_manual(values = colors,
                       name = "Covariance of successive t-stats" ) +
  
  scale_y_continuous( limits = c(0,20) )


ggsave( "2020-5-26_mle_given_previous_tstat.pdf",
        width = 15,
        height = 10)


##### Debug the Crazy Peak #####

# why the crazy peak at tstat2 = 2, tstat1 = -2?
mle2( tstat2 = 2,  
      se2 = 1, 
      tstat1 = -2, 
      se1 = 1,
      df = 50,
      t2w = 0, 
      cv = 0.8 )

log_lkl2( tstat2 = 2,  
           se2 = 1, 
           tstat1 = -2, 
           se1 = 1,
           df = 50,
           t2w = 0,
           th = 0,
           cv = .8 )$ll

# **why is optim not doing this correctly?? based on plugging numbers into the above, 
#  clearly the function is not minimized at theta = 268
optim( par = 1,  # starting point for optimization
       f = function(x) -log_lkl2(tstat2 = 2,  
                                 se2 = 1,  
                                 tstat1 = -2,
                                 se1 = 1,
                                 df = 50,
                                 t2w = 0, 
                                 cv = .8,
                                 th = x)$ll,
       method = "BFGS" )$par

# vary the p-value and calculate the yi for that p-value
dp = expand.grid( tstat1 = -2,
                  tstat2 = 2,
                  cv = .8,
                  df = 50,
                  th = seq(10,300,10))


dp = dp %>% rowwise() %>%
  mutate( ll = log_lkl2( tstat2 = tstat2,  
                             se2 = 1, 
                             tstat1 = tstat1, 
                             se1 = 1,
                             df = df,
                             t2w = 0,
                             th = th,
                             cv = .8 )$ll )


ggplot( data = dp,
        aes( x = th, 
             y = ll ) ) +

  geom_line(lwd = 1.1)



##### What Happens with Heterogeneity? #####

# seems like MLE usually gets LARGER with t2w > 0, 
#  so potentiallly another conservatism opportunity?
mle2( tstat2 = 4,  
      se2 = 1, 
      tstat1 = 1, 
      se1 = 1,
      df = 50,
      t2w = 1, 
      cv = 0 )

################################ 2020-5-22: P-VALUE VS. MLE WRONGNESS RATIO ################################ 

# ** need to look into why the regular MLE is sometimes slightly different from yi

# **seems like p<0.001 usually means wrongness ratio < 1.1
#  ** and p<0.005 usually means wrongness ratio < 1.5
#  for various choices of m, and almost invariant to SD

# **for any given study, if we have the sample size, p-value or estimate, and SD, 
#  we can calculate an adjusted point estimate assuming infinite p-hacking
#  and any amount of autocorrelation that is <1 AND no intra-study heterogeneity

# if we allow for intra-study heterogeneity, then the 1/power will be smaller
#  and the truncated density won't be t anymore, but rather integrated over the normal 
#  distribution of thi?

# ** try this next and see if the above conclusions about p-values and wrongness ratios change


##### Understand It #####

# compare regular to adjusted MLEs

# regular MLE should match yi
# ~~~ not quite
m = 15
SDi = 1
yi = .6
# pval, just for context
SEi = SDi / sqrt(m)
2 * ( 1 - pt( abs(yi/SEi), df = m-1) )

mle(m=m,
    SDi = SDi,
    yi = yi,
    regular.MLE = T)

mle(m=m,
    SDi = SDi,
    yi = yi,
    regular.MLE = F)


# and working with p-value instead of yi
m = 10
SDi = 1
pval = 0.001


( mle.reg = mle(m=m,
                SDi = SDi,
                pval = pval,
                regular.MLE = T) )

( mle.adj = mle(m=m,
                SDi = SDi,
                pval = pval,
                regular.MLE = F) )

# * ratio by which the regular MLE is off
# depends heavily on the p-value
# but less so on m, interestingly
mle.reg/mle.adj



##### Set up Plot Dataframe #####
# vary the p-value and calculate the yi for that p-value
dp = expand.grid( m = c(10, 50, 10000),
                  SDi = c(.1, 3, 50),  # ~~~ not sure the role of this...
                  pval = seq(0.0001, 0.049, .01) )

# get MLE ratio
dp = dp %>% rowwise() %>%
  mutate( mle.ratio = mle(m=m,
                          SDi = SDi,
                          pval = pval,
                          regular.MLE = T) / mle(m=m,
                                                 SDi = SDi,
                                                 pval = pval,
                                                 regular.MLE = F) )


dp$SDi.pretty = paste( "SD = ", round(dp$SDi, 4), sep = "" )
# dp$pval.pretty = factor( dp$yi.pretty,
#                          levels = c("p = 1e-04", "p = 0.02", "p = 0.049") )

dp$m.pretty = paste( "m = ", dp$m, sep = "" )

# change ordering
dp$m.pretty = factor(dp$m.pretty, levels = c("m = 10",
                                             "m = 50",
                                             "m = 10000"))



##### Make Plot #####  

colors = c("red", "orange", "black", "purple")

ggplot( data = dp,
        aes( x = pval, 
             y = mle.ratio,
             color = m.pretty ) ) +
  
  # wrongness ratio of 1.1 might be tolerable?
  geom_hline(yintercept = 1.1, 
             lty = 2,
             color = "black") +
  
  # p-val of 0.001 might get the above wrongness ratio?
  geom_vline(xintercept = .001, 
             lty = 2,
             color = "black") +
  
  # the alternative cutoff
  geom_vline(xintercept = .005, 
             lty = 2,
             color = "red") +
  
  geom_line(lwd = 1.1) +
  
  theme_bw(base_size = 18) +
  
  xlab("p-value") +
  ylab("Uncorrected MLE / corrected MLE") +
  
  facet_grid(~ SDi.pretty ) +
  
  theme_classic() +
  
  scale_color_manual(values = colors, 
                     name = "m (sample size)" ) +
  
  # zoomed in on smaller p-values
  #scale_x_continuous( limits = c(0,0.011)) +

scale_y_continuous(breaks = seq(1, 4, .25),
                   limits = c(1, 4)) +
  
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14, face = "bold") )


setwd(results.dir)
# ggsave( "pval_vs_mle_ratio_zoomed_in.pdf",
#         width = 15,
#         height = 10)

ggsave( "pval_vs_mle_ratio.pdf",
        width = 13,
        height = 5)

# version for McCormick grant

colors = c("red", "orange", "black", "purple")

ggplot( data = dp,
        aes( x = pval, 
             y = mle.ratio,
             color = m.pretty ) ) +
  
  # wrongness ratio of 1.1 might be tolerable?
  geom_hline(yintercept = 1.1, 
             lty = 2,
             color = "black") +
  
  # p-val of 0.001 might get the above wrongness ratio?
  geom_vline(xintercept = .001, 
             lty = 2,
             color = "black") +
  
  # the alternative cutoff
  geom_vline(xintercept = .005, 
             lty = 2,
             color = "red") +
  
  geom_line(lwd = 1.1) +
  
  theme_bw(base_size = 18) +
  
  xlab("p-value") +
  ylab("Maximum bias factor of uncorrected estimate") +
  
  facet_grid(~ SDi.pretty ) +
  
  theme_classic() +
  
  scale_color_manual(values = colors, 
                     name = "m (sample size)" ) +
  
  # zoomed in on smaller p-values
  #scale_x_continuous( limits = c(0,0.011)) +
  
  scale_y_continuous(breaks = seq(1, 4, .25),
                     limits = c(1, 4)) +
  
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14, face = "bold") )


setwd(results.dir)
# ggsave( "pval_vs_mle_ratio_zoomed_in.pdf",
#         width = 15,
#         height = 10)

ggsave( "pval_vs_mle_ratio.pdf",
        width = 13,
        height = 5)


################################ 2020-5-22: POINT ESTIMATE VS. LOG-LIKELIHOOD OF THI ################################ 

# ~~~ note: if both m and yi are too large, the log-lkl becomes Inf and 
#  can't be optimized


# vary the p-value and calculate the yi for that p-value
dp = expand.grid( m = c(15, 50, 100),
                  SDi = 1,  # ~~~ not sure the role of this...
                  thi = seq(0, 1.5, 0.01),
                  yi = c(.1, .5, .8) )


# calculate the observed t-stat and p-value
dp$SEi = dp$SDi / sqrt(dp$m)
dp$t = dp$yi/dp$SEi
dp$pval = 2 * ( 1 - pt( abs(dp$t), df = dp$m-1 ) )


# remove nonsignificant yis
dp = dp %>% filter( pval < 0.05 )

# get bias-corrected MLE for each 

# bm
dp = dp %>% rowwise() %>%
  mutate( ll = log_lkl1(m = m,
                        SDi = SDi,
                        thi = thi,
                        yi = yi,
                        alpha = 0.05),
          mle.adj = mle(m = m,
                        SDi = SDi,
                        yi = yi) )


# mle(m=1000,
#     SDi = 1,
#     yi = .5,
#     regular.MLE = F)
# 
# log_lkl1(m = 100,
#          SDi = 1,
#          thi = .8,
#          yi = .8,
#          alpha = 0.05)
# 


# dp$pval.pretty = paste( "p = ", round(dp$pval, 4), sep = "" )
# dp$pval.pretty = factor( dp$pval.pretty,
#                          levels = c("p = 1e-04", "p = 0.02", "p = 0.049") )

dp$yi.pretty = paste( "Point estimate = ", round(dp$yi, 4), sep = "" )
# dp$pval.pretty = factor( dp$yi.pretty,
#                          levels = c("p = 1e-04", "p = 0.02", "p = 0.049") )

dp$m.pretty = paste( "m = ", dp$m, sep = " " )



# sanity check
expect_equal( dp$yi/(dp$SDi/sqrt(dp$m)), dp$t, tol = 0.001 )


##### Make Plot #####  

colors = c("red", "orange", "black", "purple")

ggplot( data = dp,
        aes( x = thi, 
             y = ll,
             color = as.factor(m) ) ) +
  
  # regular MLE
  geom_vline( aes(xintercept = yi ),
              color = "grey",
              lty = 1,
              lwd = 1.1) +  # E[N] if tau=0 (just from false positives)
  
  # adjusted MLE
  geom_vline( aes(xintercept = mle.adj,
                  color = as.factor(m) ),
              lty = 2,
              lwd = 1.1) +  # E[N] if tau=0 (just from false positives)
  
  geom_line(lwd = 1.1) +
  
  theme_bw(base_size = 18) +
  
  xlab(bquote(theta[i])) +
  ylab("Log-likelihood") +
  
  facet_grid(~ yi.pretty ) +
  
  scale_color_manual(values = colors, 
                     name = "m (sample size)" ) +

scale_x_continuous(breaks = seq(0.2, 1, .1),
                   limits = c(0.2,1))
# 
# scale_y_continuous(breaks = seq(0, .7, .1),
#                    limits = c(0,0.7)) +


# sanity check: When k = 0 (no heterogeneity) and N = 1, 
#  we are just looking at the standard-issue lkl

setwd(results.dir)
ggsave( "yi_vs_lkl.pdf",
        width = 15,
        height = 10)


################################ 2020-5-19: P-VALUE VS. LOG-LIKELIHOOD OF THI ################################ 

# vary the p-value and calculate the yi for that p-value
dp = expand.grid( m = c(15, 50, 1000),
                  SDi = 1,  # ~~~ not sure the role of this...
                  thi = seq(0, 1.5, 0.01),
                  pval = c(1e-04, 0.02, 0.049))


dp = dp %>% rowwise() %>%
  mutate( ll = log_lkl1(m = m,
                        SDi = SDi,
                        thi = thi,
                        pval = pval,
                        alpha = 0.05) )


dp$pval.pretty = paste( "p = ", round(dp$pval, 4), sep = "" )
dp$pval.pretty = factor( dp$pval.pretty,
                         levels = c("p = 1e-04", "p = 0.02", "p = 0.049") )

dp$m.pretty = paste( "m = ", dp$m, sep = " " )

# calculate the observed yi for reference as the standard MLE
dp$t = qt(1 - dp$pval/2, df = dp$m - 1)
dp$yi = dp$t * dp$SDi/sqrt(dp$m)

# sanity check
expect_equal( dp$yi/(dp$SDi/sqrt(dp$m)), dp$t, tol = 0.001 )



##### Make Plot #####  

colors = c("red", "orange", "black")

ggplot( data = dp,
        aes( x = thi, 
            y = ll,
            color = as.factor(m) ) ) +
  
  geom_vline( aes(xintercept = yi,
                  color = as.factor(m) ),
              lty = 2 ) +  # E[N] if tau=0 (just from false positives)
  
  geom_line(lwd = 1.1) +
  
  theme_bw(base_size = 18) +
  
  xlab(bquote(theta[i])) +
  ylab("Log-likelihood") +
  
  facet_grid(~ pval.pretty ) +
  
  scale_color_manual(values = colors, 
                     name = "m (sample size)" ) 
# 
# scale_x_continuous(breaks = seq(1, 10, 2),
#                    limits = c(0,10)) +
# 
# scale_y_continuous(breaks = seq(0, .7, .1),
#                    limits = c(0,0.7)) +


# sanity check: When k = 0 (no heterogeneity) and N = 1, 
#  we are just looking at the standard-issue lkl

setwd(results.dir)
ggsave( "thi_vs_lkl.pdf",
        width = 15,
        height = 10)





################################ 2020-5-19: P-VALUE VS. E[N] ################################ 

# see helper (function EN) for explanation of k and prop
pval.vec = seq(0.001, 0.05, 0.001)
k.vec = c(0, 0.5, 2, 5)


dp = expand.grid( pval = pval.vec,
                  k = k.vec )

dp = dp %>% rowwise %>%
  mutate( Psig = Psig(thi = 0,  # true mean in this study
                      
                      k = k,
                      pval = pval,
                      
                      tails = 2 ),
          EN = EN( independent = TRUE, 
                   Psig = Psig) )


dp$k.pretty = as.character(dp$k)



##### Make Plot #####  

colors = c("lightblue", "blue", "black", "green")

ggplot( data = dp,
        aes(x = pval, 
            y = EN,
            color = k.pretty ) ) +
  
  geom_hline(yintercept = 40,
             lty = 2,
             color = "red") +  # E[N] if tau=0 (just from false positives)
  
  
  geom_line(lwd = 1.2) +
  theme_classic() +
  xlab("p-value of observed significant study") +
  ylab( bquote("Expected number of draws if" ~ theta[i] ~ "=0") ) +
  
  scale_color_manual(values = colors,
                     name = "k" ) +
  
  # scale_x_continuous(breaks = seq(0, max(mu.vec), .1)) +
  scale_y_continuous(breaks = seq(0, 20, 1),
                     limits = c(0,20)) +
  
  ggtitle("95% of true effects within k-fold of estimate, if estimate were truth")


my_ggsave( "plot_pval_vs_EN.pdf",
           width = 6, 
           height = 6)


# explore this
Psig(thi = 0,  # true mean in this study
     
     k = 2,
     pval = .0005,
     
     tails = 2 )


################################ 2020-5-19: N DRAWS VS. PAFFIRM ################################ 

dp = expand.grid( pval = c(0.0001, 0.01, 0.03, 0.05),
                  k = k.vec,
                  prop = 0.95,
                  N = seq(1, 10, 1) )



dp = dp %>% rowwise() %>%
  mutate( Psig = Psig(thi = 0,  # true mean in this study
                      k = k,
                      pval = pval,
                      tails = 2 ),
          
          PsigN = Psig_N(
            Psig = Psig,
            independent = TRUE,
            N = N) )

dp$pval.pretty = paste( "p_iN = ", dp$pval, sep = "" )


##### Make Plot #####  

colors = c("lightblue", "blue", "black", "green")

ggplot( data = dp,
        aes(x = N, 
            y = PsigN,
            color = as.factor(k) ) ) +
  
  geom_hline(yintercept = 0.025,
             lty = 2,
             color = "red") +  # E[N] if tau=0 (just from false positives)
  
  geom_line(lwd = 1.2) +
  theme_classic() +
  xlab("Number of draws") +
  ylab( bquote("P(significant) in N draws if" ~ theta[i] ~ "=0") ) +
  
  scale_color_manual(values = colors,
                     name = "k" ) +
  
  scale_x_continuous(breaks = seq(1, 10, 2),
                     limits = c(0,10)) +
  
  scale_y_continuous(breaks = seq(0, 1, .1),
                     limits = c(0,1)) +
  
  facet_grid(~ pval.pretty )



my_ggsave( "N_vs_PsigN.pdf",
           width = 6,
           height = 6)




################################ 2020-5-15: P-PVALUE: P(SIGNIF IN N DRAWS AND AS EXTREME AS OBSERVED) ################################ 

# *assume no intra-study heterogeneity* this since I think this is again conservative

ppval = function(N, 
                 pval) {
  (1-0.95^N) * pval/0.05
}

# returns p-value
ppval_inv = function(N, 
                 ppval) {
  0.05^2 / (1-0.95^N)
}


##### PPValue as Fn of N, Pval ##### 

dp = expand.grid( N = c(1, 5, 10, 20),
                  pval = seq(0.001, 0.049, 0.001) )

dp = dp %>% rowwise() %>%
  mutate( ppval = ppval(N, pval) ) 

 

colors = c("lightblue", "blue", "black", "green")

ggplot( data = dp,
        aes(x = pval, 
            y = ppval,
            color = as.factor(N) ) ) +
  
  geom_hline(yintercept = 0.05,
             lty = 2,
             color = "red") +  

  geom_line(lwd = 1.2) +
  theme_classic() +
  xlab("P(p-value)") +
  ylab( bquote("p-value") ) +
  
  scale_color_manual(values = colors,
                     name = "N draws" ) 
  
  # scale_x_continuous(breaks = seq(1, 10, 2),
  #                    limits = c(0,10)) +
  # 
  # scale_y_continuous(breaks = seq(0, .7, .1),
  #                    limits = c(0,0.7)) +
  




setwd(results.dir)
ggsave( "ppvalue.pdf",
        width = 6,
        height = 6)


##### What P-value would we need to have P-Pvalue < 0.05 as fn of N? #####


dp = expand.grid( N = seq(1,40,1) )

dp = dp %>% rowwise() %>%
  mutate( pval = ppval_inv(N, 0.05) ) 



colors = c("lightblue", "blue", "black", "green")

ggplot( data = dp,
        aes(x = N, 
            y = pval ) ) +
  
  geom_hline(yintercept = 0.05,
             lty = 2,
             color = "red") +  
  
  geom_line(lwd = 1.2) +
  theme_classic() +
  xlab("N draws") +
  ylab( bquote("p-value needed to have P(p-value) = 0.05") )


setwd(results.dir)
ggsave( "ppvalue.pdf",
        width = 6,
        height = 6)

# play with it
ppval_inv(N = 50, ppval = 0.05)

# the p<0.005 criterion is robust to 10 draws
ppval( N = 10, pval = 0.005)

################################ 2020-5-15: IS P(SIGNIF) CONSERVATIVE WRT SE? ################################ 

# seems like it is...

dp = expand.grid( thi = c(-5, -1, 0, 0.2, 5),
                  t2w = c(0, .2, 2),
                  sei = seq(0.001, 3, .001) )

dp = dp %>% rowwise() %>%
  mutate( Psig = Paffirm(thi,
                         t2w,
                         sei,
                         tails = 2) ) # look at significance rather than affirmative status )

##### Make Plot #####  

colors = c("lightblue", "blue", "black", "green")

ggplot( data = dp,
        aes(x = sei, 
            y = Psig,
            color = as.factor(t2w) ) ) +
  
  geom_hline(yintercept = 0.05,
             lty = 2,
             color = "red") +  # false positive rate with no heterogeneity
  
  geom_line(lwd = 1.2) +
  theme_classic() +
  xlab("SE") +
  ylab( bquote("P(significant)") ) +
  
  scale_color_manual(values = colors,
                     name = "t2w" ) +
  
  # scale_x_continuous(breaks = seq(1, 10, 2),
  #                    limits = c(0,10)) +
  # 
  # scale_y_continuous(breaks = seq(0, .7, .1),
  #                    limits = c(0,0.7)) +
  
  facet_grid(~ thi )


setwd(results.dir)
ggsave( "small_SE_conservative.pdf",
        width = 10,
        height = 6)


# look at just the lower tail
lower_tail = function(q = qnorm(.975),
                      sei,
                      mu,
                      t2) {
  denom = sqrt( t2 + sei^2 )
  num1 = qnorm(.975) * sei - mu
  num2 = -qnorm(.975) * sei - mu
  
  Z1 = num1 / denom
  Z2 = num2 / denom

  return( pnorm(Z2) )
}



# hope this 
mu = 0.01
lower_tail(sei = 2,
           mu = mu,
           t2 = 1)

# hope this gets BIGGER
# because we've decreased the SE
lower_tail(sei = 1,
           mu = mu,
           t2 = 1)
# obviously gets smaller for big positive mu since very low probability it would be significant in the wrong direction

# overall Psignif for these two
Paffirm(mu,
        t2 = 1,
        sei = 2,
        tails = 2)

Paffirm(mu,
        t2 = 1,
        sei = 1,
        tails = 2)






################################ TRUE MEAN VS. LIKELIHOOD, TREATING N AS KNOWN ################################ 

yi = 2
se = 1

k.vec = c(0, 0.5, 2, 5)

# look at a study with p=0.045
# 2*(1-pnorm(2))

# get yi for a vector of p-values
pvals = c(0.049, 0.02, 0.0001)
# assumes se=1, so yi is also the z-score
yi = qnorm(1-pvals/2)

dp = expand.grid( yi = yi,
                  se = 1,
                  mui = seq(0, 5, 0.01),
                  N = c(1, 3, 5),
                  k = k.vec,  
                  prop = 0.95)

dp = dp %>% rowwise() %>%
  mutate( ll = log_lkl(yi = yi,
                       se = se,
                       mui = mui,
                       N = N,
                       k = k,
                       prop = .95) )

dp$pval = 2 * (1 - pnorm(dp$yi/dp$se))

dp$pval.pretty = paste( "p = ", round(dp$pval, 4), sep = "" )
dp$pval.pretty = factor( dp$pval.pretty,
                         levels = c("p = 1e-04", "p = 0.02", "p = 0.049") )

dp$k.pretty = paste( "k = ", dp$k, sep = " " )


##### Make Plot #####  

colors = c("red", "orange", "black")

ggplot( data = dp,
        aes(x = mui, 
            y = ll,
            color = as.factor(N) ) ) +
  
  geom_vline( aes(xintercept = yi),
             lty = 2,
             color = "red") +  # E[N] if tau=0 (just from false positives)

  geom_line(lwd = 1.1) +
  
  theme_bw(base_size = 18) +

  xlab(bquote(mu[i])) +
  ylab("Log-likelihood") +
  
  facet_grid(k.pretty ~ pval.pretty ) +
  
  scale_color_manual(values = colors, 
                     name = "N" ) 
  # 
  # scale_x_continuous(breaks = seq(1, 10, 2),
  #                    limits = c(0,10)) +
  # 
  # scale_y_continuous(breaks = seq(0, .7, .1),
  #                    limits = c(0,0.7)) +
  

# sanity check: When k = 0 (no heterogeneity) and N = 1, 
#  we are just looking at the standard-issue lkl

setwd(results.dir)
ggsave( "mui_vs_lkl.pdf",
        width = 15,
        height = 10)



# next up: get the MLE of mui, omg :)

# captions on ggplots?

# 1.) p-value vs. E(N) for a given mui
# 2.) N vs. P(affirmative) for a given mui
# 3.) mui vs. likelihood for a given N
# 4.) Bias-corrected meta-analysis
 
# look at the "Problems with Vevea/Woods model" to make sure we don't have this!




################################ CORRELATED P-HACKING ################################ 

# does E[N] stay the same if the draws are correlated?
#  i.e., is expectation of number of draws until success the same if draws are correlated?)

# function that generates a binary variable correlated with previous one

# generate autocorrelated data
# https://stats.stackexchange.com/questions/29239/creating-auto-correlated-random-values-in-r
x <- diffinv(rnorm(999))
plot(x)

x2 = x<30
mean(x2)  # mean success probability



