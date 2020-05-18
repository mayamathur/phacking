
# Setting: 
# Observe a single, affirmative study. Assume that study has an underlying heterogeneous normal
#  distribution of parameters that could be estimated by different statistical models, and each model 
#  has SE equal to sei.
# 1.) On any given draw from this distribution, what's the probability of getting an affirmative estimate?
# 2.) The inverse of that probability is the expected number of draws to get the first affirmative result. 




library(testthat)
library(ggplot2)
library(dplyr)

code.dir = "~/Dropbox/Personal computer/Independent studies/*Inchoate/Sensitivity analysis for p-hacking"
results.dir = "~/Dropbox/Personal computer/Independent studies/*Inchoate/Sensitivity analysis for p-hacking/Results from R"
  
setwd(code.dir)
source("helper.R")


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

dp = expand.grid( mu = c(-5, -1, 0, 0.2, 5),
                  t2 = c(0, .2, 2),
                  sei = seq(0.001, 3, .001) )

dp = dp %>% rowwise() %>%
  mutate( Psig = Paffirm(mu,
                         t2,
                         sei,
                         tails = 2) ) # look at significance rather than affirmative status )

##### Make Plot #####  

colors = c("lightblue", "blue", "black", "green")

ggplot( data = dp,
        aes(x = sei, 
            y = Psig,
            color = as.factor(t2) ) ) +
  
  geom_hline(yintercept = 0.05,
             lty = 2,
             color = "red") +  # false positive rate with no heterogeneity
  
  geom_line(lwd = 1.2) +
  theme_classic() +
  xlab("SE") +
  ylab( bquote("P(significant)") ) +
  
  scale_color_manual(values = colors,
                     name = "t2" ) +
  
  # scale_x_continuous(breaks = seq(1, 10, 2),
  #                    limits = c(0,10)) +
  # 
  # scale_y_continuous(breaks = seq(0, .7, .1),
  #                    limits = c(0,0.7)) +
  
  facet_grid(~ mu )


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



################################ EARLIER WORK (SAVE): P-VALUE VS. E[N] ################################ 

# see helper (function EN) for explanation of k and prop
pval.vec = seq(0.001, 0.05, 0.001)

# one dataframe for each k, for plotting joy
dfs = list()

for ( i in 1:length(k.vec) ) {
  
  # EN as function of study's p-value
  EN.vec = vapply( X = pval.vec, 
                   FUN = function(.p) EN(pval = .p,
                                         k = k.vec[i],
                                         prop = 0.95),
                   FUN.VALUE = -99 )
  
  dfs[[i]] = data.frame(pval = pval.vec, 
                   k = k.vec[i],
                   EN = EN.vec )
}

dp = do.call(rbind, dfs)


##### Make Plot #####  

colors = c("lightblue", "blue", "black", "green")

ggplot( data = dp,
        aes(x = pval, 
            y = EN,
            color = as.factor(k) ) ) +
  
  geom_hline(yintercept = 40,
             lty = 2,
             color = "red") +  # E[N] if tau=0 (just from false positives)

  
  geom_line(lwd = 1.2) +
  theme_classic() +
  xlab("p-value of observed affirmative study") +
  ylab( bquote("Expected number of draws if" ~ mu[i] ~ "=0") ) +

  scale_color_manual(values = colors,
                     name = "k" ) +
  
  # scale_x_continuous(breaks = seq(0, max(mu.vec), .1)) +
  scale_y_continuous(breaks = seq(0, 40, 5),
                     limits = c(0,40)) +
  
  ggtitle("95% of true effects within k-fold of estimate, if estimate were truth")


setwd(results.dir)
ggsave( "plot_pval_vs_EN.pdf",
        width = 6, 
        height = 6)


################################ N DRAWS VS. PAFFIRM ################################ 

dp = expand.grid( pval = c(0.0001, 0.01, 0.03, 0.05),
                      k = k.vec,
                      prop = 0.95,
                      N = seq(1, 10, 1) )

dp = dp %>% rowwise() %>%
  mutate( PaffirmN = Paffirm_N(pval = pval,
                               k = k,
                               prop = 0.95,
                               N = N) )


##### Make Plot #####  

colors = c("lightblue", "blue", "black", "green")

ggplot( data = dp,
        aes(x = N, 
            y = PaffirmN,
            color = as.factor(k) ) ) +
  
  geom_hline(yintercept = 0.025,
             lty = 2,
             color = "red") +  # E[N] if tau=0 (just from false positives)
  
  geom_line(lwd = 1.2) +
  theme_classic() +
  xlab("Number of draws") +
  ylab( bquote("Expected proportion affirmative if" ~ mu[i] ~ "=0") ) +
  
  scale_color_manual(values = colors,
                     name = "k" ) +
  
  scale_x_continuous(breaks = seq(1, 10, 2),
                     limits = c(0,10)) +

  scale_y_continuous(breaks = seq(0, .7, .1),
                     limits = c(0,0.7)) +
  
  facet_grid(~ pval )


setwd(results.dir)
ggsave( "N_vs_PaffirmN.pdf",
        width = 6,
        height = 6)




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



