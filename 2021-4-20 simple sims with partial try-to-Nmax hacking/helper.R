

# DATA SIMULATION ---------------------------------------------------------------

# simulate a single study from potentially heterogeneous meta-analysis distribution 
#  and from its own heterogeneous distribution

#~~ bm: want to look at this for rho>0
#  to do that, need to draw each t-stat from Peng's conditional distribution
# e.g.: when correlated, do t-stats bunch up more?

# similar to phack_study in earlier helper script (2020-6 folder)

draw_hacked_study = function(Nmax = Inf,  # max draws to try
                       Mu,  # overall mean for meta-analysis
                       T2,  # across-study heterogeneity
                       m,  # sample size for this study
                       t2w,  # within-study heterogeneity
                       se,  # TRUE SE for this study
                       rho = 0,  # correlation of Tn with T_{n-1} 
                       hack = "affirm") {  # hack until significant or until affirmative?
  
  
  # simulate the study
  # its marginal heterogeneity reflects meta-analysis heterogeneity and within-study heterogeneity
  mui = Mu + rnorm(mean = 0,
                   sd = sqrt(T2 + t2w),
                   n = 1)
  
  success = FALSE
  N = 0  # counts draws actually made
  
  sd.y = se * sqrt(m)
  
  while ( success == FALSE ) {
    
    if ( N == Nmax ) break
    N = N + 1
    
    # draw subject-level data from this study's population effect
    y = rnorm( mean = mui,
               sd = sd.y,
               n = m)
    
    # run a one-sample t-test
    if ( rho == 0 ) {
      pval = t.test(y,
                    alternative = "two.sided")$p.value
      tstat = t.test(y,
                     alternative = "two.sided")$statistic
    } else {
      stop("rho>0 case not implemented yet")
    }
    
    
    if (hack == "signif") success = (pval < 0.05)
    if (hack == "affirm") success = (pval < 0.05 & mean(y) > 0)
  }
  
  return( data.frame(pval = pval,
                     tstat = tstat,
                     tcrit = qt(0.975, df = m-1),
                     mui = mui,
                     muHati = mean(y), 
                     N = N
  ) )
}


# draw_hacked_study( Nmax = Inf,
#              Mu = 0.1,
#              T2 = 0.1,
#              m = 50,
#              t2w = .5,
#              se = 1)

# # sanity check: should follow truncated t distribution
# d = data.frame( matrix(nrow = 500, ncol = 1))
# 
# d = d %>% rowwise() %>%
#   mutate( draw_hacked_study( Nmax = Inf,
#                                 Mu = 0.1,
#                                 T2 = 0.1,
#                                 m = 50,
#                                 t2w = .5,
#                                 se = 1) )
# 
# # compare to truncated t
# qqtrunc(x = d$tstat,
#         spec = "t",
#         df = 50-1,
#         # since I simulated iid studies here, the truncation cutoff is always the same
#         a = d$tcrit[1] )
# 
# # vs. actually drawing directly from truncated t 
# x = rtrunc( n = 500,
#                spec = "t",
#                df = 50-1,
#                a = d$tcrit[1])
# 
# plot( density(x) ) +
#   lines( density(d$tstat), col = "red")
# # looks good :)




draw_unhacked_study = function(Nmax = Inf,  # max draws to try
                             Mu,  # overall mean for meta-analysis
                             T2,  # across-study heterogeneity
                             m,  # sample size for this study
                             t2w,  # within-study heterogeneity
                             se,  # TRUE SE for this study
                             rho = 0,  # correlation of Tn with T_{n-1} 
                             hack = "affirm") {  # hack until significant or until affirmative?
  
  
  # simulate the study
  # its marginal heterogeneity reflects meta-analysis heterogeneity and within-study heterogeneity
  mui = Mu + rnorm(mean = 0,
                   sd = sqrt(T2 + t2w),
                   n = 1)
  
  success = FALSE
  N = 0  # counts draws actually made
  
  sd.y = se * sqrt(m)
  
  while ( success == FALSE ) {
    
    if ( N == Nmax ) break
    N = N + 1
    
    # draw subject-level data from this study's population effect
    y = rnorm( mean = mui,
               sd = sd.y,
               n = m)
    
    # run a one-sample t-test
    if ( rho == 0 ) {
      pval = t.test(y,
                    alternative = "two.sided")$p.value
      tstat = t.test(y,
                     alternative = "two.sided")$statistic
    } else {
      stop("rho>0 case not implemented yet")
    }
    
    
    if (hack == "signif") success = (pval < 0.05)
    if (hack == "affirm") success = (pval < 0.05 & mean(y) > 0)
  }
  
  return( data.frame(pval = pval,
                     tstat = tstat,
                     tcrit = qt(0.975, df = m-1),
                     mui = mui,
                     muHati = mean(y), 
                     N = N
  ) )
}














