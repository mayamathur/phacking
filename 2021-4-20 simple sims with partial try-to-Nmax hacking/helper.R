

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
                     success = success,
                     N = N,
                     Di = 1,  # publication indicator
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



# draws exactly Nmax results and treats the last one as the reported one
#  but returns them all
# since there's no hacking, argument "hack" is just to create indicator for whether draw was successful or not
# here the final result could be affirmative or nonaffirmative
# returns the entire study set, but only last draw is observed

draw_one_study_set = function(Nmax,  # max draws to try
                                   Mu,  # overall mean for meta-analysis
                                   T2,  # across-study heterogeneity
                                   m,  # sample size for this study
                                   t2w,  # within-study heterogeneity
                                   se,  # TRUE SE for this study
                                   rho = 0,
                                   #force.nonaffirm = FALSE, # should we force 
                                   hack ) {  # should this study set be hacked? ("no", "affirm", "signif")
  
  
  # # test only
  # Nmax = 20
  # Mu = 0.1
  # T2 = 0.1
  # m = 50
  # t2w = .5
  # se = 1
  # hack = "affirm"
  # rho=0
  
  # mean for this study set
  # doesn't have t2w because that applies to results within this study set
  mui = Mu + rnorm(mean = 0,
                   sd = sqrt(T2),
                   n = 1)
  
  sd.y = se * sqrt(m)
  
  # collect all args from outer fn, including default ones
  .args = mget(names(formals()), sys.frame(sys.nframe()))
  .args$mui = mui
  .args$sd.y = sd.y
  
 
  stop = FALSE  # indicator for whether to stop drawing results
  N = 0  # counts draws actually made
  
  
  while( stop == FALSE ) {

    newRow = do.call( draw_withinstudy_result, .args )
    
    N = N + 1
    
    # add new draw to dataset
    if ( N == 1 ) d = newRow
    if ( N > 1 ) d = rbind( d, newRow )
    
    # check if it's time to stop drawing results
    if (hack == "signif") stop = (newRow$pval < 0.05)
    if (hack == "affirm") stop = (newRow$pval < 0.05 & newRow$muHati > 0)
    # if this study set is unhacked, then success is just whether we've reached the 
    if ( hack == "no") stop = (N == Nmax)
  }
  
  
  # convenience indicators for significance and affirmative status
  d$signif = d$pval < 0.05
  d$affirm = (d$pval < 0.05 & d$muHati > 0)

  # publication indicator
  if ( hack == "signif" ) d$Di = (d$signif == TRUE)
  if (hack == "affirm") d$Di = (d$affirm == TRUE)
  # if no hacking, assume only last draw is published
  if ( hack == "no" ) d$Di = 0; d$Di[ length(d$Di) ] = 1
  
  return(d)
  
  #bm
  
}


# d = draw_one_study_set(Nmax = 20,
#                        Mu = 0.1,
#                        T2 = 0.1,
#                        m = 50,
#                        t2w = .5,
#                        se = 1,
#                        hack = "signif")
# nrow(d)
# d



# draws one unbiased result within the study
draw_withinstudy_result = function(mui,  # mean for this study set
                                   t2w,
                                   sd.y,
                                   m,
                                   #hack,
                                   rho = 0,
                                   ...) {
  
  
  # true mean for draw n (based on within-study heterogeneity)
  muin = rnorm(mean = mui,
               sd = sqrt(t2w),
               n = 1)
  
  # draw subject-level data from this study's population effect
  y = rnorm( mean = muin,
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
  
  # if (hack == "signif") success = (pval < 0.05)
  # if (hack == "affirm") success = (pval < 0.05 & mean(y) > 0)
  
  return( data.frame(pval = pval,
                     tstat = tstat,
                     tcrit = qt(0.975, df = m-1),
                     mui = mui,
                     muin = muin,
                     muHati = mean(y) ) )
  #success = success,
  #N = Nmax ) )
}

draw_withinstudy_result(mui = 0.1,
                        t2w = 0,
                        sd.y = 0.3,
                        m = 30 )


# draw_hacked_study( Nmax = Inf,
#              Mu = 0.1,
#              T2 = 0.1,
#              m = 50,
#              t2w = .5,
#              se = 1)

# # assumes rho=0; o.w. will need to be passed stats from previous draw
# draw_one_unhacked_study = function(Mu,  # overall mean for meta-analysis
#                                    T2,  # across-study heterogeneity
#                                    m,  # sample size for this study
#                                    t2w,  # within-study heterogeneity
#                                    se) {  # TRUE SE for this study
#   # draw subject-level data from this study's population effect
#   y = rnorm( mean = mui,
#              sd = sd.y,
#              n = m)
#   
#   # run a one-sample t-test
#   if ( rho == 0 ) {
#     pval = t.test(y,
#                   alternative = "two.sided")$p.value
#     tstat = t.test(y,
#                    alternative = "two.sided")$statistic
#   } else {
#     stop("rho>0 case not implemented yet")
#   }
#   
#   
#   if (hack == "signif") success = (pval < 0.05)
#   if (hack == "affirm") success = (pval < 0.05 & mean(y) > 0)
#   
#   return( data.frame(pval = pval,
#                      tstat = tstat,
#                      tcrit = qt(0.975, df = m-1),
#                      mui = mui,
#                      muHati = mean(y) ) )
# }



draw_one_unhacked_study(0, 0.1, 30, 0, .2)




