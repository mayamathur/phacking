
# NOTES ---------------------------------------------------------------

# "draw" always refers to within-study draw
# "study set" refers to all draws, observed and observed, for a given study


# DATA SIMULATION ---------------------------------------------------------------



# *note that the number of reported, hacked studies might be less than k.hacked
#  if all Nmax draws are unsuccessful

# also note that for the unhacked study sets, the single published result could be 
# affirmative OR nonaffirmative

sim_meta = function(Nmax,  # max draws to try
                    Mu,  # overall mean for meta-analysis
                    T2,  # across-study heterogeneity
                    
                    # study parameters, assumed same for all studies:
                    m,  # sample size for this study
                    t2w,  # within-study heterogeneity
                    se,  # TRUE SE
                    
                    rho = 0,
                    return.only.published = FALSE,
                    hack,  # "affirm" or "signif" only
                    
                    # args not passed to sim_one_study_set:
                    k,  # number of studies
                    k.hacked  # number of hacked studies
                    
) {
  
  # # test only
  # Nmax = 20
  # Mu = 0.1
  # T2 = 0.1
  # m = 50
  # t2w = .5
  # se = 1
  # hack = "affirm"
  # return.only.published = FALSE
  # rho=0
  # k = 30
  # k.hacked = 20
  
  #browser()
  
  # collect arguments
  .args = mget(names(formals()), sys.frame(sys.nframe()))
  # remove unnecessary args for sim_one_study_set
  .args = .args[ !names(.args) %in% c("k", "k.hacked")]
  
  
  if ( hack == "no" ) stop("hack should only be 'affirm' or 'signif' for this fn")
  
  k.unhacked = k - k.hacked
  

  
  ### Simulate the unhacked studies 
  if ( k.unhacked > 0 ) {
    for ( i in 1:(k - k.hacked) ) {
      # for unhacked studies, need to change argument "hack"

      .argsUnhacked = .args
      .argsUnhacked[ names(.args) == "hack" ] = "no"
      
      # might be multiple rows if return.only.published = FALSE
      newRows = do.call( sim_one_study_set, .argsUnhacked )
      
      # add study ID
      newRows = newRows %>% add_column( .before = 1,
                                        study = i )
      
      if ( i == 1 ) .dat = newRows else .dat = rbind( .dat, newRows )
    }
  }
  
  # #@DEBUGGING
  # browser()
  # .dat %>% group_by(Di) %>% summarise(n(), mean(mui), mean(muin), mean(yi))
  
  ### Simulate hacked studies
  if ( k.hacked > 0 ) {
    if ( exists(".dat") ) startInd = max(.dat$study) + 1 else startInd = 1
    
    for ( i in startInd:(startInd + k.hacked - 1) ) {
      # for unhacked studies, no need to change argument "hack"
      .argsHacked = .args
      
      # might be multiple rows if return.only.published = FALSE
      newRows = do.call( sim_one_study_set, .argsHacked )
      
      # add study ID
      newRows = newRows %>% add_column( .before = 1,
                                        study = i )
      
      if ( i == 1 ) .dat = newRows else .dat = rbind( .dat, newRows )
    }
  }
  
  return(.dat)
  
}



# d = sim_meta(Nmax = 20,
#              Mu = 0.1,
#              T2 = 0.1,
#              m = 50,
#              t2w = .5,
#              se = 1,
#              hack = "affirm",
#              return.only.published = FALSE,
# 
#              k = 30,
#              k.hacked = 10
# 
# )
# 
# 
# 
# 
# nrow(d)
# 
# # look at the published results only
# d %>% filter(Di == 1 ) %>%
#   group_by(hack) %>%
#   summarise( n(),
#              mean(affirm),
#              mean(mui),
#              mean(yi) )








# simulate a single study from potentially heterogeneous meta-analysis distribution 
#  and from its own heterogeneous distribution

#~~ bm: want to look at this for rho>0
#  to do that, need to draw each t-stat from Peng's conditional distribution
# e.g.: when correlated, do t-stats bunch up more?

# similar to phack_study in earlier helper script (2020-6 folder)

# for hack = "no":
# draws exactly Nmax results and treats the last one as the reported one
#  so the final result could be affirmative or nonaffirmative

# for hack = "affirm" or "signif":
# draws until affirm or signif result is obtained or Nmax is reached
# then reports the last result

sim_one_study_set = function(Nmax,  # max draws to try
                             Mu,  # overall mean for meta-analysis
                             T2,  # across-study heterogeneity
                             m,  # sample size for this study
                             t2w,  # within-study heterogeneity
                             se,  # TRUE SE for this study
                             rho = 0,
                             return.only.published = FALSE,
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
  
  # TRUE SD (not estimated)
  sd.y = se * sqrt(m)
  
  # collect all args from outer fn, including default ones
  .args = mget(names(formals()), sys.frame(sys.nframe()))
  .args$mui = mui
  .args$sd.y = sd.y
  
  
  stop = FALSE  # indicator for whether to stop drawing results
  N = 0  # counts draws actually made
  
  
  while( stop == FALSE & N < Nmax ) {
    
    newRow = do.call( make_one_draw, .args )
    
    N = N + 1
    
    # add new draw to dataset
    if ( N == 1 ) d = newRow
    if ( N > 1 ) d = rbind( d, newRow )
    
    # check if it's time to stop drawing results
    if (hack == "signif") stop = (newRow$pval < 0.05)
    if (hack == "affirm") stop = (newRow$pval < 0.05 & newRow$yi > 0)
    # if this study set is unhacked, then success is just whether we've reached the 
    if ( hack == "no") stop = (N == Nmax)
  }
  
  # number of draws made
  d$N = N
  
  d$hack = hack
  
  # convenience indicators for significance and affirmative status
  d$signif = d$pval < 0.05
  d$affirm = (d$pval < 0.05 & d$yi > 0)
  
  # publication indicator
  if ( hack == "signif" ) d$Di = (d$signif == TRUE)
  if (hack == "affirm") d$Di = (d$affirm == TRUE)
  # if no hacking, assume only last draw is published
  if ( hack == "no" ) {
    d$Di = 0
    d$Di[ length(d$Di) ] = 1
  }
  
  if ( return.only.published == TRUE ) d = d[ d$Di == 1, ]
  
  return(d)
  
}


# d = sim_one_study_set(Nmax = 2,
#                       Mu = 0.1,
#                       T2 = 0.1,
#                       m = 50,
#                       t2w = .5,
#                       se = 1,
#                       hack = "signif",
#                       return.only.published = TRUE)
# nrow(d)
# d


# # sanity check by simulation
# for ( i in 1:2000 ) {
#   newRows = sim_one_study_set(Nmax = 20,
#                          Mu = 0.1,
#                          T2 = 0.1,
#                          m = 50,
#                          t2w = .1,
#                          se = 1,
#                          hack = "no",
#                          return.only.published = FALSE )
# 
#   if ( i == 1 ) .d = newRows else .d = rbind(.d, newRows)
# 
# }
# 
# # all studies
# # note that conditional on Di == 0, variance might be off because of repeated rows
# # but means should be correct
# .d %>% group_by(Di == 1) %>%
#   summarise(n(),
#             mean(mui),
#             mean(muin),
#             var(muin),
#             mean(yi) )
# # seems fine


# # ~ Sanity check  ---------------------------------------------------------------
# #  if Nmax -> Inf and we hack until affirmative, 
# #   published results should follow truncated t distribution
# # also need to set heterogeneity to 0?
# d = data.frame( matrix(nrow = 500, ncol = 1))
# Mu = 1
# T2 = 0.5
# t2w = 0.3
# se = 1
# 
# d = d %>% rowwise() %>%
#   mutate( sim_one_study_set( Nmax = Inf,
#                              Mu = Mu,
#                              T2 = T2,
#                              m = 50,
#                              t2w = t2w,
#                              se = se,
#                              hack = "affirm",
#                              return.only.published = TRUE ) )
# 
# summary(d$N)
# 
# # calculate noncentrality parameter
# ncp = Mu / sqrt( T2 + t2w + se^2 )
# 
# # compare to truncated t
# qqtrunc(x = d$tstat,
#         spec = "t",
#         ncp = ncp,
#         df = 50-1,
#         # since I simulated iid studies here, the truncation cutoff is always the same
#         a = d$tcrit[1] )
# 
# # vs. actually drawing directly from truncated t
# x = rtrunc( n = 500,
#             spec = "t",
#             ncp = ncp,
#             df = 50-1,
#             a = d$tcrit[1])
# 
# plot( density(x) ) +
#   lines( density(d$tstat), col = "red")
# # looks great! :)




# draws one unbiased result within the study
make_one_draw = function(mui,  # mean for this study set
                         t2w,
                         sd.y,  # TRUE SD
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
    
    test = t.test(y,
                  alternative = "two.sided")
    
    pval = test$p.value
    tstat = test$statistic
    vi = test$stderr^2
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
                     yi = mean(y),
                     vi = vi) )
  #success = success,
  #N = Nmax ) )
}

# make_one_draw(mui = 0.1,
#               t2w = 0,
#               sd.y = 0.3,
#               m = 30 )

# sanity check by simulation
# for ( i in 1:5000 ) {
#   newRow = make_one_draw(mui = 0.1,
#                   t2w = 0.1,
#                   sd.y = 0.3,
#                   m = 30 )
#   
#   if ( i == 1 ) .d = newRow else .d = rbind(.d, newRow)
#   
# }
# 
# .d %>% summarise( 
#                   mean(mui),
#                   mean(muin),
#                   var(muin),
#                   mean(yi) )


#bm: this fn seems fine, but clearly something is wrong with either of the outer ones



