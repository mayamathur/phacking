
# branch v1 on git: generalizing the fns




correct_meta_phack = function(yi,
                              sei,
                              # sensitivity parameters assumed constant across studies for now
                              n,
                              t2w) {
  # 
  # # TEST ONLY
  # yi = da$yi
  # sei = da$sei
  # n = 50
  # t2w = 1
  
  zi = yi/sei
  crit = qnorm(.975)
  
  # get bias-corrected z's for each *significant* study
  # and leave the nonsignificant ones alone
  zi.adj = zi
  zi.adj[ abs(zi) > crit ] = mle4( zn = zi[ abs(zi) > crit ],  # the observed one
                                   n = n,
                                   t2w = t2w )
  
  # sanity check
  # plot(zi, zi.adj)
  
  yi.adj = zi.adj * sei
  
  library(robumeta)
  m = robu( yi.adj ~ 1, 
            d = data.frame( yi.adj, sei ),
            studynum = 1:length(yi),  # ~~~ ignore clustering for now
            var.eff.size = sei^2 )
  muhat = m$b.r
  t2 = m$mod_info$tau.sq
  mu.lo = m$reg_table$CI.L
  mu.hi = m$reg_table$CI.U
  
  cat( "CORRECTED:", round(muhat, 2), format_CI( lo = mu.lo, hi = mu.hi, digits = 2), sep = " " )
  
  # naive meta-analysis for comparison
  m = robu( yi ~ 1, 
            d = data.frame( yi, sei ),
            studynum = 1:length(yi),  # ~~~ ignore clustering for now
            var.eff.size = sei^2 )
  muhat = m$b.r
  t2 = m$mod_info$tau.sq
  mu.lo = m$reg_table$CI.L
  mu.hi = m$reg_table$CI.U
  
  cat( "\nNAIVE:", round(muhat, 2), format_CI( lo = mu.lo, hi = mu.hi, digits = 2), sep = " " )
  
}




##### 2020-6-3 #####

# the presumably conservative version:
# rho = 0, but n left as a sensitivity parameter

ll4 = Vectorize( function(Z,  # the true parameter
                          zn,  # the observed one
                          n,
                          t2w = 0,
                          include.signif.term = TRUE,
                          select.tails = 2){  # two-tailed or one-tailed selection
  
  #browser()
  
  if (n < 1) stop("n must be at least 1")
  
  # P(signif) or P(affirm) depending on tails
  crit = qnorm(.975)
  if ( abs(zn) < crit ) stop("zn must be significant")
  
  # marginal SD of the Z-statistics
  SD = t2w + 1
  
  Psig.pos = 1 - pnorm( (crit - Z) / SD )
  Psig.neg = pnorm( (-crit - Z) / SD )
  Pnonsig = 1 - Psig.pos - Psig.neg
  
  # log-prob of n-1 N.S. draws
  termA = (n-1) * log(Pnonsig)
  
  # log of inverse-power term for the single significant term
  # negative sign from the ^{-1}
  if ( select.tails == 2 ) termB = -log( Psig.neg + Psig.pos )
  if ( select.tails == 1 ) termB = -log( Psig.neg )
  
  # log-lkl of zn
  # independent of the previous N.S. draws
  termC = dnorm( x = zn,
                 mean = Z,
                 sd = SD,
                 log = TRUE )
  
  if ( include.signif.term == TRUE ) {
    return( termA + termB + termC )
  } else {
    return(termA + termC)
  }
}, vectorize.args = "zn" )

# ll4(Z = 0,  # the true parameter
#     zn = c(1.99,2),  # the observed one
#     n = 30)



mle4 = Vectorize( function( zn,  # the observed one
                            n,
                            t2w = 0,
                            include.signif.term = TRUE,
                            select.tails = 2 ){  
  
  v1 = optim( par = 0,  # starting point for optimization
              f = function(.Z) # the negative ll
                -ll4( Z = .Z,  # the true parameter
                      zn = zn,  # the observed one
                      n = n,
                      t2w = t2w,
                      include.signif.term = include.signif.term,
                      select.tails = select.tails ),
              method = "BFGS" )$par
  
  return(v1)
} )

# mle4( zn = 4,  # the observed one
#       n = 30,
#       t2w = 0 )


# mle4( zn = c(3,4),  # the observed one
#              n = 30,
#              t2w = 0 )

########################################## 2020-5-26 ########################################## 

# lkl of this t-stat given the previous one

log_lkl2 = function(tstat2,  # the current tstat
                    se2,  # se of the current tstat
                    tstat1,  # the previous tstat
                    se1,
                    df,
                    t2w, 
                    #m,  # sample size, assumed same on both draws
                    cv,# covariance of current t-stat with previous one
                    th # true mean in this study
) {
  
  
  # noncentrality parameters of each t-statistic
  # https://en.wikipedia.org/wiki/Noncentral_t-distribution#Use_in_power_analysis
  # ~~ not entirely sure about the way I incorporated t2w here:
  ncp1 = th / sqrt(se1^2 + t2w^2)
  ncp2 = th / sqrt(se2^2 + t2w^2)
  # from the "moments" part of Wiki article
  termA = sqrt(df/2) * gamma( (df-1)/2 ) / gamma(df/2)
  
  # marginal means of the two t-statistics
  # * could check these empirically
  mu1 = ncp1 * termA
  mu2 = ncp2 * termA
  
  # notation as in Ding paper ("On the Conditional Distribution of the Multivariate t Distribution")
  # regression of tstat2 on tstat1
  mu2.1 = mu2 + cv * (1/se2^2) * (tstat1-mu1)
  # = mu2 if cv=0
  
  # Mahalanobis distance of previous tstat from its mean
  d1 = (tstat1 - mu1)^2 * (1/se1^2)
  
  # covariance-related term
  sig22.1 = se2^2 - cv^2 * (1/se1^2)
  # = se2^2 if cv=0
  
  # scale parameter
  scale = ( (df + d1)/(df+1) ) * sig22.1
  # this still depends on tstat1 even if cv = 0?
  # only equals sig22.1 if d1 = 1
  
  library(metRology)
  ll = dt.scaled( x = tstat2,
                  df = df + 1,
                  mean = mu2.1,
                  sd = sqrt(scale),  # ~~ I am pretty sure we need sqrt based on the density in Ding paper
                  log = TRUE )
  
  return( data.frame(mu1 = mu1,
                     mu2 = mu2,
                     ll = ll) )
}


##### MLE for theta based on the above (conditional t) #####

# **note that sometimes the choice of starting point matters! 
#  it can come up with bullshit MLEs (e.g., crazy spikes on graph) otherwise
mle2 = function( tstat2,  # the current tstat
                 se2,  # se of the current tstat
                 tstat1,  # the previous tstat
                 se1,
                 df,
                 t2w, 
                 cv ){
  
  optim( par = 0,  # starting point for optimization
         f = function(x) -log_lkl2(tstat2 = tstat2,  
                                   se2 = se2,  
                                   tstat1 = tstat1,
                                   se1 = se1,
                                   df = df,
                                   t2w = t2w, 
                                   cv = cv,
                                   th = x)$ll,
         method = "BFGS" )$par
}


##### MLE for theta based on the conditional normal #####

# when rho = 0 and include.signif.term = 0, it's just the regular MLE
mle3 = function( z2,  # the current one
                 z1,  # the previous one
                 t2w = 0, 
                 rho,
                 include.signif.term = FALSE, # should we include the leading term, inverse prob of significance?
                 select.tails = 2 ){  
  
  # sanity check: directly maximize
  v1 = optim( par = 0,  # starting point for optimization
              f = function(.th) # the negative ll
                -ll3( th = .th,
                      z2 = z2, 
                      z1 = z1,  
                      t2w = t2w, 
                      rho = rho,
                      include.signif.term = include.signif.term,
                      select.tails = select.tails ),
              method = "BFGS" )$par
  
  # # previous version without signif term
  #   v1 = optim( par = 0,  # starting point for optimization
  #               f = function(.th) # the negative ll
  #                 # conditional lkl of z2 given z1 and autocorrelation
  #                 -dnorm( x = z2,
  #                        mean = .th + rho * (z1 - .th),
  #                        sd = sqrt(1 - rho^2),
  #                        log = TRUE ),
  #               method = "BFGS" )$par
  
  # from manual simplication (not including significance term yet)
  # note that t2w does not affect the MLE if there is no significance term
  #v2 = (z2 - rho*z1) / (1-rho)
  
  # # agree! 
  # print(v1)
  # print(v2)
  
  return(v1)
}




##### 2020-5-29: conditional lkl for Z-stat given autocorrelation and heterogeneity #####

ll3 = function(th,
               z2,  # the current one
               z1,  # the previous one
               t2w = 0, 
               rho = 0,
               include.signif.term = TRUE,
               select.tails = 2){  # two-tailed or one-tailed selection
  
  #browser()
  
  # conditional mean and variance given Z_{n-1}
  Zc = th + rho * (z1 - th)
  Vc = (1 - rho^2) * (t2w + 1)
  SDc = sqrt(Vc)
  
  # P(signif) or P(affirm) depending on tails
  crit = qnorm(.975)
  
  if ( select.tails == 2 ) termA = log( pnorm( (-crit - Zc) / SDc ) + 1 - pnorm( (crit - Zc) / SDc ) )
  if ( select.tails == 1 ) termA = log( 1 - pnorm( (crit - Zc) / SDc ) )
  
  # conditional lkl of z2 given z1 and autocorrelation
  termB = dnorm( x = z2,
                 mean = Zc,
                 sd = SDc,
                 log = TRUE )
  
  if ( include.signif.term == TRUE ) {
    # negative sign is from taking log of an inverse
    return( -termA + termB )
  } else {
    return(termB)
  }
}


########################################## BEFORE 2020-5-26: SAVE ########################################## 

my_ggsave = function(name,
                     width,
                     height,
                     .results.dir = results.dir,
                     .overleaf.dir = overleaf.dir) {
  
  setwd(.results.dir)
  ggsave( name,
          width = width, 
          height = width)
  
  setwd(.overleaf.dir)
  ggsave( name,
          width = width, 
          height = width)
}


##### For Reparameterizing Heterogeneity #####
# for reparameterizing the heterogeneity
tau = Vectorize( function(k, 
                          yi,
                          prop = 0.95) {  # percent of true effects inside caliper
  c = -qnorm( (1-prop)/2)  
  return( k*yi/c)
} )

# # sanity check: proportion below lower limit of the caliper
# yi = 3
# k = 2
# t2w = tau(k, yi)^2
# pnorm( ( (yi - k*yi) - yi ) / sqrt(t2w) ); 0.025



##### P(signif) #####
Psig = Vectorize( function(thi,  # true mean in this study
                           
                           # provide either these two...
                           t2w = NA,  # intra-study heterogeneity
                           sei = NA,  # SE 
                           
                           # or this
                           k = NA,
                           pval = NA,
                           prop = 0.95,  # proportion inside caliper
                           
                           alpha = 0.05,
                           tails = 2){  # 1 tail = affirmative
  
  crit = qnorm( 1 - alpha/2 )
  
  if ( !is.na(t2w) & !is.na(sei) ) {
    
    denom = sqrt( t2w + sei^2 )
    num1 = crit * sei - thi  # upper tail (significant and positive)
    num2 = -crit * sei - thi  # lower tail (significant and negative)
    
    Z1 = num1 / denom
    Z2 = num2 / denom
    
    if ( tails == 2 ) return(1 - pnorm(Z1) + pnorm(Z2))
    if ( tails == 1 ) return(1 - pnorm(Z1))
    
  }
  
  # the other parameterization
  if ( !is.na(k) & !is.na(pval) ) {
    if ( thi != 0 ) stop("Can't reparameterize because thi!=0")
    
    phi.term = qnorm( 1 - pval/2 )
    c = -qnorm( (1 - prop)/2)
    
    num = -crit
    denom = phi.term * sqrt( k^2/c^2 
                             + phi.term^(-2) )
    
    if ( tails == 2 ) return( 2 * pnorm(num/denom) )
    if ( tails == 1 ) return( pnorm(num/denom) )
  }
  
} )

# # sanity check: should be 0.05
# Psig(thi = 0,  # true mean in this study
#      t2w = 0,
#      sei = 5,  # shouldn't matter
#      alpha = 0.05,
#      tails = 2)

# #  when reparameterizing in terms of study's p-value and using thi=0, should be invariant to yi and sei
# #  as long as its p-value is held constant
# yi = 5
# sei = 6
# k = 1.3
# Psig(thi = 0,  # true mean in this study
#      t2w = tau(k = k,
#                yi = yi)^2,
#      sei = sei)
# 
# Psig(thi = 0,  # true mean in this study
#      pval = 2 * (1-pnorm(yi/sei)),
#      k = k)


Psig(thi = 0,  # true mean in this study
     pval = .05,
     k = 2)

##### Expected Draws #####
# allows for correlated draws
# kappa: P(p_i >= 0.05 | p_{i-1} >= 0.05) (OPPOSITE PARAMETERIZATION FROM THEORY I WROTE UP)
#  kappa could range from 1-Psig (no autocorrelation) to 1 (complete autocorrelation)
EN = Vectorize( function(kappa = NA,  
                         Psig,
                         independent = TRUE){
  
  if ( independent == TRUE & !is.na(kappa) ) stop("If independent, don't provide kappa")
  
  if ( independent == TRUE ) kappa = 1 - Psig
  
  return( 1 + (1 - Psig) / (1-kappa) )
} )

# EN( Psig = .05, kappa = 1, independent = FALSE )


##### P( at least one significant in N draws ) #####
Psig_N = function(kappa = NA,  # of the observed, affirmative study
                  Psig,
                  independent = TRUE,
                  N){
  
  # 1/EN = Paffirm for 1 draw
  1 - ( 1 - 1/EN(kappa = kappa,
                 Psig = Psig,
                 independent = independent) )^N
}



##### Likelihood for 1 Study with N->Inf and Autocorrelation < 1 #####

log_lkl1 = Vectorize( function(m, # sample size
                               
                               pval = NA,  # or could use the yi, alternatively
                               yi = NA,
                               
                               SDi,  # SD, *not* SE
                               thi, # the parameter for which we'll get the MLE
                               alpha = 0.05,
                               regular.MLE = FALSE ){  # as a sanity check, just get the regular MLE
  
  # # test only
  # m = 50
  # pval = NA
  # yi = 1
  # SDi = 1
  # thi = 0.3
  # alpha = 0.05
  
  #browser()
  
  SEi = SDi/sqrt(m)
  # ~~ not always exactly true:
  df = m - 1 
  
  # can work with either yi or pval
  if ( !is.na(yi) & !is.na(pval) ) stop("Provide either yi or pval, but not both")
  if ( is.na(yi) & is.na(pval) ) stop("Provide either yi or pval")
  if ( !is.na(yi) & is.na(pval) ) pval = 2 * ( 1 - pt( abs(yi/SEi), df = df ) )
  
  if ( pval >= 0.05 ) stop("pval needs to be less than 0.05")
  
  
  # get observed t-stat (based on H0) from pval
  # i.e., yi/SEi
  t = qt( 1 - pval/2, df = df )
  # sanity check: if providing yi, should match
  # yi/SEi; t
  
  # noncentrality parameter
  # https://en.wikipedia.org/wiki/Noncentral_t-distribution#Use_in_power_analysis
  ncp = thi/SEi
  
  # first component (normalization term for truncated lkl)
  #  is the inverse-power evaluated at thi
  crit = qt( 1 - alpha/2, df = df )  # close to 1.96 for large df
  ( power = pt( q = -crit, df = df, ncp = ncp ) + 1 - pt( q = crit, df = df, ncp = ncp ) )
  
  # truncated density
  # middle term will always be log(1)
  if ( regular.MLE == FALSE ) return( log(1/power) + log(pval < 0.05) + log( dt( x = t, df = df, ncp = ncp ) ) )
  if ( regular.MLE == TRUE ) return( log( dt( x = t, df = df, ncp = ncp ) ) )
  
  # this wold 
  
} )


##### MLE based on the above #####

mle = function( m, # sample size
                pval = NA,  # or could use the yi, alternatively
                yi = NA,
                SDi,  # SD, *not* SE
                regular.MLE = FALSE ){
  optim( par = 0,  # starting point for optimization
         f = function(x) -log_lkl1(m = m, # sample size
                                   pval = pval,  # or could use the yi, alternatively
                                   yi = yi,
                                   SDi = SDi,  # SD, *not* SE
                                   thi = x,
                                   regular.MLE = regular.MLE),
         method = "BFGS" )$par
}







# log_lkl = Vectorize( function(yi,  # same as d above
#                               se, 
#                               mui,
#                               N,
#                               k,
#                               prop,
#                               alpha = 0.05) {
# 
#   crit = qnorm(1 - alpha/2)
#   
#   if ( yi/se <= crit ) stop("yi/se indicates study is not affirmative. That won't work")
#   
#   c = -qnorm( (1 - prop) / 2 )
#   
#   # part from geometric distribution
#   # if N = 1, this term does not contribute at all to lkl because exponent = 0
#   #  this makes sense
#   # then we're just looking at the straight-up likelihood from a Normal, which makes all kinds of sense! 
#   num = crit*se - mui
#   denom = sqrt( ( (k^2 * yi^2) / c^2 ) + se^2 )
#   term1 = pnorm(num/denom)^(N-1)
#   # pnorm part (i.e., P(nonaffirm)) is decreasing in mui
#   # **so this term incentivizes mui to be large - will this get us into the same Vevea/Woods issue?
#   
#   # overall, term1 is increasing in N (lkl of getting an affirmative result, as we did, increases when
#   #  you do more p-hacking, which makes sense)
# 
#   # part from truncated normal
#   term2 = dnorm( (yi - mui) / denom )
#   
#   #browser()
#   
#   return( log(term1) + log(term2) )
# } )
# 
# 
# 
# 
# 
# 
# 
# 
# # function chooses tau such that 100*prop% of true effects are within k-fold of observed estimate,
# #  if observed estimate were true mean
# EN = function(pval,  # of the observed, affirmative study
#               k,  
#               prop ){
#   
#   
#   c = -qnorm( (1 - prop) / 2 )
#   
#   crit = qnorm(.975)  # doesn't use prop because this is referring to critical value for affirmative studies
#   num = crit * qnorm( 1 - pval/2 ) 
#   denom = sqrt( (k^2/c^2) + ( qnorm( 1 - pval/2 ) )^2 )
#   
#   Paffirm = 1 - pnorm(num/denom)
#   
#   return( 1 / Paffirm )
# }
# 
# 
# 
# 
# 
# 
# 
# ##### Just for Sanity-Checking #####
# # calculate probability affirmative for a given mean, within-study heterogeneity, and study SE
# # straight from ESC!!
# Paffirm = Vectorize( function(mu,
#                               t2,
#                               sei,
#                               tails = 1){
#   
#   denom = sqrt( t2 + sei^2 )
#   num1 = qnorm(.975) * sei - mu
#   num2 = -qnorm(.975) * sei - mu
#   
#   Z1 = num1 / denom
#   Z2 = num2 / denom
#   
#   if ( tails == 2 ) {
#     # 2-tailed power for each study
#     pwr.vec = 1 - pnorm(Z1) + pnorm(Z2)
#   }
#   
#   if ( tails == 1 ) {
#     # 1-tailed power for each study
#     pwr.vec = 1 - pnorm(Z1)
#   }
#   
#   mean(pwr.vec)
#   
# }, vectorize.args = c("t2") )
# 
# 
# 
# EN_long = function(pval,
#                    d,
#                    k,
#                    prop) {
#   
#   # e.g., for prop = 0.68, is close to 1
#   c = -qnorm( (1 - prop) / 2 )
#   ( tau = (k*d) / c )
#   
#   # sanity check
#   expect_equal( pnorm( (d - c*tau - d)/tau ),
#                 (1-prop)/2 )
#   
#   # calculate study's SE from its estimate
#   sei = d * qnorm(1 - pval/2)
#   
#   1/Paffirm( mu = 0,
#              t2 = tau^2,
#              sei = sei, 
#              tails = 1)
# }
# 
# # # sanity-check EN
# # # note that EN_long is invariant to choice of d :)
# # EN( pval = 0.01, k = 3, prop = 0.4 )
# # EN_long( pval = 0.01, k = 3, prop = 0.4, d = .2 )



########################################## FOR SIMULATION ########################################## 

# simulate a single study from potentially heterogeneous meta-analysis distribution 
#  and from its own heterogeneous distribution

#~~ bm: want to look at this for rho>0
#  to do that, need to draw each t-stat from Peng's conditional distribution
# e.g.: when correlated, do t-stats bunch up more?

phack_study = function(N.max = Inf,
                       Mu,
                       T2,
                       n,
                       t2,
                       se,
                       rho = 0,  # correlation of Tn with T_{n-1} 
                       hack = "signif") {  # hack until significant or until affirmative?
  
  
  # simulate the study
  # its marginal heterogeneity reflects meta-analysis heterogeneity and within-study heterogeneity
  mui = Mu + rnorm(mean = 0,
                   sd = sqrt(T2+t2),
                   n = 1)
  
  success = FALSE
  N = 0
  
  sd.y = se * sqrt(n)
  
  while ( success == FALSE ) {
    
    if ( N == N.max ) break
    N = N + 1
    
    y = rnorm( mean = mui,
               sd = sd.y,
               n = n)
    
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
                     ybar = mean(y),
                     N = N
  ) )
}

phack_study( N.max = Inf,
             Mu = 0,
             T2 = 0,
             n = 50,
             t2 = 1,
             se = 1)





























