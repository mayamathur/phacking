


log_lkl = Vectorize( function(yi,  # same as d above
                              se, 
                              mui,
                              N,
                              k,
                              prop,
                              alpha = 0.05) {

  crit = qnorm(1 - alpha/2)
  
  if ( yi/se <= crit ) stop("yi/se indicates study is not affirmative. That won't work")
  
  c = -qnorm( (1 - prop) / 2 )
  
  # part from geometric distribution
  # if N = 1, this term does not contribute at all to lkl because exponent = 0
  #  this makes sense
  # then we're just looking at the straight-up likelihood from a Normal, which makes all kinds of sense! 
  num = crit*se - mui
  denom = sqrt( ( (k^2 * yi^2) / c^2 ) + se^2 )
  term1 = pnorm(num/denom)^(N-1)
  # pnorm part (i.e., P(nonaffirm)) is decreasing in mui
  # **so this term incentivizes mui to be large - will this get us into the same Vevea/Woods issue?
  
  # overall, term1 is increasing in N (lkl of getting an affirmative result, as we did, increases when
  #  you do more p-hacking, which makes sense)

  # part from truncated normal
  term2 = dnorm( (yi - mui) / denom )
  
  #browser()
  
  return( log(term1) + log(term2) )
} )




# Paffirm for a single study that takes either t2, sei or k 
# then EN is just 1/this



# function chooses tau such that 100*prop% of true effects are within k-fold of observed estimate,
#  if observed estimate were true mean
EN = function(pval,  # of the observed, affirmative study
              k,  
              prop ){
  
  
  c = -qnorm( (1 - prop) / 2 )
  
  crit = qnorm(.975)  # doesn't use prop because this is referring to critical value for affirmative studies
  num = crit * qnorm( 1 - pval/2 ) 
  denom = sqrt( (k^2/c^2) + ( qnorm( 1 - pval/2 ) )^2 )
  
  Paffirm = 1 - pnorm(num/denom)
  
  return( 1 / Paffirm )
}


# P( at least one affirmative in N draws )
Paffirm_N = function(pval,  # of the observed, affirmative study
                     k,  
                     prop,
                     N){
  
  # 1/EN = Paffirm for 1 draw
  1 - ( 1 - 1/EN(pval = pval,
         k = k,
         prop = prop) )^N
}





##### Just for Sanity-Checking #####
# calculate probability affirmative for a given mean, within-study heterogeneity, and study SE
# straight from ESC!!
Paffirm = Vectorize( function(mu,
                              t2,
                              sei,
                              tails = 1){
  
  denom = sqrt( t2 + sei^2 )
  num1 = qnorm(.975) * sei - mu
  num2 = -qnorm(.975) * sei - mu
  
  Z1 = num1 / denom
  Z2 = num2 / denom
  
  if ( tails == 2 ) {
    # 2-tailed power for each study
    pwr.vec = 1 - pnorm(Z1) + pnorm(Z2)
  }
  
  if ( tails == 1 ) {
    # 1-tailed power for each study
    pwr.vec = 1 - pnorm(Z1)
  }
  
  mean(pwr.vec)
  
}, vectorize.args = c("t2") )



EN_long = function(pval,
                   d,
                   k,
                   prop) {
  
  # e.g., for prop = 0.68, is close to 1
  c = -qnorm( (1 - prop) / 2 )
  ( tau = (k*d) / c )
  
  # sanity check
  expect_equal( pnorm( (d - c*tau - d)/tau ),
                (1-prop)/2 )
  
  # calculate study's SE from its estimate
  sei = d * qnorm(1 - pval/2)
  
  1/Paffirm( mu = 0,
             t2 = tau^2,
             sei = sei, 
             tails = 1)
}

# # sanity-check EN
# # note that EN_long is invariant to choice of d :)
# EN( pval = 0.01, k = 3, prop = 0.4 )
# EN_long( pval = 0.01, k = 3, prop = 0.4, d = .2 )



########################################## FOR SIMULATION ########################################## 

# simulate a single study from potentially heterogeneous meta-analysis distribution 
#  and from its own heterogeneous distribution

phack_study = function(N.max = Inf,
                       Mu,
                       T2,
                       n,
                       t2,
                       se,
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
    
    pval = t.test(y,
                  alternative = "two.sided")$p.value
    
    if (hack == "signif") success = (pval < 0.05)
    if (hack == "affirm") success = (pval < 0.05 & mean(y) > 0)
  }
  
  return( data.frame(pval = pval,
                     N = N,
                     ybar = mean(y) ) )
}

phack_study( N.max = Inf,
             Mu = 0,
             T2 = 0,
             n = 50,
             t2 = 1,
             se = 1)





























