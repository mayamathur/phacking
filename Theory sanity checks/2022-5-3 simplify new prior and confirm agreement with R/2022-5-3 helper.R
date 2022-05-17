

simple_lli = function(.yi,
                      .sei,
                      .mu,
                      .tau,  
                      .crit ) {
  
  # # # TEST
  # .mu=0.1
  # .tau=1
  # .sei=0.5
  # .crit=1.96
  # .yi=0.5
  
  Si = sqrt(.tau^2 + .sei^2)
  alphaU = (.crit*.sei - .mu)/Si
  
  # sanity check:
  # mine = -log( Si* sqrt(2*pi) ) - 0.5 * Si^(-2) * (.yi-.mu)^2 - log( pnorm(alphaU) )
  # termA = dnorm( .yi,
  #                mean = .mu,
  #                sd = Si,
  #                log = TRUE)
  # 
  # termB = pnorm( q = alphaU,
  #                mean = 0,
  #                sd = 1,
  #                log.p = TRUE )
  # 
  # termA - termB
  # 
  # expect_equal(termA - termB, mine) 
  # 
  # # also check vs. dtruncnorm:
  # expect_equal( log( dtruncnorm(x = .yi,
  #            a = -Inf,
  #            b = .crit*.sei,
  #            mean = .mu,
  #            sd = Si) ), mine )
  
  # same as "mine" above
  -log( Si* sqrt(2*pi) ) - 0.5 * Si^(-2) * (.yi-.mu)^2 - log( pnorm(alphaU) )
}

prior = function(mu, tau, k, sei, tcrit) {
  
  if ( length(tcrit) < k ) stop("tcrit must be vector of length k")
  
  # this will be the TOTALS for all observations
  fishinfototal = matrix( 0, nrow = 2, ncol = 2 )
  
  # build a Fisher info matrix for EACH observation
  for (i in 1:k) {
    
    # for this observation
    fishinfo = matrix( NA, nrow = 2, ncol = 2 )
    
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
    
    
    # try simplifying for just the tau=0 case
    if( tau == 0 ) {
      # because e13 + e15 = 0 exactly
      # see iPad "SAPH simplify new prior"
      kss = 0
    }
    
    
    fishinfo[1,1] = -kmm
    fishinfo[1,2] = -kms
    fishinfo[2,1] = -kms
    fishinfo[2,2] = -kss
    
    #if ( is.na( det(fishinfo) ) ) browser()
    
    # add the new fisher info to the total one
    fishinfototal = fishinfototal + fishinfo
  }
  
  # # ad hoc temporary fix for numerical issues
  # fishdet = det(fishinfototal)
  # prior = max(0, fishdet)
  # return( list(Efish = fishinfototal,
  #              det = det(fishinfototal),
  #              prior = prior) )
  
  return( list(Efish = fishinfototal,
               det = det(fishinfototal),
               prior = sqrt( det(fishinfototal) ) ) )
}


prior_simp = function(mu, tau, k, sei, tcrit) {
  
  if ( length(tcrit) < k ) stop("tcrit must be vector of length k")
  
  # this will be the TOTALS for all observations
  fishinfototal = matrix( 0, nrow = 2, ncol = 2 )
  
  # build a Fisher info matrix for EACH observation
  for (i in 1:k) {
    
    # for this observation
    fishinfo = matrix( NA, nrow = 2, ncol = 2 )
    
    Si = sqrt( tau^2 + sei[i]^2 )
    cz = (sei[i] * tcrit[i] - mu) / Si
    dnor = dnorm(cz) # can't use log on this in case it's negative
    pnor = pnorm(cz)
    r = dnor/pnor
    
    kmm = Si^(-2)*(cz*r + r^2 - 1)
    kms = tau*Si^(-3)*r*( cz^2 + cz*r + 1 )
    kss = ( tau^2 * Si^(-4) ) * ( cz^3*r + cz^2*r^2 + cz*r - 2 )
    
    fishinfo[1,1] = -kmm
    fishinfo[1,2] = -kms
    fishinfo[2,1] = -kms
    fishinfo[2,2] = -kss

    # add the new fisher info to the total one
    fishinfototal = fishinfototal + fishinfo
  }
  
  return( list(Efish = fishinfototal,
               det = det(fishinfototal),
               prior = sqrt( det(fishinfototal) ) ) )
}



prior_contour_plot = function(sei) {
  
  dp = expand_grid(.mu = seq(-2, 2, .05),
                   .tau = seq(0, 2, .05))
  
  tcrit = rep( qnorm(0.975), length(sei) )
  k = length(sei)
  
  
  # make plotting dataframe by calculating log-prior for a grid of values (mu, sigma)
  dp = dp %>%
    rowwise() %>%
    mutate( prior = prior_simp(mu = .mu,
                               tau = .tau,
                               k = k,
                               sei = sei,
                               tcrit = tcrit )$prior, )
  
  
  # set up colors for contours
  get_colors = colorRampPalette( c("lemonchiffon1", "chocolate4") )
  myColors = get_colors(n=15)  # chose 11 based on errors from ggplot if it was fewer
  
  
  p = ggplot( data = dp, 
              aes(x = .mu,
                  y = .tau,
                  z = prior) ) +
    
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
  
  
  
  # information about favored tau
  cat( paste("\n\nMost favored tau: ", dp$.tau[ which.max(dp$prior) ] ) )
  
  return(p)
}




