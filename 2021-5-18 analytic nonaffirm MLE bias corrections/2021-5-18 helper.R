



# 2021-5-31: TRUNCATED NORMAL BIAS CORRECTION FNS -----------------------------

# ~ For Jeffreys prior thing -----------------------------
# my own Jeffreys prior
Jeffreys = function( params, .crit ) {
  fisher = expectFisher( params = params,
                         .crit = .crit )
  
  # this is the Jeffreys prior evaluated at these parameter values
  abs( det(fisher) )^{-1/2}
}

# joint posterior 
joint_post = function( params,
                       .xVec,
                       .crit ) {
  

  .mu = params[1]
  .sigma = params[2]
  
  # likelihood (not logged)
  joint.lkl = exp( sum(myLPDF2(mu = .mu,
                               sigma = .sigma,
                               .x = .xVec,
                               .crit = .crit) ) )
  
  # prior
  prior = Jeffreys( params = c(.mu, .sigma),
                    .crit = .crit )
  
  return( joint.lkl * prior )
  
}




# ~ For Godwin/Cordeira thing -----------------------------

# this one isn't done because I didn't get the 3rd order terms yet
godwin = function(params, .crit) {

  # TEST ONLY
  .crit = 1.96

  # inverse of expected Fisher info
  # (gives k_{ij} entries in Godwin's notation)
  fisher = expectFisher(params = params,
                        .crit = .crit)

  invFisher = solve(fisher)

  # bias of mean estimate: Godwin Eq. (6)
  p = nrow(fisher)  # number of parameters
  for ( i in 1:p ) {

    # k^{si}
    invFisher[1,i]

    for (j in 1:p) {
      for (l in 1:p) {


        H11 = Deriv(J1, "mu")


        Deriv()
      }
    }
  }

  #bm: realized I needed higher-order derivatives (work in progress on iPad)

}

godwin(params, .crit)

# 2021-5-31: GENERAL FNS FOR RIGHT-TRUNCATED NORMAL THEORY -----------------------------

# params: (mu, sigma)
# NOTE parameterization in terms of sigma rather than sigma^2
myLPDF = function(params, .x, .crit) {
  mu = params[1]
  sigma = params[2]
  -log( sqrt(2*pi) * sigma) - (.x - mu)^2 / (2 * sigma^2 ) -
    pnorm( q = .crit, mean = mu, sd = sigma, log=TRUE)
}
myLPDF( c(0, .5), 0, 1.96)

# same, but separating parameters
myLPDF2 = function(mu, sigma, .x, .crit) {
  -log( sqrt(2*pi) * sigma) - (.x - mu)^2 / (2 * sigma^2 ) -
    pnorm( q = .crit, mean = mu, sd = sigma, log=TRUE)
}


# Mills ratio for right truncation
mills = function(params, .crit) {
  mu = params[1]
  sigma = params[2]
  uStar = (.crit - mu)/sigma
  dnorm(uStar) / pnorm(uStar) 
}

# matrix (actually vector) of first derivatives of log-lkl, evaluated at .x
myJacobian = function(params, .x, .crit) {
  
  mu = params[1]
  sigma = params[2]
  
  # entry 1: derivative wrt mu
  J1 = ( (.x - mu) / sigma^2 ) +  # this part matches Deriv!
    mills(params = params, .crit = .crit) * (1/sigma)
  
  # entry 2: derivative wrt sigma
  J2 = (-1/sigma) + ( (.x - mu)^2 / sigma^3 ) +
    mills(params = params, .crit = .crit) * (.crit - mu)/sigma^2
  
  return( matrix( c(J1, J2), nrow = 1 ) )
}

# 2 x 2 matrix of second derivatives of log-lkl, evaluated at .x
myHessian = function(params, .x, .crit) {
  
  mu = params[1]
  sigma = params[2]
  mills = mills(params = params, .crit = .crit)
  uStar = (.crit - mu)/sigma
  
  # entry 1,1: dl/dmu^2
  H11 = (-1/sigma^2) + (1/sigma^2) * ( mills^2 + uStar * mills )
  
  # entry 2,2: dl/dsigma^2
  H22 = (1/sigma^2) - ( 3 * (.x - mu)^2 / sigma^4 ) +
    ( mills^2 + uStar * mills ) * ( ( .crit - mu ) / sigma^2 )^2 -
    mills * ( 2 * (.crit - mu) / sigma^3 )
  
  # entry 1,2 and 2,1: dl / dmu dsigma
  H12 = ( -2 * (.x - mu) / sigma^3 ) +
    ( mills^2 + uStar * mills ) * ( ( .crit - mu ) / sigma^3 ) -
    mills / sigma^2
  
  return( matrix( c(H11, H12, H12, H22), nrow = 2 ) )
}


# 2 x 2 matrix of -E[second derivatives of log-lkl]
expectFisher = function(params, .crit) {
  
  mu = params[1]
  sigma = params[2]
  
  # terms that will show up a lot
  mills = mills(params = params, .crit = .crit)
  uStar = (.crit - mu)/sigma
  termA = mills^2 + uStar * mills
  # termB = 2*mills + uStar
  # termC = termA * termB - mills
  
  # entry 11
  F11 = (1/sigma^2) * ( 1 - termA )
  
  # entry 22
  F22 = (1/sigma^2) * (2 - uStar*mills - 3*mills^2 - mills^2*uStar^2 - mills*uStar^3)
  
  # entry 12=21
  F12 = (mills/sigma^2) * (3 - mills*uStar - uStar^2)
  
  return( matrix( c(F11, F12, F12, F22), nrow = 2 ) )
}




# get a particular third derivative of log-lkl
# .entry: a number like "121" to say which derivative we want
#  where parameter 1 is mu and parameter 2 is sigma
# this function will either evaluate the derivative at .x
#  OR will give the negative expectation
myThirdDerivs = function(params,
                         .x = NA,
                         .giveNegExpect = FALSE,  # should we give expectation rather than evaluated at .x?
                         .crit,
                         .entry) {
  
  
  if (.giveNegExpect == FALSE & is.na(.x)) stop("Need to provide .x")
  if (.giveNegExpect == TRUE & !is.na(.x)) warning("Not using .x because you said you wanted expectation")
  
  mu = params[1]
  sigma = params[2]
  
  # terms that will show up a lot
  mills = mills(params = params, .crit = .crit)
  uStar = (.crit - mu)/sigma
  termA = mills^2 + uStar * mills
  termB = 2*mills + uStar
  termC = termA * termB - mills
  
  # for expectation, we'll use moments of truncated normal
  if (.giveNegExpect == TRUE) {
    truncNormalVar = sigma^2*(1 - uStar*mills - mills^2) 
    truncNormalMean = mills*sigma
  }
  
  ### third derivatives
  # entry 111: dl/dmu^3
  if ( .entry == 111 ) {
    
    # not stochastic, so doesn't matter if doing expectation or not
    return( (1/sigma^3) * termC )
  }
  
  # entry 221=212=122: dl/(dsigma^2 dmu)
  if ( .entry %in% c(221, 212, 122) ) {
    
    if ( .giveNegExpect == FALSE ) {
      return( (1/sigma^4) * ( 6*(.x - mu) - 2*(.crit - mu)*termA +
                                (.crit - mu)^2/sigma * termC +
                                2*mills*sigma ) )
    } else {
      return( (-1/sigma^4) * ( 6*truncNormalMean - 2*(.crit - mu)*termA +
                                 (.crit - mu)^2/sigma * termC +
                                 2*mills*sigma ) )
    }
    
  }
  
  # entry 121=211=112: dl/(dmu dsigma dmu)
  if ( .entry %in% c(121, 211, 112) ) {
    return( (1/sigma^3) * ( 2 - 2*termA + uStar*termC ) )
  }
  
  # entry 222: dl/(dsigma^3)
  if ( .entry == 222 ) {
    
    
    if ( .giveNegExpect == FALSE ) {
      return( (1/sigma^4) * ( (-2*sigma) + (12/sigma)*(.x-mu)^2 -
                                ( 4*(.crit - mu)^2/sigma )*termA +
                                ( (.crit - mu)^3/sigma^2 )*termC +
                                6*(.crit - mu)*mills + 2*( (.crit - mu)^2/sigma )*termA ) )
    } else {
      return( (-1/sigma^4) * ( (-2*sigma) + (12/sigma)*truncNormalVar -
                                 ( 4*(.crit - mu)^2/sigma )*termA +
                                 ( (.crit - mu)^3/sigma^2 )*termC +
                                 6*(.crit - mu)*mills + 2*( (.crit - mu)^2/sigma )*termA ) )
    }
    
  }
  
}









