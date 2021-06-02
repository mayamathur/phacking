
# 2021-5-31: FNS FOR RIGHT-TRUNCATED NORMAL THEORY -----------------------------

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

# get a particular third derivative of log-lkl, evaluated at .x
# .entry: a number like "121" to say which derivative we want
#  where parameter 1 is mu and parameter 2 is sigma
myThirdDerivs = function(params, .x, .crit, .entry) {
  
  mu = params[1]
  sigma = params[2]
  
  # terms that will show up a lot
  mills = mills(params = params, .crit = .crit)
  uStar = (.crit - mu)/sigma
  termA = mills^2 + uStar * mills
  termB = 2*mills + uStar
  termC = termA * termB - mills
  
  ### third erivatives wrt mu
  # entry 111: dl/dmu^3
  if ( .entry == 111 ) {
    return( (1/sigma^3) * termC )
  }
  
  # entry 221=212=122: dl/(dsigma^2 dmu)
  if ( .entry %in% c(221, 212, 122) ) {
    return( (1/sigma^4) * ( 6*(.x - mu) - 2*(.crit - mu)*termA +
                              (.crit - mu)^2/sigma * termC +
                              2*mills*sigma ) )
  }
  
  # entry 121=211=112: dl/(dmu dsigma dmu)
  if ( .entry %in% c(121, 211, 112) ) {
    return( (1/sigma^3) * ( 2 - 2*termA + uStar*termC ) )
  }
  
  # entry 222: dl/(dsigma^3)
  if ( .entry == 222 ) {
    return( (1/sigma^4) * ( (-2*sigma) + (12/sigma)*(.x-mu)^2 -
                              ( 4*(.crit - mu)^2/sigma )*termA +
                              ( (.crit - mu)^3/sigma^2 )*termC +
                              6*(.crit - mu)*mills + 2*( (.crit - mu)^2/sigma )*termA ) )
  }
  
  # H11 = (-1/sigma^2) + (1/sigma^2) * ( mills^2 + uStar * mills )
  # 
  # # entry 2,2: dl/dsigma^2
  # H22 = (1/sigma^2) - ( 3 * (.x - mu)^2 / sigma^4 ) +
  #   ( mills^2 + uStar * mills ) * ( ( .crit - mu ) / sigma^2 )^2 -
  #   mills * ( 2 * (.crit - mu) / sigma^3 )
  # 
  # # entry 1,2 and 2,1: dl / dmu dsigma
  # H12 = ( -2 * (.x - mu) / sigma^3 ) +
  #   ( mills^2 + uStar * mills ) * ( ( .crit - mu ) / sigma^3 ) -
  #   mills / sigma^2
  
}









