



# 2021-5-31: TRUNCATED NORMAL BIAS CORRECTION FNS -----------------------------

# ~ For Jeffreys prior thing -----------------------------
# my own Jeffreys prior
# uses Fisher info for ONE observation
Jeffreys = function( params, .crit ) {
  
  # the one-observation expected Fisher infoo
  fisher = expectFisher( params = params,
                         .crit = .crit )
  
  # this is the Jeffreys prior evaluated at these parameter values
  sqrt( det(fisher) )
}


# ### Sanity check:
# # in the untruncated case (very high trunc point),
# #  prior should be as for the regular normal
# # https://stats.stackexchange.com/questions/156199/jeffreys-prior-for-normal-distribution-with-unknown-mean-and-variance
# # http://www.stat.cmu.edu/~kass/papers/rules.pdf
# # should be PROPORTIONAL TO 1/sigma^2 when parameterized via sigma
# params = c(1, 2)
# mu = params[1]
# sigma = params[2]
# crit = 999  # effectively untruncated
# 
# ( prior = Jeffreys(params = params,
#                      .crit = crit) )
# expect_equal( prior, 1/(sigma^2) )
# 
# # determinant - matches
# # latter is from (2/sigma^2)*(1/sigma^2)
# myF[1,1]*myF[2,2]; 2/sigma^4
# 
# # prior vs. sqrt(determinant) - matches
# prior; sqrt(2)/sigma^2
# # MATCHES because the *2 is just a proportionality thing


# joint posterior (unnormalized; it's just prior * lkl)
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






# ~ For Godwin/Cordeiro thing -----------------------------

# # this one isn't done because I didn't get the 3rd order terms yet
# # biasedParamIndex: 1 for mu and 2 for sigma
# godwinBias = function(biasedParamIndex,
#                       params,
#                       .crit) {
#   
#   # TEST ONLY
#   #.crit = 1.96
#   
#   mu = params[1]
#   sigma = params[2]
#   
#   # inverse of expected Fisher info
#   # (gives k_{ij} entries in Godwin's notation)
#   fisher = expectFisher(params = params,
#                         .crit = .crit)
#   
#   invFisher = solve(fisher)
#   
#   # bias of mean estimate: Godwin Eq. (6)
#   p = nrow(fisher)  # number of parameters
#   
#   outerSum = 0
#   
#   for ( i in 1:p ) {
#     
#     # k^{si}
#     term1 = invFisher[biasedParamIndex,i]
#     
#     innerSum = 0
#     
#     for (j in 1:p) {
#       for (l in 1:p) {
#         
#         # k_{ij}^(l) in Godwin's notation
#         # note sign reversals throughout because Godwin's "k" terms are expectations
#         #  rather than negative expectations
#         term2 = -1 * expectFisherDerivs( params = params,
#                                          .crit = .crit,
#                                          .entry = as.numeric( paste(i, j, l, sep = "") ) )
#         
#         # k_{ijl} in Godwin's notation
#         term3 = -0.5 * myThirdDerivs(params = params,
#                                      .giveNegExpect = TRUE, 
#                                      .crit = .crit,
#                                      .entry = as.numeric( paste(i, j, l, sep = "") ) )
#         
#         # k^{jl} in Godwin's notation
#         term4 = invFisher[j,l]
#         
#         innerSum = innerSum + (term2 - term3)*term4
#         
#       }
#     }
#     
#     outerSum = outerSum + (termA*innerSum)
#   }
#   
#   return(outerSum)
# }




#godwin(params, .crit)



# ~~ Matrix version: -----------------------------------------

# gives one value of a_{ij}^{(l)} in Godwin notation (pg 1890)
# .entry: ijl
# n: sample size
get_A_entry = function( params,
                        n,
                        .crit,
                        .entry,
                        useZhou = TRUE) {
  
  # # use mine 
  # if ( useZhou == FALSE ) {
  #   
  #   
  #   # k_{ij}^(l) in Godwin's notation
  #   # note sign reversals throughout because Godwin's "k" terms are expectations
  #   #  rather than negative expectations
  #   # IMPORTANT: because expectFisherDerivs and myThirdDerivs
  #   #  are for a SINGLE observations, need to multiply them by n
  #   term1 = -1 * n * expectFisherDerivs( params = params,
  #                                        .crit = .crit,
  #                                        .entry = .entry )
  #   
  #   # k_{ijl} in Godwin's notation
  #   # again, sign reversal
  #   term2 = -0.5 * n * myThirdDerivs(params = params,
  #                                    .giveNegExpect = TRUE, 
  #                                    .crit = .crit,
  #                                    .entry = .entry )
  # }
  
  if ( useZhou == TRUE ) {
    
    # k_{ij}^(l) in Godwin's notation
    # note sign reversals because Godwin's "k" terms are expectations
    #  rather than negative expectations
    # IMPORTANT: because expectFisherDerivs and myThirdDerivs
    #  are for a SINGLE observation, need to multiply them by n
    
    mu = params[1]
    sigma = params[2]
    
    
    term1 = -1 * n * zhou_expect_fisher_derivs( mu = params[1],
                                                sigma = params[2],
                                                .crit = crit,
                                                .entry = .entry )
    
    
    
    # k_{ijl} in Godwin's notation
    # NO sign reversal this time because zhou_third_derivs
    #  didn't take the negative
    term2 = 0.5 * n * zhou_third_derivs( mu = params[1],
                                          sigma = params[2],
                                          .crit = crit,
                                          .entry = .entry )
  }
  
  return( term1 - term2 )
}


godwinBiasMatrix = function(params, n, crit) {
  
  # 2 x 2 matrices for each parameter
  A1 = matrix( c( get_A_entry( params = params, 
                               n = n,
                               .crit = crit,
                               .entry = 111),
                  get_A_entry( params = params, 
                               n = n,
                               .crit = crit,
                               .entry = 121),
                  get_A_entry( params = params, 
                               n = n,
                               .crit = crit,
                               .entry = 211),
                  get_A_entry( params = params, 
                               n = n,
                               .crit = crit,
                               .entry = 221)
  ), byrow = TRUE, nrow = 2 )
  A1
  
  A2 = matrix( c( get_A_entry( params = params,
                               n = n,
                               .crit = crit,
                               .entry = 112),
                  get_A_entry( params = params,
                               n = n,
                               .crit = crit,
                               .entry = 122),
                  get_A_entry( params = params, 
                               n = n,
                               .crit = crit,
                               .entry = 212),
                  get_A_entry( params = params, 
                               n = n,
                               .crit = crit,
                               .entry = 222)
  ), byrow = TRUE, nrow = 2 )
  A2
  
  # 2 x 4 matrix of third derivs for both parameters
  A = cbind( A1, A2 )
  A
  
  # inverse of expected Fisher info
  # (gives k_{ij} entries in Godwin's notation)
  # IMPORTANT: multiply 1-observation Fisher info by n
  # USING ZHOU
  K = n * expectFisherZhou(params = params,
                       .crit = crit)
  
  Kinv = solve(K)
  
  Kinv.vec = matrix( c(Kinv[1,1],
                       Kinv[2,1],
                       Kinv[1,2],
                       Kinv[2,2]),
                     byrow = TRUE, nrow = 4 )
  
  bias = Kinv %*% A %*% Kinv.vec
  
  # return everything
  return( llist( bias = as.matrix(bias), 
                 A = as.matrix(A),
                 Kn = as.matrix(K),  # call it "Kn" to clarify that it's for n observations
                 Kninv = as.matrix(Kinv) ) )
}

# ### Sanity check:
# # in the untruncated case (very high trunc point),
# #  Fisher info should be as for the regular MLE
# params = c(0, 1)
# mu = params[1]
# sigma = params[2]
# # NOTE: FNS BELOW WILL BREAK IF YOU CHOOSE WRONG SIGN FOR THIS
# #  BECAUSE THE NORMALIZATION TERM IN PDF AND DERIVS WILL BE O
# # CRIT NEEDS TO BE LIKE -99 FOR LEFT-TRUNC AND LIKE 99 FOR RIGHT-TRUNC
# crit = -99  # shouldn't be biased if truncation point is extreme
# n = 20
# 
# res = godwinBiasMatrix(params = params,
#                        n = n,
#                        crit = crit)
# ( theoryBias = res$bias )
# 
# expect_equal( as.numeric(res$Kinv[1,1]), sigma^2/n )
# expect_equal( as.numeric(res$Kinv[1,2]), 0 )
# expect_equal( as.numeric(res$Kinv[2,1]), 0 )
# # for this term, (2*sigma^4/n) is the variance of sigma^2;
# #  we parameterized in terms of sigma, so (0.5*sigma^(-0.5))^2
# #  is from the delta method
# expect_equal( as.numeric(res$Kinv[2,2]),
#               (2*sigma^4/n) * (0.5*sigma^(-0.5))^2 )
# # everything matches :)



# 2021-5-31: GENERAL FNS FOR RIGHT-TRUNCATED NORMAL THEORY -----------------------------

# params: (mu, sigma)
# NOTE params argument given in terms of sigma rather than sigma^2
# for ONE observation
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
# for ONE observation
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
# for ONE observation
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
# for ONE observation
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
  # untruncated case (huge crit): mills is 0 and so
  #  F11 = (1/sigma^2), which matches theory for regular normal distribution :)
  F11 = (1/sigma^2) * ( 1 - termA )
  
  # entry 22
  # untruncated case (huge crit): mills is 0 and so F22 = (2/sigma^2)
  # that seems correct because we parameterized using sigma:
  # http://confluence.marcuschiu.com/display/NOT/Fisher+Information+-+Normal+Distribution
  F22 = (1/sigma^2) * (2 - uStar*mills - mills^2*uStar^2 - mills*uStar^3)
  
  # entry 12=21
  # untruncated case (huge crit): mills is 0 and so
  #  F12 = 0, which matches theory for regular normal distribution :)
  F12 = (mills/sigma^2) * (1 - mills*uStar - uStar^2)
  
  return( matrix( c(F11, F12, F12, F22), nrow = 2 ) )
}


# ### Sanity check:
# # in the untruncated case (very high trunc point),
# #  Fisher info should be as for the regular MLE
# params = c(1, 2)
# mu = params[1]
# sigma = params[2]
# crit = 999  # effectively untruncated
# 
# ( myF = expectFisher(params = params,
#                      .crit = crit) )
# 
# expect_equal( myF[1,1], 1/(sigma^2) )
# expect_equal( myF[1,2], 0 )
# expect_equal( myF[2,1], 0 )
# expect_equal( myF[2,2], 2/(sigma^2) )


# get a particular derivative of (my) Fisher info
# .entry: a number like "121" to say which derivative we want
# i.e., "121" is d/dmu { -E[ d^2 l/ dmu dsigma ] }
#  where parameter 1 is mu and parameter 2 is sigma
# for ONE observation 
expectFisherDerivs = function(params,
                              .crit,
                              .entry) {
  
  mu = params[1]
  sigma = params[2]
  
  # terms that will show up a lot
  mills = mills(params = params, .crit = .crit)
  uStar = (.crit - mu)/sigma
  uStar2 = (.crit - mu)/sigma^2
  termA = mills^2 + uStar * mills  
  termB = 2*mills + uStar
  #termC = termA * termB - mills
  
  
  ### the derivatives
  if ( .entry == 111 ) {
    return( (-1/sigma^3) * ( (2*mills + uStar)*termA - mills ) )
  }
  
  # entry 121=211
  if ( .entry %in% c(121, 211) ) {
    
    # try checking an intermediate term
    # this one matches Deriv() and is equal to 0
    # return( mills/sigma^2*( -termA*uStar - mills*(-1/sigma) - 2*uStar*(-1/sigma) ) +
    #           termA*(1/sigma^2)*(3 - mills*uStar - uStar^2) )
    
    # # a later term
    # # matches final one and is wrong!
    # return( (1/sigma^2)*( (3 - 2*mills*uStar - uStar^2)*termA*(1/sigma) -
    #                         (mills + 2*uStar)*(-1/sigma) ) )
    
    return( (1/sigma^3) * ( (3 - 2*mills*uStar - uStar^2)*termA + mills^2 + 2*uStar*mills ) )
    
  }
  
  # entry 221
  if ( .entry == 221 ) {
    
    # handy derivatives
    d1 = termA * (1/sigma)  # derivative of termA wrt mu
    d2 = (-1/sigma) # derivative of uStar wrt mu
    
    # intermediate #1
    # return( (1/sigma^2) * ( -uStar*d1 - d2*mills - 6*mills*d1 -
    #                           mills^2*2*uStar*d2 - 2*mills*d1*uStar^2 -
    #                           mills*3*uStar^2*d2 - d1*uStar^3 ) )
    
    # intermediate #2
    # return( (1/sigma^2) * ( (-uStar - 6*mills - 2*mills*uStar^2 - uStar^3)*d1 -
    #                           (mills + mills^2*2*uStar + mills*3*uStar^2)*d2 ) )
    
    return( (1/sigma^3) * ( (-uStar - 6*mills - 2*mills*uStar^2 - uStar^3)*termA +
                              mills + (mills^2)*2*uStar + mills*3*uStar^2) )
  }
  
  
  # entry 112
  if ( .entry == 112 ) {
    
    #@TEMPORARY: just use the Deriv fn because my math is wrong
    func = Deriv(F11_func, ".sigma")
    return( func( .mu = mu,
                  .sigma = sigma,
                  .crit = .crit ) )
    
    # # #bm: I have no idea why these are wrong...
    # # intermediate
    # # matches final line below, but is wrong
    # # this term looks just like the one that checked out earlier as deriv of alpha wrt sigma:
    # #  (2*mills*uStar2 + (uStar^2/sigma) )*termA - uStar2*mills
    # return( (-2/sigma^3) - (1/sigma^2)*( (2*mills*uStar2 + (uStar^2/sigma) )*termA - uStar2*mills ) -
    #           (2/sigma^3)*termA )
    
    #wrong
    # return( (-1/sigma^2)*( (2/sigma) + ( 2*mills*uStar2 + (1/sigma)*uStar^2 + (2/sigma) )*termA -
    #                          uStar2*mills ) ) 
  }
  
  
  # entry 122=212
  if ( .entry %in% c(122, 212) ) {
    
    
    #@TEMPORARY: just use the Deriv fn because my math is wrong
    func = Deriv(F12_func, ".sigma")
    return( func( .mu = mu,
                  .sigma = sigma,
                  .crit = .crit ) )
    
    # # handy derivatives
    # d1 = 2*mills *uStar2*termA  # derivative of termA wrt sigma
    # d2 = -uStar2 # derivative of uStar wrt sigma
    # 
    # return( (mills/sigma^2)*(-d1*uStar - mills*d2 - 2*uStar*d2) +
    #           ( mills*(-2/sigma^3) + d1*(1/sigma^2) ) * (3 - mills*uStar - uStar^2) )
  }
  
  
  if ( .entry == 222 ) {
    
    #@TEMPORARY: just use the Deriv fn
    func = Deriv(F22_func, ".sigma")
    return( func( .mu = mu,
                  .sigma = sigma,
                  .crit = .crit ) )
  }
  
  
  #if ( !.entry %in% c(111, 121, 211, 221) ) stop("You haven't taken the derivatives wrt sigma yet.")
  
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



# ~~ Zhou functions ---------------------------------

# Zhou's expected Fisher info
# IMPORTANT: Theirs is for a LEFT-truncated normal, so will only
# coincide with mine when crit = mu (by symmetry)
# and I know they have mistake with variance
expectFisherZhou = function(params, .crit) {
  
  mu = params[1]
  sigma = params[2]
  uStar = (.crit - mu)/sigma
  
  # terms as in their Appendix
  LMills = dnorm(uStar) / (1 - pnorm(uStar))  # for left-truncation, but = Mills for crit = mu case
  ( A = (1/sigma^2) * (1 + uStar*LMills - LMills^2) )
  ( B = (1/sigma^2) * LMills * ( 1 - LMills*uStar + uStar^2) )
  ( C = (1/sigma^2) * ( 2 + LMills*uStar - LMills^2*uStar^2 + uStar^3*LMills ) )
  
  return( matrix( c(A, B, B, C), byrow = TRUE, nrow = 2) )
}


# same but with params separated
expectFisherZhou2 = function(mu,
                             sigma,
                             .crit) {
  
  uStar = (.crit - mu)/sigma
  
  # terms as in their Appendix
  LMills = dnorm(uStar) / (1 - pnorm(uStar))  # for left-truncation, but = Mills for crit = mu case
  A = (1/sigma^2) * (1 + uStar*LMills - LMills^2) 
  B = (1/sigma^2) * LMills * ( 1 - LMills*uStar + uStar^2) 
  C = (1/sigma^2) * ( 2 + LMills*uStar - LMills^2*uStar^2 + uStar^3*LMills )
  
  return( matrix( c(A, B, B, C), byrow = TRUE, nrow = 2) )
}


# ~~~ Zhou expected Fisher entries ------------------

# these are NEGATIVE of the the k_{ij} terms in Godwin
zhou_F11 = function(mu,
                    sigma,
                    .crit) {
  
  uStar = (.crit - mu)/sigma
  
  # terms as in their Appendix
  LMills = dnorm(uStar) / (1 - pnorm(uStar))  # for left-truncation, but = Mills for crit = mu case
  (1/sigma^2) * (1 + uStar*LMills - LMills^2) 
}

zhou_F12 = function(mu,
                    sigma,
                    .crit) {
  
  uStar = (.crit - mu)/sigma
  
  # terms as in their Appendix
  LMills = dnorm(uStar) / (1 - pnorm(uStar))  # for left-truncation, but = Mills for crit = mu case
  (1/sigma^2) * LMills * ( 1 - LMills*uStar + uStar^2) 
}


zhou_F22 = function(mu,
                    sigma,
                    .crit) {
  
  uStar = (.crit - mu)/sigma
  
  # terms as in their Appendix
  LMills = dnorm(uStar) / (1 - pnorm(uStar))  # for left-truncation, but = Mills for crit = mu case
  (1/sigma^2) * ( 2 + LMills*uStar - LMills^2*uStar^2 + uStar^3*LMills )
}



# # second derviatives without any negative sign
# # I checked this fn against expectFisherZhou in expectation
# zhou_second_derivs = function(mu,
#                              sigma, 
#                              .x,
#                              .crit, 
#                              .entry) {
#   
#   uStar = (.crit - mu)/sigma
#   uStar2 = (.crit - mu)/sigma^2
#   
#   # terms as in their Appendix
#   LMills = dnorm(uStar) / (1 - pnorm(uStar))  # for left-truncation, but = Mills for crit = mu case
#   termA = LMills^2 - uStar*LMills
#   
#   # correct in expectation
#   if ( .entry == 11 ) {
#     return( (-1/sigma^2) + (1/sigma^2)*termA )
#   }
#   
#   if ( .entry %in% c(12,21) ) {
#     return( -2*(.x - mu)/sigma^3 + termA*(.crit - mu)/sigma^3 + LMills/sigma^2 )
#   }
#   
#   if ( .entry == 22 ) {
#     return( (1/sigma^2) - (3*(.x - mu)^2 / sigma^4) + termA*uStar2^2 + LMills*2*(.crit - mu)/sigma^3 )
#   }
# } 


# ~~~ Zhou second derivatives -----------------
# these are only used as intermediates before getting the third derivatives numerically
# these have the OPPOSITE sign from the Fisher info because they they don't
#  yet take the negative
# typed in from his expressions in first column of Appendix A
zhou_J11 = function(mu,
                    sigma,
                    .x,
                    .crit) {
  
  uStar = (.crit - mu)/sigma
  uStar2 = (.crit - mu)/sigma^2
  
  # terms as in their Appendix
  LMills = dnorm(uStar) / (1 - pnorm(uStar))  # for left-truncation, but = Mills for crit = mu case
  termA = LMills^2 - uStar*LMills
  
  (-1/sigma^2) + (1/sigma^2)*termA
}


zhou_J12 = function(mu,
                    sigma,
                    .x,
                    .crit) {
  
  uStar = (.crit - mu)/sigma
  uStar2 = (.crit - mu)/sigma^2
  
  # terms as in their Appendix
  LMills = dnorm(uStar) / (1 - pnorm(uStar))  # for left-truncation, but = Mills for crit = mu case
  termA = LMills^2 - uStar*LMills
  
  -2*(.x - mu)/sigma^3 + termA*(.crit - mu)/sigma^3 + LMills/sigma^2
}

zhou_J22 = function(mu,
                    sigma,
                    .x,
                    .crit) {
  
  uStar = (.crit - mu)/sigma
  uStar2 = (.crit - mu)/sigma^2
  
  # terms as in their Appendix
  LMills = dnorm(uStar) / (1 - pnorm(uStar))  # for left-truncation, but = Mills for crit = mu case
  termA = LMills^2 - uStar*LMills
  
  (1/sigma^2) - (3*(.x - mu)^2 / sigma^4) + termA*uStar2^2 + LMills*2*(.crit - mu)/sigma^3
}



# ~~~ Zhou third derivatives -----------------
# expected third derivatives
# gotten numerically by taking derivs of JXX functions
#  and then replacing any stochastic things with their expectations

# these are Godwin's k_{ijl} terms and have the SAME sign as his 
#  because the zhou_JXX functions don't take the negative
zhou_third_derivs = function(mu,
                             sigma, 
                             .crit, 
                             .entry) {
  
  
  uStar = (.crit - mu)/sigma
  uStar2 = (.crit - mu)/sigma^2
  
  # terms as in their Appendix
  LMills = dnorm(uStar) / (1 - pnorm(uStar))  # for left-truncation, but = Mills for crit = mu case
  termA = LMills^2 - uStar*LMills
  
  # Zhou's second moment
  secondMoment = sigma^2*(1 - uStar*LMills) 
  
  
  if ( .entry == 111 ) {

    # this is already the expectation because there's nothing stochastic
    # (by looking at body of thirdDeriv and seeing there's no .x term)
    thirdDeriv = Deriv( zhou_J11, "mu" )
  }
  
  if ( .entry %in% c(112, 121, 211) ) {

    # this is already the expectation because there's nothing stochastic
    # (by looking at body of thirdDeriv and seeing there's no .x term)
    thirdDeriv = Deriv( zhou_J12, "mu" )
    
  }
  
  if ( .entry %in% c(122, 221, 212) ) {
    # from this: Deriv( zhou_J22, "mu" )
    # but making replacements to get expectations
    thirdDeriv = function(mu, sigma, .x, .crit) 
    {
      .e1 <- .crit - mu
      .e2 <- .e1/sigma
      .e3 <- dnorm(.e2)
      .e4 <- 1 - pnorm(.e2)
      .e5 <- .e3/.e4
      .e6 <- .e5 - .e2
      .e7 <- dnorm(.e2, 0, 1)
      # prior to using expectations
      # ((.e1 * (.e1 * (.e3 - (2 * .e5 - .e2) * .e7 * .e6)/sigma - 
      #            2 * (.e3 * .e6))/sigma - 2 * .e3)/.e4 + (6 * (.x - mu) - 
      #                                                       2 * (.e1 * .e7 * .e6/.e4))/sigma)/sigma^3
      # after using expectations
      ((.e1 * (.e1 * (.e3 - (2 * .e5 - .e2) * .e7 * .e6)/sigma - 
                 2 * (.e3 * .e6))/sigma - 2 * .e3)/.e4 + (6 * (0) - 
                                                            2 * (.e1 * .e7 * .e6/.e4))/sigma)/sigma^3
    }
  }
  
  if ( .entry == 222 ) {
    

    # from this: Deriv( zhou_J22, "mu" )
    # but making replacements to get expectations
    thirdDeriv = function(mu, sigma, .x, .crit) 
    {
      .e1 <- .crit - mu
      .e2 <- .e1/sigma
      .e3 <- dnorm(.e2)
      .e4 <- 1 - pnorm(.e2)
      .e5 <- .e3/.e4
      .e6 <- dnorm(.e2, 0, 1)
      .e7 <- .e5 - .e2
      # prior to using expectations
      # ((.e1 * (.e1 * (.e1 * (.e3 - (2 * .e5 - .e2) * .e6 * .e7)/sigma - 
      #                   (2 * .e6 + 4 * .e3) * .e7)/sigma - 6 * .e3)/.e4 + 12 * 
      #     ((.x - mu)^2/sigma))/sigma - 2)/sigma^3
      
      # after using expectations
      ((.e1 * (.e1 * (.e1 * (.e3 - (2 * .e5 - .e2) * .e6 * .e7)/sigma - 
                        (2 * .e6 + 4 * .e3) * .e7)/sigma - 6 * .e3)/.e4 + 12 * 
          (secondMoment/sigma))/sigma - 2)/sigma^3
    }
  }
  
  return( thirdDeriv(mu = mu,
                     sigma = sigma, 
                     .x = x, 
                     .crit = crit) )
  
}


# these are NEGATIVE of the the k_{ij}^l terms in Godwin
# because my zhou_FXX functions already return the NEGATIVE expected Fisher
zhou_expect_fisher_derivs = function(mu,
                                     sigma, 
                                     .crit, 
                                     .entry) {
  
  
  uStar = (.crit - mu)/sigma
  uStar2 = (.crit - mu)/sigma^2
  

  # terms as in their Appendix
  LMills = dnorm(uStar) / (1 - pnorm(uStar))  # for left-truncation, but = Mills for crit = mu case
  termA = LMills^2 - uStar*LMills
  
  if ( .entry == 111 ) {
    thirdDeriv = Deriv( zhou_F11, "mu" ) 
  }
  
  if ( .entry %in% c(112, 121, 211) ) {
    thirdDeriv = Deriv( zhou_F12, "mu" ) 
  }
  
  if ( .entry %in% c(122, 221, 212) ) {
    thirdDeriv = Deriv( zhou_F22, "mu" ) 
  }
  
  if ( .entry == 222 ) {
    thirdDeriv = Deriv( zhou_F22, "sigma" ) 
  }
  
  return( thirdDeriv(mu = mu,
                     sigma = sigma, 
                     .crit = crit) )
  
  
}

# zhou_expect_fisher_derivs( mu = params[1],
#                            sigma = params[2],
#                            .crit = crit,
#                            .entry = 111 )


