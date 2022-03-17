#******************************************

#Project: Heterogeneity in direct replications in psychology and its association with effect size
#Script purpose: Helper functions to prep data for tables and figures
#Code: Anton Olsson Collentine (anton@olssoncollentine.com)
#Additional comments: These functions are intended to be 'sourced' into the 
#                     tables.rmd, figures.rmd, and supplement-tables-figures.rmd files 

#******************************************



#******************************************
#Function to fit a random effects model to all empirical data----
#******************************************
#Each set of input datapoints is meta-analyzed as per the original replication project

est_heterogen_smd_raw <- function(x){
  if(any(x[, "effect_type"] == "Risk difference")){ #without the 'any' we will get warnings because we apply 'if' to a vector
    
    fit <- rma(measure = "RD", ai = outcome_t1, bi = outcome_t2, ci = outcome_c1, di = outcome_c2,  n1i = ntreatment, n2i = ncontrol,  data = x)
    
  } else if(any(x[, "effect_type"] == "Raw mean difference")){
    if(any(x[, "rp"] %in% c("RRR5", "RRR7"))){ #For these us the Knapp and Hartung adjustment of standard errors
      
      fit <- rma(measure = "MD", m1i = outcome_t1, m2i = outcome_c1, sd1i = outcome_t2, sd2i = outcome_c2, n1i = ntreatment, n2i = ncontrol, test ="knha", data = x)
    } else{
      fit <- rma(measure = "MD", m1i = outcome_t1, m2i = outcome_c1, sd1i = outcome_t2, sd2i = outcome_c2, n1i = ntreatment, n2i = ncontrol, data = x)
    }
    
  } else if(any(x[, "rp"] == "ML2")){  
    if(any(x[, "effect_type"] == "r")){
      fit <- rma(yi = effect_size, vi = ml2_ncp_variance,  data = x)
    }else{ #for ML2 cohen's q
      fit <- rma(yi = effect_size, vi = outcome_t1,  data = x) 
    }
    
  }else if(any(x[, "outcomes1_2"] == "mean _ SD")){  
    
    fit <- rma(measure = "SMD", m1i = outcome_t1, m2i = outcome_c1, sd1i = outcome_t2, sd2i = outcome_c2, n1i = ntreatment, n2i = ncontrol, data = x) 
    
  } else if(any(x[, "effect_type"] == "r")){
    if(any(x[, "rp"] == "ML1")){
      fit <- rma(measure = "COR", ri = effect_size, ni = Ntotal,  data = x)
    } else{ #for ML3
      fit <- rma(measure = "UCOR", ri = effect_size, ni = Ntotal,  data = x, vtype = "UB")
      #Note that ML3 specifies measure = "COR" with vtype = "UB", but that this is specified as above in the latest metafor 
    }
    
  } else{ #for all two-group count effects
    
    fit <- rma(measure = "OR2DL", ai = outcome_t1, bi = outcome_t2, ci = outcome_c1, di = outcome_c2, 
               n1i = ntreatment, n2i = ncontrol,  data = x)
    #standardized, ML1 used 'OR2D' which = 'OR2DL' in latest version of metafor
    
  }
  
  hetero <- confint(fit)$random[c(1, 3), ] #Gives us the tau2 and I2 estimates with confidence intervals
  data.frame(eff_size = fit$b[[1]], eff_size_p = fit$pval, #effect size (point estimate) + its p-value 
             s_I2 = hetero[2, 1], s_ci.lb = hetero[2, 2], s_ci.ub = hetero[2, 3], #I2 with Confidenc intervals
             tau2 = hetero[1, 1], tau2_ci.lb = hetero[1, 2], tau2_ci.ub = hetero[1, 3], #tau2 with confidence intervals
             H2trunc = fit$H2, H2 = fit$QE / (fit$k - 1)) #General H2  (H2trunc) and DerSimonian-Laird H2
}

#******************************************
#Functions for Figure 3 (first also used in Supplement C)----
#******************************************

MD_fit <- function(x){ #standardize mean differences in the meta-analysis
  fit <- rma(measure = "SMD", m1i = outcome_t1, m2i = outcome_c1, sd1i = outcome_t2, sd2i = outcome_c2, n1i = ntreatment, n2i = ncontrol, data = x)
  hetero <- confint(fit)$random[c(1, 3), ] 
  data.frame(eff_size = fit$b[[1]], eff_size_p = fit$pval, #effect size (point estimate) 
             s_I2 = hetero[2, 1], s_ci.lb = hetero[2, 2], s_ci.ub = hetero[2, 3],
             tau2 = hetero[1, 1], tau2_ci.lb = hetero[1, 2], tau2_ci.ub = hetero[1, 3],
             H2trunc = fit$H2, H2 = fit$QE / (fit$k - 1))
}


log_or_fit <- function(x){ #NB! if a cell is 0, 1/2 will be added to it. Fit log odds meta-analysis
  fit <- rma(measure = "OR", ai = outcome_t1, bi = outcome_t2, ci = outcome_c1, di = outcome_c2, n1i = ntreatment, n2i = ncontrol,  data = x)
  hetero <- confint(fit)$random[c(1, 3), ] 
  data.frame(eff_size = fit$b[[1]], eff_size_p = fit$pval, #effect size (point estimate) 
             s_I2 = hetero[2, 1], s_ci.lb = hetero[2, 2], s_ci.ub = hetero[2, 3],
             tau2 = hetero[1, 1], tau2_ci.lb = hetero[1, 2], tau2_ci.ub = hetero[1, 3],
             H2trunc = fit$H2, H2 = fit$QE / (fit$k - 1))
}

#labels for figure 3----
fitter <- function(df){  #function to fit correlations for annotating facet plots
  fit_tau <- summary(lm(tau ~ eff_size, data = df))
  fit_I2 <- summary(lm(s_I2 ~ eff_size, data = df))
  fit_H2 <- summary(lm(H2 ~ eff_size, data = df))
  r <- c(sqrt(fit_tau$r.squared), sqrt(fit_I2$r.squared), sqrt(fit_H2$r.squared))
  data.frame(r = format(round(r, digits = 2), nsmall = 2),
             index = c("tau", "I2", "H2")) 
}


#functions to estimate bootstrapped (pearson's) confidence intervals to add to plot correlations
pear_tau <- function(formula, data, indices){ #pearson
  d <- data[indices,]
  fit <- cor.test(formula = ~ tau + eff_size, data=d)
  return(fit$estimate)
}
pear_I2 <- function(formula, data, indices){ #pearson
  d <- data[indices,]
  fit <- cor.test(formula = ~ s_I2 + eff_size, data=d)
  return(fit$estimate)
}
pear_H2 <- function(formula, data, indices){ #pearson
  d <- data[indices,]
  fit <- cor.test(formula = ~ H2 + eff_size, data=d)
  return(fit$estimate)
}

bootfitter <- function(x){ #function that computes the bootstrapped confidence intervals using above functions to annotate plots
  bootfit <- boot(data=x, statistic=pear_tau, R=1000, formula= ~ eff_size + tau)
  tau_ci_pears <- boot.ci(bootfit, type="perc")
  bootfit <- boot(data=x, statistic=pear_I2, R=1000, formula= ~ eff_size + s_I2)
  I2_ci_pears <- boot.ci(bootfit, type="perc")
  bootfit <- boot(data=x, statistic=pear_H2, R=1000, formula= ~ eff_size + H2)
  H2_ci_pears <- boot.ci(bootfit, type="perc")
  
  ci.lb <- c(tau_ci_pears$percent[4], I2_ci_pears$percent[4], H2_ci_pears$percent[4])
  ci.ub <- c(tau_ci_pears$percent[5], I2_ci_pears$percent[5], H2_ci_pears$percent[5])
  
  ci <- paste0("[", format(round(ci.lb, 2), nsmall = 2), ", ",
               format(round(ci.ub, 2), nsmall = 2), "]")#percentile 
  data.frame(ci, index = c("tau", "I2", "H2"))
}


# spearman bootstrap CIs----
#Functions to compute confidence intervals for the spearman correlations reported

#First for I2
spear_I2 <- function(formula, data, indices){ #spearman
  d <- data[indices,]
  fit <- cor.test(formula = ~ s_I2 + eff_size, data=d, method = "spearman")
  return(fit$estimate)
}

#Second for tau
spear_tau <- function(formula, data, indices){ #spearman
  d <- data[indices,]
  fit <- cor.test(formula = ~ tau + eff_size, data=d, method = "spearman")
  return(fit$estimate)
}

#third for  H2
spear_H2 <- function(formula, data, indices){ #spearman
  d <- data[indices,]
  fit <- cor.test(formula = ~ H2 + eff_size, data=d, method = "spearman")
  return(fit$estimate)
}
#******************************************
#Functions for Supplement A----
#******************************************


#Formulas found in: 
#Borenstein, M. (2009). Effect sizes for continuous data. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), 
#The handbook of research synthesis and meta-analysis (2nd ed., pp. 221-235). New York: Russell Sage Foundation.

#Jacobs, P., & Viechtbauer, W. (2017). Estimation of the biserial correlation and its sampling variance for use
#in meta-analysis. Research Synthesis Methods, 8(2), 161-180.

#Cooper, H., Hedges, L. V., & Valentine, J. C. (2009). The handbook of research synthesis and meta-analysis 2nd
#edition. In The Hand. of Res. Synthesis and Meta-Analysis, 2nd Ed. (pp. 1-615). Russell Sage Foundation.


#For raw mean difference + SD -> Hedge's g
transform_SD <- function(m1, m2, SD1, SD2, n1, n2){ #assumes ES is raw mean difference
  sdpooled <- sqrt(((n1 - 1)*SD1^2 + (n2 - 1)*SD2^2) / (n1 + n2 - 2)) #Borenstein, M. (2009), p. 226. 
  d <- (m1 - m2) / sdpooled
  data.frame(d = d, n1 = n1, n2 = n2)
}

transform_d_to_r <- function(d, n1, n2, fisher = FALSE){ 
  m <- n1 + n2 - 2 #Jacobs and Viechtbauer (2017), p.164
  h <- m/n1 + m/n2 #p. 164
  r <- d / sqrt(d^2 + h) #point-biserial correlation, Jacobs and Viechtbauer (2017), equation 5 
  var_r <- (1-r^2)^2 / (n1 + n2 -1) #Cooper and Hedges (2009), equation 12.27
  if(fisher){ #if fisher = TRUE, do fisher transform
  r <- 0.5*log((1+r)/(1-r)) #Fishers z transformation Cooper and Hedges (2009), equation 12.28
  var_r <- 1 / (n1+n2 -3) #Fisher z variance, Cooper and Hedges (2009), equation 12.29
  }
  data.frame(r = r, vi = var_r, n = n1 + n2) 
}


#Function to apply the transformation functions to the (in the end only to data where x[, "outcomes1_2"] == "mean _ SD")
transform_MA <- function(x, fisher = FALSE){
  
  d_conversion <- transform_SD(x$outcome_t1, x$outcome_c1, x$outcome_t2, x$outcome_c2, x$ntreatment, x$ncontrol)
  transform_d_to_r(d_conversion$d, d_conversion$n1, d_conversion$n2, fisher = fisher)
  
}


#function to run random effect meta-analysis on transformed effect sizes
summarizer <- function(x){
  fitr <- rma(yi = x$r, vi = x$vi, data = x) 
  ci <- confint(fitr)$random[c(1, 3), ] #I2 and tau2 confidence intervals
  
  data.frame(r = fitr$b[[1]], tau2 = fitr$tau2, tau2.lb = ci[1, 2], tau2.ub = ci[1, 3],
             I2 = fitr$I2, I2.lb = ci[2, 2], I2.ub = ci[2, 3], H2trunc = fitr$H2, H2 = fitr$QE / (fitr$k - 1)) 
}

#Default method of metafor for calculating H2 is truncated at one. Alternative method for calculating H2 provides information also on
#excessive homogeneity, i.e, less variability than expected by chance (that is, does not have a lower limit of 1). 
#This method approximates H2 as if we were using the Dersimonian and Laird estimate of tau2, although we use REML. 
#See Higgins & Thompson, 2002 and ?print.rma.uni.

#Higgins, J., & Thompson, S. G. (2002). Quantifying heterogeneity in a meta-analysis. Statistics in medicine, 21(11), 1539-1558.



#******************************************
#Function for Supplement C ----
#******************************************

compute_precision <- function(SD1, SD2, n1, n2){
  var_pooled <- ((n1 - 1)*SD1^2 + (n2 - 1)*SD2^2) / (n1 + n2 - 2) #Borenstein, M. (2009), p. 226. 
  sampling_var <- var_pooled / (n1 + n2)
  1 / sampling_var #precision
}



#Supplement C labels
pear_sigma2 <- function(formula, data, indices){ #pearson
  d <- data[indices,]
  fit <- cor.test(formula = ~ s2 + eff_size, data=d)
  return(fit$estimate)
}

fitter_sigma2 <- function(df){  #function to fit correlations for annotating facet plots
  fit_sigma <- summary(lm(s2 ~ eff_size, data = df))
  r <- sqrt(fit_sigma$r.squared)
  data.frame(r = format(round(r, digits = 2), nsmall = 2),
             index = "sigma") 
}

bootfitter_sigma2 <- function(x){ #function that computes the bootstrapped confidence intervals using above functions to annotate plots
  bootfit <- boot(data=x, statistic=pear_sigma2, R=1000, formula= ~ eff_size + s2)
  sigma_ci_pears <- boot.ci(bootfit, type="perc")
  
  ci.lb <- sigma_ci_pears$percent[4]
  ci.ub <- sigma_ci_pears$percent[5]
  
  ci <- paste0("[", format(round(ci.lb, 2), nsmall = 2), ", ",
               format(round(ci.ub, 2), nsmall = 2), "]")#percentile 
  data.frame(ci, index = "sigma")
}
