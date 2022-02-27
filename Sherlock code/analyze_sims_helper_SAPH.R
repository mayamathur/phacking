
# ~ FNS FOR SUMMARIZING SIMULATION DATA ----------------------------------------------


# taken from TNE 2022-2-27
# fn for aggregating so we can look at different
#  iterate-level filtering rules
# .s: the iterate-level stitched data (not yet aggregated in any way)
# averagefn: fn to use when aggregating results across scenarios
# expected.sim.reps: only used for sanity checks
make_agg_data = function( .s,
                          .averagefn = "median",
                          badCoverageCutoff = 0.85,
                          expected.sim.reps = NA ){
  
  # # TEST ONLY
  # .s = s
  # .averagefn = "median",
  # badCoverageCutoff = 0.85,
  # expected.sim.reps = NA
  
  
  # make unique scenario variable, defined as scen.name AND method
  if ( !"unique.scen" %in% names(.s) ) .s$unique.scen = paste(.s$scen.name, .s$method)
  
  ##### Outcome and Parameter Variables #####
  # "outcome" variables used in analysis
  analysis.vars = c( 
    "Mhat",
    "Shat",
    
    "MLo",
    "MHi",
    
    "SLo",
    "SHi",
    
    
    ##### variables to be created in mutate below:
    
    "MhatBias",
    "ShatBias",
    
    # "MhatRelBias",
    # "VhatRelBias",
    
    "MhatCover",
    "ShatCover",
    
    "MhatWidth",
    "ShatWidth",
    
    "MhatRMSE",
    "ShatRMSE",
    
    "MhatEstSE",
    "ShatEstSE",
    
    "MhatEmpSE",
    "ShatEmpSE",
    
    # diagnostics regarding point estimation and CIs
    "MhatEstFail",
    "MhatCIFail",
    "ShatEstFail",
    "ShatCIFail",
    
    # diagnostics for main Mhat optimizer
    "OptimConverged",
    
    # diagnostics for other optimizers
    "OptimxMhatWinner",
    "OptimxPropAgreeMhatWinner",
    "OptimxPropAgreeConvergersMhatWinner",
    
    "OptimxShatWinner",
    "OptimxPropAgreeShatWinner",
    "OptimxPropAgreeConvergersShatWinner", 
    
    # Stan diagnostics, part 2
    "StanWarned",
    "MhatRhatGt1.01",
    "MhatRhatGt1.05",
    "MhatRhatGt1.10",
    "MhatRhatMax",
    
    "ShatRhatGt1.01",
    "ShatRhatGt1.05",
    "ShatRhatGt1.10",
    "ShatRhatMax"
  )
  
  
  
  
  # variables that define the scenarios
  param.vars = c("unique.scen",  
                 "method",
                 
                 "Nmax",
                 "Mu",
                 "t2a",
                 "t2w",
                 "m",
                 "true.sei.expr",
                 "hack",
                 "rho",
                 "k.pub.nonaffirm",
                 "prob.hacked",
                 "stan.adapt_delta",
                 "stan.maxtreedepth")
  
  
  # sanity check to make sure we've listed all param vars
  t = .s %>% group_by_at(param.vars) %>% summarise(n())
  if ( !is.na(expected.sim.reps) ) {
    if ( max(t$`n()`) > expected.sim.reps ) stop("param.vars in make_agg_data might be missing something because grouping that way indicated some scenarios had more than expected.sim.reps")
  }
  
  
  ##### Overwrite Analysis Variables As Their Within-Scenario Means #####
  # organize variables into 3 mutually exclusive sets: 
  # - param.vars: parameter variables for grouping
  # - toDrop: variables to drop completely
  # - firstOnly: variables that are static within a scenario, for which we
  #   should just take the first one (but not param variables)
  # - takeMean: variables for which we should take the mean within scenarios
  
  #names(.s)[ !names(.s) %in% param.vars ]  # look at names of vars that need categorizing
  
  s$V = s$t2a + s$t2w
  s$S = sqrt(s$t2a + s$t2w)
  
  toDrop = c("rep.methods",
             "get.CIs",
             "error",
             "rep.name",
             "doParallel.seconds",
             "optim.converged",
             "stan.warned",
             names_with(.dat = .s, .pattern = "optimx") )

  firstOnly = c("scen.name",
                "unique.scen",
                "V",  # calculated from scen params
                "S")
  
  ##### Add New Variables Calculated at Scenario Level #####
  
  # if you have 10K iterates, script breaks from here forward if running locally
  # "vector memory limits"
  s2 = .s %>%
    rename(
      # static within scenario
      # just renaming for clarity
      MhatEstSE = MhatSE,
      ShatEstSE = ShatSE ) %>%
    
    # take just first entry of non-parameter variables that are static within scenarios
    group_by_at(param.vars) %>%
    mutate_at( firstOnly, 
               function(x) x[1] ) %>%
    
    # make certain ad hoc variables that don't conform to below rules
    # this step creates variables that are repeated for every rep within 
    #  a scenario, which is intentional
    
    # make variables that are calculated within scenarios
    # some of the vars are defined at the iterate level (i.e., they still vary within scen), 
    #  while others are calculated at the scen level (i.e., they are static within scen)
    # after this step, we take the means within scens of all these vars, which is immaterial
    #   for the ones that are already static within scenario
    group_by_at(param.vars) %>%
    
    mutate( sim.reps.actual = n(),
            
            # varies within scenario
            MhatBias = Mhat - Mu,
            ShatBias = Shat - S,
            
            # varies within scenario
            MhatCover = covers(truth = Mu, lo = MLo, hi = MHi),
            ShatCover = covers(truth = S, lo = SLo, hi = SHi),
            
            # varies within scenario
            MhatWidth = MHi - MLo,
            ShatWidth = SHi - SLo,
            
            # static within scenario
            MhatRMSE = sqrt( meanNA( (Mhat - Mu)^2 ) ),
            ShatRMSE = sqrt( meanNA( ( Shat - S )^2 ) ),
            
            # static within scenario
            MhatEstFail = mean(is.na(Mhat)),
            MhatCIFail = mean(is.na(MLo)),
            ShatEstFail = mean(is.na(Shat)),
            ShatCIFail = mean(is.na(SLo)),
            
            # static within scenario
            MhatEmpSE = sd(Mhat, na.rm = TRUE),
            ShatEmpSE = sd(Shat, na.rm = TRUE),
            
    
            # varies within scenario
            # how much smaller is estimated SE compared to empirical one?
            MhatSEBias = MhatEstSE - MhatEmpSE,
            ShatSEBias = ShatEstSE - ShatEmpSE,
            
            # varies within scenario
            MhatSERelBias = (MhatEstSE - MhatEmpSE) / MhatEmpSE, 
            ShatSERelBias = (ShatEstSE - ShatEmpSE) / ShatEmpSE,
            
            # static within scenario
            OptimConverged = meanNA(optim.converged),
            
            OptimxMhatWinner = meanNA(optimx.Mhat.winner),
            OptimxPropAgreeMhatWinner = meanNA(optimx.Pagree.Mhat.winner),
            OptimxPropAgreeConvergersMhatWinner = meanNA(optimx.Pagree.of.convergers.Mhat.winner),
            
            OptimxShatWinner = meanNA(optimx.Shat.winner),
            OptimxPropAgreeShatWinner = meanNA(optimx.Pagree.Shat.winner),
            OptimxPropAgreeConvergersShatWinner = meanNA(optimx.Pagree.of.convergers.Shat.winner),
            
            # # proportions of reps within scenarios for which 2 optimizers disagreed by 
            # #  more than various thresholds:
            # MhatOptimDisagree0.001 = meanNA( abs(Mhat.opt.diff) > 0.001 ),
            # MhatOptimDisagree0.01 = meanNA( abs(Mhat.opt.diff) > 0.01 ),
            # MhatOptimDisagree0.05 = meanNA( abs(Mhat.opt.diff) > 0.05 ),
            
            # static within scenario
            StanWarned = meanNA(stan.warned),
            MhatRhatGt1.01 = meanNA(MhatRhat > 1.01),
            MhatRhatGt1.05 = meanNA(MhatRhat > 1.05),
            MhatRhatGt1.10 = meanNA(MhatRhat > 1.10),
            MhatRhatMax = max(MhatRhat, na.rm = TRUE),
            
            ShatRhatGt1.01 = meanNA(ShatRhat > 1.01),
            ShatRhatGt1.05 = meanNA(ShatRhat > 1.05),
            ShatRhatGt1.10 = meanNA(ShatRhat > 1.10),
            ShatRhatMax = max(ShatRhat, na.rm = TRUE)
    ) 
  
  
  # now look for which variables should have their means taken
  # this step must happen here, after we've started making s2, 
  #  so that the takeMean vars are actually in s2
  ( takeMean = names(s2)[ !names(s2) %in% c(param.vars, toDrop, firstOnly) ] )
  # sanity check: have all variables been sorted into these categories?
  expect_equal( TRUE,
                all( names(s2) %in% c(param.vars, toDrop, firstOnly, takeMean) ) )
  
  
  ##### Aggregate to Scenario Level #####
  
  # calculate scenario-level averages, but keep dataset at the rep level
  #  for now to facilitate sanity checks
  # IMPORTANT: this uses meanNA regardless of the passed .avgfun 
  #  because right now we are only aggregating WITHIN scens
  #  so we should never use median 
  
  # don't try to drop vars that don't exist
  toDrop = toDrop[ toDrop %in% names(s2) ]
  
  s3 = s2 %>%
    # take averages of numeric variables
    group_by_at(param.vars) %>%
    mutate_at( takeMean,
               function(x) meanNA(x) ) %>%
    select( -all_of(toDrop) )
  
  
  
  # sanity check: name mismatches
  if ( length( analysis.vars[ !analysis.vars %in% names(s2) ] ) > 0 ) {
    stop("Might have name mismatches; edit analysis.vars in make_agg_data")
  }
  
  # sanity check: SDs of all analysis variables should be 0 within unique scenarios
  t = data.frame( s3 %>% group_by(unique.scen) %>%
                    summarise_at( analysis.vars, sd ) )
  
  t = t %>% select(-unique.scen)
  expect_equal( FALSE,
                any( !as.matrix( t[, 2:(ncol(t)) ] ) %in% c(0, NA, NaN) ) )
  # end sanity checks
  
  
  ##### Aggregate Data at Scenario Level #####
  # make aggregated data by keeping only first row for each 
  #  combination of scenario name and calib.method
  agg = s3[ !duplicated(s3$unique.scen), ]
  
  ##### create Variables That Are Defined At Scenario Rather Than Iterate Level #####
  agg = agg %>% mutate( BadMhatCover = MhatCover < badCoverageCutoff,
                        BadShatCover = ShatCover < badCoverageCutoff )
  
  # # absolute bias is now just the absolute value of bias
  # #  and this is only used for the regressions in Supplement
  # agg$PhatAbsBias = abs(agg$PhatBias)
  # agg$DiffAbsBias = abs(agg$DiffBias)
  
  ##### Make New Variables At Scenario Level ##### 
  # label methods more intelligently for use in plots
  agg$method.pretty.est = NA
  agg$method.pretty.est[ agg$method == c("naive") ] = "Naive"
  agg$method.pretty.est[ agg$method == c("gold-std") ] = "Gold standard"
  agg$method.pretty.est[ agg$method == c("maon") ] = "MAON"
  agg$method.pretty.est[ agg$method == c("2PSM") ] = "2PSM"
  agg$method.pretty.est[ agg$method %in% c("jeffreys-sd") ] = "Jeffreys-SD mode"
  agg$method.pretty.est[ agg$method %in% c("jeffreys-var") ] = "Jeffreys-var mode"
  agg$method.pretty.est[ agg$method %in% c("jeffreys-mcmc-pmed") ] = "Jeffreys-SD median"
  agg$method.pretty.est[ agg$method %in% c("jeffreys-mcmc-pmean") ] = "Jeffreys-SD mean"
  agg$method.pretty.est[ agg$method %in% c("jeffreys-mcmc-max-lp-iterate") ] = "Jeffreys-SD maxLP"
  agg$method.pretty.est[ agg$method == c("mle-sd") ] = "MLE-sd"
  agg$method.pretty.est[ agg$method == c("mle-var") ] = "MLE-var"
  table(agg$method, agg$method.pretty.est)
  
  
  agg$method.pretty.inf = NA
  agg$method.pretty.inf[ agg$method == c("naive") ] = "Naive"
  agg$method.pretty.inf[ agg$method == c("gold-std") ] = "Gold standard"
  agg$method.pretty.inf[ agg$method == c("maon") ] = "MAON"
  agg$method.pretty.inf[ agg$method == c("2PSM") ] = "2PSM"
  agg$method.pretty.inf[ agg$method %in% c("jeffreys-sd") ] = "Jeffreys-SD mode Wald"
  agg$method.pretty.inf[ agg$method %in% c("jeffreys-var") ] = "Jeffreys-var mode Wald"
  agg$method.pretty.inf[ agg$method %in% c("jeffreys-mcmc-pmed", "jeffreys-mcmc-pmean", "jeffreys-mcmc-max-lp-iterate") ] = "Jeffreys posterior quantiles"
  agg$method.pretty.inf[ agg$method == c("mle-sd") ] = "MLE-sd Wald"
  agg$method.pretty.inf[ agg$method == c("mle-var") ] = "MLE-var Wald"

  table(agg$method, agg$method.pretty.inf)

  return(agg %>% ungroup() )
}



# GENERIC SMALL HELPERS -------------------------------------------------------------

# quick mean with NAs removed
meanNA = function(x){
  mean(x, na.rm = TRUE)
}

# quick median with NAs removed
medNA = function(x){
  median(x, na.rm = TRUE)
}


# quick median with NAs removed and 10th and 90th percentiles
medNA_pctiles = function(x){
  paste( round( median(x, na.rm = TRUE), 2 ),
         " (",
         round( quantile(x, probs = 0.10, na.rm = TRUE), 2 ),
         ", ",
         round( quantile(x, probs = 0.90, na.rm = TRUE), 2 ),
         ")",
         sep = "" )
}



names_with = function(.dat, .pattern) {
  names(.dat)[ grepl(pattern = .pattern, x = names(.dat) ) ]
}


# # STRAIGHT FROM MRM:
# # summarize performance metrics given a dataset (dat) that is already scenario-aggregated
# #  looks for all variables with "Bias" or "Cover" in their names and takes their means
# # description: description of the row to make a nice table
# # selectVars: "Phat", "Diff" (by default looks for a global variable by this name)
# my_summarise = function(dat,
#                         description = NA,
#                         .selectVars = selectVars,
#                         badCoverageCutoff = 0.85,
#                         badWidthCutoff = 0.90,
#                         averagefn = "mean"
# ){
#   
#   # variables whose average should be taken
#   meanVars = c( namesWith(pattern = "Bias", dat = dat), 
#                 namesWith(pattern = "AbsErr", dat = dat),
#                 namesWith(pattern = "Cover", dat = dat),
#                 namesWith(pattern = "Width", dat = dat) )
#   
#   if (.selectVars == "Phat") meanVars = meanVars[ !grepl(pattern = "Diff", x = meanVars) ]
#   
#   if (.selectVars == "Diff") meanVars = meanVars[ !grepl(pattern = "Phat", x = meanVars) ]
#   
#   
#   if (averagefn == "mean") avgfun = function(x) meanNA(x)
#   if (averagefn == "median") avgfun = function(x) medNA(x)
#   if (averagefn == "median.pctiles") avgfun = function(x) medNA_pctiles(x)
#   
#   # make a one-row summary
#   tab = dat %>% 
#     summarise_at( .vars = meanVars, 
#                   function(x) avgfun(x) )
#   
#   # proportion of SCENARIOS with bad MEAN coverage
#   tab = tab %>% add_column(BadPhatCover = mean(dat$CoverPhat < badCoverageCutoff),
#                            BadDiffCover = mean(dat$CoverDiff < badCoverageCutoff))
#   if (.selectVars == "Phat") tab = tab %>% select(-BadDiffCover)
#   if (.selectVars == "Diff") tab = tab %>% select(-BadPhatCover)
#   
#   # proportion of SCENARIOS with bad median CI width
#   tab = tab %>% add_column(BadPhatWidth = mean(dat$PhatCIWidth > badWidthCutoff),
#                            BadDiffWidth = mean(dat$DiffCIWidth > badWidthCutoff) )
#   if (.selectVars == "Phat") tab = tab %>% select(-BadDiffWidth)
#   if (.selectVars == "Diff") tab = tab %>% select(-BadPhatWidth)
#   
#   # only round selected columns if we have strings
#   #  (e.g., median with percentiles)
#   if ( averagefn == "median.pctiles" ){
#     
#     if (.selectVars == "Phat"){
#       tab$BadPhatCover = round(tab$BadPhatCover, 2)
#       tab$BadPhatWidth = round(tab$BadPhatWidth, 2)
#     }
#     
#     if (.selectVars == "Diff"){
#       tab$BadDiffCover = round(tab$BadDiffCover, 2)
#       tab$BadDiffWidth = round(tab$BadDiffWidth, 2)
#     } 
#     
#   } else {  # otherwise round all columns
#     tab = round( tab, 2 )
#   }
#   
#   tab = tab %>% add_column(n.scens = nrow(dat), .before = 1 )
#   
#   if ( !is.na(description) ) tab = tab %>% add_column(Scenarios = description, .before = 1)
#   return(tab)
# }

 

