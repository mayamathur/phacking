
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
  
  
  # make unique scenario variable, defined as scen.name AND method
  if ( !"unique.scen" %in% names(.s) ) .s$unique.scen = paste(.s$scen.name, .s$method)
  
  ##### Outcome and Parameter Variables #####
  # "outcome" variables used in analysis
  analysis.vars = c( 
    "Mhat",
    "Vhat",
    "Shat",
    
    "MLo",
    "MHi",
    
    "VLo",
    "VHi",
    
    "SLo",
    "SHi",
    
    
    ##### variables to be created in mutate below:
    
    "MhatBias",
    "VhatBias",
    "ShatBias",
    
    # "MhatRelBias",
    # "VhatRelBias",
    
    "MhatCover",
    "VhatCover",
    "ShatCover",
    
    "MhatWidth",
    "VhatWidth",
    "ShatWidth",
    
    "MhatRMSE",
    "VhatRMSE",
    "ShatRMSE",
    
    "MhatEstSE",
    "VhatEstSE",
    "ShatEstSE",
    
    "MhatEmpSE",
    "VhatEmpSE",
    "ShatEmpSE",
    
    # diagnostics regarding point estimation and CIs
    "MhatEstFail",
    "MhatCIFail",
    "ShatEstFail",
    "ShatCIFail",
    
    # diagnostics regarding bootstraps
    "BtPropResamplesFail",
    "BtMhatCIFail",
    "BtVhatCIFail",
    "BtShatCIFail",
    
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
                 "boot.reps",
                 "stan.iter",
                 "stan.adapt_delta",
                 "stan.maxtreedepth",
                 "trunc.type",
                 "prop.retained",
                 "mu",
                 "V",
                 "n")
  
  
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
  #   should just take the first one
  # - takeMean: variables for which we should take the mean within scenarios
  
  names(.s)[ !names(.s) %in% param.vars ]  # look at names of vars that need categorizing
  toDrop = c("rep.methods",
             "get.CIs",
             "error",
             "bt.prop.resamples.failed",
             "sim.reps",  # this is the INTENDED sim.reps, so confusing to retain it
             "rep.name",
             "doParallel.seconds",
             "optim.converged",
             "stan.warned",
             names_with(.dat = .s, .pattern = "optimx") )
  
  firstOnly = c("scen.name",
                "unique.scen",
                "Za",  # calculated from theory, so fixed within scen params
                "Zb"
  )
  
  ##### Add New Variables Calculated at Scenario Level #####
  
  # prevent errors for non-bootstrap methods
  if ( ! "bt.prop.resamples.failed" %in% names(.s) ) .s$bt.prop.resamples.failed = NA
  
  # if you have 10K iterates, script breaks from here forward if running locally
  # "vector memory limits"
  s2 = .s %>%
    rename(
      # static within scenario
      # just renaming for clarity
      MhatEstSE = MhatSE,
      VhatEstSE = VhatSE,
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
    
    mutate( sim.reps = n(),
            
            # varies within scenario
            MhatBias = Mhat - mu,
            VhatBias = Vhat - V,
            ShatBias = Shat - sqrt(V),
            
            # varies within scenario
            MhatCover = covers(truth = mu, lo = MLo, hi = MHi),
            VhatCover = covers(truth = V, lo = VLo, hi = VHi),
            ShatCover = covers(truth = sqrt(V), lo = SLo, hi = SHi),
            
            # varies within scenario
            MhatWidth = MHi - MLo,
            VhatWidth = VHi - VLo,
            ShatWidth = SHi - SLo,
            
            # static within scenario
            MhatRMSE = sqrt( meanNA( (Mhat - mu)^2 ) ),
            VhatRMSE = sqrt( meanNA( (Vhat - V)^2 ) ),
            ShatRMSE = sqrt( meanNA( ( Shat - sqrt(V) )^2 ) ),
            
            # static within scenario
            MhatEstFail = mean(is.na(Mhat)),
            MhatCIFail = mean(is.na(MLo)),
            ShatEstFail = mean(is.na(Shat)),
            ShatCIFail = mean(is.na(SLo)),
            
            # static within scenario
            BtPropResamplesFail = mean(bt.prop.resamples.failed),
            BtMhatCIFail = mean( is.na(MLo) ),
            BtVhatCIFail = mean( is.na(VLo) ),
            BtShatCIFail = mean( is.na(SLo) ),
            
            # static within scenario
            MhatEmpSE = sd(Mhat, na.rm = TRUE),
            VhatEmpSE = sd(Vhat, na.rm = TRUE),
            ShatEmpSE = sd(Shat, na.rm = TRUE),
            
            
            
            # varies within scenario
            # how much smaller is estimated SE compared to empirical one?
            MhatSEBias = MhatEstSE - MhatEmpSE,
            VhatSEBias = VhatEstSE - VhatEmpSE,
            ShatSEBias = ShatEstSE - ShatEmpSE,
            
            # varies within scenario
            MhatSERelBias = (MhatEstSE - MhatEmpSE) / MhatEmpSE, 
            VhatSERelBias = (VhatEstSE - VhatEmpSE) / VhatEmpSE,
            ShatSERelBias = (ShatEstSE - ShatEmpSE) / ShatEmpSE,
            
            # static within scenario
            #bm
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
  agg$method.pretty.est[ agg$method %in% c("mle-wald", "mle-profile") ] = "MLE"
  agg$method.pretty.est[ agg$method %in% c("boot-mle-wald") ] = "MLE + boot"
  agg$method.pretty.est[ agg$method %in% c("jeffreys-wald-map", "jeffreys-profile-map") ] = "Jeffreys mode"
  agg$method.pretty.est[ agg$method %in% c("boot-jeffreys-wald-map") ] = "Jeffreys mode + boot"
  agg$method.pretty.est[ agg$method %in% c("jeffreys-mcmc-pmed") ] = "Jeffreys median"
  agg$method.pretty.est[ agg$method %in% c("jeffreys-mcmc-pmean") ] = "Jeffreys mean"
  table(agg$method, agg$method.pretty.est)
  
  agg$method.pretty.inf = NA
  agg$method.pretty.inf[ agg$method %in% c("mle-wald") ] = "MLE Wald"
  agg$method.pretty.inf[ agg$method %in% c("mle-profile") ] = "MLE profile"
  agg$method.pretty.inf[ agg$method %in% c("boot-mle-wald") ] = "MLE BCa"
  agg$method.pretty.inf[ agg$method %in% c("jeffreys-wald-map") ] = "Jeffreys mode Wald"
  agg$method.pretty.inf[ agg$method %in% c("jeffreys-profile-map") ] = "Jeffreys profile"
  agg$method.pretty.inf[ agg$method %in% c("boot-jeffreys-wald-map") ] = "Jeffreys mode BCa"
  agg$method.pretty.inf[ agg$method %in% c("jeffreys-mcmc-pmed", "jeffreys-mcmc-pmean") ] = "Jeffreys posterior quantiles"
  table(agg$method, agg$method.pretty.inf)
  
  agg$`Truncation type` = agg$trunc.type
  agg$`Truncation type`[ agg$`Truncation type` == "single" ] = "Single"
  agg$`Truncation type`[ agg$`Truncation type` == "double-symm" ] = "Symmetric double"
  agg$`Truncation type`[ agg$`Truncation type` == "double-asymm" ] = "Asymmetric double"
  
  agg$MethodUsesBoot = grepl( "boot", agg$method )
  
  
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

 

