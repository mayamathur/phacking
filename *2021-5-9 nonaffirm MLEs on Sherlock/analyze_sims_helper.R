
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


# STRAIGHT FROM MRM:
# summarize performance metrics given a dataset (dat) that is already scenario-aggregated
#  looks for all variables with "Bias" or "Cover" in their names and takes their means
# description: description of the row to make a nice table
# selectVars: "Phat", "Diff" (by default looks for a global variable by this name)
my_summarise = function(dat,
                        description = NA,
                        .selectVars = selectVars,
                        badCoverageCutoff = 0.85,
                        badWidthCutoff = 0.90,
                        averagefn = "mean"
){
  
  # variables whose average should be taken
  meanVars = c( namesWith(pattern = "Bias", dat = dat), 
                namesWith(pattern = "AbsErr", dat = dat),
                namesWith(pattern = "Cover", dat = dat),
                namesWith(pattern = "Width", dat = dat) )
  
  if (.selectVars == "Phat") meanVars = meanVars[ !grepl(pattern = "Diff", x = meanVars) ]
  
  if (.selectVars == "Diff") meanVars = meanVars[ !grepl(pattern = "Phat", x = meanVars) ]
  
  
  if (averagefn == "mean") avgfun = function(x) meanNA(x)
  if (averagefn == "median") avgfun = function(x) medNA(x)
  if (averagefn == "median.pctiles") avgfun = function(x) medNA_pctiles(x)
  
  # make a one-row summary
  tab = dat %>% 
    summarise_at( .vars = meanVars, 
                  function(x) avgfun(x) )
  
  # proportion of SCENARIOS with bad MEAN coverage
  tab = tab %>% add_column(BadPhatCover = mean(dat$CoverPhat < badCoverageCutoff),
                           BadDiffCover = mean(dat$CoverDiff < badCoverageCutoff))
  if (.selectVars == "Phat") tab = tab %>% select(-BadDiffCover)
  if (.selectVars == "Diff") tab = tab %>% select(-BadPhatCover)
  
  # proportion of SCENARIOS with bad median CI width
  tab = tab %>% add_column(BadPhatWidth = mean(dat$PhatCIWidth > badWidthCutoff),
                           BadDiffWidth = mean(dat$DiffCIWidth > badWidthCutoff) )
  if (.selectVars == "Phat") tab = tab %>% select(-BadDiffWidth)
  if (.selectVars == "Diff") tab = tab %>% select(-BadPhatWidth)
  
  # only round selected columns if we have strings
  #  (e.g., median with percentiles)
  if ( averagefn == "median.pctiles" ){
    
    if (.selectVars == "Phat"){
      tab$BadPhatCover = round(tab$BadPhatCover, 2)
      tab$BadPhatWidth = round(tab$BadPhatWidth, 2)
    }
    
    if (.selectVars == "Diff"){
      tab$BadDiffCover = round(tab$BadDiffCover, 2)
      tab$BadDiffWidth = round(tab$BadDiffWidth, 2)
    } 
    
  } else {  # otherwise round all columns
    tab = round( tab, 2 )
  }
  
  tab = tab %>% add_column(n.scens = nrow(dat), .before = 1 )
  
  if ( !is.na(description) ) tab = tab %>% add_column(Scenarios = description, .before = 1)
  return(tab)
}



