

# PRELIMINARIES ---------------------------------------------

#@SPECIFIC TO HAGGER
root.dir = setwd("~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/Applied examples/Lodder")


# ~~ Load Data and Packages -----------------------------------------------


toLoad = c("crayon",
           "dplyr",
           "foreach",
           "doParallel",
           "boot",
           "metafor",
           "robumeta",
           "data.table",
           "purrr",
           "metRology",
           "fansi",
           "MetaUtility",
           "ICC",
           "cfdecomp",
           "tidyr",
           "truncdist",
           "tibble",
           "tmvtnorm",
           "testthat",
           "truncreg",
           "truncnorm",
           "rstan",
           "optimx",
           "weightr",
           "here")

lapply( toLoad,
        require,
        character.only = TRUE)


# helper fns
code.dir = here("Sherlock code")
setwd(code.dir)
source("helper_SAPH.R")


# dataset
setwd(root.dir)
# for codebook, see the second sheet of xlsx file
# from their script b/c format is weird:
dp <-read.xlsx("MoneyPrimingMetaAnalysis_Dataset_FINAL_correction2.xlsx",sheetName="values",keepformulas=FALSE,startRow=2,header=T)

# from their analysis script:
# Only select studies from dataset that meet inclusion criteria
#d<-subset(d,d$Included.==1)
expect_equal( all( as.numeric(dp$Included.) == 1 ), TRUE ) 

# PREP DATASET ------------------------------

# Add standard error and p-value of Hedges' g
metdat<-cbind(metdat,gse=sqrt(metdat$var.of.g))
metdat<-cbind(metdat,gp=round(1-pnorm(metdat$g/metdat$gse),4))

dp$yi = dp$g
dp$vi = dp$var.of.g
dp$sei = sqrt(dp$vi)
dp$Zi = dp$yi/dp$sei
dp$pval.two = 2 * ( 1 - pnorm( abs(dp$Zi) ) ) 
dp$affirm = (dp$pval.two <= 0.05) & (dp$yi > 0)
expect_equal( dp$affirm, dp$Zi > qnorm(0.975) ) 

# explore
#@TOTAL COUNT OF STUDIES AND OF PREREG STUDIES DOESN'T AGREE WITH THEIR FIGURE 1
# RETURN TO THIS LATER!
# IS THIS RELATED TO THE CORRECTIONS? I DON'T THINK THAT'S ENOUGH TO EXPLAIN THE DIFFERENCE
# E.G., I HAVE TOTAL K = 287, but they report k = 246 in Fig 1
table(dp$Preregistered)
table(dp$Preregistered, dp$replication)
table(dp$dep)  # DVs

# exclude missing data because estimate_jeffreys_RTMA, etc., aren't designed to handle that
dp = dp %>% filter( !is.na(yi) & !is.na(vi) )



setwd(results.dir)
fwrite(dp, "lodder_prepped.csv")




