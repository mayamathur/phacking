

# PRELIMINARIES ---------------------------------------------

# ~ Load Data and Packages -----------------------------------------------

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
           "here",
           "xlsx")

lapply( toLoad,
        require,
        character.only = TRUE)


# helper fns
general.code.dir = "~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Linked to OSF (SAPH)/Code (git)/Sherlock code"

setwd(general.code.dir)
source("helper_SAPH.R")


# dataset
raw.data.dir = here("Dataset from their repo")
setwd(raw.data.dir)
# codebook: second sheet of xlsx file
# read in per authors' analysis script b/c format is weird:
dp = suppressWarnings( read.xlsx("MoneyPrimingMetaAnalysis_Dataset_FINAL_correction2.xlsx",
               sheetName="values",
               keepformulas=FALSE,
               startRow=2,
               header=T) )

# given the strange format, also save raw data as csv
setwd(raw.data.dir)
fwrite(dp, "lodder_raw.csv")

# from their analysis script:
# Only select studies from dataset that meet inclusion criteria
#d<-subset(d,d$Included.==1)
expect_equal( all( as.numeric(dp$Included.) == 1 ), TRUE ) 


# their old (uncorrected) dataset for sanity checks
dp.old = suppressWarnings( read.xlsx("MoneyPrimingMetaAnalysis_Dataset_FINAL.xlsx",
                                 sheetName="values",
                                 keepformulas=FALSE,
                                 startRow=2,
                                 header=T) )

setwd(raw.data.dir)
fwrite(dp, "lodder_raw_uncorrected.csv")


# SANITY CHECKS FOR NUMBERS OF STUDIES ------------------------------

# Important: The numbers of studies reported in Lodder article
#  are slightly off becuase they corrected their dataset on OSF
#  after publishing the paper, but never published an erratum.
# For this reason, I am doing sanity checks using their OLD
#  (uncorrected) dataset and then running the same code on their
#  new dataset, since they didn't report stats for the corrected dataset. 

# ~ Appendix A: including all interaction rows -----------------
# this is NOT the main analysis
expect_equal( nrow(dp.old), 289 )
expect_equal( sum(dp.old$Preregistered == 1), 51 )


# ~ Fig 1 (main analysis): including only 1 interaction row per study -----------------
# this is the main analysis
# from their R script (section "all studies combined"):
# # Subset of studies based on interaction effects
# metdat.intno<-subset(metdat,metdat$Interaction.==0) # Create subset of studies without interaction effects
# metdat.intyes<-subset(metdat,metdat$Interaction.==1) # Create subset of studies with interaction effects
# metdat.inthigh<-subset(metdat.intyes,metdat.intyes$Interaction.Identification..1.Largest.predicted.effect.==1) # Create subset of studies with largest predicted interaction effect                
# metdat.high<-rbind(metdat.inthigh,metdat.intno) # Combine studies with largest predicted interaction effect with studies without interaction effects

dp.old = dp.old %>% filter( Interaction. == 0 |
                       (Interaction. == 1 &
                          Interaction.Identification..1.Largest.predicted.effect.==1 ) )

# main analysis (Fig 1)
expect_equal( nrow(dp.old), 246 )
expect_equal( sum(dp.old$Preregistered == 1), 47 )
# matches :)


# FILTER DATASET TO REPRODUCE MAIN ANALYSES ------------------------------

# filter the new, corrected dataset in same way as above

dp = dp %>% filter( Interaction. == 0 |
                              (Interaction. == 1 &
                                 Interaction.Identification..1.Largest.predicted.effect.==1 ) )


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

# counts of study types
table(dp$Preregistered)
table(dp$Preregistered, dp$replication)
table(dp$dep)  # DVs

# exclude missing data because estimate_jeffreys_RTMA, etc., aren't designed to handle that
dp = dp %>% filter( !is.na(yi) & !is.na(vi) )



setwd( here("Prepped data") )
fwrite(dp, "lodder_prepped.csv")




