

# IMPORTANT NOTES -----------------------------

# This script goes with 2022-3-11_applied_doParallel_SAPH.R.
#  That script is in the Sherlock folder b/c it's designed to be run via sbatch file. 

# SAPB-E documentation for these datasets:
# codebook_for_b2_data_prepped_step2.csv

# PRELIMINARIES ----------------------------------------

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
           "weightr")
lapply( toLoad,
        require,
        character.only = TRUE)


# helper fns
code.dir = here("Sherlock code")
setwd(code.dir)
source("helper_SAPH.R")



# ~ PREP DATA -----------------------------

# these datasets are directly copied from this dir:
# "~/Dropbox/Personal computer/Independent studies/Sensitivity analysis for publication bias (SAPB)/Linked to OSF (SAPB)/Empirical benchmarks/Data collection/Prepped data for analysis"
setwd("~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/Applied examples/2022-3-11 prep SAPBE dataset for SAPH/Datasets from SAPBE")
b2 = fread("**b2_data_prepped_step2.csv")
f2 = fread("**f2_data_aggregated_step2.csv")

# exclude duplicated Metalab sensitivity estimates
unique( b2$meta.name[ grepl(pattern = "[.]sens", b2$meta.name) == TRUE ] )
b2 = b2 %>% filter( grepl(pattern = "[.]sens", meta.name) == FALSE ) 

expect_equal(63, nuni(b2$meta.name))

# how much clustering?
# average number of studies per cluster
t = b2 %>% group_by(meta.name) %>%
  summarise( n.randomly.chosen = sum(randomly.chosen),
             n.total = n(),
             k.per.cluster = round(n.total/n.randomly.chosen, 3),
             # to see if dataset is public vs. scraped
             Data.source = Data.source[1] )
View(t)

# only analyze ones with relatively little clustering
meta.keepers = t$meta.name[ t$k.per.cluster < 1.2 ]

b2 = b2 %>% filter(meta.name %in% meta.keepers)
f2 = f2 %>% filter(meta.name %in% meta.keepers) 

nuni(b2$meta.name)

#**quick look at these metas' characteristics
t = f2 %>% select(meta.name, group, discipline,
              k.all, k.rc, k.nonaffirm.rc,
              Mhat, Mhat.Lo, Mhat.Hi,
              Mhat.Worst, Mhat.Worst.Lo,
              LogEta) %>%
  mutate_if(is.numeric, function(x) round(x,2))

View(t)

setwd("~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/Applied examples/2022-3-11 prep SAPBE dataset for SAPH/Datasets from SAPBE")

fwrite(t, "quick_summary_sapbe_meta_subset.csv")

# ~ WRITE PREPPED DATA AND ALSO SEND TO CLUSTER -----------------------------

setwd("~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/Applied examples/2022-3-11 prep SAPBE dataset for SAPH/Datasets prepped for SAPH")

fwrite(b2, "b2_long_prepped.csv")
fwrite(f2, "f2_short_prepped.csv")

# send to cluster
scp /Users/mmathur/Dropbox/Personal\ computer/Independent\ studies/2021/Sensitivity\ analysis\ for\ p-hacking\ \(SAPH\)/Code\ \(git\)/Applied\ examples/2022-3-11\ prep\ SAPBE\ dataset\ for\ SAPH/Datasets\ prepped\ for\ SAPH/* mmathur@login.sherlock.stanford.edu:/home/groups/manishad/SAPH/applied_examples/data/sapbe













