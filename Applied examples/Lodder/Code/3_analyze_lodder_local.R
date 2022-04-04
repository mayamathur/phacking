
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
           "RColorBrewer")

lapply( toLoad,
        require,
        character.only = TRUE)

# helper fns
general.code.dir = "~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Linked to OSF (SAPH)/Code (git)/Sherlock code"
setwd(general.code.dir)
source("helper_SAPH.R")
source("analyze_sims_helper_SAPH.R")

# get prepped data
setwd( here("Prepped data") )
dp = fread("lodder_prepped.csv")

# get results from Sherlock 
#  from 2_analyze_lodder_sherlock.R
setwd( here("Results from Sherlock") )
rs = fread("results_lodder.csv")

# for this script's own results
results.dir = here("Results from local analysis")
overleaf.dir.figs = "/Users/mmathur/Dropbox/Apps/Overleaf/P-hacking (SAPH)/figures_SAPH/lodder"
overleaf.dir.nums = "/Users/mmathur/Dropbox/Apps/Overleaf/P-hacking (SAPH)/results_from_R_SAPH"

# Z-SCORE DENSITY  ------------------------------


# these are WITHIN-study Z-scores
xmin = floor(min(dp$Zi))
xmax = ceiling(max(dp$Zi))

p = ggplot(data = data.frame(x = c(xmin, 3)),
            aes(x)) +
  
  geom_vline(xintercept = 0,
             lwd = 1,
             color = "gray") +
  
  geom_vline(xintercept = qnorm(0.975),
             lty = 2,
             color = "red") +
  
  # estimated density of estimates
  geom_density( data = dp,
                aes(x = Zi),
                adjust = .3 ) +
  
  ylab("") +
  scale_x_continuous( breaks = seq(-2, 8, 1),
                      limits = c(-2, 8)) +
  
  # put the Z=1.96 label slightly off from that location so it's visible with the line
  annotate(geom = "text", x = 1.75, y = 0.08, label = "Z = 1.96", color = "red",
           angle = 90) +
  
  xlab("Z-score") +
  theme_minimal() +
  scale_y_continuous(breaks = NULL) +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=16))



setwd(results.dir)
my_ggsave( name = paste( "lodder_z_density.pdf", sep="_" ),
           .width = 7, 
           .height = 6,
           .results.dir = results.dir,
           .overleaf.dir = overleaf.dir.figs )




# FOREST PLOT OF ESTIMATES BY METHOD ------------------------------

# ~ Recode method ------------------------

rs$method.pretty = NA
rs$method.pretty[ rs$method == c("naive") ] = "Uncorrected"
rs$method.pretty[ rs$method == c("gold-std") ] = "Gold standard"
rs$method.pretty[ rs$method == c("maon") ] = "MAN"
rs$method.pretty[ rs$method == c("2psm") ] = "SM"
rs$method.pretty[ rs$method == c("csm-dataset-2psm") ] = "SMKH" # "known hacking"
rs$method.pretty[ rs$method == c("prereg-naive") ] = "Preregistered only"
rs$method.pretty[ rs$method %in% c("jeffreys-mcmc-pmed") ] = "RTMA"
table(rs$method, rs$method.pretty)




# ~ Plotting dataframe ----------------------
rsp = rs %>% filter( !is.na(method.pretty) &
                       # I accidentally included 2 rows for naive
                       !duplicated(method.pretty) ) %>%
  droplevels()
table(rsp$method.pretty)


# force ordering of methods
correct.order = c("Uncorrected",
                  "SM",
                  "Preregistered only",
                  
                  "RTMA",
                  "MAN",
                  "SMKH")

rsp$method.pretty = factor(rsp$method.pretty, levels = rev(correct.order))
levels(rsp$method.pretty)

# ~ Make plot ----------------------

my.shapes = c(16, 2)

# SAVE: this is to set colors automatically
# # have colors and shapes match sim study
# # taken from inside analyze_sims_helper::sim_plot_multiple_outcomes
# n.colors.needed = length(unique(rsp$method.pretty))
# .colors = brewer.pal(n = n.colors.needed, name = "Dark2")
# if( length(.colors) > n.colors.needed ) .colors = .colors[1:n.colors.needed]
# # this maps the colors onto levels of the factor
# names(.colors) = levels( factor(rsp$method.pretty) )
# 
# # highlight certain methods
# .colors[ names(.colors) == "RTMA" ] = "red"

# hard-code colors to match simulations
.colors = c(SMKH = "#1B9E77",
            MAN = "#ff9900",
            RTMA = "red",
            `Preregistered only` = "#3399ff",
            SM = "#00cc00",
            Uncorrected = "black")


# ~~ Set ggplot linetype scale ----
# by default, dotted lines
# but use solid lines for new proposed methods
.lty = rep("dashed", nuni(rsp$method.pretty))
names(.lty) = names(.colors)

newMethods = c("SMKH",
               "MAN",
               "RTMA")

.lty[ names(.lty) %in% newMethods ] = "solid"



# ~ Make plot ----------------------

p = ggplot( data = rsp,
            aes( x = Mhat,
                y = method.pretty, 
                color = method.pretty,
                shape = (method.pretty == "Uncorrected"),
                linetype = method.pretty ) ) +
  
  geom_vline(xintercept = 0,
             lwd = .8,
             color = "gray") +
  
  geom_point(size = 3) +
  geom_errorbarh( aes(xmax = MHi,
                      xmin = MLo),
                  height = 0,
                  lwd = 0.8 ) +
  
  # manually provided colors
scale_colour_manual(values = .colors ) +
  
  # manually provided linetypes
  scale_linetype_manual(values = .lty,
                        guide = "none") +

  scale_shape_manual(values = my.shapes,
                     guide = "none") +
  
  #coord_cartesian( xlim = c(xmin, xmax) )+
  scale_x_continuous(breaks= seq(-0.05, 0.35, 0.05 ) ) +
  
  xlab( "Pooled estimate with 95% CI" ) +
  ylab("") +
  labs(color  = "Method", linetype = "Method", shape = "Method") +
  # manually set shapes in color legend to match the second legend
  guides(colour = guide_legend( override.aes = list( shape = c(2,16,16,16,16,16) ),
                                reverse = TRUE) ) +
  
  theme_bw() +
  theme( text = element_text(face = "bold"),
         panel.grid.major.y = element_blank(),
         panel.grid.minor.y = element_blank(),
         legend.position = "right" )
  

setwd(results.dir)
my_ggsave( name = paste( "lodder_forest.pdf", sep="_" ),
           .width = 7, 
           .height = 3,
           .results.dir = results.dir,
           .overleaf.dir = overleaf.dir.figs )


# LOOK AT SEI DISTRIBUTION TO INFORM SIMS ------------------------------

# seems like I already have a folder to this effect?
















