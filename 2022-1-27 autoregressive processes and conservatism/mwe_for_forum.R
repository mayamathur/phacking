


library(dplyr)
library(tidyverse)
library(simts)

# number of time series to simulate
k = 1000

# number of draws in each series
draws = 100

# simulate Gaussian AR(1)'s with autocorrelation = 0.9
for ( i in 1:k ) {
  .d = data.frame( yi = as.numeric( gen_gts( draws, AR1(phi = 0.9, sigma2 = 0.5) ) ) )
  .d$iterate = i
  .d$draw.index = 1:nrow(.d)
  if ( i == 1 ) d = .d else d = bind_rows(d, .d)
}

agg = d %>% group_by(draw.index) %>%
  summarise( Mean = mean(yi),
             SD = sd(yi) )

plot(agg$draw.index, agg$SD, type="l")

