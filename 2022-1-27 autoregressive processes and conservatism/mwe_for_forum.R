

# QUESTION 1 ---------------------------------------------------------------

# https://stats.stackexchange.com/questions/562548/why-empirically-does-strict-stationarity-only-hold-asymptotically-for-an-ar/562564?noredirect=1#comment1035828_562564

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



# QUESTION 2 ---------------------------------------------------------------

# similar to above, but more iterates
library(dplyr)
library(tidyverse)
library(simts)
library(ggplot2)

# number of time series to simulate
k = 5000

# number of draws in each series
draws = 100

# simulate Gaussian AR(1)'s with autocorrelation = 0.9 and errors ~ N(0, 0.1)
for ( i in 1:k ) {
  .d = data.frame( yi = as.numeric( gen_gts( draws, AR1(phi = 0.9, sigma2 = 0.1) ) ) )
  .d$iterate = i
  .d$draw.index = 1:nrow(.d)
  if ( i == 1 ) d = .d else d = bind_rows(d, .d)
}



# avoid any asymptotic issues by discarding warmup draws
d2 = d %>% filter(draw.index >= 50) %>%
  group_by(iterate) %>%
  mutate(max.yi = max(yi))

# classify half of time series as "always small" based on median of max yi
threshold = median(d2$max.yi)
d2$always.small = d2$max.yi < threshold
# vs. whether the CURRENT draw is small
d2$draw.small = d2$yi < threshold

# marginal stationarity
agg = d2 %>% group_by(draw.index) %>%
  summarise( Mean = mean(yi),
             SD = sd(yi) )
plot(agg$draw.index, agg$Mean, type="l")
plot(agg$draw.index, agg$SD, type="l")



# conditional stationarity
agg2 = d2 %>% group_by(always.small, draw.index) %>%
  summarise( Mean = mean(yi),
             SD = sd(yi) )


ggplot( data = agg2,
        aes(x = draw.index,
            y = Mean,
            color = always.small) ) + 
  geom_line() +
  theme_bw()




# individual draws

#** version that compares to individual draws that are small
agg3 = d2 %>% group_by(draw.small, draw.index) %>%
  summarise( Mean = mean(yi),
             SD = sd(yi) )


ggplot( data = agg3,
        aes(x = draw.index,
            y = Mean,
            color = draw.small) ) + 
  
  geom_line() +
  
  theme_bw()



#***both in same plot
#*# important because I might be able to argue that even though E[X_t | ALL nonaffirm] isn't stationary, it's less than E[X_t | X_t nonaffirm] for EVERY value of t

ggplot() + 
  
  geom_line( data = agg2,
             aes(x = draw.index,
                 y = Mean,
                 color = always.small) ) +
  
  # c.f. all individual draws that were always small
  geom_line( data = agg3,
             aes(x = draw.index,
                 y = Mean,
                 color = draw.small,
                 ),
             lty = 2) +
  
  theme_bw()















ggplot( data = d2,
        aes(x = as.factor(draw.index),
            y = yi,
          color = always.small) ) + 
  geom_boxplot() +
  theme_bw()




















