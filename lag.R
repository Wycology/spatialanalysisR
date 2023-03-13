######################################
# Task: Calculating spatial lag in R #
# Author: Wyclife Agumba Oluoch      #
# See: https://github.com/Wycology   #
# Last edited: 12th March 2023       #
######################################

library(sf)           # Version 1.0.1
library(sfdep)        # Version 0.2.3
library(tidyverse)    # Version 2.0.0
library(patchwork)    # Version 1.1.2

pkgs <- c("sf", "sfdep", "tidyverse", "patchwork")

for (pkg in pkgs) {
  print(paste0(pkg, " ", packageVersion(pkg)))
}

data(guerry, 
     package = "sfdep")

st_geometry(guerry) %>% 
  plot()

nb <- st_contiguity(guerry)
attributes(nb)

nb[1:3]

wt <- st_weights(nb)
wt[1:3]

x <- guerry$crime_pers

st_lag(x, nb, wt)

ij <- nb[[1]]
wij <- wt[[1]]
xij <- x[ij]
xij * wij
sum(xij * wij)

gg_crime_obs <- ggplot(guerry, 
                       aes(fill = crime_pers)) +
  geom_sf(color = "black", 
          lwd = 0.15) +
  scale_fill_viridis_c(limits = 
                         range(guerry$crime_pers)) +
  theme_void()

crime_lags <- guerry %>% 
  mutate(
    nb = st_contiguity(geometry), 
    wt = st_weights(nb),
    crime_lag = st_lag(crime_pers,
                       nb, 
                       wt)
  )

gg_crime_lag <- ggplot(crime_lags, 
                       aes(fill = crime_lag)) +
  geom_sf(color = "black", 
          lwd = 0.15) +
  scale_fill_viridis_c(limits = range(guerry$crime_pers)) +
  theme_void() 

wrap_plots(gg_crime_obs, 
           gg_crime_lag)

# Kenya data -------------------------------------------------------------------------------------

kenya <- st_read("D:/FROM HP/DATA/DATABASE/Kenya_admin_2014_WGS84.shp")
head(kenya)
# maize <- readxl::read_xlsx("C:/Users/ZEF/Downloads/kenya-maize-production-by-counties.xlsx")

st_geometry(kenya) %>% plot()

nb <- st_contiguity(kenya)
attributes(nb)

nb[1:3]

wt <- st_weights(nb)
wt[1:3]

kenya_bound <- kenya %>% 
  left_join(maize, 
            by = "COUNTY")

x <- kenya_bound$yield_2016

st_lag(x, nb, wt)

ij <- nb[[1]]
wij <- wt[[1]]
xij <- x[ij]
xij * wij
sum(xij * wij)

gg_crime_obs <- ggplot(kenya_bound, 
                       aes(fill = rand)) +
  geom_sf(color = "black", 
          lwd = 0.15) +
  scale_fill_viridis_c(limits = range(kenya_bound$yield_2016)) +
  theme_void()

crime_lags <- kenya_bound %>% 
  mutate(
    nb = st_contiguity(geometry), 
    wt = st_weights(nb),
    crime_lag = st_lag(yield_2016, 
                       nb, 
                       wt)
  )

gg_crime_lag <- ggplot(crime_lags, 
                       aes(fill = crime_lag)) +
  geom_sf(color = "black", 
          lwd = 0.15) +
  scale_fill_viridis_c(limits = range(kenya_bound$yield_2016)) +
  theme_void() 

wrap_plots(gg_crime_obs, 
           gg_crime_lag)

CPCOLS <- c("#010305",
                     "#33A02C", 
                     "#D5E61C")

ggplot(iris, aes(Petal.Length, 
                 Petal.Width)) + 
  geom_point(aes(color = Species), 
             size = 3) + 
  scale_color_manual(values = CPCOLS)