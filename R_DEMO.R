library(raster)
library(dplyr)
library(sdm)

set.seed(2014) # For reproducibility reasons
r <- raster(nrow = 10, ncol = 10) # Creating raster of 100 cells
r[] <- rnorm(1:ncell(r)) # Filling the raster with random values

names(r) <- "min_temp" # Re-naming the raster layer
df <- data.frame(r = as.data.frame(r), coordinates(r))|> dplyr::arrange(min_temp) # Dataframe the raster 
 
df <- df |> mutate(pa = case_when(min_temp <= -1 ~ 1,
                                  min_temp > -1 ~ 0,
                                  TRUE ~ 2))
head(df)

d <- sdmData(pa ~ min_temp, train = df)
m <- sdm(pa ~ min_temp, data = d, methods = 'rf', replica = 'boot', n = 4)

coordinates(df) <- ~x+y
gridded(df) <- TRUE
rea <- raster(df) 
ens <- predict(m, rea, overwrite = T)
plot(rea)
plot(ens)