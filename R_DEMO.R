library(raster) # version 3.5.21
library(dplyr)  # version 1.0.9
library(sdm)    # version 1.1.8

set.seed(2014) # Setting seed for reproducibility reasons.
r <- raster(nrow = 10, ncol = 10) # Creating raster of 100 (10 rows and 10 columns) cells.
r[] <- rnorm(1:ncell(r)) # Filling the raster cells with random values of length equal to cell numbers.

names(r) <- "min_temp" # Re-naming the raster layer as min_temp
df <- data.frame(r = as.data.frame(r), coordinates(r)) |> dplyr::arrange(min_temp) # Dataframe the raster 
 
df <- df |> mutate(pa = case_when(min_temp <= -1 ~ 1,
                                  min_temp > -1 ~ 0,
                                  TRUE ~ 2)) # Adding classed column to the dataframe
head(df) # Checking first few rows of the df


d <- sdmData(pa ~ min_temp, train = df) # Building sdm data object
m <- sdm(pa ~ min_temp, data = d, methods = 'rf', replica = 'boot', n = 4) # Building sdm model

coordinates(df) <- ~x+y
gridded(df) <- TRUE
rea <- raster(df) 
ens <- predict(m, rea, overwrite = T)
plot(rea)
plot(ens)