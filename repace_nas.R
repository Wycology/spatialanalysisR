# Replacing NA values with median values in every column of a tibble object in R

library(dplyr)      # Version 1.1.4
library(tidyr)      # Version 1.3.0
library(tibble)     # Version 3.2.1
library(foreach)    # Version 1.5.2
library(doParallel) # Version 1.0.17

pkgs  <- c("dplyr", "tidyr", "tibble", "foreach", "doParallel")

foreach(i = pkgs, .combine = c) %do% {packageVersion(i)}

my_tibble <- tibble(column1 = c(4, 9, 10, NA, 5, 12, NA, 7, 11, 8), 
                    column2 = c(2, 1, 4, 5, 6, 7, 4, 3, 8, 12), 
                    column3 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                    column4 = c(NA, NA, NA, NA, NA, NA, NA, NA, 11, 8), 
                    column5 = c(NA, NA, NA, NA, NA, NA, NA, NA, 11, 8), 
                    column6 = c(NA, 9, 10, NA, 5, 12, NA, 7, 11, 8),
                    column7 = c(NA, NA, -4, NA, NA, NA, NA, NA, 0, -9999))

my_tibble |> mutate(across(where(is.numeric), ~ replace_na(., median(., na.rm = TRUE))))
