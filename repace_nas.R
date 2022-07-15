# Replacing the NAs with Median

library(dplyr)
library(tidyr)

column1 <- c(4, 9, 10, NA, 5, 12, NA, 7, 11, 8)
column2 <- c(rep(NA, 10))
column3 <- c(rep(NA, 8), 11, 8)
column4 <- c(NA, 9, 10, 12, 5, 12, NA, 7, 11, 8)
column5 <- c(NA, NA, NA, NA, NA, NA, NA, NA, 0, -9999)

my_df <- data.frame(column1, column2, column3, column4, column5)

my_df |> mutate(across(where(is.numeric), ~ replace_na(., median(., na.rm = TRUE))))