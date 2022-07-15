# Replacing the NAs with Median values per column

library(dplyr) # Version 1.0.9
library(tidyr) # Version 1.2.0
library(tibble) # Version 3.1.7

my_tibble <- tribble(
  ~column1, ~column2, ~column3, ~column4, ~column5,
  4, 9, 10, NA, 5, 12, NA, 7, 11, 8,
  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
  NA, NA, NA, NA, NA, NA, NA, NA, 11, 8,
  NA, 9, 10, 12, 5, 12, NA, 7, 11, 8,
  NA, NA, NA, NA, NA, NA, NA, NA, 0, -9999
)

my_tibble |> mutate(across(where(is.numeric), ~ replace_na(., median(., na.rm = TRUE))))