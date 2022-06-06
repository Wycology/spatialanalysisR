library(dplyr)
library(tidyr)

column1 <- c(4, 9, 10, NA, 5, 12, NA, 7, 11, 8)
column2 <- c(4, NA, 10, 7, 5, 12, 12, 7, NA, 8)
column3 <- c(4, 9, NA, NA, 5, 12, 6, 7, 11, 8)
column4 <- c(NA, 9, 10, 12, 5, 12, NA, 7, 11, 8)
column5 <- c(NA, NA, NA, NA, NA, NA, NA, NA, 0, -9999)

my_df <- data.frame(column1, column2, column3, column4, column5)

my_df


my_df |> mutate(across(where(is.numeric), ~ replace_na(., median(., na.rm = TRUE))))
