# Preparation

library(haven)
library(tidyverse)

# Data


ias2009febmay <- read_dta("./data-raw/stata-dataset.dta") |>
  filter(yyyyqq >= 200901 & yyyyqq <= 200902) |>
  select(q2, sex, age, work, educ, tenure) |>
  filter(q2 != 9) |>
  as_factor() |>
  mutate(
    sex  = fct_rev(fct_drop(sex)),
    age  = fct_drop(age),
    work = fct_rev(work),
    q2   = fct_drop(q2)
  ) |>
  na.omit()
usethis::use_data(ias2009febmay, overwrite = TRUE)
