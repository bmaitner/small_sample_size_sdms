library(questionr)
library(tidyverse)
library(renv)

packages_used <- questionr::qscan(load = FALSE) %>%
  unlist() %>%
  as.vector() %>%
  unique()

data.frame(packages = packages_used)

#renv::init()
renv::snapshot()
