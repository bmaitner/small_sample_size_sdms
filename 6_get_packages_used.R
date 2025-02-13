library(questionr)
library(tidyverse)

packages_used <- questionr::qscan(load = FALSE) %>%
  unlist() %>%
  as.vector() %>%
  unique()

data.frame(packages = packages_used)
