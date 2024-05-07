library(disdat)
library(np)
library(AUC)
library(rvinecopulib)
library(pROC)
library(lemon)
library(pbsdm)
library(tidyverse)
library(sf)
library(DescTools)
library(foreach)
library(doParallel)
source("R/evaluate_disdat.R")
source("R/rarified_eval_disdat.R")

# Load the full model output

  stop("this needs updating after models finish running")

model_vector = c("maxnet",
                 "rangebagging/none",
                 "kde/kde",
                 "gaussian/gaussian",
                 "gaussian/none",
                 "kde/none",
                 "rulsif",
                 "ulsif",
                 "gaussian/kde",
                 "vine / vine",
                 "kde / vine"
                 #,"lobagoc/none" # need debugging AUC calc
                 #,"vine/none"    # need debugging AUC calc
                 ) #need to add more selected models to this


#Lower priority stuff to add                     
  # "rangebagging / rangebagging"
  # "gaussian / rangebagging"
  # "kde / rangebagging"
  # "kde / gaussian"
  # "gaussian / vine"
  # "vine / rangebagging"
  # "vine / gaussian"
  # "vine / kde"                 


rarified_eval_disdat(presence_vector = (1:10)^2,
                     n_reps = 3,
                     model_vector,
                     quantile = 0.05,
                     temp_full_RDS = "outputs/temp_rarified_full.RDS",
                     temp_fold_RDS = "outputs/temp_rarified_fold.RDS",
                     verbose = TRUE,
                     ncl = 5,
                     seed = 2005)


# Setting levels: control = 0, case = 1
# Setting direction: controls < cases
