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
                 "kde/kde") #need to add more selected models to this


rarified_eval_disdat(presence_vector = (2:10)^2,
                     n_reps = 3,
                     model_vector,
                     quantile = 0.05,
                     temp_full_RDS = "outputs/temp_rarified_full.RDS",
                     temp_fold_RDS = "outputs/temp_rarified_fold.RDS",
                     verbose = TRUE,
                     ncl = 5,
                     seed = 2005)


