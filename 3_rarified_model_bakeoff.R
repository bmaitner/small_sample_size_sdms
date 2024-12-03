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

  stop("add code to remove any spaces in the model vector. E.g., 'kde / kde' should be converted to 'kde/kde'")
  stop("alternatively, just throw an error telling me to remove the spaces")


model_vector = c("maxnet",
                 "rangebagging/none",
                 "kde/kde",
                 "gaussian/gaussian",
                 "gaussian/none",
                 "kde/none",
                 "rulsif",
                 "ulsif",
                 "gaussian/kde",
                 "gaussian/vine",
                 "gaussian/rangebagging",
                 "vine/vine",
                 "kde/vine",
                 "kde/gaussian"
                 #"lobagoc/none" # need debugging AUC calc
                 ,"vine/none"    # need debugging AUC calc
                 ) #need to add more selected models to this


#Lower priority stuff to add                     
  # "rangebagging/rangebagging"
  # "kde/rangebagging"
  # "vine/rangebagging"
  # "vine/gaussian"
  # "vine/kde"


if(file.exists("outputs/temp_rarified_full.RDS")){
  
  file.copy(from = "outputs/temp_rarified_full.RDS",
            to = "outputs/temp_rarified_full_backup.RDS",
            overwrite = TRUE)
}

if(file.exists("outputs/temp_rarified_fold.RDS")){
  
  file.copy(from = "outputs/temp_rarified_fold.RDS",
            to = "outputs/temp_rarified_fold_backup.RDS",
            overwrite = TRUE)
}

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
