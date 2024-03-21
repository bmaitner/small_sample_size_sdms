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

  fold_model_outputs <- readRDS("outputs/bake_off_pnp_fold_model_outputs.RDS")
  full_model_outputs <- readRDS("outputs/bake_off_pnp_full_model_outputs.RDS")

  
  full_model_output_all <- readRDS(file = "outputs/full_model_output_all.RDS")
  fold_model_output_all <- readRDS(file = "outputs/fold_model_output_all.RDS")
  
  
  presences <- unique(full_model_outputs[c("species","n_presence")])

(1:10)^2

stop("Coding still in progress")

model_vector <- unique(full_model_output_all$method)
model_vector <- model_vector[1]


rarified_eval_disdat(presence_vector = (2:10)^2,
                     n_reps = 3,
                     model_vector,
                     quantile = 0.05,
                     temp_full_RDS = "outputs/temp_rarified_full.RDS",
                     temp_fold_RDS = "outputs/temp_rarified_fold.RDS")


rarified_eval_disdat(presence_vector = (2:3)^2,
                     n_reps = 2,
                     model_vector,
                     quantile = 0.05,
                     temp_full_RDS = "outputs/temp_rarified_full.RDS",
                     temp_fold_RDS = "outputs/temp_rarified_fold.RDS")

