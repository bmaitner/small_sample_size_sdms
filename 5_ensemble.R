# The goal of this script is to evaluate the performance of an ensemble of models for small sample size species (which will be compared to the model bakeoff and rarified bakeoff)

library(disdat)
library(np)
library(AUC)
library(rvinecopulib)
library(pROC)
library(lemon)
library(pbsdm)
library(kernlab)
library(tidyverse)
library(sf)
library(DescTools)
library(foreach)
library(doParallel)
source("R/evaluate_disdat_ensemble.R")


# Select models to use in ensemble

  model_vector <- c("maxnet","rulsif","kde/kde")

# Get ensemble performance data

  ensemble_performance <- evaluate_ensemble_disdat(model_vector = model_vector,
                         quantile = 0.05,
                         verbose = TRUE,
                         ncl = 5)
