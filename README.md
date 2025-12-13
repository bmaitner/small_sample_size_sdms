# Flexible methods for species distribution modeling with small samples

This dataset contains the code, intermediate data products, and figures associated with the paper "Flexible Methods for Species Distribution Modeling with Small Samples" by Maitner *et al* (2025). This paper focuses on how best to model the distributions of species with small sample sizes (e.g., rare or poorly-sample species), and compares the performances of several different approaches and algorithms.  Data used in this study have been previously published, but intermediate data products and results are included in this repository in order to increase transparency and save any interested parties from having to re-run these time-consuming analyses.


## Description of the data and file structure

The key scripts for the analyses begin with a number and are meant to be run in the corresponding order.

The **R folder** contains R functions written for use in this manuscript, with documentation provided in [Roxygen](https://roxygen2.r-lib.org/) syntax.

The **outputs** folder contains data on model performances. Files beginning with "temp" are temporary files used to resume interrupted workflows. In file names, "dr" refers to density-ratio models while "pnp" refers to plug-and-play models plus environmental-range models and "all" includes density-ratio, plug-and-play models, and environmental range models.  "Full" refers to models fitted with the full set of training data while "fold" refers to data pertaining to k-fold cross-validation of model.  The folder **model_predictions** contains parquet files which describe, for each species and model algorithm combination, which locations were predicted to be within a species' range (1) vs. outside of a species' range (0).

The .RDS files in the outputs folders can be opened in R using the function readRDS.  The parquet files in the model_predictions folder can be opened with the [arrow package](https://cran.r-project.org/package=arrow) for R, using either the function read_parquet to open a single file or open_dataset to access multiple files at once.  

The **figures** and **tables** folders contain figures and tables (respectively) that are generated from the numbered scripts.

## Sharing/Access information

Links to other publicly accessible locations of the data:

 - [Access this dataset on Dryad](https://doi.org/10.5061/dryad.0vt4b8hc2)

 - [Access this dataset on Zenodo](https://doi.org/10.5281/zenodo.14902719)

 - [Access this dataset on Github](https://github.com/bmaitner/small_sample_size_sdms)


Data were derived from the following sources:

 - [Elith, J., Graham, C., Valavi, R., Abegg, M., Bruce, C., Ford, A., Guisan, A., Hijmans, R. J., Huettmann, F., Lohmann, L., & Others. (2020). Presence-only and presence-absence data for comparing species distribution modeling methods. Biodiversity Informatics, 15(2), 69–80. https://doi.org/10.17161/bi.v15i2.13384](https://doi.org/10.17161/bi.v15i2.13384)
 
 - [Maitner, B. S., Boyle, B., Casler, N., Condit, R., Donoghue, J., Durán, S. M., Guaderrama, D., Hinchliff, C. E., Jørgensen, P. M., Kraft, N. J. B., McGill, B., Merow, C., Morueta-Holme, N., Peet, R. K., Sandel, B., Schildhauer, M., Smith, S. A., Svenning, J.-C., Thiers, B., … Enquist, B. J. (2017). The BIEN R package: A tool to access the Botanical Information and Ecology Network (BIEN) database. Methods in Ecology and Evolution / British Ecological Society. https://doi.org/10.1111/2041-210X.12861](https://doi.org/10.1111/2041-210X.12861)
 
 - [Fick, S. E., & Hijmans, R. J. (2017). WorldClim 2: new 1‐km spatial resolution climate surfaces for global land areas: NEW CLIMATE SURFACES FOR GLOBAL LAND AREAS. International Journal of Climatology: A Journal of the Royal Meteorological Society, 37(12), 4302–4315. https://doi.org/10.1002/joc.5086](https://doi.org/10.1002/joc.5086)


## Code/Software

The key scripts for the analyses begin with a number and are meant to be run in the corresponding order. Files 1 and 2 contain analyses which evaluate the performances of different species distribution modelling algorithms. 3 is focused on understanding how well model performances evaluated with fitting or cross-validation data compare with independent, presence-absence performance. 4 evaluates the performance of model ensembles. 5 evaluates how range model predicitons converge as sample size increases and provides an example visualization of an ensemble of three algorithms. 6 creates a figure depicting the counter-intuitive shape of relative occurrence rate funcitons when combining gaussian and range-bagging algorithms. 7 contains code used to track the packages used in analyses. 8 contains code underlying figures that were added at the suggestions of reviewers in the final update of the manuscript.

[![DOI](https://zenodo.org/badge/381449808.svg)](https://doi.org/10.5281/zenodo.14902719)
