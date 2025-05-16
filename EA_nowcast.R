# ==============================================================================
#             Dynamic Matrix Factor Models and the EM Algorithm: 
#         A Nowcasting Framework for Mixed Frequency Data in Euro Area
#                               - Main Script - 
# ==============================================================================
# This script estimates a Dynamic Matrix Factor Model (DMFM) by the EM-algorithm 
# and performs a real-time nowcasting exercise for GDP growth in the Euro Area's 
# main economies:
# Germany (DE), France (FR), Italy (IT), and Spain (ES). 
# It includes:
#   1. Data preparation and harmonization across countries
#   2. Construction of a 3D data tensor
#   3. Estimation of the number of latent factors and model lags
#   4. Model initialization and parameter estimation using the EM algorithm on 
#      the whole dataset
#   5. Pseudo real-time rolling nowcasting of GDP
#   6. Saving results for further use
# DMFM and DFM sections are modular and clearly separated for experimentation.
#
# ==============================================================================
# Author      : Davide Delfino
# Institution : Alma Mater Studiorum - University of Bologna
# Dataset     : Barigozzi & Lissona (2024), EA-MD-QD
#
# Disclaimer  : DMFM functions adapted from Barigozzi & Trapin (2024)
#
# ==============================================================================
# Notes       : Other EA countries such as BE, NL, AT, PT, IR, GR and different 
#               set of Monthly and Quarterly variables can be added
# ==============================================================================
# Script Type : Data Preparation / Estimation / Forecasting
# ==============================================================================

# ==============================================================================
# SET WORKING DIRECTORY
# ==============================================================================
path <- "C:/Users/david/Desktop/University/Master's Thesis/Nowcasting_DMFM"
setwd(path)


# ==============================================================================
# LOAD REQUIRED LIBRARIES
# ==============================================================================

# Data Handling
library(tidyverse)
library(lubridate)
library(abind)

# Time Series Analysis
library(tseries)
library(zoo)
library(fBasics)
library(vars)

# Linear Algebra & Utilities
library(Matrix)
library(RSpectra)
library(MASS)
library(pracma)

# I/O
library(writexl)
library(readxl)
library(xtable)

# Plotting
library(ggplot2)
library(reshape2)
library(patchwork)


# ==============================================================================
# SOURCE USER-DEFINED FUNCTIONS
# ==============================================================================

# Data preparation & transformation
source("Functions/R/dmfm.preparation.R")        # ---> ✅

# Utilities
source("Functions/R/matrix.utils.R")            # ---> ✅

# DMFM Estimation Pipeline
source("Functions/R/dmfm.initialization.R")     # ---> ✅
source("Functions/R/kalman.R")                  # ---> ✅
source("Functions/R/dmfm.estimation.R")         # ---> ✅
source("Functions/R/dmfm.results.R")            # ---> ✅

# Visualization utilities
source("Functions/R/dmfm.visualization.R")      # ---> ✅


# ==============================================================================
# DEFINE PARAMETERS
# ==============================================================================

countries <- c("DE", "FR", "IT", "ES")  # Countries to include in the model

# Configuration parameters
P <- list(
  modelM = "Large",                     # Monthly model size
  modelQ = "Small",                     # Quarterly model size
  covid_start = as.Date("2020-01-01"),  # covid_start = as.Date("2020-03-01"),
  covid_end   = as.Date("2021-07-01"),  # covid_end   = as.Date("2020-12-01"),
  startEV     = as.Date("2017-01-01"),
  endEV       = as.Date("2025-01-01"),
  Tmax        = 300                     # Max time dimension
)

kmax <- c(3, 7)  # Initial upper bounds for number of factors


# ==============================================================================
# STEP 1: DATA PREPARATION AND ALIGNMENT
# ==============================================================================

# Step 1.1: Select base variable names most correlated with GDP per country
base_names <- select_base_variable_names(countries, P)

# Step 1.2: Identify common variables across all countries (monthly & quarterly)
common_vars <- filter_common_variables_across_countries(
  countries, base_names$monthly, base_names$quarterly
)

# Step 1.3: Extract base names (same across countries)
all_vars <- names(common_vars$monthly[[1]])

# Step 1.4: Create compatible data for each country using common variables
country_list <- list()
for (cc in countries) {
  res <- prepare_country_data(
    cc, P,
    selected_vars_m = common_vars$monthly,
    selected_vars_q = common_vars$quarterly,
    Covid = TRUE
  )
  country_list[[cc]] <- list(
    Data = res$Data,
    Dates = res$DatesM,
    Series = res$Series
  )
}


# ==============================================================================
# STEP 2: BUILD TENSOR DATASET (T × p1 × p2)
# ==============================================================================

tensor <- tensor_data(countries, P)
Y <- tensor$Y  # Data tensor
W <- tensor$W  # Mask (missing indicator)

# Standardize and demean
std <- standardize_Y(Y)
Y_std <- std$Y_scaled

# Check percentage of missing values
nan_percent_Y(Y_std)



################################# !is.na GDP ###################################


ea_index <- which(dimnames(tensor$Y)[[2]] == "IT")
gdp_idx <- which(dimnames(tensor$Y)[[3]] == "GDP")

# Estrai la serie del GDP (standardizzata)
gdp_ts <- tensor$Y[, ea_index, gdp_idx]


# Serie temporale con Date
gdp_df <- data.frame(
  Date = res$DatesM[1:P$Tmax],              # oppure converti in trimestre con get_quarter()
  GDP  = gdp_ts
)

ggplot(gdp_df, aes(x = Date, y = GDP)) +
  geom_line(color = "#1f77b4", size = 1) +
  geom_point(data = subset(gdp_df, !is.na(GDP)), color = "black", size = 1.5) +
  labs(
    title = "Observed GDP Series from Tensor",
    subtitle = "Standardized GDP values for Italy (example)",
    x = "Date", y = "GDP (Standardized)"
  ) +
  theme_minimal(base_size = 13)


# ==============================================================================
# STEP 3: FACTOR ESTIMATION
# ==============================================================================

## === DMFM Pipeline === ##
k_hat <- mfm.f(Y_std, W, kmax)  # Estimate number of factors
f.lag_mar <- mar_model_selection_auto(Y_std, W, k_hat)  # Select lag order

## === DFM Alternative === ##
# k_hat <- dfm.f.vec(Y_std, W, kmax)
# f.lag_mar <- mar_model_selection_vec(Y_std, W, k_hat)


# ==============================================================================
# STEP 4: INITIALIZE MODEL INPUTS
# ==============================================================================

## === DMFM Initialization === ##
imp <- mfm.cl(Y_std, W, k_hat)
inputs <- dmfm.na.2.sv(imp$Y_imputed, k_hat, W, t = "dmfm")

## === DFM Alternative Initialization === ##
# imp <- mfm.cl.vec(Y_std, W, k_hat)
# inputs <- dmfm.na.2.sv.vec(imp$Y_imputed, k_hat, W, t = "dmfm")


# ==============================================================================
# STEP 5: FIT DMFM USING FULL DATASET
# ==============================================================================

out <- dmfm.na.em(
  . = inputs,
  Y = Y_std,
  k = k_hat,
  W = W,
  t = "dmfm"
)


# ==============================================================================
# STEP 6: ROLLING NOWCASTING (PSEUDO REAL-TIME)
# ==============================================================================

roll_nowcast <- rolling_nowcast_dmfm_by_release(
  Y_full = Y,
  W_full = W,
  k = k_hat,
  dates = res$DatesM,
  group = res$Group,
  gdp_idx = res$quarterly_start_idx,
  startEV = P$startEV,
  endEV = P$endEV,
  model_type = "dmfm",
  max_iter = 1000,
  eps = 1e-3
)


# ==============================================================================
# STEP 7: SAVE OUTPUT
# ==============================================================================

save(
  res, tensor, k_hat, W, inputs, std, out, roll_nowcast, 
  file = "Results/dmfm_rollnowcastLS_11_204_40varCovidpluse-03.RData"
)

# (Optional) Load previously saved results
# load("Results/dmfm_rollnowcastLS_11_204_40var.RData")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================