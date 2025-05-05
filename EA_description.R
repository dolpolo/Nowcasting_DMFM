# **************************************************************************** #
# ************************ Dynamic Matrix Factor Model *********************** #
# **************************************************************************** #  

# This code aims to Estimate a Dynamic Matrix Factor Model and providing Forecasts
# on the GDP growth for EA and specific countries such as IT,FR,ES,DE.

################################################################################

# ==============================================================================
# SET DIRECTORY
# ==============================================================================
getwd()
path <- "C:/Users/david/Desktop/Thesis"
setwd(path)


# ==============================================================================
# CALL LIBRARIES
# ==============================================================================

# Data Manipulation
library(tidyverse)
library(lubridate)
library(tidyr)
library(dplyr)

# Time Series
library(tseries)
library(zoo)
library(fBasics)

# Linear Algebra
library(Matrix)
library(RSpectra)
library(MASS)
library(pracma)

# Input/Output
library(writexl)
library(readxl)
library(xtable)

# Visualization
library(ggplot2)
library(reshape2)
library(scales)

# Matlab reading
library(R.matlab)

# GIS 
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# LaTeX
library(kableExtra)
library(knitr)



# ==============================================================================
# CALL FUNCTIONS
# ==============================================================================

# Data Preparation 
source("C:/Users/david/Desktop/Thesis/Functions/R/DataPreparation.R")


# ==============================================================================
# GDP CONTRIBUTION FOR SINGLE COUNTRY TO EA
# ==============================================================================

GDP_com <- GDP_communality(c("DE", "FR", "IT", "ES"))

# ==============================================================================
# COUNTRY NATIONS IN THE DATASET
# ==============================================================================

plot_euro_area_gdp(file_path = "Data/GDP_millionsEA.xlsx", quarter = "2024-Q4")

