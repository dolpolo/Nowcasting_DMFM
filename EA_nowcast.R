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
path <- "C:/Users/david/Desktop/University/Master's Thesis/Nowcasting_DMFM"
setwd(path)


# ==============================================================================
# CALL LIBRARIES
# ==============================================================================

# Data Manipulation
library(tidyverse)
library(lubridate)
library(abind)

# Time Series
library(tseries)
library(zoo)
library(fBasics)
library(vars)

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
library(patchwork)


# ==============================================================================
# CALL FUNCTIONS
# ==============================================================================

# Data Preparation 
source("C:/Users/david/Desktop/University/Master's Thesis/Nowcasting_DMFM/Functions/R/DataPreparation.R")

# General utilities
source("C:/Users/david/Desktop/University/Master's Thesis/Nowcasting_DMFM/Functions/R/matrix.utils.R")

# DMFM
source("C:/Users/david/Desktop/University/Master's Thesis/Nowcasting_DMFM/Functions/R/dmfm.initialization.R")
source("C:/Users/david/Desktop/University/Master's Thesis/Nowcasting_DMFM/Functions/R/kalman.R")
source("C:/Users/david/Desktop/University/Master's Thesis/Nowcasting_DMFM/Functions/R/dmfm.estimation.R")
source("C:/Users/david/Desktop/University/Master's Thesis/Nowcasting_DMFM/Functions/R/dmfm.results.R")

# ==============================================================================
# SET PARAMETERS
# ==============================================================================

# No data for NL and BE in 2024 and afterwards
countries <- c("DE","FR","IT","ES")

# Parameters for the data preparation
P <- list()

P$modelM <- "Large"
P$modelQ <- "Small"
P$covid_start <- as.Date("2020-03-01")
P$covid_end <- as.Date("2020-12-01")
P$startEV <- as.Date("2017-01-01")
P$endEV <- as.Date("2025-01-01")

P$Tmax <- 300

# Initial guess for the number of factors
kmax <- c(3, 7)


# ==============================================================================
# STEP 1: PREPARE COUNTRY-SPECIFIC 2D DATASETS
# ==============================================================================

# Step 1: selezione dei nomi base delle variabili (più correlate al GDP per paese)
base_names <- select_base_variable_names(countries, P)

# Step 2: filtra le variabili base comuni a tutti i paesi (con suffissi)
common_vars <- filter_common_variables_across_countries(countries, base_names$monthly, base_names$quarterly)

# Step 3: estrai i nomi base comuni (sono gli stessi per tutti i paesi)
all_vars <- names(common_vars$monthly[[1]])  # nomi base comuni, senza suffisso

# Step 4: prepara i dati per ogni paese usando solo le variabili comuni
country_list <- list()
for (cc in countries) {
  res <- prepare_country_data(cc, P, 
                              selected_vars_m = common_vars$monthly, 
                              selected_vars_q = common_vars$quarterly)
  country_list[[cc]] <- list(
    Data = res$Data,
    Dates = res$DatesM,
    Series = res$Series
  )
}

# ==============================================================================
# STEP 2: BUILD TENSOR Y (T × p1 × p2) AND MASK W
# ==============================================================================

tensor <- tensor_data(countries, P)
Y <- tensor$Y
W <- tensor$W
summary(Y)


# Demean and Standardize
std <- standardize_Y(Y)
Y_std <- std$Y_scaled
summary(Y_std)

# Verify the percentage of NaN
nan_percent_Y(Y_std)

# ==============================================================================
# STEP 3: ANALYSIS ON FACTORS
# ==============================================================================

# ---- DMFM ----
# Number of factors: ER by Yu Yu et al. (JoE, 2021)
k_hat <- mfm.f(Y_std, W, kmax)

# Number of lags in MAR equation
# f.lag_mar <- mar_model_selection(Y_std, W, k_hat)
f.lag_mar <- mar_model_selection_auto(Y_std, W, k_hat)

# ---- DFM ----
# Number of factors: ER by Yu Yu et al. (JoE, 2021)
# k_hat <- dfm.f.vec(Y_std, W, kmax)

# Number of lags in MAR equation
# f.lag_mar <- factorize_bic_mar_vec(Y_std, W, k_hat)
# f.lag_mar <- mar_model_selection_vec(Y_std, W, k_hat)


# ==============================================================================
# STEP 3: LIST OF INPUTS --> . 
# ==============================================================================

# Just to try
# k_hat <- c(3,3)

# ---- DMFM ----
# Can and Lam imputation
imp <- mfm.cl(Y_std, W, k_hat)
# Trapin Method (Barigozzi's variant) for initialization
inputs <- dmfm.na.2.sv(imp$Y_imputed, k_hat, W, t = "dmfm")


# ---- DFM ----
# Can and Lam imputation
# imp <- mfm.cl.vec(Y_std, W, k_hat)
# Trapin Method (Barigozzi's variant) for initialization
# inputs <- dmfm.na.2.sv.vec(imp$Y_imputed, k_hat, W, t = "dmfm")


# ==============================================================================
# STEP 4: INTERPRETATION OF FACTORS 
# ==============================================================================

# Analyse the loading structure and plot the factor 
l.R_str <- inputs$R
l.C_str <- inputs$C

plot_factors_tensor(inputs$fs)

# ==============================================================================
# STEP 5: DMFM ESTIMATION (FULL DATASET)
# ==============================================================================

# Modello con input Dolp
out <- dmfm.na.em(
  . = inputs,
  Y = Y_std,
  k = k_hat,
  W = W,
  t = "dmfm"
)


# ==============================================================================
# STEP 6: ROLLING NOWCAST
# ==============================================================================

# ---- DMFM ----

# Rolling Nowcast
Y_nowcast <- rolling_nowcast_dmfm(Y_std, W, k = k_hat, t_vec = 276:300)
nowcasts <- Y_nowcast$nowcasts

# ---- DFM ----

# Rolling Nowcast
# Y_nowcast <- rolling_nowcast_dfm(Y_std, W, k = k_hat, t_vec = 276:300)
# nowcasts <- Y_nowcast$nowcasts


# ==============================================================================
# STEP 5: SAVE RESULTS
# ==============================================================================

# save(res,tensor,inputs, std, roll_nowcast, k_hat, W, file = "dmfm_rollnowcastLS_11_201_40var.RData")

# load("dmfm_shiftnowcastLS_11_201_40var.RData")


shift_release <- function(Y, W, quarterly_var_indices) {
  T <- dim(Y)[1]
  p1 <- dim(Y)[2]
  p2 <- dim(Y)[3]
  
  for (idx in quarterly_var_indices) {
    for (j in 1:p1) {
      tmp <- Y[, j, idx]
      tmp_shifted <- c(NA, tmp[1:(T - 1)])  # shift di 1 mese in avanti
      Y[, j, idx] <- tmp_shifted
      W[, j, idx] <- as.numeric(!is.na(tmp_shifted))
    }
  }
  return(list(Y_shifted = Y, W_shifted = W))
}

nQ <- res$quarterly_start_idx

shift <- shift_release(Y_std,W,nQ)
Y_s <- shift$Y_shifted
W_s <- shift$W_shifted

t_vec <- which(res$DatesM >= P$startEV & res$DatesM <= P$endEV)

# Rolling Nowcast
Y_nowcast <- rolling_nowcast_dmfm(Y_s, W_s, k = k_hat, t_vec = t_vec)




################################################################################
######### AGGIORNAMENTO NOWCAST MENSILE CONSIDERANDO RELEASES BOZZA ###########
Y_full <- Y_std
W_full <- W
k <- k_hat
dates <- res$DatesM
group <- res$Group
quarterly_var_indices <- res$quarterly_start_idx

rolling_nowcast_dmfm_release <- function(Y_full, W_full, k, dates, group,
                                         quarterly_var_indices = NULL,
                                         model_type = "dmfm",
                                         max_iter = 1000,
                                         eps = 1e-3,
                                         target_idx = NULL) {
  T <- dim(Y_full)[1]
  p1 <- dim(Y_full)[2]
  p2 <- dim(Y_full)[3]
  
  t_vec <- which(dates >= P$startEV & dates <= P$endEV)
  
  nowcast_list <- list()
  factors_list <- list()
  
  for (t in t_vec) {
    cat(sprintf("Rolling nowcast at vintage t = %d (%s)\n", t, as.character(dates[t])))
    
    #' t = 288 # mese di riferimento per Q(0)
    #' t = 289
    #' t = 290
    #' t = 291 # mese di release per Q(0)
    
    
    Y_t <- Y_full[1:t,,, drop = FALSE]
    W_t <- W_full[1:t,,, drop = FALSE]
    
    for (j in 1:p2) {
      d_j <- group[j]
      last_obs_time <- t - d_j  # ultimo mese il cui valore è disponibile al tempo t
      if (last_obs_time < 1) {
        # Nessuna osservazione disponibile per questa variabile
        Y_t[, , j] <- NA
        W_t[, , j] <- 0
      } else {
        # Maschera i valori da (last_obs_time + 1) in poi
        Y_t[(last_obs_time + 1):t, , j] <- NA
        W_t[(last_obs_time + 1):t, , j] <- 0
      }
    }
    
    
    # prova1 <- Y_t[,1,]
    # prova2 <- Y_t[,1,]
    # prova3 <- Y_t[,1,]
    # prova4 <- Y_t[,1,]
    
    # --- Inizializzazione EM ---
    imp <- mfm.cl(Y_t, W_t, k)
    inputs <- dmfm.na.2.sv(imp$Y_imputed, k, W_full, t = model_type)
    
    fit <- dmfm.na.em(inputs, Y = Y_t, k = k, W = W_t, t = model_type,
                      max.iter = max_iter, eps = eps)
    model <- fit$model
    
    f_t <- model$ff[t,,]
    Y_hat <- model$R %*% f_t %*% t(model$C)
    
    nowcast_value <- if (!is.null(target_idx)) {
      Y_hat[target_idx[[1]], target_idx[[2]]]
    } else {
      Y_hat
    }
    
    nowcast_list[[as.character(t)]] <- nowcast_value
    factors_list[[as.character(t)]] <- f_t
  }
  
  return(list(nowcasts = nowcast_list, factors = factors_list))
}




################################################################################
########################## AGGIORNAMENTO NOWCAST MENSILE CONSIDERANDO RELEASES##
Y_full <- Y_std
W_full <- W
k <- k_hat
dates <- res$DatesM
group <- res$Group
gdp_idx <- res$quarterly_start_idx



rolling_nowcast_dmfm_by_release <- function(Y_full, W_full, k, dates, group, gdp_idx, 
                                            model_type = "dmfm", max_iter = 1000, eps = 1e-3,
                                            t_start, t_end = NULL) {
  T <- dim(Y_full)[1]
  if (is.null(t_end)) t_end <- T
  p1 <- dim(Y_full)[2]
  p2 <- dim(Y_full)[3]
  
  # Output lists for each month in the quarter
  nowcast_M1 <- list()
  nowcast_M2 <- list()
  nowcast_M3 <- list()
  
  for (t in seq(t_start, t_end)) {
    
    # t = 288
    # t = 289
    # t = 290
    # t = 291 
    # t = 292
    
    cat(sprintf(">>> Rolling nowcast at t = %d (%s)\n", t, as.character(dates[t])))
    
    # Subset up to time t
    Y_t <- Y_full[1:t, , , drop = FALSE]
    W_t <- W_full[1:t, , , drop = FALSE]
    
    # Demean and Standardize
    std <- standardize_Y(Y_t)
    Y_t_std <- std$Y_scaled
    
    # pre_prova1 <- Y_t_std[,1,]
    # pre_prova2 <- Y_t_std[,1,]
    # pre_prova3 <- Y_t_std[,1,]
    # pre_prova4 <- Y_t_std[,1,]
    # pre_prova5 <- Y_t_std[,1,]
    
    for (j in 1:p2) {
      freq_j <- group[j]
      delay <- freq_j - 1
      release_time <- t - delay
      
      if (release_time >= 1) {
        # Cancelliamo solo le osservazioni future non ancora rilasciate
        if (release_time + 1 <= t) {
          Y_t_std[(release_time + 1):t, , j] <- NA
          W_t[(release_time + 1):t, , j] <- 0
        }
        # Le osservazioni precedenti al release_time rimangono
      } else {
        # Se la release_time è prima dell'inizio dei dati, maschera tutto
        Y_t_std[, , j] <- NA
        W_t[, , j] <- 0
      }
    }
    
    # prova1 <- Y_t_std[ ,1,]
    # prova2 <- Y_t_std[ ,1,]
    # prova3 <- Y_t_std[ ,1,]
    # prova4 <- Y_t_std[ ,1,]
    # prova5 <- Y_t_std[ ,1,]
    
    
    # EM Estimation
    imp <- mfm.cl(Y_t_std, W_t, k)
    inputs <- dmfm.na.2.sv(imp$Y_imputed, k, W_full, t = model_type)
    fit <- dmfm.na.em(inputs, Y = Y_t_std, k = k, W = W_t, t = model_type, max.iter = max_iter, eps = eps)
    model <- fit$model
    
    # Reconstruct nowcast from filtered factors
    f_t <- model$ff[t, , ]
    Y_hat <- model$R %*% f_t %*% t(model$C)
    Y_hat_true <- Y_hat * std$sd + std$mean
    gdp_forecast <- Y_hat[, gdp_idx]
    
    # Identify month in the quarter: 1, 2, 3
    m_trimestre <- ((as.integer(format(dates[t], "%m")) - 1) %% 3) + 1
    
    # Save nowcast by quarter month
    date_str <- as.character(dates[t])
    if (m_trimestre == 1) {
      nowcast_M1[[date_str]] <- gdp_forecast
    } else if (m_trimestre == 2) {
      nowcast_M2[[date_str]] <- gdp_forecast
    } else if (m_trimestre == 3) {
      nowcast_M3[[date_str]] <- gdp_forecast
    }
  }
  
  return(list(M1 = nowcast_M1, M2 = nowcast_M2, M3 = nowcast_M3))
}

roll_nowcast <- rolling_nowcast_dmfm_by_release(Y_full = Y, # lo standardizzo ogni volta?
                                                W_full = W,
                                                k = k_hat,
                                                dates = res$DatesM,
                                                group = res$Group,
                                                gdp_idx = res$quarterly_start_idx,
                                                model_type = "dmfm",
                                                max_iter = 1000,
                                                eps = 1e-3,
                                                t_start = 288,
                                                t_end = NULL)



# ==============================================================================
# NOWCAST (FINAL VERSION) Pseudo real time exercise
# ==============================================================================

rolling_nowcast_dmfm_by_release <- function(Y_full, W_full, k, dates, group, gdp_idx, 
                                            startEV, endEV,
                                            model_type = "dmfm", max_iter = 1000, eps = 1e-3) {
  T <- dim(Y_full)[1]
  p1 <- dim(Y_full)[2]
  p2 <- dim(Y_full)[3]
  
  # Calcola gli indici temporali corrispondenti alle date
  t_start <- which(dates == startEV)
  t_end <- which(dates == endEV)
  if (length(t_start) == 0 || length(t_end) == 0) {
    stop("startEV o endEV non trovate nel vettore dates.")
  }
  
  # Output lists per i tre mesi del trimestre
  nowcast_M1 <- list()
  nowcast_M2 <- list()
  nowcast_M3 <- list()
  
  for (t in seq(t_start, t_end)) {
    cat(sprintf(">>> Rolling nowcast at t = %d (%s)\n", t, as.character(dates[t])))
    
    # Subset e maschere release
    Y_t <- Y_full[1:t, , , drop = FALSE]
    W_t <- W_full[1:t, , , drop = FALSE]
    
    # Standardizza il dataset fino a t
    std <- standardize_Y(Y_t)
    Y_t_std <- std$Y_scaled
    
    # Applica regole di rilascio
    for (j in 1:p2) {
      freq_j <- group[j]
      delay <- freq_j - 1
      release_time <- t - delay
      if (release_time >= 1) {
        if (release_time + 1 <= t) {
          Y_t_std[(release_time + 1):t, , j] <- NA
          W_t[(release_time + 1):t, , j] <- 0
        }
      } else {
        Y_t_std[, , j] <- NA
        W_t[, , j] <- 0
      }
    }
    
    # EM Estimation
    imp <- mfm.cl(Y_t_std, W_t, k)
    inputs <- dmfm.na.2.sv(imp$Y_imputed, k, W_full, t = model_type)
    fit <- dmfm.na.em(inputs, Y = Y_t_std, k = k, W = W_t, t = model_type, max.iter = max_iter, eps = eps)
    model <- fit$model
    
    # Ricostruisci il nowcast
    f_t <- model$ff[t, , ]
    Y_hat <- model$R %*% f_t %*% t(model$C)
    Y_hat_true <- Y_hat * std$sd + std$mean
    gdp_forecast <- Y_hat_true[, gdp_idx]
    
    # Salva per il mese del trimestre corrispondente
    m_trimestre <- ((as.integer(format(dates[t], "%m")) - 1) %% 3) + 1
    date_str <- as.character(dates[t])
    if (m_trimestre == 1) {
      nowcast_M1[[date_str]] <- gdp_forecast
    } else if (m_trimestre == 2) {
      nowcast_M2[[date_str]] <- gdp_forecast
    } else if (m_trimestre == 3) {
      nowcast_M3[[date_str]] <- gdp_forecast
    }
  }
  
  return(list(M1 = nowcast_M1, M2 = nowcast_M2, M3 = nowcast_M3))
}

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

save(res,tensor,inputs, std, roll_nowcast, k_hat, W, file = "dmfm_rollnowcastLS_11_201_40var.RData")

load("dmfm_rollnowcastLS_11_201_40var.RData")
# ==============================================================================
# RMSFE FOR ANY MONTH NOWCAST 
# ==============================================================================
DatesM <- res$DatesM[1:300]
Y_std <- std$Y_scaled


compute_rmsfe <- function(nowcast_list, Y_true, gdp_idx, DatesM) {
  # nowcast_list: lista con nomi-date (es. "2024-01-01"), ciascun elemento è un vettore di p1 paesi
  # Y_true: array [T, p1, p2] destandardizzato
  # gdp_idx: indice colonna GDP
  # DatesM: vettore delle date mensili corrispondenti a Y_true (es. res$DatesM)
  
  dates_vec <- as.Date(names(nowcast_list))  # vettore di date
  n <- length(dates_vec)
  p1 <- dim(Y_true)[2]
  rmsfe_mat <- matrix(NA, nrow = p1, ncol = n)
  
  for (i in seq_along(dates_vec)) {
    date_t <- dates_vec[i]
    t <- which(DatesM == date_t)
    if (length(t) == 0) next
    
    # Calcola il mese all’interno del trimestre
    m_trimestre <- (as.integer(format(date_t, "%m")) - 1) %% 3 + 1
    t_gdp_true <- t + (3 - m_trimestre) + 1 # Fine trimestre ( wirdoooo)
    
    if (t_gdp_true > dim(Y_true)[1]) next  # fuori dai dati
    for (j in 1:p1) {
      forecast <- nowcast_list[[i]][j]
      actual <- Y_true[t_gdp_true, j, gdp_idx]
      if (!is.na(actual) && !is.na(forecast)) {
        rmsfe_mat[j, i] <- (forecast - actual)^2
      }
    }
  }
  
  # Media per paese
  rmsfe_vec <- sqrt(rowMeans(rmsfe_mat, na.rm = TRUE))
  names(rmsfe_vec) <- paste0("Country_", 1:p1)
  return(rmsfe_vec)
}

# Destandardizza i dati veri
Y_true <- inverse_standardize_Y(Y_std, std$mean, std$sd)

# Calcola RMSFE per ogni mese del trimestre
rmsfe_M1 <- compute_rmsfe(roll_nowcast$M1, Y_true, gdp_idx = res$quarterly_start_idx, DatesM)
rmsfe_M2 <- compute_rmsfe(roll_nowcast$M2, Y_true, gdp_idx = res$quarterly_start_idx, DatesM)
rmsfe_M3 <- compute_rmsfe(roll_nowcast$M3, Y_true, gdp_idx = res$quarterly_start_idx, DatesM)


# Crea dataframe
rmsfe_df <- data.frame(
  Country = names(rmsfe_M1),
  M1 = rmsfe_M1,
  M2 = rmsfe_M2,
  M3 = rmsfe_M3
)

# Plot con ggplot2
library(tidyr)
library(ggplot2)

rmsfe_long <- pivot_longer(rmsfe_df, cols = c(M1, M2, M3), names_to = "Month", values_to = "RMSFE")

ggplot(rmsfe_long, aes(x = Country, y = RMSFE, fill = Month)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "RMSFE of GDP Nowcast - Month of Quarter",
    x = "Country", y = "RMSFE"
  ) +
  scale_fill_manual(values = c("M1" = "steelblue", "M2" = "orange", "M3" = "darkred")) +
  theme_minimal(base_size = 14)



# =============================================================================
# Boxplot RMSFME
# =============================================================================

compute_rmsfe_long <- function(nowcast_list, Y_true, gdp_idx, DatesM, month_label = "M1") {
  dates_vec <- as.Date(names(nowcast_list))
  p1 <- dim(Y_true)[2]
  result_df <- data.frame()
  
  for (i in seq_along(dates_vec)) {
    date_t <- dates_vec[i]
    t <- which(DatesM == date_t)
    if (length(t) == 0) next
    
    m_trimestre <- (as.integer(format(date_t, "%m")) - 1) %% 3 + 1
    t_gdp_true <- t + (3 - m_trimestre) + 1
    if (t_gdp_true > dim(Y_true)[1]) next
    
    for (j in 1:p1) {
      forecast <- nowcast_list[[i]][j]
      actual <- Y_true[t_gdp_true, j, gdp_idx]
      if (!is.na(actual) && !is.na(forecast)) {
        err <- forecast - actual
        result_df <- rbind(result_df, data.frame(
          Country = paste0("Country_", j),
          Date = date_t,
          Month = month_label,
          Error = err,
          RMSFE = sqrt(err^2)
        ))
      }
    }
  }
  return(result_df)
}

# Calcola errori per ogni mese del trimestre
rmsfe_long_M1 <- compute_rmsfe_long(roll_nowcast$M1, Y_true, gdp_idx = res$quarterly_start_idx, DatesM, "M1")
rmsfe_long_M2 <- compute_rmsfe_long(roll_nowcast$M2, Y_true, gdp_idx = res$quarterly_start_idx, DatesM, "M2")
rmsfe_long_M3 <- compute_rmsfe_long(roll_nowcast$M3, Y_true, gdp_idx = res$quarterly_start_idx, DatesM, "M3")

# Combina in un unico dataframe
rmsfe_all_long <- bind_rows(rmsfe_long_M1, rmsfe_long_M2, rmsfe_long_M3)

# Boxplot
ggplot(rmsfe_all_long, aes(x = Month, y = RMSFE, fill = Month)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, color = "gray40", alpha = 0.5) +
  facet_wrap(~ Country, scales = "free_y") +
  labs(
    title = "Distribuzione RMSFE per Paese nei Mesi del Trimestre",
    x = "Mese del Trimestre", y = "RMSFE"
  ) +
  scale_fill_manual(values = c("M1" = "steelblue", "M2" = "orange", "M3" = "darkred")) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")


# =============================================================================
# Plot aggiornamenti del nowcast 
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

plot_nowcast_vs_true_gdp <- function(roll_nowcast, Y_true, gdp_idx, DatesM) {
  p1 <- dim(Y_true)[2]  # numero paesi
  
  # Nomi coerenti: "Country_1", ..., "Country_p1"
  country_names <- paste0("Country_", 1:p1)
  
  # Funzione ausiliaria per costruire il data.frame da ciascuna lista di nowcast
  build_df <- function(nc_list, label) {
    df <- data.frame(
      Date = as.Date(names(nc_list)),
      do.call(rbind, nc_list)
    )
    colnames(df)[-1] <- country_names
    df$Month <- label
    return(df)
  }
  
  # Costruzione tabelle per ciascun mese
  df_M1 <- build_df(roll_nowcast$M1, "M1")
  df_M2 <- build_df(roll_nowcast$M2, "M2")
  df_M3 <- build_df(roll_nowcast$M3, "M3")
  
  # Unione di tutti i nowcast
  df_all <- bind_rows(df_M1, df_M2, df_M3) %>%
    pivot_longer(cols = all_of(country_names), names_to = "Country", values_to = "Nowcast") %>%
    mutate(Country = factor(Country, levels = country_names))
  
  # Serie reale del GDP (solo nei mesi in cui esce: ogni 3 mesi)
  df_true <- data.frame(
    Date = DatesM,
    Y_true_gdp = Y_true[ , , gdp_idx]
  )
  colnames(df_true)[-1] <- country_names
  df_true <- df_true %>%
    slice(seq(3, n(), by = 3)) %>%
    pivot_longer(cols = all_of(country_names), names_to = "Country", values_to = "True")
  
  # Merge tra nowcast e serie reale
  plot_df <- left_join(df_all, df_true, by = c("Date", "Country"))
  
  # Plot
  ggplot(plot_df, aes(x = Date)) +
    geom_line(aes(y = Nowcast, color = Month), linewidth = 1) +
    geom_point(aes(y = True), color = "black", size = 2, shape = 16, na.rm = TRUE) +
    facet_wrap(~ Country, scales = "free_y") +
    labs(
      title = "Nowcast del GDP vs Valori Reali per Paese",
      x = "Data", y = "GDP",
      color = "Mese del Trimestre"
    ) +
    scale_color_manual(values = c("M1" = "steelblue", "M2" = "orange", "M3" = "darkgreen")) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")
}


plot_nowcast_vs_true_gdp(roll_nowcast, Y_true, gdp_idx = res$quarterly_start_idx, DatesM)



# =============================================================================
# Plot aggiornamenti del nowcast per ogni mese
# =============================================================================
plot_monthly_nowcast_series <- function(nowcast_list, Y_true, gdp_idx, DatesM, month_label = "M1") {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  
  p1 <- dim(Y_true)[2]  # numero di paesi
  country_names <- paste0("Country_", 1:p1)
  
  # Costruisci tabella dei nowcast
  df_nowcast <- data.frame(
    Date = as.Date(names(nowcast_list)),
    do.call(rbind, nowcast_list)
  )
  colnames(df_nowcast)[-1] <- country_names
  
  df_nowcast <- df_nowcast %>%
    pivot_longer(-Date, names_to = "Country", values_to = "Nowcast") %>%
    mutate(Country = factor(Country, levels = country_names)) %>%
    mutate(
      MonthNum = as.integer(format(Date, "%m")),
      Shift = 3 - ((MonthNum - 1) %% 3),
      Date_adjusted = Date + months(Shift)
    )
  
  # Serie reale del GDP (ogni 3 mesi)
  df_true <- data.frame(
    Date = DatesM,
    Y_true_gdp = Y_true[, , gdp_idx]
  )
  colnames(df_true)[-1] <- country_names
  df_true <- df_true %>%
    slice(seq(3, n(), by = 3)) %>%
    pivot_longer(-Date, names_to = "Country", values_to = "True")
  
  # Merge tra Date_adjusted e Date vera
  df_plot <- df_nowcast %>%
    left_join(df_true, by = c("Date_adjusted" = "Date", "Country"))
  
  # Plot
  ggplot(df_plot, aes(x = Date)) +
    geom_line(aes(y = Nowcast), color = "blue", linewidth = 1) +
    geom_point(aes(y = True), color = "red", size = 2, shape = 16, na.rm = TRUE) +
    facet_wrap(~ Country, scales = "free_y") +
    labs(
      title = paste("Nowcast Serie Storica da", month_label, "vs True GDP"),
      x = "Data (mese previsione)", y = "GDP"
    ) +
    theme_minimal(base_size = 13)
}

plot_monthly_nowcast_series(roll_nowcast$M1, Y_true, res$quarterly_start_idx, DatesM, month_label = "M1")
plot_monthly_nowcast_series(roll_nowcast$M2, Y_true, res$quarterly_start_idx, DatesM, month_label = "M2")
plot_monthly_nowcast_series(roll_nowcast$M3, Y_true, res$quarterly_start_idx, DatesM, month_label = "M3")

