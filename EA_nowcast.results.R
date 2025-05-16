# ==============================================================================
#               Dynamic Matrix Factor Models and the EM Algorithm: 
#          A Nowcasting Framework for Mixed Frequency Data in Euro Area
#                               - Results Script - 
# ==============================================================================
# This script analyzes the outputs from a Dynamic Matrix Factor Model (DMFM)
# estimated via the EM algorithm. It performs a comprehensive evaluation of the
# nowcasting performance for GDP growth in major Euro Area economies:
# Germany (DE), France (FR), Italy (IT), and Spain (ES).
#
# Main Components:
#   1. Load results from estimation
#   2. Interpret latent factors and loadings
#   3. Compute RMSFE for each forecast month in the quarter
#   4. Evaluate pre-COVID vs post-COVID nowcasting performance
#   5. Visualize results (factors, errors, nowcasts)
#
# ==============================================================================
# Author   : Davide Delfino
# University: Alma Mater Studiorum - University of Bologna
# Dataset  : Barigozzi & Lissona (2024), EA-MD-QD
#
# ==============================================================================
# Script Type : Descriptive / Visualization
# ==============================================================================

# ==============================================================================
# SET WORKING DIRECTORY
# ==============================================================================
path <- "C:/Users/david/Desktop/University/Master's Thesis/Nowcasting_DMFM"
setwd(path)

# ==============================================================================
# LOAD REQUIRED LIBRARIES
# ==============================================================================
library(tidyverse)
library(lubridate)
library(purrr)
library(abind)
library(writexl)
library(readxl)
library(xtable)
library(ggplot2)
library(reshape2)
library(patchwork)

# Matlab reading
library(R.matlab)

# ==============================================================================
# 0. LOAD DEPENDENCIES & USER FUNCTIONS
# ==============================================================================
library(ggplot2)
library(tidyr)
library(patchwork)

source("Functions/R/dmfm.preparation.R")     # Destandardization, utilities
source("Functions/R/dmfm.visualization.R")   # Plotting tools

# ==============================================================================
# 1. PARAMETERS
# ==============================================================================
countries    <- c("DE", "FR", "IT", "ES")
month_labels <- paste0("M", 1:3)
n_months     <- length(month_labels)
month_index  <- 1  # <== USER SELECTION: 1 = M1, 2 = M2, 3 = M3
country_code <- "IT"


P <- list(
  modelM      = "Large",
  modelQ      = "Small",
  covid_start = as.Date("2020-03-01"),
  covid_end   = as.Date("2020-12-01"),
  startEV     = as.Date("2017-01-01"),
  endEV       = as.Date("2025-01-01"),
  Tmax        = 300
)


dmfm_file <- "Results/dmfm_rollnowcastLS_11_204_40var.RData"
# dmfm_file <- "Results/dmfm_rollnowcastLS_12_204_12varCovidpluse-04.RData"
# dmfm_file <- "Results/dmfm_rollnowcastLS_11_204_40varCovidpluse-03.RData"


# dfm_file  <- "Results/DFM/fcst_DE_11.mat"
# dfm_file  <- "Results/DFM/fcst_FR_11.mat"
dfm_file  <- "Results/DFM/fcst_IT_11.mat"
# dfm_file  <- "Results/DFM/fcst_ES_11.mat"


# ==============================================================================
# 2. LOAD DMFM DATA
# ==============================================================================
load(dmfm_file)
DatesM <- res$DatesM[1:P$Tmax]
Y_std  <- std$Y_scaled
Y_true <- inverse_standardize_Y(Y_std, std$mean, std$sd)

# ==============================================================================
# 3. LOAD DFM DATA
# ==============================================================================
mat_data <- readMat(dfm_file)
DateQQ_dt <- as.Date(as.numeric(mat_data$DateQQ), origin = "0000-01-01") - days(366)
FcstQQ <- mat_data$FcstQQ
TrueQQ <- mat_data$TrueQQ

# Rinomina colonne FcstQQ
n_quarters <- ncol(FcstQQ) / n_months
Q_vals     <- seq(from = -1, length.out = n_quarters)
colnames(FcstQQ) <- unlist(lapply(Q_vals, function(q) paste0(month_labels, "(Q(", q, "))")))

# ==============================================================================
# 4. PREPARE GDP ARRAY (DMFM)
# ==============================================================================
get_quarter <- function(date) {
  paste0(format(as.Date(date), "%Y"), "Q", ceiling(as.numeric(format(as.Date(date), "%m")) / 3))
}

for (m in month_labels) {
  names(roll_nowcast[[m]]) <- sapply(names(roll_nowcast[[m]]), get_quarter)
}

common_quarters <- sort(Reduce(intersect, lapply(roll_nowcast[month_labels], names)))
n_qtr       <- length(common_quarters)
n_countries <- length(roll_nowcast$M1[[1]])

gdp_array <- array(NA, dim = c(n_qtr, n_months, n_countries),
                   dimnames = list(common_quarters, month_labels, paste0("Country", 1:n_countries)))

for (i in seq_along(common_quarters)) {
  qtr <- common_quarters[i]
  for (m in 1:n_months) {
    gdp_array[i, paste0("M", m), ] <- roll_nowcast[[paste0("M", m)]][[qtr]]
  }
}

# ==============================================================================
# 5. SELECT COUNTRY AND PREPARE COMPARISON DATA
# ==============================================================================
ea_index     <- which(dimnames(tensor$Y)[[2]] == country_code)
gdp_array    <- gdp_array[1:30, , ea_index]
gdp_true_ea  <- TrueQQ[1:30]
gdp_array_dfm <- FcstQQ[, 4:6]  # M1, M2, M3

# Set dates and quarter names
DatesQ       <- unique(sapply(DatesM, get_quarter))
DatesQ_30    <- DatesQ[69:98]
DatesQ_real  <- DatesM[match(DatesQ_30, sapply(DatesM, get_quarter))]
dimnames(gdp_array)[[1]] <- DatesQ_30
rownames(gdp_array_dfm)  <- DatesQ_30

# ==============================================================================
# 6. COMPUTE RMSE
# ==============================================================================
rmse_dmfm <- sapply(1:3, function(m) {
  sqrt(mean((gdp_array[, m] - gdp_true_ea)^2, na.rm = TRUE))
})
names(rmse_dmfm) <- month_labels
print(round(rmse_dmfm, 4))

# ==============================================================================
# 7. PREPARE DATA FOR PLOT
# ==============================================================================
df_plot <- data.frame(
  Date        = DatesQ_real,
  Quarter     = factor(DatesQ_30, levels = DatesQ_30),
  TrueGDP     = gdp_true_ea,
  NowcastGDP  = gdp_array[, month_index]
)

df_long <- pivot_longer(df_plot, cols = c(TrueGDP, NowcastGDP),
                        names_to = "Type", values_to = "Value")

# ==============================================================================
# 8. PLOT FULL TIME SERIES WITH COVID HIGHLIGHT
# ==============================================================================
covid_start_q <- get_quarter(P$covid_start)
covid_end_q   <- get_quarter(P$covid_end)

g_full <- ggplot() +
  annotate("rect", xmin = covid_start_q, xmax = covid_end_q,
           ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.3) +
  geom_col(data = df_long, aes(x = Quarter, y = Value, fill = Type),
           position = "dodge", width = 0.9, alpha = 0.5, show.legend = FALSE) +
  geom_line(data = df_plot, aes(x = Quarter, y = TrueGDP, color = "True GDP", group = 1), size = 1.1) +
  geom_line(data = df_plot, aes(x = Quarter, y = NowcastGDP, color = "Nowcast GDP", group = 1),
            size = 1.1, linetype = "dashed") +
  scale_color_manual(values = c("True GDP" = "#1f77b4", "Nowcast GDP" = "#ff7f0e")) +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  labs(
    title = paste("GDP Nowcasting for", country_code, "-", month_labels[month_index]),
    subtitle = "True vs Nowcasted GDP with COVID period highlighted",
    x = "Quarter", y = "GDP Growth", color = "Model"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ==============================================================================
# 9. SPLIT PRE/POST COVID
# ==============================================================================
df_pre  <- df_plot[df_plot$Date < P$covid_start, ]
df_post <- df_plot[df_plot$Date > P$covid_end, ]
df_long_pre  <- df_long[df_long$Date < P$covid_start, ]
df_long_post <- df_long[df_long$Date > P$covid_end, ]

plot_gdp_panel <- function(df_line, df_bar, title) {
  ggplot() +
    geom_col(data = df_bar, aes(x = Quarter, y = Value, fill = Type),
             position = "dodge", width = 0.9, alpha = 0.5, show.legend = FALSE) +
    geom_line(data = df_line, aes(x = Quarter, y = TrueGDP, group = 1, color = "True GDP"), size = 1.1) +
    geom_line(data = df_line, aes(x = Quarter, y = NowcastGDP, group = 1, color = "Nowcast GDP"),
              size = 1.1, linetype = "dashed") +
    scale_color_manual(values = c("True GDP" = "#1f77b4", "Nowcast GDP" = "#ff7f0e")) +
    coord_cartesian(ylim = c(-1.5, 1.5)) +
    labs(title = title, x = "Quarter", y = "GDP Growth", color = "Model") +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

g_pre  <- plot_gdp_panel(df_pre, df_long_pre, paste("DMFM Nowcasting -", month_labels[month_index], "(Pre-COVID)"))
g_post <- plot_gdp_panel(df_post, df_long_post, paste("DMFM Nowcasting -", month_labels[month_index], "(Post-COVID)"))

# ==============================================================================
# 10. VISUALIZE
# ==============================================================================
g_full
g_pre / g_post

# ==============================================================================
# 11. BOXPLOT — Forecast Errors per Modello e Periodo
# ==============================================================================

df_errors_dmfm <- do.call(rbind, lapply(1:n_months, function(m) {
  data.frame(
    Quarter = dimnames(gdp_array)[[1]],
    Error   = gdp_array[, m] - gdp_true_ea,
    Model   = paste0("DMFM_M", m)
  )
}))

# Aggiungi date e classifica periodo COVID
df_errors_dmfm$Date <- DatesM[match(df_errors_dmfm$Quarter, sapply(DatesM, get_quarter))]
df_errors_dmfm$CovidPeriod <- ifelse(df_errors_dmfm$Date < P$covid_start, "Pre-COVID",
                                     ifelse(df_errors_dmfm$Date > P$covid_end, "Post-COVID", NA))
df_errors_dmfm <- na.omit(df_errors_dmfm)

# BOXPLOT ERRORI
g_errors_dmfm <- ggplot(df_errors_dmfm, aes(x = Model, y = Error, fill = CovidPeriod)) +
  geom_boxplot(alpha = 0.85, outlier.shape = 16, outlier.size = 1.5) +
  scale_fill_manual(values = c("Pre-COVID" = "#1f77b4", "Post-COVID" = "#ff7f0e")) +
  labs(
    title = "Forecast Errors (DMFM): M1, M2, M3",
    subtitle = "Separated by COVID Period",
    x = "Model", y = "Error (Nowcast - True)",
    fill = "Period"
  ) +
  theme_minimal(base_size = 13)

# ==============================================================================
# 12. BOXPLOT — RMSFE distribuiti per Periodo (2 boxplot: Pre/Post)
# ==============================================================================

df_rmsfe_dmfm <- do.call(rbind, lapply(1:n_months, function(m) {
  data.frame(
    Quarter = dimnames(gdp_array)[[1]],
    RMSFE   = sqrt((gdp_array[, m] - gdp_true_ea)^2),
    Model   = paste0("M", m)
  )
}))

df_rmsfe_dmfm$Date <- DatesM[match(df_rmsfe_dmfm$Quarter, sapply(DatesM, get_quarter))]
df_rmsfe_dmfm$Period <- ifelse(df_rmsfe_dmfm$Date < P$covid_start, "Pre-COVID",
                               ifelse(df_rmsfe_dmfm$Date > P$covid_end, "Post-COVID", NA))
df_rmsfe_dmfm <- na.omit(df_rmsfe_dmfm)

# BOXPLOT RMSFE
g_rmsfe_dmfm <- ggplot(df_rmsfe_dmfm, aes(x = Period, y = RMSFE, fill = Period)) +
  geom_boxplot(alpha = 0.85, outlier.shape = 16, outlier.size = 1.5) +
  geom_jitter(width = 0.1, size = 2.5, aes(color = Model), show.legend = TRUE) +
  scale_fill_manual(values = c("Pre-COVID" = "#1f77b4", "Post-COVID" = "#ff7f0e")) +
  scale_color_manual(values = c("M1" = "#1f77b4", "M2" = "#2ca02c", "M3" = "#d62728")) +
  labs(
    title = "RMSFE Distribution (DMFM)",
    subtitle = "Models M1, M2, M3 across Pre- and Post-COVID",
    x = "Period", y = "RMSFE"
  ) +
  theme_minimal(base_size = 13)

# ==============================================================================
# 13. VISUALIZE
# ==============================================================================
g_errors_dmfm
g_rmsfe_dmfm

# ==============================================================================
# 14. STAMPA RMSFE per M1, M2, M3 — Pre e Post COVID (DMFM)
# ==============================================================================

# Funzione di calcolo RMSFE per un modello e un periodo
compute_rmsfe <- function(errors, date_vec, period = "Pre-COVID") {
  idx <- if (period == "Pre-COVID") {
    which(date_vec < P$covid_start)
  } else if (period == "Post-COVID") {
    which(date_vec > P$covid_end)
  } else {
    stop("Periodo non valido")
  }
  sqrt(mean(errors[idx]^2, na.rm = TRUE))
}

# Calcolo RMSFE per ogni modello e periodo
cat("=== RMSFE (DMFM) by Model and Period ===\n")
for (m in 1:n_months) {
  pred <- gdp_array[, m]
  error <- pred - gdp_true_ea
  date_vec <- DatesM[match(dimnames(gdp_array)[[1]], sapply(DatesM, get_quarter))]
  
  rmsfe_pre  <- compute_rmsfe(error, date_vec, period = "Pre-COVID")
  rmsfe_post <- compute_rmsfe(error, date_vec, period = "Post-COVID")
  
  cat(paste0("M", m, ":  Pre-COVID RMSFE = ", round(rmsfe_pre, 4),
             "   |   Post-COVID RMSFE = ", round(rmsfe_post, 4), "\n"))
}

# ==============================================================================
# 15. DFM — Preparazione dati
# ==============================================================================

# Ritaglia le stime per 30 trimestri e crea il dataframe
DatesQ      <- unique(sapply(DatesM, get_quarter))
DatesQ_30   <- DatesQ[69:98]
DateQQ_30   <- DateQQ_dt[1:30]

dfm_df <- data.frame(
  Quarter = factor(DatesQ_30, levels = DatesQ_30),
  Date    = DateQQ_30,
  TrueGDP = TrueQQ[1:30],
  M1 = FcstQQ[1:30, 4],
  M2 = FcstQQ[1:30, 5],
  M3 = FcstQQ[1:30, 6]
)

# ==============================================================================
# 16. DFM — Calcolo RMSFE globale per M1, M2, M3
# ==============================================================================
rmsfe_dfm <- sapply(1:3, function(m) {
  pred <- dfm_df[[paste0("M", m)]]
  sqrt(mean((pred - dfm_df$TrueGDP)^2, na.rm = TRUE))
})
names(rmsfe_dfm) <- paste0("M", 1:3)
print(round(rmsfe_dfm, 4))

# ==============================================================================
# 17. DFM — Grafico completo con COVID evidenziato
# ==============================================================================
df_plot_dfm <- data.frame(
  Quarter    = dfm_df$Quarter,
  Date       = dfm_df$Date,
  TrueGDP    = dfm_df$TrueGDP,
  NowcastGDP = dfm_df[[paste0("M", month_index)]]
)

df_long_dfm <- pivot_longer(df_plot_dfm, cols = c("TrueGDP", "NowcastGDP"),
                            names_to = "Type", values_to = "Value")

g_full_dfm <- ggplot() +
  annotate("rect", xmin = get_quarter(P$covid_start), xmax = get_quarter(P$covid_end),
           ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.3) +
  geom_col(data = df_long_dfm, aes(x = Quarter, y = Value, fill = Type),
           position = "dodge", width = 0.9, alpha = 0.5, show.legend = FALSE) +
  geom_line(data = df_plot_dfm, aes(x = Quarter, y = TrueGDP, color = "True GDP", group = 1), size = 1.1) +
  geom_line(data = df_plot_dfm, aes(x = Quarter, y = NowcastGDP, color = "Nowcast GDP", group = 1),
            size = 1.1, linetype = "dashed") +
  scale_color_manual(values = c("True GDP" = "#1f77b4", "Nowcast GDP" = "#ff7f0e")) +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  labs(
    title = paste("DFM GDP Nowcasting -", month_labels[month_index]),
    subtitle = "True vs Nowcasted GDP with COVID period highlighted",
    x = "Quarter", y = "GDP Growth", color = "Model"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ==============================================================================
# 18. DFM — Grafici separati pre/post COVID
# ==============================================================================
df_pre_dfm  <- df_plot_dfm[df_plot_dfm$Date < P$covid_start, ]
df_post_dfm <- df_plot_dfm[df_plot_dfm$Date > P$covid_end, ]
df_long_pre_dfm  <- df_long_dfm[df_long_dfm$Date < P$covid_start, ]
df_long_post_dfm <- df_long_dfm[df_long_dfm$Date > P$covid_end, ]

plot_dfm_panel <- function(df_line, df_bar, title) {
  ggplot() +
    geom_col(data = df_bar, aes(x = Quarter, y = Value, fill = Type),
             position = "dodge", width = 0.9, alpha = 0.5, show.legend = FALSE) +
    geom_line(data = df_line, aes(x = Quarter, y = TrueGDP, group = 1, color = "True GDP"), size = 1.1) +
    geom_line(data = df_line, aes(x = Quarter, y = NowcastGDP, group = 1, color = "Nowcast GDP"),
              size = 1.1, linetype = "dashed") +
    scale_color_manual(values = c("True GDP" = "#1f77b4", "Nowcast GDP" = "#ff7f0e")) +
    coord_cartesian(ylim = c(-1.5, 1.5)) +
    labs(title = title, x = "Quarter", y = "GDP Growth", color = "Model") +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

g_pre_dfm  <- plot_dfm_panel(df_pre_dfm, df_long_pre_dfm,  paste("DFM Nowcasting -", month_labels[month_index], "(Pre-COVID)"))
g_post_dfm <- plot_dfm_panel(df_post_dfm, df_long_post_dfm, paste("DFM Nowcasting -", month_labels[month_index], "(Post-COVID)"))

# ==============================================================================
# 19. DFM — Boxplot Errori M1, M2, M3
# ==============================================================================
df_errors_dfm <- do.call(rbind, lapply(1:3, function(m) {
  data.frame(
    Quarter = dfm_df$Quarter,
    Date    = dfm_df$Date,
    Error   = dfm_df[[paste0("M", m)]] - dfm_df$TrueGDP,
    Model   = paste0("DFM_M", m)
  )
}))
df_errors_dfm$Period <- ifelse(df_errors_dfm$Date < P$covid_start, "Pre-COVID",
                               ifelse(df_errors_dfm$Date > P$covid_end, "Post-COVID", NA))
df_errors_dfm <- na.omit(df_errors_dfm)

g_errors_dfm <- ggplot(df_errors_dfm, aes(x = Model, y = Error, fill = Period)) +
  geom_boxplot(alpha = 0.85, outlier.shape = 16, outlier.size = 1.5) +
  scale_fill_manual(values = c("Pre-COVID" = "#1f77b4", "Post-COVID" = "#ff7f0e")) +
  labs(title = "Forecast Errors (DFM)", x = "Model", y = "Forecast Error", fill = "Period") +
  theme_minimal(base_size = 13)

# ==============================================================================
# 20. DFM — Boxplot RMSFE (2 box: Pre e Post COVID)
# ==============================================================================
df_rmsfe_dfm <- do.call(rbind, lapply(1:3, function(m) {
  data.frame(
    Quarter = dfm_df$Quarter,
    Date    = dfm_df$Date,
    RMSFE   = sqrt((dfm_df[[paste0("M", m)]] - dfm_df$TrueGDP)^2),
    Model   = paste0("M", m)
  )
}))
df_rmsfe_dfm$Period <- ifelse(df_rmsfe_dfm$Date < P$covid_start, "Pre-COVID",
                              ifelse(df_rmsfe_dfm$Date > P$covid_end, "Post-COVID", NA))
df_rmsfe_dfm <- na.omit(df_rmsfe_dfm)

g_rmsfe_dfm <- ggplot(df_rmsfe_dfm, aes(x = Period, y = RMSFE, fill = Period)) +
  geom_boxplot(width = 0.5, alpha = 0.85, outlier.shape = 16) +
  geom_jitter(width = 0.1, size = 3, aes(color = Model)) +
  scale_fill_manual(values = c("Pre-COVID" = "#1f77b4", "Post-COVID" = "#ff7f0e")) +
  scale_color_manual(values = c("M1" = "#1f77b4", "M2" = "#2ca02c", "M3" = "#d62728")) +
  labs(
    title = "RMSFE Distribution (DFM)",
    subtitle = "Models M1, M2, M3 — Pre/Post COVID",
    x = "Period", y = "RMSFE"
  ) +
  theme_minimal(base_size = 13)

# ==============================================================================
# 21. VISUALIZE
# ==============================================================================
g_full_dfm
g_pre_dfm / g_post_dfm
g_errors_dfm
g_rmsfe_dfm

# ==============================================================================
# 22. DMFM vs DFM — Serie storica + Boxplot Errori
# ==============================================================================

# STEP 1 — Prepara la struttura dati combinata per confronto
month_label <- paste0("M", month_index)

df_compare <- data.frame(
  Quarter   = factor(DatesQ_30, levels = DatesQ_30),
  Date      = DatesM[match(DatesQ_30, sapply(DatesM, get_quarter))],
  TrueGDP   = gdp_true_ea,
  DMFM      = gdp_array[, month_index],
  DFM       = FcstQQ[1:30, 3 + month_index]
)

# STEP 2 — Grafico True vs DMFM vs DFM (linee sovrapposte)
df_long_all <- pivot_longer(df_compare, cols = c("TrueGDP", "DMFM", "DFM"),
                            names_to = "Model", values_to = "Value")

g_dmfm_dfm_series <- ggplot(df_long_all, aes(x = Quarter, y = Value, color = Model, group = Model)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c(
    "TrueGDP" = "black",
    "DMFM"    = "#1f77b4",  # blu
    "DFM"     = "#ff7f0e"   # arancio
  )) +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  labs(
    title = paste("Nowcasting Comparison: DMFM vs DFM -", month_label),
    subtitle = "True GDP vs Nowcasted GDP (30 Quarters)",
    x = "Quarter", y = "GDP Growth", color = "Series"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# STEP 3 — Calcolo errori + classificazione periodo
df_errors_compare <- data.frame(
  Quarter    = df_compare$Quarter,
  Date       = df_compare$Date,
  Error_DMFM = df_compare$DMFM - df_compare$TrueGDP,
  Error_DFM  = df_compare$DFM - df_compare$TrueGDP
)

df_errors_compare$Period <- ifelse(df_errors_compare$Date < P$covid_start, "Pre-COVID",
                                   ifelse(df_errors_compare$Date > P$covid_end, "Post-COVID", NA))

df_long_errors <- pivot_longer(df_errors_compare,
                               cols = c("Error_DMFM", "Error_DFM"),
                               names_to = "Model", values_to = "Error")
df_long_errors <- na.omit(df_long_errors)

# Rinomina modelli per legenda
df_long_errors$Model <- factor(df_long_errors$Model,
                               levels = c("Error_DMFM", "Error_DFM"),
                               labels = c("DMFM", "DFM"))

# STEP 4 — Boxplot degli errori: DMFM vs DFM per Pre/Post COVID
g_dmfm_dfm_errors <- ggplot(df_long_errors, aes(x = Model, y = Error, fill = Period)) +
  geom_boxplot(alpha = 0.85, outlier.shape = 16, outlier.size = 1.5, width = 0.6) +
  scale_fill_manual(values = c("Pre-COVID" = "#1f77b4", "Post-COVID" = "#ff7f0e")) +
  labs(
    title = paste("Forecast Errors (", month_label, "): DMFM vs DFM"),
    subtitle = "Pre- and Post-COVID Comparison",
    x = "Model", y = "Forecast Error (Nowcast - True)",
    fill = "Period"
  ) +
  theme_minimal(base_size = 13)

# ==============================================================================
# 23. VISUALIZZAZIONE COMBINATA (opzionale)
# ==============================================================================
g_dmfm_dfm_series       # Serie storica a linee
g_dmfm_dfm_errors       # Boxplot confronto errori




# ==============================================================================
# 24. Serie storica: True GDP vs Nowcast M1-M3 (DMFM vs DFM)
# ==============================================================================

# Costruzione base: 30 trimestri
n_qtr <- length(DatesQ_30)
true_gdp <- TrueQQ[1:n_qtr]
quarter_vec <- factor(DatesQ_30, levels = DatesQ_30)

# ---- DMFM ----
df_dmfm_long <- do.call(rbind, lapply(1:3, function(m) {
  data.frame(
    Quarter = quarter_vec,
    Value   = gdp_array[, m],
    Month   = paste0("M", m),
    Model   = "DMFM"
  )
}))

# ---- DFM ----
df_dfm_long <- do.call(rbind, lapply(1:3, function(m) {
  data.frame(
    Quarter = quarter_vec,
    Value   = FcstQQ[1:n_qtr, 3 + m],
    Month   = paste0("M", m),
    Model   = "DFM"
  )
}))

# ---- TRUE GDP (uguale per entrambi) ----
df_true <- data.frame(
  Quarter = rep(quarter_vec, 2),
  Value   = rep(true_gdp, 2),
  Month   = "True",
  Model   = rep(c("DMFM", "DFM"), each = n_qtr)
)

# ---- Combina tutto ----
df_all <- rbind(df_dmfm_long, df_dfm_long, df_true)
df_all$Month <- factor(df_all$Month, levels = c("True", "M1", "M2", "M3"))



library(ggplot2)

ggplot(df_all, aes(x = Quarter, y = Value, color = Month, group = Month)) +
  geom_line(size = 1.2) +
  facet_wrap(~ Model, ncol = 1) +
  scale_color_manual(values = c(
    "True" = "black",
    "M1"   = "#1f77b4",
    "M2"   = "#2ca02c",
    "M3"   = "#d62728"
  )) +
  labs(
    title = "True GDP vs Nowcast Updates (M1–M3) — DMFM vs DFM",
    subtitle = "30 Quarters – Comparison by Month of Forecast",
    x = "Quarter", y = "GDP Growth",
    color = "Series"
  ) +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ==============================================================================
# 25. Serie continua Nowcast: M1 → M2 → M3 concatenati
# ==============================================================================

# 1. Vettore temporale esteso per i mesi (1 punto per ciascun M1, M2, M3)
# Genera le "date" artificiali per ogni mese all'interno del trimestre
Dates_extended <- rep(DatesQ_30, each = 3)
Month_extended <- rep(c("M1", "M2", "M3"), times = length(DatesQ_30))

# 2. Serie nowcast concatenate
nowcast_dmfm <- as.vector(t(gdp_array[, 1:3]))         # riga per riga: M1, M2, M3, ...
nowcast_dfm  <- as.vector(t(FcstQQ[1:30, 4:6]))

# 3. Crea dataframe completo long
df_nowcast_long <- data.frame(
  Quarter     = rep(DatesQ_30, each = 3),
  Month       = Month_extended,
  SeriesMonth = paste0(rep(DatesQ_30, each = 3), "_", Month_extended),
  DMFM        = nowcast_dmfm,
  DFM         = nowcast_dfm
)

# 4. Asse temporale "artificiale" sequenziale
df_nowcast_long$TimeIndex <- 1:nrow(df_nowcast_long)

# Serie GDP vera ogni 3 mesi (da posizioni 2, 5, 8, ...)
true_gdp_series <- data.frame(
  TimeIndex = seq(2, by = 3, length.out = length(TrueQQ[1:30])),
  GDP       = TrueQQ[1:30]
)

library(ggplot2)

# Long format
df_plot_long <- pivot_longer(df_nowcast_long, cols = c("DMFM", "DFM"),
                             names_to = "Model", values_to = "Nowcast")

ggplot() +
  geom_line(data = df_plot_long, aes(x = TimeIndex, y = Nowcast, color = Model), size = 1.1) +
  geom_point(data = true_gdp_series, aes(x = TimeIndex, y = GDP), color = "black", size = 2.5, shape = 16) +
  geom_line(data = true_gdp_series, aes(x = TimeIndex, y = GDP), color = "black", size = 1, linetype = "dotted") +
  scale_color_manual(values = c("DMFM" = "#1f77b4", "DFM" = "#ff7f0e")) +
  labs(
    title = "Nowcast Evolution (M1 → M2 → M3) vs True Quarterly GDP",
    subtitle = "DMFM and DFM nowcasts joined in monthly sequence",
    x = "Nowcast Month Index", y = "GDP Growth",
    color = "Model"
  ) +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  theme_minimal(base_size = 13)

# ==============================================================================
# 26. Last trimester convergence
# ==============================================================================

# Trova l'ultimo trimestre per cui abbiamo il valore reale del GDP
last_q_idx <- max(which(!is.na(TrueQQ)))
last_qtr   <- DatesQ_30[last_q_idx]
true_val   <- TrueQQ[last_q_idx]

dmfm_vals <- gdp_array[last_q_idx, 1:3]
dfm_vals  <- FcstQQ[last_q_idx, 4:6]

df_progress <- data.frame(
  Month  = factor(c("M1", "M2", "M3", "Release"), levels = c("M1", "M2", "M3", "Release")),
  DMFM   = c(dmfm_vals, true_val),
  DFM    = c(dfm_vals, true_val)
)

df_long <- tidyr::pivot_longer(df_progress, cols = c("DMFM", "DFM"),
                               names_to = "Model", values_to = "Value")

ggplot(df_long, aes(x = Month, y = Value, group = Model, color = Model)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("DMFM" = "#1f77b4", "DFM" = "#ff7f0e")) +
  geom_hline(yintercept = true_val, linetype = "dotted", color = "black") +
  labs(
    title = paste("Nowcast Evolution for Last Quarter:", last_qtr),
    subtitle = "From M1 to M3 and Actual Release",
    x = "Month", y = "GDP Growth", color = "Model"
  ) +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  theme_minimal(base_size = 13)

# Costruzione dati nowcast + release
df_progress <- data.frame(
  Month  = factor(c("M1", "M2", "M3", "Release"), levels = c("M1", "M2", "M3", "Release")),
  DMFM   = c(dmfm_vals, true_val),
  DFM    = c(dfm_vals, true_val),
  LineGroup = c("Forecast", "Forecast", "Forecast", NA)  # Release = NA → non collegata
)

# Long format con gruppo per le linee
df_long <- tidyr::pivot_longer(df_progress, cols = c("DMFM", "DFM"),
                               names_to = "Model", values_to = "Value")

df_long$Group <- paste0(df_long$Model, "_", df_progress$LineGroup)

# Plot
ggplot(df_long, aes(x = Month, y = Value, group = Group, color = Model)) +
  geom_line(data = subset(df_long, Month != "Release"), size = 1.2) +  # linee solo su M1–M3
  geom_point(size = 3) +
  geom_hline(yintercept = true_val, linetype = "dotted", color = "black") +
  scale_color_manual(values = c("DMFM" = "#1f77b4", "DFM" = "#ff7f0e")) +
  labs(
    title = paste("Nowcast Evolution for Last Quarter:", last_qtr),
    subtitle = "M1–M3 connected, Release shown as point only",
    x = "Month", y = "GDP Growth", color = "Model"
  ) +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  theme_minimal(base_size = 13)

# Errore nel tempo 
df_errors <- data.frame(
  Month = factor(c("M1", "M2", "M3"), levels = c("M1", "M2", "M3")),
  DMFM_Error = dmfm_vals - true_val,
  DFM_Error  = dfm_vals - true_val
)

df_errors_long <- tidyr::pivot_longer(df_errors, cols = c("DMFM_Error", "DFM_Error"),
                                      names_to = "Model", values_to = "Error")

ggplot(df_errors_long, aes(x = Month, y = Error, group = Model, color = Model)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("DMFM_Error" = "#1f77b4", "DFM_Error" = "#ff7f0e")) +
  labs(
    title = paste("Forecast Error Convergence –", last_qtr),
    subtitle = "Nowcast error = forecast - true GDP",
    x = "Month", y = "Forecast Error", color = "Model"
  ) +
  theme_minimal(base_size = 13)



# ==============================================================================
# 27. FACTORS AND LOADINGS
# ==============================================================================

l.R_str <- inputs$R             # Raw loadings
l.C_str <- inputs$C             # Column loadings

plot_factors_tensor(inputs$fs)


# ==============================================================================
# 28. Last Run
# ==============================================================================
# Plot convergence of log-likelihood
plot_em_llk(out$history)


# Smoothing

country_code <- "IT"  # Paese da visualizzare
GDP_varname  <- "GDP"

# ==============================================================================
# STEP 1 — SELEZIONE VARIABILI E COSTRUZIONE DATASET
# ==============================================================================

# 1.1 Variabili più correlate al GDP per paese
base_names <- select_base_variable_names(countries, P)

# 1.2 Filtra variabili comuni a tutti i paesi
common_vars <- filter_common_variables_across_countries(
  countries, base_names$monthly, base_names$quarterly
)

# 1.3 Base name delle variabili comuni
all_vars <- names(common_vars$monthly[[1]])

# 1.4 Prepara i dati per ogni paese (COVID = FALSE per evitare NaN)
country_list <- list()
for (cc in countries) {
  res <- prepare_country_data(
    cc, P,
    selected_vars_m = common_vars$monthly,
    selected_vars_q = common_vars$quarterly,
    Covid = FALSE
  )
  country_list[[cc]] <- list(
    Data   = res$Data,
    Dates  = res$DatesM,
    Series = res$Series
  )
}

# ==============================================================================
# STEP 2 — COSTRUZIONE TENSORE CON COVID INCLUSO
# ==============================================================================

tensor_covid <- tensor_data(countries, P)  # contiene i NaN dovuti al COVID
Y  <- tensor_covid$Y
W  <- tensor_covid$W

# Standardizzazione (media e sd per la destandardizzazione)
std <- standardize_Y(Y)
Y_std <- std$Y_scaled

# ==============================================================================
# STEP 3 — SMOOTHING E PREPARAZIONE SERIE GDP
# ==============================================================================

# 3.1 Stima smussata (fatta altrove → qui la carichi)
Y_smoothed <- inverse_standardize_Y(Y_scaled = out$model$Y, mean_Y = std$mean, sd_Y = std$sd)

# 3.2 Serie reale con COVID incluso
Y_true_real <- Y  # già destandardizzato

# 3.3 Estrai date e indici
Dates     <- dimnames(Y)[[1]]
Dates_fixed <- as.Date(paste0(Dates, "-01"))  # "YYYY-MM" → Date valida

ea_index  <- which(dimnames(Y)[[2]] == country_code)
GDP_index <- which(dimnames(Y)[[3]] == GDP_varname)

# 3.4 Estrai le due serie
gdp_true     <- Y_true_real[, ea_index, GDP_index]
gdp_smoothed <- Y_smoothed[, ea_index, GDP_index]

# ==============================================================================
# STEP 4 — GRAFICO: GDP OSSERVATO vs STIMA SMOOTHED
# ==============================================================================

library(tidyr)
library(ggplot2)

# Costruisci dataframe
df_comp <- data.frame(
  Date         = Dates_fixed,
  GDP_True     = gdp_true,
  GDP_Smoothed = gdp_smoothed
)

# Long format
df_long <- pivot_longer(df_comp, cols = c("GDP_True", "GDP_Smoothed"),
                        names_to = "Series", values_to = "Value")

df_long_clean <- na.omit(df_long)

ggplot(df_long_clean, aes(x = Date, y = Value, color = Series)) +
  geom_line(size = 1) +
  labs(
    title = paste("GDP: True vs Smoothed –", country_code),
    x = "Date", y = "GDP Growth", color = "Series"
  ) +
  theme_minimal(base_size = 13)


# ==============================================================================
# END OF SCRIPT
# ==============================================================================
