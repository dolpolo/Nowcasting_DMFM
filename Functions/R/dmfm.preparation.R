################################################################################
########################### DATASET PREPARATION ################################
################################################################################

# ==============================================================================
# STANDARDIZATION AND DE-STANDARDIZATION
# ==============================================================================


standardize_Y <- function(Y) {
  
  #' Standardize a 3D array along the time dimension
  #'
  #' Applies z-score standardization (mean 0, sd 1) independently to each variable
  #' across time in a 3D array Y of dimensions (T x p1 x p2). Missing values are ignored.
  #'
  #' @param Y A 3D array (T x p1 x p2) of time-series data with possible NA values.
  #' @return A list with:
  #'   \item{Y_scaled}{Standardized 3D array}
  #'   \item{mean}{(p1 x p2) matrix of means used}
  #'   \item{sd}{(p1 x p2) matrix of standard deviations used}

  T      <- dim(Y)[1]
  p1     <- dim(Y)[2]
  p2     <- dim(Y)[3]
  Y_std  <- array(NA, dim = dim(Y), dimnames = dimnames(Y))
  mean_Y <- matrix(0, p1, p2)
  sd_Y   <- matrix(1, p1, p2)
  
  for (i in 1:p1) {
    for (j in 1:p2) {
      y_ij <- Y[, i, j]
      mu   <- mean(y_ij, na.rm = TRUE)
      sigma <- sd(y_ij, na.rm = TRUE)
      if (is.na(sigma) || sigma == 0) sigma <- 1
      Y_std[, i, j] <- (y_ij - mu) / sigma
      mean_Y[i, j]  <- mu
      sd_Y[i, j]    <- sigma
    }
  }
  
  return(list(Y_scaled = Y_std, mean = mean_Y, sd = sd_Y))
}


# ******************************************************************************

inverse_standardize_Y <- function(Y_scaled, mean_Y, sd_Y) {
  
  #' Revert standardization of a 3D array
  #'
  #' Restores the original scale of a standardized 3D array using stored means and sds.
  #'
  #' @param Y_scaled Standardized 3D array
  #' @param mean_Y (p1 x p2) matrix of means
  #' @param sd_Y (p1 x p2) matrix of standard deviations
  #' @return Original scale 3D array

  T <- dim(Y_scaled)[1]
  p1 <- dim(Y_scaled)[2]
  p2 <- dim(Y_scaled)[3]
  
  Y_original <- array(NA, dim = dim(Y_scaled), dimnames = dimnames(Y_scaled))
  
  for (t in 1:T) {
    Y_original[t, , ] <- Y_scaled[t, , ] * sd_Y + mean_Y
  }
  
  return(Y_original)
}

# ******************************************************************************


center_Y <- function(Y) {
  
  #' Center a 3D array by subtracting the mean of each variable
  #'
  #' Removes the mean from each (i, j) variable across time (T) in a 3D array.
  #'
  #' @param Y A 3D array (T x p1 x p2)
  #' @return A list with:
  #'   \item{Y_centered}{Centered 3D array}
  #'   \item{mean}{(p1 x p2) matrix of means}
  
  T         <- dim(Y)[1]
  p1        <- dim(Y)[2]
  p2        <- dim(Y)[3]
  Y_centered <- array(NA, dim = dim(Y), dimnames = dimnames(Y))
  mean_Y     <- matrix(0, p1, p2)
  
  for (i in 1:p1) {
    for (j in 1:p2) {
      y_ij <- Y[, i, j]
      mu   <- mean(y_ij, na.rm = TRUE)
      Y_centered[, i, j] <- y_ij - mu
      mean_Y[i, j]       <- mu
    }
  }
  
  return(list(Y_centered = Y_centered, mean = mean_Y))
}

# ******************************************************************************


decenter_Y <- function(Y_centered, mean_Y) {
  
  #' Restore mean to a previously centered 3D array
  #'
  #' Adds back the mean to a centered 3D array to recover original data.
  #'
  #' @param Y_centered Centered 3D array
  #' @param mean_Y (p1 x p2) matrix of means
  #' @return 3D array restored to original mean

  T <- dim(Y_centered)[1]
  p1 <- dim(Y_centered)[2]
  p2 <- dim(Y_centered)[3]
  
  Y_original <- array(NA, dim = dim(Y_centered), dimnames = dimnames(Y_centered))
  
  for (t in 1:T) {
    Y_original[t, , ] <- Y_centered[t, , ] + mean_Y
  }
  
  return(Y_original)
}


# ==============================================================================
# COVID TREATMENT
# ==============================================================================

# COVID MASK FUNCTION
CovidNaN <- function(data, dates, class, start_covid, end_covid) {
  real_idx <- which(grepl("r", tolower(class)))
  mask <- dates >= as.Date(as.yearmon(start_covid)) & dates <= as.Date(as.yearmon(end_covid))
  covid_matrix <- matrix(FALSE, nrow = nrow(data), ncol = ncol(data))
  covid_matrix[mask, real_idx] <- TRUE
  data[covid_matrix] <- NA
  list(data, covid_matrix)
}

# ******************************************************************************

nan_percent_Y <- function(Y) {
  
  #' Compute percentage of missing values in a 3D array
  #'
  #' Calculates the proportion of NA values across the entire array.
  #'
  #' @param Y A 3D array with NA values
  #' @return Numeric: percentage of missing values
  
  total_values <- length(Y)
  total_NAs    <- sum(is.na(Y))
  percent      <- total_NAs / total_values * 100
  return(percent)
}

# ==============================================================================
# GENERAL PREPARATION OF DATA 
# ==============================================================================


select_base_variable_names <- function(countries, P) {
  
  # Selects the most GDP-correlated monthly and quarterly variables for each country.
  # Inputs:
  #   - countries: vector of country codes (e.g., c("IT", "DE", "FR"))
  #   - P: list of parameters including model size (P$modelM, P$modelQ)
  # Output:
  #   - A list with two elements:
  #       $monthly: vector of base names of selected monthly variables (common across countries)
  #       $quarterly: vector of base names of selected quarterly variables (common across countries)
  
  modelM <- P$modelM
  modelQ <- P$modelQ
  
  model_sizes <- list(
    small = list(Q = 1, M = 5),  # small: solo GDP tra le trimestrali
    medium = list(Q = 5, M = 10),
    large = list(Q = 10, M = 50)
  )
  
  n_m <- model_sizes[[tolower(modelM)]]$M
  n_q <- model_sizes[[tolower(modelQ)]]$Q
  
  selected_m_names <- c()
  selected_q_names <- c()
  
  for (cc in countries) {
    DataFile <- paste0("Data/", cc, "/Processed/Data_", cc, ".xlsx")
    monthly_data <- readxl::read_excel(DataFile, sheet = "MonthlyLong")
    quarterly_data <- readxl::read_excel(DataFile, sheet = "QuarterlyLong")
    
    DataM <- monthly_data[, -1]
    DataQ <- quarterly_data[, -1]
    
    gdp_q_idx <- which(grepl("^GDP_", colnames(DataQ), ignore.case = TRUE))[1]
    gdp_q <- DataQ[[gdp_q_idx]]
    gdp_rep <- rep(gdp_q, each = 3)[1:nrow(DataM)]
    
    # Correlazioni mensili
    cors_m <- apply(DataM, 2, function(x) suppressWarnings(cor(x, gdp_rep, use = "pairwise.complete.obs")))
    top_m <- names(cors_m)[order(abs(cors_m), decreasing = TRUE)[1:min(n_m, length(cors_m))]]
    
    # Correlazioni trimestrali (solo se medium/large)
    top_q <- NULL
    if (n_q > 1) {
      cors_q <- apply(DataQ[, -gdp_q_idx], 2, function(x) suppressWarnings(cor(x, gdp_q, use = "pairwise.complete.obs")))
      top_q <- names(DataQ)[-gdp_q_idx][order(abs(cors_q), decreasing = TRUE)[1:min(n_q, length(cors_q))]]
    }
    
    # Rimozione suffisso paese
    base_m <- gsub(paste0("_", cc, "$"), "", top_m)
    base_q <- if (!is.null(top_q)) gsub(paste0("_", cc, "$"), "", top_q) else character()
    
    selected_m_names <- union(selected_m_names, base_m)
    selected_q_names <- union(selected_q_names, base_q)
  }
  
  # Per modello small: trimestrali = solo GDP
  if (n_q == 1) {
    selected_q_names <- "GDP"
  }
  
  return(list(monthly = selected_m_names, quarterly = selected_q_names))
}

# ******************************************************************************

filter_common_variables_across_countries <- function(countries, selected_base_m,
                                                     selected_base_q) {
  
  # Filters selected variable names, retaining only those present in all countries.
  # Inputs:
  #   - countries: vector of country codes
  #   - selected_base_m: vector of selected base monthly variable names
  #   - selected_base_q: vector of selected base quarterly variable names
  # Output:
  #   - A list with two elements:
  #       $monthly: list of named vectors (per country) with full monthly variable names
  #       $quarterly: list of named vectors (per country) with full quarterly variable names
  
  valid_m_vars <- list()
  valid_q_vars <- list()
  
  for (cc in countries) {
    DataFile <- paste0("Data/", cc, "/Processed/Data_", cc, ".xlsx")
    monthly_data <- readxl::read_excel(DataFile, sheet = "MonthlyLong")
    quarterly_data <- readxl::read_excel(DataFile, sheet = "QuarterlyLong")
    
    DataM <- monthly_data[, -1]
    DataQ <- quarterly_data[, -1]
    
    col_m <- colnames(DataM)
    col_q <- colnames(DataQ)
    
    # Rimuove suffisso country e costruisce mappa
    base_m <- gsub(paste0("_", cc, "$"), "", col_m)
    base_q <- gsub(paste0("_", cc, "$"), "", col_q)
    
    matched_m <- selected_base_m[selected_base_m %in% base_m]
    matched_q <- selected_base_q[selected_base_q %in% base_q]
    
    # Mapping: base → full
    matched_full_m <- setNames(paste0(matched_m, "_", cc), matched_m)
    matched_full_q <- setNames(paste0(matched_q, "_", cc), matched_q)
    
    valid_m_vars[[cc]] <- matched_full_m
    valid_q_vars[[cc]] <- matched_full_q
  }
  
  # Trova solo i base_name presenti in tutti i paesi
  common_base_m <- Reduce(intersect, lapply(valid_m_vars, names))
  common_base_q <- Reduce(intersect, lapply(valid_q_vars, names))
  
  # Per ogni paese, ritorna solo le variabili comuni
  filtered_vars_m <- lapply(valid_m_vars, function(x) x[common_base_m])
  filtered_vars_q <- lapply(valid_q_vars, function(x) x[common_base_q])
  
  return(list(monthly = filtered_vars_m, quarterly = filtered_vars_q))
}

# ******************************************************************************

prepare_country_data <- function(country_code, P, selected_vars_m, 
                                 selected_vars_q, Covid = FALSE) {
  
  # Loads, transforms, and harmonizes monthly and quarterly macro data for a single country.
  # Inputs:
  #   - country_code: string with the country code (e.g., "IT")
  #   - P: list of parameters for transformation and filtering
  #   - selected_vars_m: named vector of monthly variable names to keep (for that country)
  #   - selected_vars_q: named vector of quarterly variable names to keep (for that country)
  #   - Covid: if TRUE, real variables affected by the pandemic shock are masked (defult = FALSE)
  # Output:
  #   - A list with two matrices:
  #       $DataM: transformed monthly data (T_months x n_m)
  #       $DataQ: transformed quarterly data aligned to monthly frequency (T_months x n_q)
  
  
  # ---- PATHS ----
  DataFile <- paste0("Data/", country_code, "/Processed/Data_", country_code, ".xlsx")
  LegendFile <- paste0("Data/", country_code, "/Original/Legend_", country_code, ".xlsx")
  
  # ---- MONTHLY DATA ----
  monthly_data <- readxl::read_excel(DataFile, sheet = "MonthlyLong", col_names = TRUE)
  DatesM <- as.Date(monthly_data[[1]])
  DataM <- monthly_data[, -1]
  vars_m_country <- selected_vars_m[[country_code]]
  DataM <- DataM[, vars_m_country]
  
  # Monthly Legend
  legend_m <- readxl::read_excel(LegendFile, sheet = "MDescriptionFull")
  idx_m <- match(vars_m_country, legend_m[[1]])
  ClassM <- legend_m$Class[idx_m]
  TransfM <- floor(legend_m[idx_m, 8])
  GroupM <- legend_m$M1[idx_m]
  
  # Annualization for monthly
  DataMTrf <- DataM
  annualize_idx <- which(TransfM %in% c(2, 4))
  DataMTrf[, annualize_idx] <- DataMTrf[, annualize_idx] * 100
  
  if (Covid) {
    resM <- CovidNaN(DataMTrf, DatesM, ClassM, P$covid_start, P$covid_end)
    DataMTrf <- resM[[1]]
    Covid_obsM <- resM[[2]]
  } else {
    Covid_obsM <- matrix(FALSE, nrow = nrow(DataMTrf), ncol = ncol(DataMTrf))
  }
  T_M <- nrow(DataMTrf)
  
  # ---- QUARTERLY DATA ----
  quarterly_data <- readxl::read_excel(DataFile, sheet = "QuarterlyLong", col_names = TRUE)
  DatesQ <- as.Date(quarterly_data[[1]])
  DataQ <- quarterly_data[, -1]
  vars_q_country <- selected_vars_q[[country_code]]
  DataQ <- DataQ[, vars_q_country]
  
  # Quarterly Legend
  legend_q <- readxl::read_excel(LegendFile, sheet = "QDescriptionFull")
  idx_q <- match(vars_q_country, legend_q[[1]])
  ClassQ <- legend_q$Class[idx_q]
  TransfQ <- floor(legend_q[idx_q, 8])
  GroupQ <- legend_q$M1[idx_q]
  
  DataQTrf <- DataQ
  annualize_idx <- which(TransfQ %in% c(2, 4))
  DataQTrf[, annualize_idx] <- DataQTrf[, annualize_idx] * 100
  
  if (Covid) {
    resQ <- CovidNaN(DataQTrf, DatesQ, ClassQ, P$covid_start, P$covid_end)
    DataQTrf <- resQ[[1]]
    Covid_obsQ <- resQ[[2]]
  } else {
    Covid_obsQ <- matrix(FALSE, nrow = nrow(DataQTrf), ncol = ncol(DataQTrf))
  }
  
  # ---- ESPANSIONE TRIMESTRALE A MENSILE ----
  DataQMTrf <- do.call(rbind, lapply(1:nrow(DataQTrf), function(i) {
    rbind(rep(NA, ncol(DataQTrf)),
          rep(NA, ncol(DataQTrf)),
          DataQTrf[i, ])
  }))
  T_Q <- nrow(DataQMTrf)
  
  # ---- UNIFICA LUNGHEZZA ----
  T_final <- max(T_M, T_Q)
  if (T_M < T_final) {
    paddingM <- data.frame(matrix(NA, T_final - T_M, ncol(DataMTrf)))
    colnames(paddingM) <- colnames(DataMTrf)
    DataMTrf <- rbind(DataMTrf, paddingM)
    Covid_obsM <- rbind(Covid_obsM, matrix(FALSE, T_final - T_M, ncol(Covid_obsM)))
  }
  if (T_Q < T_final) {
    paddingQ <- data.frame(matrix(NA, T_final - T_Q, ncol(DataQMTrf)))
    colnames(paddingQ) <- colnames(DataQMTrf)
    DataQMTrf <- rbind(DataQMTrf, paddingQ)
    Covid_obsQ <- rbind(Covid_obsQ, matrix(FALSE, T_final - T_Q, ncol(Covid_obsQ)))
  }
  
  # ---- OUTPUT ----
  Data <- cbind(DataMTrf, DataQMTrf)
  Series <- c(vars_m_country, vars_q_country)
  Group <- c(GroupM, GroupQ)
  
  result <- list(
    Data = Data,
    DatesM = DatesM,
    DatesQ = DatesQ,
    Series = Series,
    UnbalancedPattern = NULL,
    Covid_obsM = Covid_obsM,
    Covid_obsQ = Covid_obsQ,
    Group = Group,
    Name = paste0("Data_", country_code),
    quarterly_start_idx = ncol(DataMTrf) + 1
  )
  
  assign(result$Name, result$Data, envir = .GlobalEnv)
  return(result)
}

# ==============================================================================
# REMOVE COUNTRY CODES FROM VARIABLES NAME
# ==============================================================================

remove_country_code <- function(names, ccodes) {
  
  # Removes country code suffixes (e.g., "_IT") from a vector of variable names.
  # Inputs:
  #   - names: vector of variable names (e.g., "GDP_IT", "CPI_DE")
  #   - ccodes: vector of country codes to be stripped (e.g., c("IT", "DE"))
  # Output:
  #   - A character vector with suffixes removed (e.g., "GDP", "CPI")
  
  pattern <- paste0("_(", paste(ccodes, collapse = "|"), ")$")
  gsub(pattern, "", names)
}

# ==============================================================================
# BUILD METRIX-VARIATE DATA
# ==============================================================================

tensor_data <- function(countries, P) {
  
  # Constructs a 3D tensor (array) of harmonized macro data across countries.
  # Inputs:
  #   - countries: vector of country codes
  #   - P: list of parameters including model size, transformation options, etc.
  # Output:
  #   - A 3D array (T x N x V) where:
  #       T = number of months,
  #       N = number of countries,
  #       V = number of harmonized variables
  #   - Missing values are preserved as NA
  
  T_max <- P$Tmax
  
  # Carica i dataset per ciascun paese
  data_list <- lapply(countries, function(cc) get(paste0("Data_", cc)))
  names(data_list) <- countries
  
  # Rimuove i codici paese dalle colonne
  cleaned_data_list <- list()
  for (i in seq_along(countries)) {
    cc <- countries[i]
    data <- data_list[[cc]]
    
    # Tronca a Tmax osservazioni se serve
    if (nrow(data) > T_max) {
      data <- data[1:T_max, ]
    }
    
    colnames(data) <- remove_country_code(colnames(data), countries)
    cleaned_data_list[[cc]] <- data
  }
  
  # Trova tutte le variabili comuni
  all_vars <- unique(unlist(lapply(cleaned_data_list, colnames)))
  
  # Crea tensor Y
  Y <- array(NA, dim = c(T_max, length(countries), length(all_vars)),
             dimnames = list(
               format(res$DatesM[1:T_max], "%Y-%m"),
               countries,
               all_vars
             ))
  
  # Allinea le variabili per ciascun paese
  for (i in seq_along(countries)) {
    cc <- countries[i]
    data <- cleaned_data_list[[cc]]
    
    # Nuova matrice con tutte le variabili, NA dove mancano
    aligned_data <- matrix(NA, nrow = T_max, ncol = length(all_vars))
    colnames(aligned_data) <- all_vars
    matched_vars <- intersect(colnames(data), all_vars)
    aligned_data[1:nrow(data), matched_vars] <- as.matrix(data[, matched_vars])
    
    Y[, i, ] <- aligned_data
  }
  
  # Crea matrice W (1 se osservato, 0 se NA)
  W <- ifelse(is.na(Y), 0, 1)
  
  return(list(Y = Y, W = W))
}
