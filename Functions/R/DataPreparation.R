################################################################################
############################ DESCRIPTIVE ANALYSIS ##############################
################################################################################

# ==============================================================================
# GDP COMMUNALITIES IN EA 
# ==============================================================================

GDP_communality <- function(country_codes) {
  
  all_gdp <- list()
  
  for (cc in country_codes) {
    data_path <- paste0("Data/", cc, "/Original/", cc, "DataQ_LT.xlsx")
    gdp_label <- paste0("GDP_", cc)
    
    df <- read_excel(data_path)
    
    # Special case for EA: compute PACF, ACF, correlation heatmap
    if (cc == "EA") {
      EAdataQ <- df
      
      # ---- Autocorrelation analysis for GDP variables ----
      target_variables <- dplyr::select(EAdataQ, dplyr::starts_with("GDP_"))
      
      for (var in names(target_variables)) {
        ts_var <- target_variables[[var]]
        if (is.numeric(ts_var)) {
          cat(paste("Autocorrelation for:", var, "\n"))
          acf(ts_var, main = paste("ACF for", var))
          pacf(ts_var, main = paste("PACF for", var))
        }
      }
      
      # ---- Correlation heatmap between variables (no labels inside cells) ----
      EA_corr_data <- dplyr::select(EAdataQ, -Time)
      correlation_matrix <- cor(EA_corr_data, use = "pairwise.complete.obs")
      melted_corr <- melt(correlation_matrix)
      
      p_heat <- ggplot(data = melted_corr, aes(x = Var1, y = Var2, fill = value)) +
        geom_tile(color = "white") +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                             midpoint = 0, limits = c(-1, 1))+
        theme_minimal() +
        labs(title = "EA Variables Correlation Heatmap", x = NULL, y = NULL) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid = element_blank())
      print(p_heat)
    }
    
    # ---- Extract GDP series and add to list ----
    gdp_idx <- which(tolower(colnames(df)) == tolower(gdp_label))
    
    gdp_series <- df[[gdp_idx]]*100
    date_series <- df[[1]]
    
    gdp_df <- data.frame(
      date = date_series,
      gdp = gdp_series,
      country = cc
    )
    
    all_gdp[[cc]] <- gdp_df
  }
  
  # ---- Combine all GDP series ----
  full_gdp <- bind_rows(all_gdp)
  
  # === Line Plot ===
  p1 <- ggplot(full_gdp, aes(x = date, y = gdp, color = country)) +
    geom_line(linewidth = 1) +
    labs(title = "GDP per Country Over Time", x = "Date", y = "GDP") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # === Facet Plot ===
  p2 <- ggplot(full_gdp, aes(x = date, y = gdp)) +
    geom_line(color = "steelblue") +
    facet_wrap(~ country, scales = "free_y") +
    labs(title = "GDP Faceted by Country", x = "Date", y = "GDP") +
    theme_minimal()
  
  print(p1)
  print(p2)
  
  # === Comovement Analysis ===
  wide_gdp <- full_gdp %>%
    filter(!is.na(gdp)) %>%
    pivot_wider(names_from = country, values_from = gdp)
  
  gdp_matrix <- as.matrix(wide_gdp[,-1])
  rownames(gdp_matrix) <- wide_gdp$date
  
  comov_matrix <- cor(gdp_matrix, use = "pairwise.complete.obs")
  comov_melt <- melt(comov_matrix)
  colnames(comov_melt) <- c("Country1", "Country2", "Correlation")
  
  p3 <- ggplot(comov_melt, aes(x = Country1, y = Country2, fill = Correlation)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0.5, limits = c(0, 1), name = "Correlation") +
    labs(title = "GDP Comovement Between Countries") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
  
  print(p3)
  
  assign("gdp_results", list(
    data = all_gdp,
    combined = full_gdp,
    comovement = comov_matrix
  ), envir = .GlobalEnv)
}

# ==============================================================================
# GDP BY COUNTRU CONTRIBUTION
# ==============================================================================


plot_euro_area_gdp <- function(file_path = "Data/GDP_millionsEA.xlsx",
                               sheet_name = "GDP",
                               quarter = "2024-Q4",
                               special_countries = c("Italy", "Germany", "France", "Spain", "Ireland",
                                                     "Greece", "Netherlands", "Belgium", "Austria"),
                               show_table = TRUE) {
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(sf)
  library(rnaturalearth)
  library(knitr)
  library(kableExtra)
  
  # 1. Importa e pulisci dati
  gdp_data_raw <- read_excel(file_path, sheet = sheet_name)
  
  gdp_data <- gdp_data_raw %>%
    dplyr::select(Country = 1, value = !!sym(quarter)) %>%
    mutate(
      Country = case_when(
        Country == "Germany (until 1990 former territory of the FRG)" ~ "Germany",
        Country == "Euro area" ~ NA_character_,
        TRUE ~ Country
      ),
      value = suppressWarnings(as.numeric(value))
    ) %>%
    filter(!is.na(Country), !Country %in% c(
      "European Union - 27 countries (from 2020)",
      "Euro area – 20 countries (from 2023)")
    )
  
  # 2. Mappa base
  world <- ne_countries(scale = "large", returnclass = "sf")
  europe <- world %>% filter(region_un == "Europe")
  
  euro_area_countries <- c(
    "Austria", "Belgium", "Croatia", "Cyprus", "Estonia", "Finland",
    "France", "Germany", "Greece", "Ireland", "Italy", "Latvia",
    "Lithuania", "Luxembourg", "Malta", "Netherlands", "Portugal",
    "Slovakia", "Slovenia", "Spain"
  )
  
  # 3. Calcolo share
  gdp_total_ea <- gdp_data %>%
    filter(Country %in% euro_area_countries) %>%
    summarise(total = sum(value, na.rm = TRUE)) %>%
    pull(total)
  
  gdp_data <- gdp_data %>%
    mutate(GDP_Share_EA = (value / gdp_total_ea) * 100)
  
  # 4. Join mappa + dati
  europe_gdp <- europe %>%
    left_join(gdp_data, by = c("name" = "Country")) %>%
    mutate(
      EuroArea = ifelse(name %in% euro_area_countries, "Euro Area", "Rest of Europe"),
      SpecialCountry = name %in% special_countries
    )
  
  # 5. Plot
  gdp_plot <- ggplot() +
    geom_sf(data = europe, fill = "grey95", color = "black", size = 0.2) +  # contorni neri
    geom_sf(data = europe_gdp %>% filter(EuroArea == "Euro Area"),
            aes(fill = value), color = "black", size = 0.2) +  # anche qui contorno nero sottile
    geom_sf(data = europe_gdp %>% filter(SpecialCountry),
            fill = NA, color = "red", size = 0.9) +
    geom_sf_text(data = europe_gdp %>% filter(SpecialCountry),
                 aes(label = paste0(round(GDP_Share_EA, 1), "%")),
                 size = 3.5, color = "black", fontface = "bold") +
    scale_fill_gradientn(
      colours = c("white", "#B3DFF2", "#4A90E2", "#003366"),
      na.value = "grey90",
      name = paste0("GDP (M€)\n", quarter)
    ) +
    coord_sf(xlim = c(-25, 45), ylim = c(33, 72), expand = FALSE) +
    theme_minimal() +
    labs(
      title = paste("GDP Contribution by Country in Euro Area -", quarter),
      caption = "Source: EUROSTAT, Natural Earth"
    ) +
    theme(
      panel.background = element_rect(fill = "aliceblue", color = NA),
      legend.position = "right",
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )
  
  print(gdp_plot)
  
  # 6. Tabella
  gdp_table <- gdp_data %>%
    filter(Country %in% euro_area_countries) %>%
    arrange(desc(GDP_Share_EA)) %>%
    dplyr::select(Country, GDP_Million_Euros = value, Contribution_Percentage = GDP_Share_EA)
  
  if (show_table) {
    return(kable(gdp_table, format = "latex", booktabs = TRUE, digits = 2,
                 caption = paste("GDP Contribution by Country in Euro Area -", quarter),
                 col.names = c("Country", "GDP (Million €)", "Contribution (%)")) %>%
             kable_styling(latex_options = c("hold_position", "striped", "scale_down")))
  } else {
    return(list(map = gdp_plot, table = gdp_table))
  }
}


################################################################################
############################## DATSET PREPARATION ##############################
################################################################################

# ==============================================================================
# STANDARDIZATION AND DESTANDARDIZATION
# ==============================================================================

standardize_Y <- function(Y) {
  T <- dim(Y)[1]
  p1 <- dim(Y)[2]
  p2 <- dim(Y)[3]
  
  Y_std <- array(NA, dim = dim(Y))
  dimnames(Y_std) <- dimnames(Y) 
  mean_Y <- matrix(0, p1, p2)
  sd_Y <- matrix(1, p1, p2)
  
  for (i in 1:p1) {
    for (j in 1:p2) {
      y_ij <- Y[, i, j]
      mu <- mean(y_ij, na.rm = TRUE)
      sigma <- sd(y_ij, na.rm = TRUE)
      if (is.na(sigma) || sigma == 0) sigma <- 1
      Y_std[, i, j] <- (y_ij - mu) / sigma
      mean_Y[i, j] <- mu
      sd_Y[i, j] <- sigma
    }
  }
  
  return(list(Y_scaled = Y_std, mean = mean_Y, sd = sd_Y))
}

inverse_standardize_Y <- function(Y_scaled, mean_Y, sd_Y) {
  T <- dim(Y_scaled)[1]
  p1 <- dim(Y_scaled)[2]
  p2 <- dim(Y_scaled)[3]
  
  Y_original <- array(NA, dim = dim(Y_scaled))
  dimnames(Y_original) <- dimnames(Y_centered)
  
  for (t in 1:T) {
    Y_original[t,,] <- Y_scaled[t,,] * sd_Y + mean_Y
  }
  
  return(Y_original)
}

center_Y <- function(Y) {
  T <- dim(Y)[1]
  p1 <- dim(Y)[2]
  p2 <- dim(Y)[3]
  
  Y_centered <- array(NA, dim = dim(Y))
  dimnames(Y_centered) <- dimnames(Y) 
  mean_Y <- matrix(0, p1, p2)
  
  for (i in 1:p1) {
    for (j in 1:p2) {
      y_ij <- Y[, i, j]
      mu <- mean(y_ij, na.rm = TRUE)
      Y_centered[, i, j] <- y_ij - mu
      mean_Y[i, j] <- mu
    }
  }
  
  return(list(Y_centered = Y_centered, mean = mean_Y))
}

decenter_Y <- function(Y_centered, mean_Y) {
  T <- dim(Y_centered)[1]
  p1 <- dim(Y_centered)[2]
  p2 <- dim(Y_centered)[3]
  
  Y_original <- array(NA, dim = dim(Y_centered))
  dimnames(Y_original) <- dimnames(Y_centered)
  
  for (t in 1:T) {
    Y_original[t,,] <- Y_centered[t,,] + mean_Y
  }
  
  return(Y_original)
}


# ==============================================================================
# NOT AVAILABLE PERCENTAGE
# ==============================================================================

nan_percent_Y <- function(Y) {
  total_values <- length(Y)
  total_NAs <- sum(is.na(Y))
  percent <- total_NAs / total_values * 100
  return(percent)
}


# ==============================================================================
# GENERAL PREPARATION OF DATA 
# ==============================================================================
prepare_country_data <- function(country_code, P) {
  modelM <- P$modelM
  modelQ <- P$modelQ
  
  # ---- SET DYNAMIC PATHS ----
  DataFile <- paste0("Data/", country_code, "/Processed/Data_", country_code, ".xlsx")
  LegendFile <- paste0("Data/", country_code, "/Original/Legend_", country_code, ".xlsx")
  
  # ---- LOAD MONTHLY DATA ----
  monthly_data <- readxl::read_excel(DataFile, sheet = "MonthlyLong", col_names = TRUE)
  DataM <- monthly_data[, -1]
  DatesM <- as.Date(monthly_data[[1]])
  
  legend_m <- readxl::read_excel(LegendFile, sheet = "MDescriptionFull")
  model_col_m <- which(colnames(legend_m) == modelM)
  ListM <- which((legend_m[, model_col_m]) == 1)
  
  DataM <- DataM[, ListM]
  GroupM <- legend_m[ListM, 16]
  SeriesM <- legend_m[ListM, 1]
  TransfM <- floor(legend_m[ListM, 8])
  UnbM <- legend_m[ListM, 12:14]
  ClassM <- legend_m[ListM, 11]
  
  DataMTrf <- DataM
  annualize_idx <- which(TransfM[[1]] %in% c(2, 4))
  DataMTrf[, annualize_idx] <- DataMTrf[, annualize_idx] * 100
  
  CovidNaN <- function(data, dates, class, start_covid, end_covid) {
    real_idx <- which(grepl("r", tolower(class)))
    mask <- dates >= as.Date(as.yearmon(start_covid)) & dates <= as.Date(as.yearmon(end_covid))
    covid_matrix <- matrix(FALSE, nrow = nrow(data), ncol = ncol(data))
    covid_matrix[mask, real_idx] <- TRUE
    data[covid_matrix] <- NA
    list(data, covid_matrix)
  }
  
  resM <- CovidNaN(DataMTrf, DatesM, ClassM[[1]], P$covid_start, P$covid_end)
  DataMTrf <- resM[[1]]
  Covid_obsM <- resM[[2]]
  T_M <- nrow(DataMTrf)
  
  # ---- LOAD QUARTERLY DATA ----
  quarterly_data <- readxl::read_excel(DataFile, sheet = "QuarterlyLong", col_names = TRUE)
  DataQ <- quarterly_data[, -1]
  DatesQ <- as.Date(quarterly_data[[1]])
  
  legend_q <- readxl::read_excel(LegendFile, sheet = "QDescriptionFull")
  model_col_q <- which(colnames(legend_q) == modelQ)
  ListQ <- which((legend_q[, model_col_q]) == 1)
  
  DataQ <- DataQ[, ListQ]
  GroupQ <- legend_q[ListQ, 16]
  SeriesQ <- legend_q[ListQ, 1]
  TransfQ <- floor(legend_q[ListQ, 8])
  UnbQ <- legend_q[ListQ, 12:14]
  ClassQ <- legend_q[ListQ, 11]
  
  DataQTrf <- DataQ
  annualize_q <- which(TransfQ[[1]] %in% c(2, 4))
  DataQTrf[, annualize_q] <- DataQTrf[, annualize_q] * 100
  
  resQ <- CovidNaN(DataQTrf, DatesQ, ClassQ[[1]], P$covid_start, P$covid_end)
  DataQTrf <- resQ[[1]]
  Covid_obsQ <- resQ[[2]]
  
  DataQMTrf <- do.call(rbind, lapply(1:nrow(DataQTrf), function(i) {
    rbind(rep(NA, ncol(DataQTrf)),
          rep(NA, ncol(DataQTrf)),
          DataQTrf[i, ])
  }))
  T_Q <- nrow(DataQMTrf)
  
  T_final <- max(T_M, T_Q)
  if (T_M < T_final) {
    paddingM <- data.frame(matrix(NA, T_final - T_M, ncol(DataMTrf)))
    colnames(paddingM) <- colnames(DataMTrf)
    DataMTrf <- rbind(DataMTrf, paddingM)
  }
  if (T_Q < T_final) {
    paddingQ <- data.frame(matrix(NA, T_final - T_Q, ncol(DataQMTrf)))
    colnames(paddingQ) <- colnames(DataQMTrf)
    DataQMTrf <- rbind(DataQMTrf, paddingQ)
  }
  
  Data <- cbind(DataMTrf, DataQMTrf)
  Series <- c(SeriesM, SeriesQ)
  Group <- c(GroupM, GroupQ)
  UnbPatt <- rbind(UnbM, UnbQ)
  
  result <- list(
    Data = Data,
    DatesM = DatesM,
    DatesQ = DatesQ,
    Series = Series,
    Group = Group,
    UnbalancedPattern = UnbPatt,
    Covid_obsM = Covid_obsM,
    Covid_obsQ = Covid_obsQ,
    Name = paste0("Data_", country_code)
  )
  
  assign(result$Name, result$Data, envir = .GlobalEnv)
  result$quarterly_start_idx <- ncol(DataMTrf) + 1
  return(result)
}

# ==============================================================================
# REMOVE COUNTRY CODES FROM VARIABLES NAME
# ==============================================================================

# Funzione per rimuovere suffissi paese
remove_country_code <- function(names, ccodes) {
  pattern <- paste0("_(", paste(ccodes, collapse = "|"), ")$")
  gsub(pattern, "", names)
}

# ==============================================================================
# BUILD METRIX-VARIATE DATA
# ==============================================================================

tensor_data <- function(countries, P) {
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
               as.character(1:T_max),
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
