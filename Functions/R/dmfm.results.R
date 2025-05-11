# ==============================================================================
# FACTORS VISUALIZATION
# ==============================================================================

plot_factors_tensor <- function(Factors, title = "Dynamic Factors", free_y = TRUE) {
  T <- dim(Factors)[1]
  k1 <- dim(Factors)[2]
  k2 <- dim(Factors)[3]
  
  df <- data.frame()
  for (i in 1:k1) {
    for (j in 1:k2) {
      df_temp <- data.frame(
        time = 1:T,
        value = Factors[, i, j],
        row_factor = paste0("F_row_", i),
        col_factor = paste0("F_col_", j),
        component = paste0("F(", i, ",", j, ")")
      )
      df <- rbind(df, df_temp)
    }
  }
  
  p <- ggplot(df, aes(x = time, y = value)) +
    geom_line(color = "steelblue", linewidth = 0.8) +
    facet_wrap(~ component, scales = if (free_y) "free_y" else "fixed") +
    theme_minimal(base_size = 14) +
    labs(
      title = title,
      x = "Tempo",
      y = "Valore"
    )
  
  return(p)
}


# ==============================================================================
# NOWCAST
# ==============================================================================


rolling_nowcast_dmfm <- function(Y, W, k, t_vec, model_type = "dmfm", max_iter = 1000, eps = 1e-03, target_idx = NULL) {
  # Y: array [T, p1, p2]
  # W: array [T, p1, p2]
  # k: vector c(k1, k2)
  # t_vec: vector of vintages (e.g. 20:300)
  # model_type: "dmfm", "fm", "fv"
  # target_idx: optional list(row, col) to extract specific nowcast
  
  T <- dim(Y)[1]
  p1 <- dim(Y)[2]
  p2 <- dim(Y)[3]
  
  nowcast_list <- list()
  factors_list <- list()
  
  for (t in t_vec) {
    cat(sprintf("Rolling nowcast at vintage t = %d\n", t))
    
    Y_t <- Y[1:t,, , drop = FALSE]
    W_t <- W[1:t,, , drop = FALSE]
    
    # Initial values using imputation
    imp <- mfm.cl(Y_t, W_t, k)
    
    # Initial parameters for MFM and MAR
    inputs <- dmfm.na.2.sv(imp$Y_imputed, k, W, t = "dmfm")
    
    # EM estimation
    fit <- dmfm.na.em(inputs, Y_t, k, W_t, t = model_type, max.iter = max_iter, eps = eps)
    model <- fit$model
    
    # Nowcast: filtered factors at time t
    f_t <- model$ff[t,,]
    
    # Reconstruct Y_t_hat
    Y_hat <- model$R %*% f_t %*% t(model$C)
    
    # Store full matrix or only selected entry
    if (!is.null(target_idx)) {
      nowcast_value <- Y_hat[target_idx[[1]], target_idx[[2]]]
    } else {
      nowcast_value <- Y_hat
    }
    
    nowcast_list[[as.character(t)]] <- nowcast_value
    factors_list[[as.character(t)]] <- f_t
  }
  
  return(list(nowcasts = nowcast_list, factors = factors_list))
}

# ---- 
rolling_nowcast_dfm <- function(Y, W, k, t_vec, model_type = "dfm", max_iter = 1000, eps = 1e-1, target_idx = NULL) {
  # Y: array [T, 1, p2]
  # W: array [T, 1, p2]
  # k: vector c(1, k2)
  # t_vec: vector of vintages
  # model_type: "dfm"
  # target_idx: optional column index for nowcast extraction
  
  T <- dim(Y)[1]
  p2 <- dim(Y)[3]
  
  nowcast_list <- list()
  factors_list <- list()
  
  for (t in t_vec) {
    cat(sprintf("Rolling nowcast at vintage t = %d\n", t))
    
    Y_t <- Y[1:t,, , drop = FALSE]
    W_t <- W[1:t,, , drop = FALSE]
    
    # Initial values
    imp <- mfm.cl.vec(Y_t, W_t, k)
    inputs <- dmfm.na.2.sv.vec(imp$Y_imputed, k, W_t, t = model_type)
    
    # EM estimation
    fit <- dmfm.na.em.vec(inputs, Y_t, k, W_t, t = model_type, max.iter = max_iter, eps = eps)
    model <- fit$model
    
    # Nowcast: filtered factor f_t
    f_t <- model$ff[t,,]  # dimensione: [1 x k2]
    
    # Reconstruct Y_hat: Y_t = f_t %*% t(C)
    Y_hat <- f_t %*% t(model$C)  # dimensione: [1 x p2]
    
    # Store result
    if (!is.null(target_idx)) {
      nowcast_value <- Y_hat[1, target_idx]
    } else {
      nowcast_value <- Y_hat
    }
    
    nowcast_list[[as.character(t)]] <- nowcast_value
    factors_list[[as.character(t)]] <- f_t
  }
  
  return(list(nowcasts = nowcast_list, factors = factors_list))
}


# ==============================================================================
# PREDICTION
# ==============================================================================


# dmfm.predict----
dmfm.predict <- function(., Y, k, W, t){
  
  n <- dim(Y)[1]
  p <- dim(Y)[-1]
  
  O <- lapply(split(W,seq(nrow(W))), function(x){which(vec(x)==1)})
  if (t=="fm"){
    kout  <- kalman.na(list(t=.$B%x%.$A,
                            l=.$C%x%.$R,
                            q=.$Q%x%.$P,
                            h=.$K%x%.$H),
                       apply(Y,1,vec), prod(k),
                       list(f=rep(0, prod(k)),
                            P=diag(prod(k))),
                       O)
  } else if (t=="fv"){
    kout  <- kalman.na(list(t=.$BA,
                            l=.$CR,
                            q=.$QP,
                            h=.$KH),
                       apply(Y,1,vec), prod(k),
                       list(f=rep(0, prod(k)),
                            P=diag(prod(k))),
                       O)
  }else{
    kout  <- kalman.na(list(t=.$BA,
                            l=.$C%x%.$R,
                            q=.$QP,
                            h=.$K%x%.$H),
                       apply(Y,1,vec), prod(k),
                       list(f=rep(0, prod(k)),
                            P=diag(prod(k))),
                       O)
  }
  
  .$ff <- .$fs <- array(0, c(n, k[1], k[2]))
  .$Pf <- .$Ps <- .$Cs <- array(0, c(n, prod(k), prod(k)))
  .$Y  <- array(0, c(n, p[1], p[2]))
  for (i in 1:n){
    .$ff[i,,] <- matrix(kout$ff[,i],k[1],k[2])
    .$fs[i,,] <- matrix(kout$fs[,i],k[1],k[2])
    .$Pf[i,,] <- kout$Pf[,,i]
    .$Ps[i,,] <- kout$Ps[,,i]
    .$Cs[i,,] <- kout$Cs[,,i]
    if (t=="fv"){
      .$Y[i,,] <- matrix(.$CR%*%kout$ff[,i], p[1], p[2])
    }else{
      .$Y[i,,] <- .$R%*%.$ff[i,,]%*%t(.$C)
    }
  }
  
  return(.)
  
}