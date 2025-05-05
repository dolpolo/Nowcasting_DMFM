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