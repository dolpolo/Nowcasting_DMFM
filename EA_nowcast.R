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
source("C:/Users/david/Desktop/University/Master's Thesis/Nowcasting_DMFM/Functions/R/dmfm.forecast.R")

# ==============================================================================
# SET PARAMETERS
# ==============================================================================

# No data for NL and BE in 2024 and afterwards
countries <- c("DE", "FR", "IT", "ES")

# Parameters for the data preparation
P <- list()

P$modelM <- "Large"
P$modelQ <- "Small"
P$covid_start <- as.Date("2020-03-01")
P$covid_end <- as.Date("2020-12-01")
P$EV_dfm <- as.Date("2017-01-01")
P$Tmax <- 300 

# Initial guess for the number of factors
kmax <- c(4, 15)


# ==============================================================================
# STEP 1: PREPARE COUNTRY-SPECIFIC 2D DATASETS
# ==============================================================================

country_list <- list()

for (cc in countries) {
  res <- prepare_country_data(cc, P)
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

# Demean and Standardize
std <- standardize_Y(Y)
Y_std <- std$Y_scaled

# Verify the percentage of NaN
nan_percent_Y(Y_std)

# ==============================================================================
# STEP 3: ANALYSIS ON FACTORS
# ==============================================================================

# Number of factors: ER by Yu Yu et al. (JoE, 2021)
k_hat <- mfm.f(Y_std, W, kmax)

# Number of lags in MAR equation
f.lag_mar <- factorize_bic_matrix(Y_std, W, k_hat)
f.lag_mar <- mar_model_selection(Y_std, W, k_hat)

# ==============================================================================
# STEP 3: LIST OF INPUTS --> . 
# ==============================================================================

# Just to try
# k_hat <- c(3,3)

# Can and Lam imputation
imp <- mfm.cl(Y_std, W, k_hat)

# Trapin Method (Barigozzi's variant) for initialization
inputs <- dmfm.na.2.sv(imp$Y_imputed, k_hat, W, t = "dmfm")


# ==============================================================================
# STEP 4: INTERPRETATION OF FACTORS 
# ==============================================================================

# Analyse the loading structure and plot the factor 
l.R_str <- inputs$R
l.C_str <- inputs$C

plot_factors_tensor(inputs$fs)



# ==============================================================================
# STEP 4: ESTIMATION
# ==============================================================================

# dmfm.em.llk----
# Prediction error gaussian log-likelihood used to control the convergence
# of the EM algorithm. Coincides with Eq.(17) in Barigozzi and Luciani (2022)
dmfm.em.llk <- function(., Y, t){
  
  n <- dim(Y)[1]
  p <- dim(Y)[-1]
  
  if (t=="fv"){
    V   <- array(0, c(n, prod(p), prod(p)))
    m <- matrix(0, n, prod(p))
    llk <- rep(0, n)
    for (i in 1:n){
      V[i,,] <- .$CR%*%.$Pf[i,,]%*%t(.$CR) + .$KH
      m[i,]  <- .$CR%*%vec(.$ff[i,,])
      llk[i] <- log(det(V[i,,])) + t(vec(Y[i,,])-m[i,])%*%solve(V[i,,])%*%(vec(Y[i,,])-m[i,])
    }
  } else {
    V   <- array(0, c(n, prod(p), prod(p)))
    m   <- array(0, c(n, p))
    d   <- rep(0,n)
    llk <- rep(0, n)
    for (i in 1:n){
      V[i,,] <- (.$C%x%.$R)%*%.$Pf[i,,]%*%t(.$C%x%.$R) + .$K%x%.$H
      m[i,,] <- .$R%*%.$ff[i,,]%*%t(.$C)
      d[i]   <- log(det(solve(.$Pf[i,,]) + (t(.$C)%*%solve(.$K)%*%.$C)%x%(t(.$R)%*%solve(.$H)%*%.$R))) + log(det(as.matrix(.$Pf[i,,]))) + p[1]*log(det(.$K)) + p[2]*log(det(.$H))
      llk[i] <- d[i] + t(vec(Y[i,,]-m[i,,]))%*%solve(V[i,,])%*%(vec(Y[i,,]-m[i,,]))
    }
  }
  
  out <- -0.5*sum(llk) 
  
  return(out)
  
}

# dmfm.na.em----
dmfm.na.em <- function(., Y, k, W, t, max.iter = 1000, eps = 1e-04){
  
  # Dimensions
  n <- dim(Y)[1]
  p <- dim(Y)[-1]
  
  # Required matrices: Apply a kronaker on Factor structure
  M <- matrix(0, k[1]*k[2], k[1]*k[2])
  for (i in 1:k[1]){
    for (j in 1:k[2]){
      i.k1 <- matrix(0, k[1], 1)
      i.k2 <- matrix(0, k[2], 1)
      i.k1[i] <- 1
      i.k2[j] <- 1
      M <- M + (i.k1%*%t(i.k2))%x%t(i.k1%*%t(i.k2))
    }
  }
  
  I1 <- diag(k[1])
  I2 <- diag(k[2])
  
  # Selection matrix vectorized 
  D.w  <- array(0, c(n, p[1]*p[2], p[1]*p[2]))
  D.wt <- array(0, c(n, p[1]*p[2], p[1]*p[2]))
  for (i in 1:n){
    i <- 1
    D.w[i,,] <- diag(as.vector(vec(W[i,,])))
    D.wt[i,,] <- diag(as.vector(vec(t(W[i,,]))))
  }
  
  criterion <- T
  iter <- 0
  while(criterion && iter < max.iter){
    cat(sprintf(">>> Iterazione EM: %d\n", iter))
    
    # Save parameters previous iteration
    .. <- .
    
    ## MAXIMIZATION STEP ##
    
    ## Loadings (R - C) ##
    R.1 <- matrix(0, p[1]*k[1], p[1]*k[1])
    R.2 <- matrix(0, p[1]*k[1], 1)
    C.1 <- matrix(0, p[2]*k[2], p[2]*k[2])
    C.2 <- matrix(0, p[2]*k[2], 1)
    # for (i in 1:n){
    #   Y[i,,][is.na(Y[i,,])] <- 0
    #   for (r in 1:p[1]){
    #     for (q in 1:p[1]){
    #       R.1 <- R.1 + star(t(..$C)%*%D.w[i,(0:(p[2]-1))*p[1]+r,(0:(p[2]-1))*p[1]+q]%*%solve(..$K)%*%..$C, vec(.$fs[i,,])%*%t(vec(.$fs[i,,])) + .$Ps[i,,])%x%(E.basis(r,q,p[1],p[1])%*%solve(..$H))
    #     }
    #   }
    #   R.2 <- R.2 + vec((W[i,,]*(solve(..$H)%*%Y[i,,]%*%solve(..$K)))%*%..$C%*%t(matrix(.$fs[i,,],k[1],k[2])))
    #   
    #   for (r in 1:p[2]){
    #     for (q in 1:p[2]){
    #       C.1 <- C.1 + star(t(..$R)%*%D.wt[i,(0:(p[1]-1))*p[2]+r,(0:(p[1]-1))*p[2]+q]%*%solve(..$H)%*%..$R, M%*%(vec(.$fs[i,,])%*%t(vec(.$fs[i,,])) + .$Ps[i,,])%*%t(M))%x%(E.basis(r,q,p[2],p[2])%*%solve(..$K))
    #     }
    #   }
    #   C.2 <- C.2 + vec(t(W[i,,]*(solve(..$H)%*%Y[i,,]%*%solve(..$K)))%*%..$R%*%matrix(.$fs[i,,],k[1],k[2]))
    # }
    # .$R <- matrix(solve(R.1)%*%R.2, p[1], k[1])
    # .$C <- matrix(solve(C.1)%*%C.2, p[2], k[2])
    Y[is.na(Y)] <- 0
    for (i in 1:n){
      for (r in 1:p[1]){
        for (q in 1:p[1]){
          R.1 <- R.1 + star(t(..$C)%*%D.w[i,(0:(p[2]-1))*p[1]+r,(0:(p[2]-1))*p[1]+q]%*%solve(..$K)%*%..$C, vec(.$fs[i,,])%*%t(vec(.$fs[i,,])) + .$Ps[i,,])%x%(E.basis(r,q,p[1],p[1])%*%solve(..$H))
        }
      }
      R.2 <- R.2 + vec((W[i,,]*(solve(..$H)%*%Y[i,,]%*%solve(..$K)))%*%..$C%*%t(matrix(.$fs[i,,],k[1],k[2])))
    }
    .$R <- matrix(solve(R.1)%*%R.2, p[1], k[1])
    for (i in 1:n){
      for (r in 1:p[2]){
        for (q in 1:p[2]){
          C.1 <- C.1 + star(t(.$R)%*%D.wt[i,(0:(p[1]-1))*p[2]+r,(0:(p[1]-1))*p[2]+q]%*%solve(..$H)%*%.$R, M%*%(vec(.$fs[i,,])%*%t(vec(.$fs[i,,])) + .$Ps[i,,])%*%t(M))%x%(E.basis(r,q,p[2],p[2])%*%solve(..$K))
        }
      }
      C.2 <- C.2 + vec(t(W[i,,]*(solve(..$H)%*%Y[i,,]%*%solve(..$K)))%*%.$R%*%matrix(.$fs[i,,],k[1],k[2]))
    }
    .$C <- matrix(solve(C.1)%*%C.2, p[2], k[2])
    
    ## Measurement variance (H - K) ##
    H. <- matrix(0, p[1], p[1])
    K. <- matrix(0, p[2], p[2])
    # for (i in 1:n){
    #   H. <- H. + t((W[i,,]*Y[i,,])%*%solve(..$K)%*%t(W[i,,]*Y[i,,]))
    #   H. <- H. - t((W[i,,]*Y[i,,])%*%solve(..$K)%*%t(W[i,,]*(.$R%*%.$fs[i,,]%*%t(.$C))))
    #   H. <- H. - t((W[i,,]*(.$R%*%.$fs[i,,]%*%t(.$C)))%*%solve(..$K)%*%t(W[i,,]*Y[i,,]))
    #   H. <- H. + t(star(solve(..$K), D.w[i,,]%*%(.$C%x%.$R)%*%(vec(.$fs[i,,])%*%t(vec(.$fs[i,,]))+.$Ps[i,,])%*%t(.$C%x%.$R)%*%t(D.w[i,,])))
    #   H. <- H. + (..$H%*%matrix(1,p[1],p[2])%*%..$K)%*%solve(..$K)%*%t(1-W[i,,])
    #   
    #   K. <- K. + t(t(W[i,,]*Y[i,,])%*%solve(..$H)%*%(W[i,,]*Y[i,,]))
    #   K. <- K. - t(t(W[i,,]*Y[i,,])%*%solve(..$H)%*%(W[i,,]*(.$R%*%.$fs[i,,]%*%t(.$C))))
    #   K. <- K. - t(t(W[i,,]*(.$R%*%.$fs[i,,]%*%t(.$C)))%*%solve(..$H)%*%(W[i,,]*Y[i,,]))
    #   K. <- K. + t(star(solve(..$H), D.wt[i,,]%*%(.$R%x%.$C)%*%(M%*%(vec(.$fs[i,,])%*%t(vec(.$fs[i,,]))+.$Ps[i,,])%*%t(M))%*%t(.$R%x%.$C)%*%t(D.wt[i,,])))
    #   K. <- K. + t(1-W[i,,])%*%solve(..$H)%*%(..$H%*%matrix(1,p[1],p[2])%*%..$K)
    # }
    # .$H <- diag(diag(H.)/(p[2]*n), p[1], p[1])
    # .$K <- diag(diag(K.)/(p[1]*n), p[2], p[2])
    for (i in 1:n){
      H. <- H. + t((W[i,,]*Y[i,,])%*%solve(..$K)%*%t(W[i,,]*Y[i,,]))
      H. <- H. - t((W[i,,]*Y[i,,])%*%solve(..$K)%*%t(W[i,,]*(.$R%*%.$fs[i,,]%*%t(.$C))))
      H. <- H. - t((W[i,,]*(.$R%*%.$fs[i,,]%*%t(.$C)))%*%solve(..$K)%*%t(W[i,,]*Y[i,,]))
      H. <- H. + t(star(solve(..$K), D.w[i,,]%*%(.$C%x%.$R)%*%(vec(.$fs[i,,])%*%t(vec(.$fs[i,,]))+.$Ps[i,,])%*%t(.$C%x%.$R)%*%t(D.w[i,,])))
      H. <- H. + (..$H%*%matrix(1,p[1],p[2])%*%..$K)%*%solve(..$K)%*%t(1-W[i,,])
    }
    .$H <- diag(diag(H.)/(p[2]*n), p[1], p[1])
    for (i in 1:n){
      K. <- K. + t(t(W[i,,]*Y[i,,])%*%solve(.$H)%*%(W[i,,]*Y[i,,]))
      K. <- K. - t(t(W[i,,]*Y[i,,])%*%solve(.$H)%*%(W[i,,]*(.$R%*%.$fs[i,,]%*%t(.$C))))
      K. <- K. - t(t(W[i,,]*(.$R%*%.$fs[i,,]%*%t(.$C)))%*%solve(.$H)%*%(W[i,,]*Y[i,,]))
      K. <- K. + t(star(solve(.$H), D.wt[i,,]%*%(.$R%x%.$C)%*%(M%*%(vec(.$fs[i,,])%*%t(vec(.$fs[i,,]))+.$Ps[i,,])%*%t(M))%*%t(.$R%x%.$C)%*%t(D.wt[i,,])))
      K. <- K. + t(1-W[i,,])%*%solve(.$H)%*%(.$H%*%matrix(1,p[1],p[2])%*%..$K)
    }
    .$K <- diag(diag(K.)/(p[1]*n), p[2], p[2])
    
    ## VAR parameters (A - B) ##
    if (t=="fm"){
      A.1 <- matrix(0, k[1], k[1])
      A.2 <- matrix(0, k[1], k[1])
      B.1 <- matrix(0, k[2], k[2])
      B.2 <- matrix(0, k[2], k[2])
      # for (i in 2:n){
      #   A.1 <- A.1 + star(solve(..$Q)%*%..$B, vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,])
      #   A.2 <- A.2 + star(t(..$B)%*%solve(..$Q)%*%..$B, vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,])) + .$Ps[i-1,,])
      #   B.1 <- B.1 + star(solve(..$P)%*%..$A, M%*%(vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,])%*%t(M))
      #   B.2 <- B.2 + star(t(..$A)%*%solve(..$P)%*%..$A, M%*%(vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,]))+.$Ps[i-1,,])%*%t(M))
      # }
      # .$A <- A.1%*%solve(A.2)
      # .$B <- B.1%*%solve(B.2)
      for (i in 2:n){
        A.1 <- A.1 + star(solve(..$Q)%*%..$B, vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,])
        A.2 <- A.2 + star(t(..$B)%*%solve(..$Q)%*%..$B, vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,])) + .$Ps[i-1,,])
      }
      .$A <- A.1%*%solve(A.2)
      for (i in 2:n){
        B.1 <- B.1 + star(solve(..$P)%*%.$A, M%*%(vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,])%*%t(M))
        B.2 <- B.2 + star(t(.$A)%*%solve(..$P)%*%.$A, M%*%(vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,]))+.$Ps[i-1,,])%*%t(M))
      }
      .$B <- B.1%*%solve(B.2)
    }else{
      BA.1 <- matrix(0, prod(k), prod(k))
      BA.2 <- matrix(0, prod(k), prod(k))
      for (i in 2:n){
        BA.1 <- BA.1 + vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,]
        BA.2 <- BA.2 + vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,])) + .$Ps[i-1,,]
      }
      .$BA <- BA.1%*%solve(BA.2)
    }
    
    ## State variance (P - Q) ##
    if (t=="fm"){
      P. <- matrix(0, k[1], k[1])
      Q. <- matrix(0, k[2], k[2])
      # for (i in 2:n){
      #   P.1 <- star(solve(..$Q), (vec(.$fs[i,,])%*%t(vec(.$fs[i,,])) + .$Ps[i,,]) )
      #   P.2 <- star(solve(..$Q)%*%.$B, (vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,]))%*%t(.$A)
      #   P.3 <- .$A%*%t(star(solve(..$Q)%*%.$B, (vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,])))
      #   P.4 <- .$A%*%star(t(.$B)%*%solve(..$Q)%*%.$B, (vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,])) + .$Ps[i-1,,]))%*%t(.$A)
      #   P.  <- P. + P.1 - P.2 - P.3 + P.4
      #   Q.1 <- star(solve(..$P), M%*%(vec(.$fs[i,,])%*%t(vec(.$fs[i,,]))+.$Ps[i,,])%*%t(M))
      #   Q.2 <- star(solve(..$P)%*%.$A, M%*%((vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,]))+.$Cs[i-1,,])%*%t(M)))%*%t(.$B)
      #   Q.3 <- .$B%*%t(star(solve(..$P)%*%.$A, M%*%((vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,]))+.$Cs[i-1,,])%*%t(M))))
      #   Q.4 <- .$B%*%star(t(.$A)%*%solve(..$P)%*%.$A, M%*%(vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,])) + .$Ps[i-1,,])%*%t(M))%*%t(.$B)
      #   Q.  <- Q. + Q.1 - Q.2 - Q.3 + Q.4
      # }
      # .$P <- P./((n-1)*k[2])
      # .$Q <- Q./((n-1)*k[1])
      for (i in 2:n){
        P.1 <- star(solve(..$Q), (vec(.$fs[i,,])%*%t(vec(.$fs[i,,])) + .$Ps[i,,]) )
        P.2 <- star(solve(..$Q)%*%.$B, (vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,]))%*%t(.$A)
        P.3 <- .$A%*%t(star(solve(..$Q)%*%.$B, (vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,])))
        P.4 <- .$A%*%star(t(.$B)%*%solve(..$Q)%*%.$B, (vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,])) + .$Ps[i-1,,]))%*%t(.$A)
        P.  <- P. + P.1 - P.2 - P.3 + P.4
      }
      .$P <- P./((n-1)*k[2])
      for (i in 2:n){
        Q.1 <- star(solve(.$P), M%*%(vec(.$fs[i,,])%*%t(vec(.$fs[i,,]))+.$Ps[i,,])%*%t(M))
        Q.2 <- star(solve(.$P)%*%.$A, M%*%((vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,]))+.$Cs[i-1,,])%*%t(M)))%*%t(.$B)
        Q.3 <- .$B%*%t(star(solve(.$P)%*%.$A, M%*%((vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,]))+.$Cs[i-1,,])%*%t(M))))
        Q.4 <- .$B%*%star(t(.$A)%*%solve(.$P)%*%.$A, M%*%(vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,])) + .$Ps[i-1,,])%*%t(M))%*%t(.$B)
        Q.  <- Q. + Q.1 - Q.2 - Q.3 + Q.4
      }
      .$Q <- Q./((n-1)*k[1])
    }else{
      QP <- matrix(0, prod(k), prod(k))
      for (i in 2:n){
        QP <- QP + vec(.$fs[i,,])%*%t(vec(.$fs[i,,])) + .$Ps[i,,] - (vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,])%*%t(.$BA)
      }
      .$QP <- QP/(n-1)
    }
    
    ## Expectation step ##
    O <- lapply(split(W,seq(nrow(W))), function(x){which(vec(x)==1)})
    if (t=="fm"){
      kout  <- kalman.na(list(t=.$B%x%.$A,
                              l=.$C%x%.$R,
                              q=.$Q%x%.$P,
                              h=.$K%x%.$H),
                         apply(Y,1,vec), prod(k),
                         list(f=vec(.$fs[1,,]),
                              P=diag(prod(k))),
                         O)
    }else{
      kout  <- kalman.na(list(t=.$BA,
                              l=.$C%x%.$R,
                              q=.$QP,
                              h=.$K%x%.$H),
                         apply(Y,1,vec), prod(k),
                         list(f=vec(.$fs[1,,]),
                              P=diag(prod(k))),
                         O)
    }
    .$ff <- .$fs <- array(0, c(n, k[1], k[2]))
    .$Pf <- .$Ps <- .$Cs <- array(0, c(n, prod(k), prod(k)))
    for (i in 1:n){
      .$ff[i,,] <- matrix(kout$ff[,i],k[1],k[2])
      .$fs[i,,] <- matrix(kout$fs[,i],k[1],k[2])
      .$Pf[i,,] <- kout$Pf[,,i]
      .$Ps[i,,] <- kout$Ps[,,i]
      .$Cs[i,,] <- kout$Cs[,,i]
    }
    
    ## Check convergence ##
    llk.old <- dmfm.em.llk(.., Y, t)
    llk.new <- dmfm.em.llk(., Y, t)
    delta <- 2*abs(llk.new-llk.old)/abs(llk.new+llk.old)
    
    cat(sprintf("Log-likelihood old: %.4f, new: %.4f, delta: %.6f\n", llk.old, llk.new, delta))
    
    print(c(llk.old,llk.new,delta))
    if (delta < eps) {
      criterion <- F
      .$Y <- array(0, c(n, p[1], p[2]))
      for (i in 1:n){
        .$Y[i,,] <- .$R%*%.$fs[i,,]%*%t(.$C)
      }
    }
    
    iter <- iter + 1
    
  }
  
  return(.)
  
}

apply(W[,,41], 2, sum)  # oppure:
sum(W[,,41])

# ==============================================================================
# STEP 4: ROLLING NOWCAST
# ==============================================================================

rolling_nowcast_dmfm <- function(Y, W, k, t_vec, model_type = "dmfm", max_iter = 1000, eps = 1e-1, target_idx = NULL) {
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
    inputs <- dmfm.na.2.sv(imp$Y_imputed, k_hat, W, t = "dmfm")
    
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

Y_nowcast <- rolling_nowcast_dmfm(Y_std, W, k = k_hat, t_vec = 290:300)

# Per ottenere i nowcast:
nowcasts <- Y_nowcast$nowcasts


save(res,tensor,inputs, std, Y_nowcast, k_hat, W, file = "dmfm_nowcastMS_11.RData")









gdp_col <- 41  # esempio: colonna in cui si trova il GDP
nowcast_gdp <- sapply(Y_nowcast$nowcasts, function(Yhat) Yhat[1, gdp_col])

plot(as.numeric(names(nowcast_gdp)), nowcast_gdp, type = "l", col = "blue",
     xlab = "Vintage", ylab = "Nowcast GDP", main = "Nowcast GDP nel tempo")
# Supponiamo che Y_true siano i dati veri per i mesi 290–300
Y_real <- Y_std[290:300,,]

# Estrai i nowcast della variabile gdp per i mesi corrispondenti
vintages <- as.numeric(names(Y_nowcast$nowcasts))
gdp_nowcast <- sapply(Y_nowcast$nowcasts, function(Yhat) Yhat[1, gdp_col])
gdp_true <- Y_real[, 1, gdp_col]

# Calcola RMSFE
rmsfe_gdp <- sqrt(mean((gdp_nowcast - gdp_true)^2, na.rm = TRUE))
cat("RMSFE GDP =", rmsfe_gdp, "\n")

# Imposta l’indice della variabile e del paese
row_idx <- 1       # es: primo paese
col_idx <- 41      # es: colonna in cui si trova GDP (modifica con la tua)

# Vettore dei vintages usati
vintages <- as.numeric(names(Y_nowcast$nowcasts))

# Estrai nowcast per ogni vintage
nowcast_series <- sapply(Y_nowcast$nowcasts, function(Yhat) Yhat[row_idx, col_idx])

# Estrai veri valori dallo stesso periodo
true_series <- Y_std[vintages, row_idx, col_idx]

# Supponiamo che la variabile di interesse sia alla riga 1, colonna 1
idx_row <- 1
idx_col <- 41

# Estrai le date dai nomi
vintages <- as.numeric(names(Y_nowcast$nowcasts))

# Estrai nowcast e veri valori
nowcast_series <- sapply(Y_nowcast$nowcasts, function(y) y[idx_row, idx_col])
true_series <- Y_std[vintages, idx_row, idx_col]

df <- data.frame(
  Date = vintages,  # oppure una sequenza di mesi: seq(as.Date("2000-01-01"), ..., by = "month")
  Nowcast = nowcast_series,
  RealGDP = true_series
)

# Plot: Nowcast vs Verità
plot(vintages, true_series, type = "l", col = "black", lwd = 2,
     ylab = "Valore (standardizzato)", xlab = "Vintage (t)", 
     main = "Nowcast vs. Valore Osservato", ylim = range(c(true_series, nowcast_series), na.rm = TRUE))
lines(vintages, nowcast_series, col = "blue", lwd = 2, lty = 2)
legend("topleft", legend = c("Vero", "Nowcast"), col = c("black", "blue"), lty = c(1, 2), lwd = 2)

library(ggplot2)

ggplot(df, aes(x = Date)) +
  geom_line(aes(y = Nowcast), color = "blue", linewidth = 1) +
  geom_point(aes(y = RealGDP), color = "red", size = 2) +
  labs(title = "Confronto tra Nowcast e PIL Reale",
       x = "Data",
       y = "Valore") +
  theme_minimal()

