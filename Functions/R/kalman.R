kalman.na <- function(par, y, k, start, O){

  # Size
  n <- ncol(y)
  d <- nrow(y)

  # Parameters
  t <- par$t # VAR parameters
  Z <- par$l # Factor loadings
  Q <- par$q # Q is the state variance
  H <- par$h # H is the measurement variance

  # Initialize
  f.f <- f.u <- f.s <- matrix(0, k, n)
  v <- matrix(0, d, n)
  P.f <- P.u <- P.s <- array(0, c(k, k, n))
  S.1 <- array(0, c(d, d, n)) # inverse of F in Durbin (2012)
  L <- array(0, c(k, k, n))
  qK  <- array(0, c(k, d, n)) # quasi Kalman gain
  C.f <- C.s <- array(0, c(k, k, n))

  # kalman filter
  f.f[,1]  <- f.u[,1]  <- start$f
  P.f[,,1] <- P.u[,,1] <- start$P
  v[O[[1]],1]  <- y[O[[1]],1] - Z[O[[1]],]%*%f.f[,1,drop=F]
  S.1[O[[1]],O[[1]],1] <- solve(Z[O[[1]],]%*%as.matrix(P.f[,,1])%*%t(Z[O[[1]],]) + H[O[[1]],O[[1]]])
  for (i in 2:n){
    f.f[,i]              <- t%*%f.u[,i-1, drop=F]
    P.f[,,i]             <- t%*%P.u[,,i-1]%*%t(t) + Q
    S.1[O[[i]],O[[i]],i] <- solve(Z[O[[i]],]%*%as.matrix(P.f[,,i])%*%t(Z[O[[i]],]) + H[O[[i]],O[[i]]])
    qK[,O[[i]],i]        <- as.matrix(P.f[,,i])%*%t(Z[O[[i]],])%*%S.1[O[[i]],O[[i]],i]
    v[O[[i]],i]          <- y[O[[i]],i] - Z[O[[i]],]%*%f.f[,i,drop=F]
    f.u[,i]              <- f.f[,i,drop=F] + qK[,O[[i]],i]%*%v[O[[i]],i]
    P.u[,,i]             <- as.matrix(P.f[,,i]) - qK[,O[[i]],i]%*%Z[O[[i]],]%*%as.matrix(P.f[,,i])
    L[,,i]               <- t - t%*%qK[,O[[i]],i]%*%Z[O[[i]],]
    C.f[,,i]             <- as.matrix(P.f[,,i])%*%t(L[,,i])
  }

  # Kalman smoother
  C.s[,,n] <- t(C.f[,,n])
  r <- matrix(0, k, 1)
  N <- matrix(0, k ,k)
  for (i in (n):1){
    if (i<n){C.s[,,i] <- t(C.f[,,i]%*%(diag(k) - N%*%P.f[,,i+1]))}
    r <- t(Z[O[[i]],])%*%S.1[O[[i]],O[[i]],i]%*%v[O[[i]],i] + t(L[,,i])%*%r
    N <-  t(Z[O[[i]],])%*%S.1[O[[i]],O[[i]],i]%*%Z[O[[i]],] + t(L[,,i])%*%N%*%L[,,i]
    f.s[,i]  <- f.f[,i,drop=F] + as.matrix(P.f[,,i])%*%r
    P.s[,,i] <- as.matrix(P.f[,,i]) - as.matrix(P.f[,,i])%*%N%*%as.matrix(P.f[,,i])
  }

  return(list(ff=f.f, Pf=P.f, fs=f.s, Ps=P.s, Cs=C.s))

}

kalman <- function(par, y, k, start){

  # Size
  n <- ncol(y)
  d <- nrow(y)

  # Parameters
  t <- par$t # VAR parameters
  Z <- par$l # Factor loadings
  Q <- par$q # Q is the state variance
  H <- par$h # H is the measurement variance

  # Initialize
  f.f <- f.u <- f.s <- matrix(0, k, n)
  v   <- matrix(0, d, n)
  P.f <- P.u <- P.s <- array(0, c(k, k, n))
  S.1 <- array(0, c(d, d, n)) # inverse of F in Durbin (2012)
  L   <- array(0, c(k, k, n))
  qK  <- array(0, c(k, d, n)) # quasi Kalman gain
  C.f <- C.s <- array(0, c(k, k, n))

  # kalman filter
  f.f[,1]  <- f.u[,1]  <- start$f
  P.f[,,1] <- P.u[,,1] <- start$P
  v[,1]    <- y[,1] - Z%*%f.f[,1]
  S.1[,,1] <- solve(Z%*%P.f[,,1]%*%t(Z) + H)
  for (i in 2:n){
    f.f[,i]  <- t%*%f.u[,i-1, drop=F]
    P.f[,,i] <- t%*%P.u[,,i-1]%*%t(t) + Q
    S.1[,,i] <- solve(Z%*%P.f[,,i]%*%t(Z) + H)
    qK[,,i]  <- P.f[,,i]%*%t(Z)%*%S.1[,,i]
    v[,i]    <- y[,i] - Z%*%f.f[,i]
    f.u[,i]  <- f.f[,i] + qK[,,i]%*%v[,i]
    P.u[,,i] <- P.f[,,i] - qK[,,i]%*%Z%*%P.f[,,i]
    L[,,i]   <- t - t%*%qK[,,i]%*%Z
    C.f[,,i] <- P.f[,,i]%*%t(L[,,i])
  }

  # Kalman smoother
  C.s[,,n] <- t(C.f[,,n])
  r <- matrix(0, k, 1)
  N <- matrix(0, k ,k)
  for (i in (n):1){
    if (i<n){C.s[,,i] <- t(C.f[,,i]%*%(diag(k) - N%*%P.f[,,i+1]))}
    r <- t(Z)%*%S.1[,,i]%*%v[,i] + t(L[,,i])%*%r
    N <-  t(Z)%*%S.1[,,i]%*%Z + t(L[,,i])%*%N%*%L[,,i]
    f.s[,i]  <- f.f[,i] + P.f[,,i]%*%r
    P.s[,,i] <- P.f[,,i] - P.f[,,i]%*%N%*%P.f[,,i]
  }

  return(list(ff=f.f, Pf=P.f, fs=f.s, Ps=P.s, Cs=C.s))

}
