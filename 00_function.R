#------------------------------------------------------------------------------#
# Paper:  Kin selection and population dynamics
# Title:  Functions for stochastic sensitivity
# Data:   HMD & HFD
# Author: Rahul Mondal
# Date:   27/03/2026
#------------------------------------------------------------------------------#

# make list of matrices
leslie_matrices <- function(Umat, Fmat, SRB = 1.05){
  na <- dim(Umat)[1]
  nt <- dim(Umat)[2]
  Ftilde <- matrix(0, na, nt)
  for(t in 1:nt) for(a in 1:(na-1)){
    Ftilde[a,t] <- Umat[1,t]*(1/(1+SRB))*0.5*(Fmat[a,t]+Fmat[a+1,t]*Umat[a+1,t])
  }
  A_list <- vector("list", nt)
  for(t in seq_len(nt)){
    A <- matrix(0, na, na)
    A[1,] <- Ftilde[,t]
    for(a in 1:(na-1)){
      A[a+1, a] <- Umat[a,t]
    }
    A_list[[t]] <- A
  }
  return(A_list)
}

# extract kin
extract_kin <- function(kin_out){
  kin_df <- kin_out$kin_full %>% 
    group_by(year, age_focal) %>% 
    summarise(
      gm = sum(living[kin == "gm"]),
      os = sum(living[age_kin <=10 & kin == "os"]),
      oso = sum(living[age_kin > 10 & kin == "os"]),
      m = sum(living[kin == "m"]),
      .groups = "drop"
    ) %>% 
    # considering a Poisson distribution for sisters
    mutate(os = 1-exp(-os),
           oso = 1-exp(-oso)) %>%
    filter(age_focal<5)
  to_matrix <- function(df, col){
    df %>% 
      select(age_focal, year, value = {{col}}) %>% 
      pivot_wider(names_from = year, values_from = value) %>% 
      arrange(age_focal) %>% 
      column_to_rownames("age_focal") %>% 
      as.matrix()
  }
  kin_list <- list(
    gm  = to_matrix(kin_df, gm),
    os  = to_matrix(kin_df, os),
    oso = to_matrix(kin_df, oso),
    m   = to_matrix(kin_df, m)
  )
  # kin_list <- lapply(kin_list, function(m) m - mean(m))
  
  kin_list
}

# infant survival as function of kins
compute_qx <- function(lambda, delta, beta, B, X, years){
  na <- dim(lambda)[1]
  nt <- length(years)
  qx <- matrix(0, nrow = na, ncol = nt)
  for(a in 1:na) for(t in 1:nt){
    mu_logit = sum(B[years[t],]*lambda[a,]) + 
      X$gm[a,t]*beta["gm"] + X$os[a,t]*beta["os"] +
      X$oso[a,t]*beta["oso"] + X$m[a,t]*beta["m"] +
      delta[years[t]]
    qx[a,t] <- plogis(mu_logit)
  }
  return(1-qx)
}

# stochastic growth rate
lambda_s <- function(matrices, niter){
  matrices <- matrix(unlist(matrices), ncol = length(matrices))
  a <- sqrt(dim(matrices)[1]) # number of ages
  n <- dim(matrices)[2] # number of matrices
  r <- numeric(niter)
  n0 <- rep(1, times = a) # inital age distirbution
  for(t in 1:niter){
    A <- matrix(matrices[,t], nrow = a)
    n0 <- A %*% n0
    N <- sum(n0)
    r[t] <- log(N)
    n0 <- n0/N
  }
  loglsim <- mean(r)
  dse <- 1.96 * sqrt(var(r) / niter)
  CI <- c(loglsim - dse, loglsim + dse)
  stoch <- exp(loglsim)
  
  return(stoch)
}

# compute d A/ d theta i
perturb_A <- function(lambda, delta, beta, B, perturbation, which.perturb, 
                      Usim, Fsim, kins, years){
  beta.P1 <- beta
  beta.P2 <- beta
  beta.P1[[which.perturb]] <- beta[[which.perturb]] + perturbation
  beta.P2[[which.perturb]] <- beta[[which.perturb]] - perturbation
  # P1
  U5.P1 <- compute_qx(lambda = lambda, delta = delta, beta = beta.P1, B = B,
                      X = kins, years = years)
  Uhat.P1 <- Usim
  Uhat.P1[1:5,] <- U5.P1
  A.P1 <- leslie_matrices(Umat = Uhat.P1, Fmat = Fsim)
  
  # P2
  U5.P2 <- compute_qx(lambda = lambda, delta = delta, beta = beta.P2, B = B,
                      X = kins, years = years)
  Uhat.P2 <- Usim
  Uhat.P2[1:5,] <- U5.P2
  A.P2 <- leslie_matrices(Umat = Uhat.P2, Fmat = Fsim)
  
  dAt.dthetai <- list()
  for(t in 1:length(A.P1)){
    dAt.dthetai[[t]] <- (A.P1[[t]]-A.P2[[t]])/(2*perturbation)
  }
  return(dAt.dthetai)
}

# stochastic sensitivity
stoch_sens <- function(A, dAt.dthetai, tlimit){
  k <- ncol(A[[1]])
  
  wvec <- rep(1/k, k)
  w <- cbind(wvec)
  r <- rep(0, tlimit)
  for(i in 1:tlimit){
    a <- A[[i]]
    wvec <- a %*% wvec
    r[i] <- sum(wvec)
    wvec <- wvec/r[i]
    w <- cbind(w, wvec)
  }
  
  vvec <- rep(1/k, k)
  v <- cbind(vvec)
  for (i in rev(1:tlimit)) {
    a <- A[[i]]
    vvec <- vvec %*% a
    v <- cbind(t(vvec), v)
  }
  
  sensmat <- matrix(0, nrow = k, ncol = k)
  beta_sens <- 0
  for(i in 1:tlimit){
    denom <- as.numeric(r[i] *t(v[, i + 1]) %*% w[, i + 1])
    sens_i <- (v[, i + 1] %*% t(w[, i]))/denom
    sensmat <- sensmat + sens_i
    beta_sens <- beta_sens + sum(sens_i*dAt.dthetai[[i]])
  }
  sensmat <- sensmat/tlimit
  beta_sens <- beta_sens/tlimit
  
  list(sensmat = sensmat, beta_sens = beta_sens)
}

# helper for reshaping inputs for nonlinear equation solver
kins2vec <- function(kins){
  c(as.vector(kins$gm), as.vector(kins$os), as.vector(kins$oso), as.vector(kins$m))
}

vec2kins <- function(vec,na,nt){
  n <- na*nt
  list(
    gm = matrix(vec[1:n], nrow = na, ncol = nt),
    os = matrix(vec[(n+1):(2*n)], nrow = na, ncol = nt),
    oso = matrix(vec[(2*n+1):(3*n)], nrow = na, ncol = nt),
    m = matrix(vec[(3*n+1):(4*n)], nrow = na, ncol = nt)
  )
}

# function to compute difference between guessed kin and the estimate
solveXkin <- function(Umat, Fmat, lambda, delta, beta, B, Xvec, years){
  nt <- length(years)
  na <- dim(lambda)[1]
  Xguess <- vec2kins(Xvec, na = na, nt = nt)
  # compute child survival as a function of kins
  U5 <- compute_qx(lambda = lambda, delta = delta, beta = beta, 
                   B = B, X = Xguess, years = years)
  Uhat <- Umat
  Uhat[1:5,] <- U5
  # construct kins
  kin_out <- kin(p = Uhat, f = Fmat, time_invariant = FALSE)
  Xvalue <- extract_kin(kin_out)
  # compute difference between estimated and guessed kins
  diff_gm <- Xguess$gm - Xvalue$gm
  diff_os <- Xguess$os - Xvalue$os
  diff_oso <- Xguess$oso - Xvalue$oso
  diff_m <- Xguess$m - Xvalue$m
  
  return(c(diff_gm, diff_os, diff_oso, diff_m))
}

# final function to run the model
kin_dependent_mat <- function(Fmat, Umat, lambda, delta, beta, B, niter){
  a <- dim(Fmat)[1]
  n <- dim(Fmat)[2]
  # simulate environments (IID)
  prob <- rep(1/n,n)
  col <- sample(1:n, niter, replace = TRUE, prob = prob)
  Fsim <- Fmat[,col]
  Usim <- Umat[,col]
  # construct kinship over simulated matrix sequences
  years <- colnames(Usim)
  colnames(Fsim) <- colnames(Usim) <- as.character(1:niter)
  kin_out <- kin(p = Usim, f = Fsim, time_invariant = FALSE)
  kins <- extract_kin(kin_out)
  # solve the nonlinear feedback to get K-hat
  na <- dim(lambda)[1]
  X0 <- kins2vec(kins)
  sol <- nleqslv(
    x = X0,
    fn = function(Xvec) solveXkin(
      Umat = Usim,
      Fmat = Fsim,
      lambda = lambda,
      delta = delta,
      beta = beta,
      B = B, 
      Xvec = Xvec,
      years = years
    ),
    method = "Newton",
    control = list(ftol = 1e-2, xtol = 1e-2, maxit = 200)
  )
  Khat <- vec2kins(sol$x, na = na, nt = niter)
  # compute child survival as a function of kins
  U5 <- compute_qx(lambda = lambda, delta = delta, beta = beta, 
                   B = B, X = Khat, years = years)
  Uhat <- Usim
  Uhat[1:5,] <- U5
  # construct Leslie matrices 
  A <- leslie_matrices(Umat = Uhat, Fmat = Fsim)
  # compute stochastic growth rate
  lambdas <- lambda_s(matrices = A, niter = niter)
  
  list(
    Fsim = Fsim, Uhat = Uhat, A = A, lambdas = lambdas, Khat = Khat, 
    years = years, niter = niter
  )
}




