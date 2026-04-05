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

extract_kin <- function(kin){
  kin_df <- kin %>% 
    filter(kin %in% c("gm","m","os") & age_focal == 0) %>% 
    dplyr::select(-c(age_focal, cohort, count_cum_dead, mean_age_lost)) %>% 
    pivot_wider(
      id_cols = year,
      names_from = kin,
      values_from = c(count_living, mean_age, sd_age, count_dead)
    )
  gm <- kin_df$count_living_gm
  os <- kin_df$count_living_os
  age_m <- kin_df$mean_age_m
  list(gm = gm, os = os, age_m = age_m)
}

# infant survival as function of kins
U1_fn <- function(params, kins, years){
  nt <- length(years)
  qx = rep(0, nt)
  for(t in 1:nt){
    logit_qx = params$intercept + params$deltas + params$kappas[years[t]] + 
      params$betas["gm"]*kins$gm[t]+ params$betas["os"]*kins$os[t]+
      params$betas["age_m"]*kins$age_m[t]
    qx[t] <- plogis(logit_qx)
  }
  U1 <- 1-qx
  return(U1)
}

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

perturb_A <- function(params, perturbation, which.perturb, Usim, Fsim, kins, years){
  parameters.P1 <- params
  parameters.P2 <- params
  parameters.P1[["betas"]][[which.perturb]]<-params[["betas"]][[which.perturb]]+perturbation
  parameters.P2[["betas"]][[which.perturb]]<-params[["betas"]][[which.perturb]]-perturbation
  # P1
  U1.P1 <- U1_fn(params = parameters.P1, kins = kins, years = years)
  Uhat.P1 <- Usim
  Uhat.P1[1,] <- U1.P1
  A.P1 <- leslie_matrices(Umat = Uhat.P1, Fmat = Fsim)
  # P2
  U1.P2 <- U1_fn(params = parameters.P2, kins = kins, years = years)
  Uhat.P2 <- Usim
  Uhat.P2[1,] <- U1.P2
  A.P2 <- leslie_matrices(Umat = Uhat.P2, Fmat = Fsim)
  
  dAt.dthetai <- list()
  for(t in 1:length(A.P1)){
    dAt.dthetai[[t]] <- (A.P1[[t]]-A.P2[[t]])/(2*perturbation)
  }
  return(dAt.dthetai)
}

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

kin_dependent_mat <- function(Fmat, Umat, niter, params, perturbation, pert_param){
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
  kins <- extract_kin(kin_out$kin_summary)
  # compute infant survival as a function of kin
  U1 <- U1_fn(params = params, kins = kins, years = years)
  Uhat <- Usim
  Uhat[1,] <- U1
  # construct Leslie matrices 
  A <- leslie_matrices(Umat = Uhat, Fmat = Fsim)
  # compute stochastic growth rate
  lambdas <- lambda_s(matrices = A, niter = niter)
  # compute sensitivity
  dA.dthetai <- perturb_A(params = params, perturbation = perturbation, 
                          which.perturb = pert_param, Usim = Usim, Fsim = Fsim,
                          kins = kins, years = years)
  sens <- stoch_sens(A = A, dAt.dthetai = dA.dthetai, tlimit = niter)
  # output
  list(kin = kins, lambdas = lambdas, sens = sens)
}


