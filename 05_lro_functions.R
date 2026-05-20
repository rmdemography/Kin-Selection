#------------------------------------------------------------------------------#
# Paper:  Kin selection and population dynamics
# Title:  LRO functions
# Data:   HMD & HFD
# Author: Rahul Mondal
# Date:   02/04/2026
#------------------------------------------------------------------------------#
pacman::p_load(Matrix)

mat2list <- function(mat){
  veclist <- lapply(seq_len(ncol(mat)),function(j) mat[,j])
  names(veclist) <- colnames(mat)
  return(veclist)
}

Ui.fn <- function(veclist){
  Ui <- lapply(veclist, function(s){
    na <- length(s)
    U <- matrix(0, na, na)
    for(i in 1:(na-1)){
      U[i+1,i] <- s[i]
    }
    U
  })
  return(Ui)
}

# iid environment
Dblock.fn <- function(g,s){
  prob <- rep(1/g,g) # iid environment
  Dj <- Matrix(prob %o% rep(1,g), sparse = TRUE)
  D <- kronecker(Diagonal(s),Dj)
  return(D)
}

# deterministic but aperiodic environment
Dblock.fn <- function(g,s){
  Dj <- Matrix(0, nrow = g, ncol = g, sparse = TRUE)
  for(i in 1:(g-1)){
    Dj[i+1,i] <- 1
  }
  Dj[g,g] <- 1
  D <- kronecker(Diagonal(s),Dj)
  return(D)
}

vecperm.fn <- function(g,s){
  gs <- g*s
  from_id <- integer(gs)
  to_id <- integer(gs)
  k <- 1L
  for(i in 1:g)for(j in 1:s){
    from_id[k] <- (i-1L)*s+j
    to_id[k] <- (j-1L)*g+i
    k <- k+1L
  }
  K <- sparseMatrix(
    i = to_id,
    j = from_id,
    x = rep(1,gs),
    dims = c(gs,gs)
  )
  return(K)
}

Utilde.fn <- function(D, K, Ui){
  Kt <- t(K)
  Ublock <- bdiag(Ui)
  step1 <- Ublock %*% Kt
  step2 <- K %*% step1
  Utilde <- D %*% step2
  return(Utilde)
}

ftilde.fn <- function(fmat, K){
  fvec <- as.vector(fmat)
  ftilde <- as.vector(K%*%fvec)
  return(ftilde)
}

R1tilde.fn <- function(ftilde){
  gs <- length(ftilde)
  ftildeZ <- c(ftilde,0)
  ones <- rep(1, gs+1)
  R1tilde <- Matrix(outer(ones, ftildeZ), sparse = TRUE)
  return(R1tilde)
}

Ptilde.fn <- function(Utilde){
  gs <- dim(Utilde)[1]
  dvec <- 1-colSums(Utilde)
  Ptilde <- rbind(cbind(Utilde, Matrix(0, gs, 1)), Matrix(c(dvec,1),1, gs+1))
  return(Ptilde)
}

rho1.fn <- function(Ptilde, ftilde, Utilde){
  gs <- dim(Utilde)[1]
  ftildeZ <- c(ftilde,0)
  PR <- Ptilde %*% Diagonal(x=ftildeZ)
  v <- colSums(PR)[1:gs]
  N <- Diagonal(gs)-t(Utilde)
  rho1 <- solve(N,v)
  return(as.vector(rho1))
}

#------------------------------------------------------------------------------#
# Sensitivity
#------------------------------------------------------------------------------#
u1.fn <- function(params, kins, nt){
  q1 <- rep(0,nt)
  for(t in 1:nt){
    logit_q1 = params$intercept + params$deltas + params$kappas[t] + 
      params$betas["gm"]*kins$gm[t]+ params$betas["os"]*kins$os[t]+
      params$betas["oso"]*kins$oso[t] + params$betas["age_m"]*kins$age_m[t]
    q1[t] <- plogis(logit_q1)
  }
  u1 <- 1-q1
  return(u1)
}

# block-construction matrices
Li.fn <- function(i,g,s){
  Li <- Matrix(0,g*s,s,sparse = TRUE)
  Li[((i-1)*s+1):(i*s),] <- diag(s)
  return(Li)
}
Qi.fn <- function(i,g,s){
  Qi <- Matrix(0, s, g*s, sparse = TRUE)
  Qi[,((i-1)*s+1):(i*s)] <- diag(s)
  return(Qi)
}

dUtildedthetai.fn <- function(dUidthetai,Dblock,K,i,g,s){
  Li <- Li.fn(i,g,s)
  Qi <- Qi.fn(i,g,s)
  dUblock <- Li %*% dUidthetai[[i]]%*% Qi
  DK <- Dblock %*% K
  dUtilde.dthetai <- DK %*% dUblock %*% t(K)
  # dUtilde.dthetai <- as.vector(dUtilde.dthetai)
  return(dUtilde.dthetai)
}

dPtildedthetai.fn <- function(dUtildedthetai){
  gs <- dim(dUtildedthetai)[1]
  C1_term <- rbind(cbind(dUtildedthetai, matrix(0,gs,1)),
                         matrix(0,1,gs+1))
  cs_dUtdthetai <- matrix(colSums(dUtildedthetai), nrow = 1)
  C2_term <- rbind(matrix(0, gs, gs+1), cbind(-cs_dUtdthetai,0))
  dPtildedthetai <- C1_term+C2_term
  return(dPtildedthetai)
}

term1.fn <- function(dPtildedthetai, R1tilde){
  gs <- dim(dPtildedthetai)[1]-1
  Z <- cbind(Diagonal(gs), matrix(0, gs, 1))
  ones <- matrix(1, gs+1, 1)
  term1 <- Z%*%t(R1tilde*dPtildedthetai)%*%ones
  return(term1)
}

term2.fn <- function(dUtildedthetai, rho1){
  term2 <- t(dUtildedthetai)%*%rho1
  return(term2)
}

drho1dthetai.fn <- function(term1, term2, Utilde){
  gs <- length(term1)
  Ntilde <- solve(Diagonal(gs)-Utilde)
  drho1dthetai <- t(Ntilde)%*%(term1+term2)
  return(drho1dthetai)
}


dUidthetai.fn <- function(lambda, delta, beta, B, perturbation, which.perturb, Umat, Khat){
  nt <- dim(Umat)[2]
  years <- colnames(Umat)
  
  beta.P1 <- beta
  beta.P2 <- beta
  beta.P1[[which.perturb]] <- beta[[which.perturb]] + perturbation
  beta.P2[[which.perturb]] <- beta[[which.perturb]] - perturbation
  # P1
  U5.P1 <- compute_qx(lambda = lambda, delta = delta, beta = beta.P1, B = B,
                      X = Khat, years = years)
  Umat.P1 <- Umat
  Umat.P1[1:5,] <- U5.P1
  Uvec.P1 <- mat2list(Umat.P1)
  Uhati.p1 <- Ui.fn(Uvec.P1)
  # P2
  U5.P2 <- compute_qx(lambda = lambda, delta = delta, beta = beta.P2, B = B,
                      X = Khat, years = years)
  Umat.P2 <- Umat
  Umat.P2[1:5,] <- U5.P2
  Uvec.P2 <- mat2list(Umat.P2)
  Uhati.p2 <- Ui.fn(Uvec.P2)
  
  dUidthetai <- list()
  for(t in 1:nt){
    dUidthetai[[t]] <- (Uhati.p1[[t]]-Uhati.p2[[t]])/(2*perturbation)
    # dUidthetai[[t]] <- as.vector(dUidthetai[[t]])
  }
  return(dUidthetai)
}

kin_dependent_Ui <- function(Smat, Fmat, lambda, delta, beta, B){
  years <- colnames(Smat)
  
  # compute kinship 
  kin_out <- kin(p = Smat, f = Fmat, time_invariant = FALSE)
  kins <- extract_kin(kin_out)
  
  # solve the nonlinear feedback to get K-hat
  na <- dim(lambda)[1]
  X0 <- kins2vec(kins)
  sol <- nleqslv(
    x = X0,
    fn = function(Xvec) solveXkin(
      Umat = Smat,
      Fmat = Fmat,
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
  Khat <- vec2kins(sol$x, na = na, nt = length(years))
  
  # compute child survival as a function of kins
  U5 <- compute_qx(lambda = lambda, delta = delta, beta = beta, 
                   B = B, X = Khat, years = years)
  Uhat <- Smat
  Uhat[1:5,] <- U5
  
  list(Khat = Khat, Uhat = Uhat)
}

lro_sens <- function(Uhat, Fmat, lambda, delta, beta, B, Khat, pert_param){
  s <- dim(Uhat)[1] # ages
  g <- dim(Uhat)[2] # years/environments
  years <- colnames(Uhat)
  
  # compute mean LRO
  Svec   <- mat2list(Uhat)
  Ui     <- Ui.fn(Svec)
  Dblock <- Dblock.fn(g, s)
  K      <- vecperm.fn(g, s)
  Utilde <- Utilde.fn(Dblock, K, Ui)
  ftilde <- ftilde.fn(Fmat, K)
  R1tilde <- R1tilde.fn(ftilde)
  Ptilde  <- Ptilde.fn(Utilde)
  rho1    <- rho1.fn(Ptilde, ftilde, Utilde)
  
  # perturb Ui
  dUidthetai <- dUidthetai.fn(
    lambda = lambda,
    delta = delta,
    beta = beta,
    B = B,
    perturbation = 0.01,
    which.perturb = pert_param,
    Umat = Uhat,
    Khat = Khat
  )
  
  # compute sensitivity
  sensmat <- matrix(0, g * s, g)
  colnames(sensmat) <- years
  
  for (i in seq_len(g)) {
    dUtildedthetai <- dUtildedthetai.fn(dUidthetai, Dblock, K, i, g, s)
    dPtildedthetai <- dPtildedthetai.fn(dUtildedthetai)
    term1          <- term1.fn(dPtildedthetai, R1tilde)
    term2          <- term2.fn(dUtildedthetai, rho1)
    drho1dthetai   <- drho1dthetai.fn(term1, term2, Utilde)
    sensmat[, i]   <- as.vector(drho1dthetai)
  }
  
  return(sensmat)
}


sensall2df <- function(sensall, g, s, ages, years){
  countries <- names(sensall)
  betas <- names(sensall[[1]])
  gs <- g*s
  
  dflist <- list()
  for(ctry in countries)for(beta in betas){
    sensmat <- sensall[[ctry]][[beta]]
    k <- 1:gs
    start_age <- floor((k-1)/g)
    start_env <- ((k-1)%%g)+1
    dfwide <- as.data.frame(sensmat)
    colnames(dfwide) <- as.character(years)
    dflong <- dfwide %>% 
      mutate(
        start_age = start_age,
        start_year = years[start_env],
        country = ctry,
        beta = beta
      ) %>% 
      pivot_longer(
        cols = as.character(years),
        names_to = "pert_year",
        values_to = "sensitivity"
      ) %>% 
      mutate(pert_year = as.integer(pert_year))
    dflist[[paste(ctry, beta, sep = "_")]] <- dflong
  }
  dffull <- bind_rows(dflist)
  return(dffull)
}

