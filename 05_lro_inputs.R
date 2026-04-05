#------------------------------------------------------------------------------#
# Paper:  Kin selection and population dynamics
# Title:  Construct LRO inputs
# Data:   HMD & HFD
# Author: Rahul Mondal
# Date:   02/04/2026
#------------------------------------------------------------------------------#
pacman::p_load(Matrix)
# survival matrix: U_t
construct_U <- function(df){
  nt <- length(unique(df$Year))
  na <- length(unique(df$Age))
  Uvec <- lapply(split(df, df$Year), function(d){
    1-d$qx[order(d$Age)]
  })
  qvec <- lapply(split(df, df$Year), function(d){
    d$qx[order(d$Age)]
  })
  Ut <- list()
  for(yr in names(Uvec)){
    s <- Uvec[[yr]]
    U <- matrix(0, na, na)
    for(i in 1:(na-1)){
      U[i+1,i] <- s[i]
    }
    Ut[[yr]] <- U
  }
  list(Ut = Ut, qx = qvec)
}

# vec permutation matrix K
vec_permutation_mat <- function(g,s){
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

# IID environmental transition matrix: D[g,g]
construct_D <- function(g,s){
  prob <- rep(1/g,g)
  gs <- g*s
  D_j <- Matrix(prob %o% rep(1,g), sparse = TRUE)
  D <- kronecker(Diagonal(s),D_j)
  return(D)
}

construct_Utilde <- function(D, K, U){
  Kt <- t(K)
  Ublock <- bdiag(U)
  step1 <- Ublock %*% Kt
  step2 <- K %*% step1
  Utilde <- D %*% step2
  return(Utilde)
}

construct_ftilde <- function(fx, K, g, s){
  fvec <- as.vector(fx)
  ftilde <- as.vector(K%*%fvec)
  return(ftilde)
}

construct_R1tilde <- function(ftilde, g,s){
  gs <- g*s
  ftilde_Z <- c(ftilde,0)
  ones <- rep(1, gs+1)
  R1tilde <- Matrix(outer(ones, ftilde_Z), sparse = TRUE)
  return(R1tilde)
}

construct_Ptilde <- function(Utilde, g, s){
  gs <- g*s
  dvec <- 1-colSums(Utilde)
  Ptilde <- rbind(cbind(Utilde, Matrix(0, gs, 1)), Matrix(c(dvec,1),1, gs+1))
  return(Ptilde)
}

# LRO mean
compute_rho1 <- function(Ptilde, ftilde, Utilde, g,s){
  gs <- g*s
  ftilde_full <- c(ftilde,0)
  PR <- Ptilde %*% Diagonal(x=ftilde_full)
  v <- colSums(PR)[1:gs]
  N <- Diagonal(gs)-t(Utilde)
  rho1 <- solve(N,v)
  return(as.vector(rho1))
}

# sensitivity function (van Daalen & Caswell 2020)
compute_sensitivity_DC2020 <- function(Utilde, ftilde, g, s, rho1, i){
  gs <- g*s
  alpha <- 1L
  K <- vec_permutation_mat(g,s)
  N <- Diagonal(gs)-Utilde
  N_factor <- Matrix::lu(t(N))
  Ptilde <- construct_Ptilde(Utilde, g, s)
  ftilde_Z <- c(ftilde,0)
  ones <- rep(1, gs+1)
  R1vec <- rep(ftilde_Z, each = gs+1)
  # R1vec <- as.vector(outer(ones, ftilde_Z))
  K1 <- vec_permutation_mat(gs+1, gs+1)
  K2 <- vec_permutation_mat(gs, gs)
  embed_alive <- rbind(Diagonal(gs), matrix(0, alpha, gs))
  embed_dead <- rbind(Matrix(0, gs, alpha), Diagonal(alpha))
  C1 <- kronecker(embed_alive, embed_alive)
  C2 <- kronecker(embed_alive, embed_dead)
  I1 <- kronecker(Diagonal(gs), Matrix(1,1,gs))
  dPdU <- C1-C2 %*% I1
  Zmat <- cbind(Diagonal(gs), Matrix(0,gs,1))
  Z1 <- kronecker(Matrix(1,1,gs+1), Zmat)
  DvecR1_dPdU <- .sparseDiagonal(x=R1vec)%*%dPdU
  inner1 <- Z1 %*% K1 %*% DvecR1_dPdU
  rho1I <- kronecker(Matrix(rho1, 1, gs), Diagonal(gs))
  inner2 <- rho1I %*% K2
  bracket <- inner1+inner2
  sens_vecU <- solve(t(N), bracket)
  dvecU_dqi <- dvecUtilde_dqi(K,g,s,i)
  drho1_dqi <- sens_vecU %*% dvecU_dqi
  return(drho1_dqi)
}

dvecUtilde_dqi <- function(K, g, s, i){
  gs <- g*s
  D <- construct_D(g,s)
  rid <- 2:s
  cid <- 1:(s-1)
  vec_id <- (cid-1)*s+rid
  dvecUi_dqi <- sparseMatrix(
    i = vec_id,
    j = 1:(s-1),
    x = rep(-1, s-1),
    dims = c(s^2,s)
  )
  Li <- sparseMatrix(
    i = ((i-1)*s+1):(i*s),
    j = 1:s,
    x = rep(1,s),
    dims = c(gs, s)
  )
  Qi <- t(Li)
  QL <- kronecker(t(Qi), Li)
  KDK <- kronecker(K, D %*% K)
  dvecUtilde_dqi <- KDK %*% QL %*% dvecUi_dqi
  return(dvecUtilde_dqi)
}

# reformulated sensitivity
compute_sensitivity <- function(Utilde, ftilde, g, s, rho1, ei){
  gs <- g*s
  alpha <- 1L
  pi <- rep(1/g, g)
  ftilde_z <- c(ftilde,0)
  Nt <- t(Diagonal(gs) - Utilde)
  Nfactor <- Matrix::lu(Nt)
  drho1_dqi <- Matrix(0, gs, s, sparse = TRUE)
  for(a in 0:(s-2)){
    colU <- a*g+ei
    rowU <- (a+1)*g+(1:g)
    # term 1
    t1 <- numeric(gs)
    t1[colU] <- -ftilde_z[colU]
    # term 2
    t2 <- numeric(gs)
    t2[colU] <- -sum(pi*rho1[rowU])
    rhs <- t1+t2
    drho1_dqi[,a+1] <- solve(Nfactor, rhs)
  }
  return(drho1_dqi)
}

# sensitivity to kin coefficients beta
sensitivity_beta <- function(Utilde, ftilde, g, s, rho1, ei, qvec, qpert, eps){
  gs <- g*s
  pi <- rep(1/g,g)
  drho1_dqi <- compute_sensitivity(Utilde, ftilde, g, s, rho1, ei)
  dq_dtheta <- numeric(s)
  for(a in 1:s){
    dq_dtheta[a] <- (qpert[a]-qvec[a])/eps
  }
  drho1_dthetai <- Matrix(0, gs, s, sparse = TRUE)
  for(a in 1:s){
    drho1_dthetai[,a] <- drho1_dqi[,a]*dq_dtheta[a]
  }
  return(drho1_dthetai)
}

qpert_fn <- function(kin, params, which.perturb, qx, eps){
  nt <- length(qx)
  kins <- extract_kin(kin)
  parameters <- params
  parameters[["betas"]][[which.perturb]]<-params[["betas"]][[which.perturb]]+eps
  q0 <- rep(0, nt)
  for(t in 1:nt){
    logit_q0 = parameters$intercept + parameters$deltas + 
      parameters$kappas[t] + parameters$betas["gm"]*kins$gm[t]+ 
      parameters$betas["os"]*kins$os[t]+ parameters$betas["age_m"]*kins$age_m[t]
    q0[t] <- plogis(logit_q0)
    qx[[t]][1] <- q0[t]
  }
  return(qx)
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
