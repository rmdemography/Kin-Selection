#------------------------------------------------------------------------------#
# Paper:  Kin selection and population dynamics
# Title:  Stochastic sensitivity
# Data:   HMD & HFD
# Author: Rahul Mondal
# Date:   27/03/2026
#------------------------------------------------------------------------------#
pacman::p_load(future, furrr, purrr)
set.seed(1234)
source("03_reconstruct_qx.R")
source("00_function.R")
# source("02_fit_mortality.R")

plan(multisession, workers = 5)

# compute stochastic sensitivity
countries <- unique(df1$country)
beta_names <- names(beta)
fx_list <- matlist$fx
Sx_list <- matlist$Sx

results_par <- future_map(countries, function(country){
  
  lambda_c <- lambda[,country,]
  delta_c <- delta[,country]
  
  message(sprintf("[%s] Solving nonlinear equations for \n", country))
  sens_inputs <- kin_dependent_mat(
    Fmat = fx_list[[country]],
    Umat = Sx_list[[country]],
    lambda = lambda_c,
    delta = delta_c,
    beta = beta,
    B = B,
    niter = 1000
  )
  
  map(beta_names, function(beta_i){
    message(sprintf("country=%s beta=%s\n", country, beta_i))
    
    dAt.dthetai = perturb_A(
      lambda = lambda_c, delta = delta_c, beta = beta, B = B,
      perturbation = 0.01, which.perturb = beta_i, 
      Usim = sens_inputs$Uhat, Fsim = sens_inputs$Fsim,
      kins = sens_inputs$Khat, years = sens_inputs$years
    )
    sens <- stoch_sens(
      A = sens_inputs$A, 
      dAt.dthetai = dAt.dthetai, 
      tlimit = sens_inputs$niter
    )
    list(
      summary = data.frame(
        country   = country,
        beta      = beta_i,
        beta_val  = beta[beta_i],
        lambda_s  = sens_inputs$lambdas,
        beta_sens = sens$beta_sens
      ),
      sensmat = data.frame(
        country = country,
        beta    = beta_i,
        row     = rep(1:nrow(sens$sensmat), each=ncol(sens$sensmat)),
        col     = rep(1:ncol(sens$sensmat), times=nrow(sens$sensmat)),
        value   = as.vector(sens$sensmat)
      )
    )
  })
}, .options = furrr_options(seed=TRUE))

plan(sequential)

# flatten results
results_flat <- unlist(results_par, recursive=FALSE)
lambdasens   <- do.call(rbind, lapply(results_flat, `[[`, "summary"))
sensmat_df   <- do.call(rbind, lapply(results_flat, `[[`, "sensmat"))
# save results
save(lambdasens, sensmat_df, file = "result/lambdasens.RData")
