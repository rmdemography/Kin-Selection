#------------------------------------------------------------------------------#
# Paper:  Kin selection and population dynamics
# Title:  Stochastic sensitivity
# Data:   HMD & HFD
# Author: Rahul Mondal
# Date:   27/03/2026
#------------------------------------------------------------------------------#
pacman::p_load(future, furrr, purrr, progressr)

set.seed(1234)
source("00_function.R")
source("03_reconstruct_q0.R")

# plan(multisession, workers = availableCores() - 1)
plan(multisession, workers = 50)
handlers(global = TRUE)
handlers("progress")

# compute stochastic sensitivity
countries <- c("SWE", "JPN", "FRATNP")
beta_names <- names(params[["betas"]])
fx_list <- lapply(result, function(x) x$fx)
Sx_list <- lapply(result, function(x) x$Sx)

with_progress({
  p <- progressor(steps = length(countries) * length(beta_names))
  
  results_par <- future_map(
    countries,
    function(country){
      
      params_c <- params
      params_c$deltas <- params_c$deltas[country]
      
      map(beta_names, function(beta_i){
        p(sprintf("country=%s beta=%s", country, beta_i))
        
        res <- kin_dependent_mat(
          Fmat         = fx_list[[country]],
          Umat         = Sx_list[[country]],
          niter        = 1000,
          params       = params_c,
          perturbation = 0.01,
          pert_param   = beta_i
        )
        
        list(
          summary = data.frame(
            country   = country,
            beta      = beta_i,
            beta_val  = params[["betas"]][beta_i],
            lambda_s  = res$lambdas,
            beta_sens = res$sens$beta_sens
          ),
          sensmat = data.frame(
            country = country,
            beta    = beta_i,
            row     = rep(1:nrow(res$sens$sensmat), each = ncol(res$sens$sensmat)),
            col     = rep(1:ncol(res$sens$sensmat), times = nrow(res$sens$sensmat)),
            value   = as.vector(res$sens$sensmat)
          )
        )
      })
    },
    .options = furrr_options(seed = TRUE)
  )
})

plan(sequential)

results_flat <- unlist(results_par, recursive = FALSE)
summary_df   <- do.call(rbind, lapply(results_flat, `[[`, "summary"))
sensmat_df   <- do.call(rbind, lapply(results_flat, `[[`, "sensmat"))

