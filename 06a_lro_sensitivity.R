#------------------------------------------------------------------------------#
# Paper:  Kin selection and population dynamics
# Title:  Compute LRO sensitivity — parallelised with future
# Data:   HMD & HFD
# Author: Rahul Mondal
# Date:   02/04/2026
#------------------------------------------------------------------------------#

pacman::p_load(future, future.apply, furrr, purrr)

# source("02_fit_mortality.R")
source("03_reconstruct_qx.R")
source("00_function.R")
source("05_lro_functions.R")

plan(multisession, workers = 5)

countries <- unique(df1$country)
beta_names <- names(beta)
fx_list <- matlist$fx
Sx_list <- matlist$Sx

results_flat <- future_lapply(countries, function(country){
  Smat <- Sx_list[[country]]
  Fmat <- fx_list[[country]]
  
  lambda_c <- lambda[,country,]
  delta_c <- delta[,country]
  
  # compute Khat and Uhat solving nonlienar equation
  message(sprintf("[%s] Solving nonlinear equations for \n", country))
  res <- kin_dependent_Ui(
    Smat = Smat, 
    Fmat = Fmat, 
    lambda = lambda_c,
    delta = delta_c,
    beta = beta,
    B = B
  )
  Khat <- res$Khat
  Uhat <- res$Uhat
  
  # compute sensitivty for each beta
  sens_country <- vector("list", length(beta))
  names(sens_country) <- beta_names
  for(beta_i in beta_names){
    sens_country[[beta_i]] <- lro_sens(
      Uhat = Uhat,
      Fmat = Fmat,
      lambda = lambda_c,
      delta = delta_c,
      beta = beta,
      B = B,
      Khat = Khat,
      pert_param = beta_i
    )
  }
  list(country = country, sens = sens_country)
  
}, future.seed = TRUE)

# Reconstruct the nested list
lrosens <- vector("list", length(countries))
names(lrosens) <- countries
for(res in results_flat){
  lrosens[[res$country]] <- res$sens
}

# format the result
dim_Umat <- dim(matlist$Sx[[1]])
ages <- 0:(dim_Umat[1]-1)
s <- dim_Umat[1]
g <- dim_Umat[2]
years <- colnames(matlist$Sx[[1]])

sens_beta <- sensall2df(lrosens,g,s,ages,years)

sens_diagonal <- sens_beta %>%
  filter(start_age == 0, start_year == pert_year) %>%
  select(country, beta, year = pert_year, sensitivity)
save(sens_diagonal, file = "result/sens_lro.RData")

sens_mean <- sens_diagonal %>%
  group_by(country, beta) %>%
  summarise(sensitivity = mean(sensitivity), .groups = "drop")
save(sens_mean, file = "result/lrosens_mean.RData")

plan(sequential)

