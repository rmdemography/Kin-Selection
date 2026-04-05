#------------------------------------------------------------------------------#
# Paper:  Kin selection and population dynamics
# Title:  Compute LRO sensitivity
# Data:   HMD & HFD
# Author: Rahul Mondal
# Date:   02/04/2026
#------------------------------------------------------------------------------#
source("00_function.R")
source("03_reconstruct_q0.R")
source("05_lro_inputs.R")
# run full loop
g <- 74
s <- 111
gs <- g*s
years <- 1950:2023
countries <- names(result)
betas <- names(params$betas)
nc <- length(countries)
np <- length(betas)
sensall <- vector("list", nc)
names(sensall) <- countries

for(ctry in countries){
  cat("Country:", ctry, "\n")
  
  #input
  df <- result[[ctry]]$df
  fx <- result[[ctry]]$fx
  kin <- result[[ctry]]$kin$kin_summary
  
  Ui <- construct_U(df)$Ut
  qvec <- construct_U(df)$qx
  K <- vec_permutation_mat(g,s)
  D <- construct_D(g,s)
  Utilde <- construct_Utilde(D, K, Ui)
  ftilde <- construct_ftilde(fx, K, g, s)
  Ptilde <- construct_Ptilde(Utilde, g, s)
  rho1 <- compute_rho1(Ptilde, ftilde, Utilde, g, s)
  
  sensall[[ctry]] <- vector("list", np)
  names(sensall[[ctry]]) <- betas
  for(beta in betas){
    cat("\n Beta:", beta, "\n")
    params_c <- params
    params_c$deltas <- params_c$deltas[ctry]
    qpert <- qpert_fn(kin = kin, params = params_c, which.perturb = beta,
                      qx = qvec, eps = 0.01)
    sensmat <- matrix(0, gs, g)
    colnames(sensmat) <- years
    for(ei in 1:g){
      cat("Environment", ei, "of", g, "(year", years[ei], ")\r")
      sensbetai <- sensitivity_beta(Utilde = Utilde, ftilde = ftilde, g = g, s = s, 
                                    rho1 = rho1, ei = ei, qvec = qvec[[ei]], 
                                    qpert = qpert[[ei]], eps = 0.01)
      sensmat[,ei] <- as.vector(sensbetai[,1])
    }
    sensall[[ctry]][[beta]] <- sensmat
  }
}
ages <- 0:110
sens_beta <- sensall2df(sensall, g, s, ages, years)
sens_diagonal <- sens_beta %>%
  filter(
    start_age  == 0,
    start_year == pert_year   # diagonal condition
  ) %>%
  select(country, beta, year = pert_year, sensitivity)

sens_mean <- sens_diagonal %>% 
  group_by(country, beta) %>% 
  summarise(sensitivity = mean(sensitivity), .groups = "drop")
save(sens_mean, file = "sens_lro.RData")
