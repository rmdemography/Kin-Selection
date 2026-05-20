#------------------------------------------------------------------------------#
# Paper:  Kin selection and population dynamics
# Title:  Reconstruct qx from posteriors of the parameters
# Data:   WPP 2024
# Author: Rahul Mondal
# Date:   27/03/2026
#------------------------------------------------------------------------------#

# posterior draw
lambda_d <- fit$draws("lambda")
delta_d <- fit$draws("delta")
beta_d <- fit$draws("beta")

# median draw
lambda_med <- apply(as_draws_matrix(lambda_d), 2, median)
delta_med <- apply(as_draws_matrix(delta_d), 2, median)
beta_med <- apply(as_draws_matrix(beta_d), 2, median)

# reshape
A <- stan_data$A
L <- stan_data$L
Ks <- stan_data$K
lambda <- array(lambda_med, dim = c(A,L,Ks))
delta <- matrix(delta_med, nrow = T, ncol = L)
beta <- as.numeric(beta_med)


xid <- df1$xid
tid <- df1$tid
lid <- df1$lid
X <- stan_data$X

# reconstruct qx
compute_qx_fitted <- function(lambda, delta, beta, xid, tid, lid, B, X){
  N <- length(xid)
  mu_logit <- numeric(N)
  for(n in 1:N){
    mu_logit[n] <- sum(B[tid[n],]*lambda[xid[n], lid[n],]) + sum(X[n,]*beta) + delta[tid[n], lid[n]]
  }
  return(plogis(mu_logit))
}

df1$qx_med <- compute_qx_fitted(lambda = lambda, delta = delta, beta = beta,
                        xid = xid, tid = tid, lid = lid, B = B, X = X)

ggplot(df1)+
  geom_line(aes(x=year, y = qx_median), color="red", lty = "dashed")+
  geom_line(aes(x=year, y=qx_med), color = "blue")+
  facet_grid(age~country)

# add dimension names
countries <- unique(df1$country)
years <- unique(df1$year)
rownames(B) <- years
rownames(delta) <- years
colnames(delta) <- countries
dimnames(lambda) <- list(NULL, countries, NULL)

cov_names <- c("gm", "os", "oso", "m")
beta <- setNames(as.numeric(beta_med), cov_names)



