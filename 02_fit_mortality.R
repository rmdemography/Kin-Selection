#------------------------------------------------------------------------------#
# Paper:  Kin selection and population dynamics
# Title:  Fit Bayesian model for child mortality
# Data:   WPP 2024
# Author: Rahul Mondal
# Date:   22/04/2026
#------------------------------------------------------------------------------#
pacman::p_load(splines, cmdstanr, posterior, bayesplot, loo, ggplot2, RColorBrewer)
source("01a_compute_kin.R")
df1 <- kin_df %>%
  rename(age = age_focal) %>%
  mutate(os = 1-exp(-os),
         oso = 1-exp(-oso)) %>%
  left_join(df, by = c("country", "year", "age")) %>%
  filter(age < 5)

# indices
A <- length(unique(df1$age))
T <- length(unique(df1$year))
L <- length(unique(df1$country))

df1 <- df1 %>% 
  mutate(
    xid = age + 1L,
    tid = as.integer(factor(year)),
    lid = as.integer(factor(country))
  )
# B-spline
years <- unique(df1$year)
knots <- seq(min(years), max(years), by = 2.5)
B <- bs(years, knots = knots[-c(1,length(knots))], 
        degree = 3,  intercept = TRUE)
K <- ncol(B)

# Second-order difference matrix D_K
make_diff_matrix <- function(K) {
  D <- matrix(0, nrow = K - 2, ncol = K)
  for (i in 1:(K - 2)) {
    D[i, i]     <-  1
    D[i, i + 1] <- -2
    D[i, i + 2] <-  1
  }
  return(D)
}

D  <- make_diff_matrix(K)
Q  <- K - 2
Pmat <- t(D) %*% solve(D%*%t(D)) # Projection matrix: D'(DD')^-1

# center covariates
cov_vars <- c("gm", "os", "oso", "m")
# df1 <- df1 %>%
#   mutate(across(all_of(cov_vars), ~. - mean(., na.rm = TRUE)))
# df1 %>% summarise(across(all_of(cov_vars), mean))

df1 <- df1 %>% arrange(xid, tid, lid)

# informative prior on beta
# Jamison et al. (2002)
m1 <- c(0.655, 0.377, 1.138)
gm1 <- c(0.884, 0.434, 1.802)
os <- c(1.084, 0.816, 1.440)
# Sear and Mace (2009)
m2 <- c(beta = -1.66, se = 0.61)
gm2 <- c(beta = -0.55, se = 0.27)
oso <- c(beta = -0.48, se = 0.24)

or2beta <- function(or, lci, uci){
  beta <- log(or)
  se <- (log(uci)-log(lci))/(2*1.96)
  return(c(beta = beta, se = se))
}

combine_beta <- function(...){
  estimates <- list(...)
  betas <- sapply(estimates, function(x) x["beta"])
  ses <- sapply(estimates, function(x) x["se"])
  weights <- 1/ses^2
  beta <- sum(weights*betas)/sum(weights)
  se <- sqrt(1/sum(weights))
  return(c(beta = beta, se = se))
}

# compute beta and SE from OR
m1 <- or2beta(0.655, 0.377, 1.138)
gm1 <- or2beta(0.884, 0.434, 1.802)
os <- or2beta(1.084, 0.816, 1.440)

# combine estimates across studies by inverse-variance weights
m <- combine_beta(m1, m2)
gm <- combine_beta(gm1, gm2)

# define the priors
mu_beta <- c(gm["beta"], os["beta"], oso["beta"], m["beta"])
sigma_beta <- c(gm["se"], os["se"], oso["se"], m["se"])

# horseshoe hyperparameter
p0 <- 2
n_delta <- T*L
scale_global <- (p0/(n_delta-p0))*(1/sqrt(A))

# stan data
stan_data <- list(
  L = L,
  A = A,
  T = T,
  K = K,
  Q = Q,
  N = nrow(df1),
  P = length(cov_vars),
  B = B,
  Pmat = Pmat,
  logit_qx = qlogis(df1$qx),
  X = as.matrix(df1[, cov_vars]),
  xid = df1$xid,
  tid = df1$tid,
  lid = df1$lid,
  mu_beta = mu_beta,
  sigma_beta = sigma_beta,
  scale_global = scale_global,
  nu_global = 2,
  nu_local = 1,
  slab_scale = 1,
  slab_df = 4
)


# Compile
mod <- cmdstan_model("model1.stan")

# Fit
fit <- mod$sample(
  data            = stan_data,
  seed            = 123,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 2000,
  iter_sampling   = 2000,
  adapt_delta     = 0.9,
  max_treedepth   = 12
)
fit$cmdstan_diagnose()

# Posterior Predictive Checks
# density
qx_rep <- fit$draws("qx_rep", format = "matrix")
ppc_dens_overlay(y = df1$qx, yrep = qx_rep)
ggsave("result/A1a_density.png", width = 10, height = 8, units = "cm", dpi = 300)

ppc_dens_overlay_grouped(y = df1$qx, yrep = qx_rep, group = factor(df1$age))
ggsave("result/A1b_density_age.png", width = 16, height = 14, units = "cm", dpi = 300)

ppc_dens_overlay_grouped(y = df1$qx, yrep = qx_rep, group = factor(df1$country))

# residuals
qx_rep_mean <- colMeans(qx_rep)
logit_res <- (df1$qx)-(qx_rep_mean)
df1 <- df1 %>% 
  mutate(qx_pred = qx_rep_mean,
         res = logit_res)
ggplot(df1, aes(x = (qx_pred), y = res)) +
  geom_point(alpha = 0.3)+
  geom_hline(yintercept = 0, colour = "red")+
  labs(x = "Fitted q", y = "q Residual")+
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 12, color = "black", face = "bold",
                                margin = margin(t = 0.5, unit = "cm")),
    axis.title.y = element_text(size = 12, color = "black", face = "bold",
                                margin = margin(r = 0.5, unit = "cm")),
    plot.margin = margin(1,1,1,1, "cm")
  )
ggsave("result/A2_residual.png", width = 14, height = 14, units = "cm", dpi = 300)

# posterior median and credible interval
qx_rep_median <- apply(qx_rep, 2, median)
qx_rep_lower <- apply(qx_rep, 2, quantile, 0.025)
qx_rep_upper <- apply(qx_rep, 2, quantile, 0.975)
df1 <- df1 %>% 
  mutate(
    qx_median = qx_rep_median,
    qx_lower = qx_rep_lower, 
    qx_upper = qx_rep_upper
  )

display.brewer.all(colorblindFriendly = TRUE)
col1 <- brewer.pal(9, "Blues")[8]
col2 <- brewer.pal(9, "Blues")[5]
col3 <- brewer.pal(9, "Greys")[9]
ggplot(df1, aes(x = year))+
  geom_ribbon(aes(ymin = qx_lower, ymax = qx_upper), 
              alpha = 0.9, fill = col2)+
  geom_line(aes(y = qx_median , linetype = "Fitted"), colour = col1, lwd = 0.5)+
  geom_line(aes(y = qx, linetype = "Observed"), colour = col3, lwd = 0.5)+
  facet_grid(country~age)+
  scale_x_continuous(breaks = c(1960, 1990, 2020), 
                     label = c("1960", "'90", "2020"))+
  scale_y_continuous(n.breaks = 3)+
  scale_linetype_manual(values = c("Observed" = "dashed", "Fitted" = "solid"))+
  labs(x = "Year", y = expression(q[x]), linetype = NULL)+
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 14, color = "black", face = "bold", 
                                margin = margin(t=0.5, unit = "cm")),
    axis.title.y = element_text(size = 14, color = "black", face = "bold", 
                                margin = margin(r = 0.5, unit = "cm")),
    axis.text = element_text(size = 11, color = "black"),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_blank(),
    legend.text = element_text(size = 11, color = "black"),
    legend.position = "bottom",
    plot.margin = margin(1,1,1,1,"cm")
  )
ggsave("result/A3_fit.png", width = 20, height = 18, units = "cm", dpi = 300)  

# posterior predictive test statistics
ppc_stat(y = df1$qx, yrep = qx_rep, stat = "mean")
ppc_stat(y = df1$qx, yrep = qx_rep, stat = "sd")

ppc_stat_grouped(y = df1$qx, yrep = qx_rep, group = factor(df1$age), stat = "mean")+
  scale_x_continuous(labels = scales::label_number(accuracy = 0.001))
ggsave("result/A4a_ppc_mean.png", width = 20, height = 18, units = "cm", dpi = 300)

ppc_stat_grouped(y = df1$qx, yrep = qx_rep, group = factor(df1$age), stat = "sd")+
  scale_x_continuous(labels = scales::label_number(accuracy = 0.001))
ggsave("result/A4b_ppc_sd.png", width = 20, height = 18, units = "cm", dpi = 300)

# RMSE
rmse <- sqrt(mean((df1$qx - qx_rep_mean)^2))
cat("RMSE:", round(rmse, 4), "\n")

# MAE
mae <- mean(abs(df1$qx - qx_rep_mean))
cat("MAE:", round(mae, 4), "\n")

# Coverage of 95% credible interval
coverage <- mean(df1$qx >= qx_rep_lower & df1$qx <= qx_rep_upper)
cat("95% CI coverage:", round(coverage * 100, 1), "%\n")

# LOO
log_lik    <- fit$draws("log_lik", format = "matrix")
loo_result <- loo(log_lik)
print(loo_result)

# trace plots
draws <- fit$draws(format = "array")
betas <- paste0("beta[", 1:length(cov_vars), "]")
mcmc_trace(draws, pars = betas, 
           facet_args = list(
             labeller = as_labeller(c(
               "beta[1]" = "gamma[g]",
               "beta[2]" = "gamma[m<=10]",
               "beta[3]" = "gamma[m>10]",
               "beta[4]" = "gamma[d]"
             ), label_parsed)
           ))+
  theme(
    strip.text = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )
ggsave("result/A5_traceplot.png", width = 20, height = 18, units = "cm", dpi = 300)
