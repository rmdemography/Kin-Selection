#------------------------------------------------------------------------------#
# Paper:  Kin selection and population dynamics
# Title:  Reconstruct q0 from regression coefficients
# Data:   HMD & HFD
# Author: Rahul Mondal
# Date:   27/03/2026
#------------------------------------------------------------------------------#
source("02_grandmother_effect.R")

# extract fixed effects via lm()
fe3_lm <- lm(logit_qx ~ gm + os + age_m + factor(country) + factor(Year),
             data = df1)
b_lm <- coef(fe3_lm)[c("gm", "os", "age_m")] # coefficients 

# country-specific fixed effects
df_fit <- df1 %>% mutate(country = factor(country))
country_levels     <- levels(factor(df1$country))
ref_country        <- country_levels[1]
country_coef_names <- paste0("factor(country)", country_levels[-1])
country_coef_vals  <- coef(fe3_lm)[country_coef_names]
names(country_coef_vals) <- country_levels[-1]

# year-specific fixed effects
year_levels     <- levels(factor(df1$Year))
ref_year        <- year_levels[1]
year_coef_names <- paste0("factor(Year)", year_levels[-1])
year_coef_vals  <- coef(fe3_lm)[year_coef_names]
names(year_coef_vals) <- year_levels[-1]

df_fit <- df_fit %>%
  mutate(
    # lm fitted values
    fitted_logit_lm = fitted(fe3_lm),
    fitted_qx_lm    = 1 / (1 + exp(-fitted_logit_lm)),
    resid_logit_lm  = logit_qx - fitted_logit_lm,
    resid_qx_lm     = qx - fitted_qx_lm,
    # intercept goes into country_fe
    country_fe = coef(fe3_lm)["(Intercept)"] +
      ifelse(
        as.character(country) == ref_country, 0,
        country_coef_vals[as.character(country)]
      ),
    year_fe = ifelse(
      as.character(Year) == ref_year, 0,
      year_coef_vals[as.character(Year)]
    ),
    fitted_logit_manual = country_fe + year_fe +
      b_lm["gm"]*gm + b_lm["os"]*os + b_lm["age_m"]*age_m,
    fitted_qx_manual    = 1 / (1 + exp(-fitted_logit_manual))
  )
# check model fit
diff_lm <- df_fit$fitted_logit_lm - df_fit$fitted_logit_manual
cat("\n── Verification (lm fitted vs manual) ──\n")
cat(sprintf("Max absolute difference:  %.2e\n", max(abs(diff_lm))))

cat("\nSanity check — first 10 rows:\n")
print(dplyr::select(df_fit, country, Year, qx,
                    fitted_logit_manual, fitted_qx_manual, resid_qx_lm) %>%
        head(10), digits = 5)

cat("\nOverall fit (qx scale):\n")
cat(sprintf("  RMSE: %.6f\n", sqrt(mean(df_fit$resid_qx_lm^2))))
cat(sprintf("  MAE:  %.6f\n", mean(abs(df_fit$resid_qx_lm))))
cat(sprintf("  Cor:  %.6f\n", cor(df_fit$qx, df_fit$fitted_qx_lm)))

#------------------------------------------------------------------------------#
# construct list of parameters
#------------------------------------------------------------------------------#

# intercept
intercept <- coef(fe3_lm)["(Intercept)"]

# betas
betas <- coef(fe3_lm)[c("gm", "os", "age_m")]

# country fixed effects
country_level <- levels(factor(df1$country))
ref_country <- country_level[1]
country_fe <- c(
  setNames(0, ref_country),
  coef(fe3_lm)[paste0("factor(country)", country_level[-1])]
)
names(country_fe) <- country_level

# year fixed effects
year_level <- levels(factor(df1$Year))
ref_year <- year_level[1]
year_fe <- c(
  set_names(0, ref_year),
  coef(fe3_lm)[paste0("factor(Year)", year_level[-1])]
)
names(year_fe) <- year_level

# parameter list
params <- list(
  intercept = intercept,
  betas = betas,
  deltas = country_fe,
  kappas = year_fe
)
