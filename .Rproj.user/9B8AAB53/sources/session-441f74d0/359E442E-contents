#------------------------------------------------------------------------------#
# Paper:  Kin selection and population dynamics
# Title:  Model for grandmother's effect on infant survival
# Data:   HMD & HFD
# Author: Rahul Mondal
# Date:   27/03/2026
#------------------------------------------------------------------------------#
pacman::p_load(plm, lmtest, sandwich, patchwork, modelsummary, ggplot2)
source("01_compute_kin.R")

df1 <- df1 %>% 
  rename(
    "gm" = "count_living_gm", 
    "os" = "count_living_os", 
    "age_m" = "mean_age_m"
    ) %>% 
  mutate(logit_qx = log(qx/(1-qx)),
         trend   = as.integer(Year) - 1950L) %>% 
  arrange(country, Year)

# make panel data object
pdata <- pdata.frame(df1, index = c("country", "Year"))

# Model 1: grandmothers only
fe1 <- plm(logit_qx ~ gm,
           data = pdata, model = "within", effect = "twoways")

# Model 2: + siblings
fe2 <- plm(logit_qx ~ gm + os,
           data = pdata, model = "within", effect = "twoways")

# Model 3: + mean age of mother (full model)
fe3 <- plm(logit_qx ~ gm + os + age_m,
           data = pdata, model = "within", effect = "twoways")

# Clustered standard errors (by country)
cl_se <- function(model) coeftest(model, vcov = vcovHC(model, type = "HC1", cluster = "group"))

cat("\n=== Fixed Effects Model 1 (grandmothers only) ===\n"); print(cl_se(fe1))
cat("\n=== Fixed Effects Model 2 (+ os) ===\n");        print(cl_se(fe2))
cat("\n=== Fixed Effects Model 3 (full) ===\n");              print(cl_se(fe3))

# Diagnostic test
# 1. Hausman test
re3 <- plm(logit_qx ~ gm + os + age_m, data = pdata, model = "random",
           random.method = "amemiya")
hausman <- phtest(fe3, re3)
cat("Decision:", ifelse(hausman$p.value < 0.05,
                        "Reject RE → Fixed Effects appropriate.",
                        "Fail to reject RE → RE may be preferred."), "\n")
# 2. Serial Correlation (Wooldridge/Breusch-Godfrey)
sc_test <- pbgtest(fe3)
print(sc_test)
cat("Decision:", ifelse(sc_test$p.value < 0.05,
                        "Serial correlation present → clustered SEs justified.",
                        "No evidence of serial correlation."), "\n")
# 3. Heteroskedasticity (Breusch-Pagan)
bp_test <- bptest(fe3)
print(bp_test)
cat("Decision:", ifelse(bp_test$p.value < 0.05,
                        "Heteroskedasticity present → robust SEs justified.",
                        "No evidence of heteroskedasticity."), "\n")
# 4. Cross-Sectional Dependence (Pesaran CD)
cd_test <- pcdtest(fe3, test = "cd")
print(cd_test)
cat("Decision:", ifelse(cd_test$p.value < 0.05,
                        "Cross-sectional dependence detected → Driscoll-Kraay SEs used in R1.",
                        "No significant cross-sectional dependence."), "\n")
# 5. Unit Root Tests (ADF per country per variable)
vars_to_test <- c("logit_qx", "gm", "os", "age_m")
countries_v  <- levels(df1$country)
adf_rows     <- list()

for (ctry in countries_v) {
  for (vname in vars_to_test) {
    x      <- df1 %>% filter(country == ctry) %>% arrange(Year) %>% pull(all_of(vname))
    result <- tryCatch(adf.test(x), error = function(e) NULL)
    adf_rows[[length(adf_rows) + 1]] <- data.frame(
      country    = ctry,
      variable   = vname,
      adf_stat   = if (!is.null(result)) round(result$statistic[[1]], 3) else NA,
      p_value    = if (!is.null(result)) round(result$p.value, 4) else NA,
      stationary = if (!is.null(result)) ifelse(result$p.value < 0.05, "Yes", "No (unit root)") else "Error",
      stringsAsFactors = FALSE
    )
  }
}

adf_results <- do.call(rbind, adf_rows)
print(adf_results, row.names = FALSE)
cat("\nNote: If unit roots present in levels, first-difference model (R3) addresses this.\n")

# 6. F-tests for FE significance
fe_ind  <- plm(logit_qx ~ gm + os + age_m, data = pdata, model = "within", effect = "individual")
fe_time <- plm(logit_qx ~ gm + os + age_m, data = pdata, model = "within", effect = "time")
pool    <- plm(logit_qx ~ gm + os + age_m, data = pdata, model = "pooling")

cat("Country FE:"); print(pFtest(fe_ind,  pool))
cat("Year FE:");    print(pFtest(fe_time, pool))

# Robustness checks
# 1. Driscoll-Kraay / SCC SEs
vcov_dk    <- vcovSCC(fe3, type = "HC1", cluster = "group")
fe_dk_cl   <- coeftest(fe3, vcov = vcov_dk)
print(fe_dk_cl)
cat("(Spatial correlation consistent SEs, robust to cross-sectional dependence)\n")

# 2. Country-Specific Linear Time Trends
fe_trends <- plm(logit_qx ~ gm + os + age_m + trend,
                 data = pdata, model = "within", effect = "twoways")
fe_trends2 <- lm(logit_qx ~ gm + os + age_m +
                   country:trend +          # country-specific slopes
                   factor(country) +        # country FE
                   factor(Year),            # year FE
                 data = df1)

vcov_cl2     <- vcovCL(fe_trends2, cluster = ~country)
fe_trends_cl <- coeftest(fe_trends2, vcov = vcov_cl2)
cat("Key predictors:\n")
rows_keep <- rownames(fe_trends_cl)[rownames(fe_trends_cl) %in% c("gm", "os", "age_m")]
print(fe_trends_cl[rows_keep, ])

# 3. First-Differences Model
fe_fd    <- plm(logit_qx ~ gm + os + age_m, data = pdata, model = "fd")
fe_fd_cl <- coeftest(fe_fd, vcov = vcovHC(fe_fd, type = "HC1", cluster = "group"))
print(fe_fd_cl)

# 4. Omitted Variable Check (gm excluded)
fe_no_gm    <- plm(logit_qx ~ os + age_m, data = pdata, model = "within", effect = "twoways")
fe_no_gm_cl <- coeftest(fe_no_gm, vcov = vcovHC(fe_no_gm, type = "HC1", cluster = "group"))
print(fe_no_gm_cl)
cat("Note: os coefficient inflates when gm is excluded → gm absorbs part of sibling effect.\n")

# 5. Sub-Period Analysis
for (period in list(c(1950L, 1989L), c(1990L, 2023L))) {
  df_sub  <- df1 %>% filter(Year >= period[1], Year <= period[2])
  pd_sub  <- pdata.frame(df_sub, index = c("country", "Year"))
  m_sub   <- plm(logit_qx ~ gm + os + age_m, data = pd_sub,
                 model = "within", effect = "twoways")
  cl_sub  <- coeftest(m_sub, vcov = vcovHC(m_sub, type = "HC1", cluster = "group"))
  cat(sprintf("\nPeriod %d–%d:\n", period[1], period[2]))
  print(cl_sub[c("gm", "os", "age_m"), ])
}

