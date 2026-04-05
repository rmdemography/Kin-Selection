#------------------------------------------------------------------------------#
# Paper:  Kin selection and population dynamics
# Title:  Program for kinship network construction
# Data:   HMD & HFD
# Author: Rahul Mondal
# Date:   27/03/2026
#------------------------------------------------------------------------------#
rm(list = ls())
pacman::p_load(tidyverse, HMDHFDplus, DemoKin)

# credentials
cat("Enter HMD username:\n"); hmd_us <- userInput(silent = TRUE)
cat("Enter HMD password:\n"); hmd_pw <- userInput(silent = TRUE)
cat("Enter HFD username:\n"); hfd_us <- userInput(silent = TRUE)
cat("Enter HFD password:\n"); hfd_pw <- userInput(silent = TRUE)

# countries
countries <- c("SWE", "JPN", "FRATNP")

# construct kinship network by country
consturct_kin <- function(country, hmd_us, hmd_pw, hfd_us, hfd_pw){
  
  message("Loading life table (HMD):", country)
  lt <- readHMDweb(CNTRY = country, item = "fltper_1x1", 
                   username = hmd_us, password = hmd_pw)
  lt <- lt %>% 
    filter(Year>=1950 & Year<=2023) %>% 
    mutate(Sx = 1-qx) %>% 
    arrange(Year, Age)
  
  message("Loading ASFR (HFD):", country)
  asfr <- readHFDweb(CNTRY = country, item = "asfrRR",
                     username = hfd_us, password = hfd_pw)
  # overlapping years
  common_years <- intersect(unique(lt$Year), unique(asfr$Year))
  message("Common years: ", min(common_years), "-",  max(common_years),
          "(n = ", length(common_years), ")")
  df <- asfr %>% 
    right_join(lt, by = c("Year", "Age")) %>% 
    filter(Year %in% common_years) %>% 
    mutate(ASFR = ifelse(is.na(ASFR), 0, ASFR)) %>% 
    arrange(Year, Age)
  
  na <- length(unique(df$Age))
  nt <- length(unique(df$Year))
  years <- sort(unique(df$Year))
  # vital rate matrices for DemoKin input: Sx[a,t] and fx[a,t]
  Sx <- matrix(df$Sx, nrow = na, ncol = nt, dimnames = list(NULL, years))
  fx <- matrix(df$ASFR, nrow = na, ncol = nt, dimnames = list(NULL, years))
  
  message("Constructing kin for:", country)
  kin_out <- kin(p = Sx, f = fx, time_invariant = FALSE)
  list(kin = kin_out, df = df, Sx = Sx, fx = fx)
}

# output for countries
result <- list()
for(cty in countries){
  result[[cty]] <- consturct_kin(
    country = cty,
    hmd_us = hmd_us,
    hmd_pw = hmd_pw,
    hfd_us = hfd_us,
    hfd_pw = hfd_pw
  )
}

kin_df <- imap_dfr(result, function(res, cty) {
  res$kin$kin_summary %>% mutate(country = cty)
})

# prepare data for infant survival model
kin_df <- kin_df %>% 
  filter(kin %in% c("gm","m","os") & age_focal == 0) %>% 
  dplyr::select(-c(age_focal, cohort, count_cum_dead, mean_age_lost)) %>% 
  pivot_wider(
    id_cols = c(country, year),
    names_from = kin,
    values_from = c(count_living, mean_age, sd_age, count_dead)
  )

# combined life table and asfr dataframe
df <- imap_dfr(result, function(res,cty){
  res$df %>% mutate(country = cty)
})

# input data for grandmother model
df1 <- df %>% 
  filter(Age==0) %>% 
  full_join(kin_df, by = c("country", "Year" = "year")) %>% 
  dplyr::select(country, Year, qx, count_living_gm, count_living_os, mean_age_m)
