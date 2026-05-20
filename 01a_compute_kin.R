#------------------------------------------------------------------------------#
# Paper:  Kin selection and population dynamics
# Title:  Program for kinship network construction
# Data:   WPP 2024
# Author: Rahul Mondal
# Date:   08/04/2026
#------------------------------------------------------------------------------#
rm(list = ls())
pacman::p_load(tidyverse, readxl, httr, purrr, DemoKin, ggplot2, viridis)
options(timeout = 900)

# select countries
countries <- c("Bangladesh", "Ethiopia", "Gambia", "India", "Japan")

# asfr
f.url <- "https://population.un.org/wpp/assets/Excel%20Files/1_Indicator%20(Standard)/EXCEL_FILES/3_Fertility/WPP2024_FERT_F01_FERTILITY_RATES_BY_SINGLE_AGE_OF_MOTHER.xlsx"
temp.file <- tempfile(fileext = ".xlsx")
on.exit(unlink(temp.file))
GET(f.url,write_disk(temp.file, overwrite = TRUE))
excel_sheets(temp.file)
asfr <- read_xlsx(temp.file)
colnames(asfr) <- asfr[12,]
asfr <- asfr[-c(1,2,4:10)]
colnames(asfr)[1:2] <- c("country", "year")
asfr <- asfr %>% 
  filter(country %in% countries) %>%
  pivot_longer(!c(1:2), names_to = "age", values_to = "fx")
asfr[2:4] <- lapply(asfr[2:4], as.numeric)

# life table
lt.url <- "https://population.un.org/wpp/assets/Excel%20Files/1_Indicator%20(Standard)/EXCEL_FILES/4_Mortality/WPP2024_MORT_F06_3_SINGLE_AGE_LIFE_TABLE_ESTIMATES_FEMALE.xlsx"
temp.file <- tempfile(fileext = ".xlsx")
on.exit(unlink(temp.file))
GET(lt.url,write_disk(temp.file, overwrite = TRUE))
sheets <- excel_sheets(temp.file)
sheets <- sheets[-length(sheets)]
lt <- map_dfr(sheets, ~ read_xlsx(temp.file, sheet = .x))
colnames(lt) <- lt[12,]
lt <- lt[c(3,11:12,15)]
colnames(lt) <- c("country", "year", "age", "qx")
lt <- lt %>% filter(country %in% countries)
lt[2:4] <- lapply(lt[2:4], as.numeric)

df <- lt %>% 
  left_join(asfr, by = c("country", "year", "age")) %>% 
  mutate(Sx = 1-qx,
         fx = ifelse(is.na(fx), 0, fx/1000)) %>% 
  arrange(country, year, age)

df2matlist <- function(df, vital.rate = c("fx", "Sx")){
  lapply(setNames(vital.rate, vital.rate), function(v){
    mats <- tapply(df[[v]], list(df$age, df$year, df$country), identity)
    asplit(mats, 3)
  })
}
matlist <- df2matlist(df)

kin_df <- lapply(names(matlist$Sx), function(ctry){
  kin_out <- kin(p = matlist$Sx[[ctry]], f = matlist$fx[[ctry]], 
                 time_invariant = FALSE)
  ma <- kin_out$kin_summary %>% 
    filter(kin=="m") %>% 
    select(age_focal, year, mean_age) %>% 
    rename("ma" = "mean_age")
  gm.os.m <- kin_out$kin_full %>% 
    group_by(year, age_focal) %>% 
    summarise(
      gm = sum(living[kin == "gm"]),
      os = sum(living[age_kin <=10 & kin == "os"]),
      oso = sum(living[age_kin > 10 & kin == "os"]),
      m = sum(living[kin == "m"]),
      .groups = "drop"
    ) %>% 
    left_join(ma, by = c("age_focal", "year")) %>% 
    mutate(country = ctry)
}) %>% 
  bind_rows() %>% 
  select(country, everything())


