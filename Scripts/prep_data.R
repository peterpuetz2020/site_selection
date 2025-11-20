# clear environment
rm(list = ls())

# in case not installed, install pacman package
if (!require("pacman"))
  install.packages("pacman")
# install (if not done yet) and load required packages and read in self-written
# functions
pacman::p_load(here)
source(here("Scripts", "setup.R"))

# read in data
df <- readRDS(here(read_data_here, "all_data.RDS")) %>% 
  # count share of values below lod
  group_by(uwwtd_code) %>% 
  mutate(mean_blod = mean(below_lod)) %>% 
  ungroup()

df <- df %>% 
  # use only values below limit of detection
  filter(below_lod == 0) %>% 
  mutate(logvalue = log10(value),
         Standort = uwwtd_code) %>% 
  # reformat date variable
  mutate(Anfang_Probenahme = as.Date(format(
    as.POSIXct(Anfang_Probenahme, format =
                 '%Y-%m-%d %H:%M:%S'),
    format = '%Y-%m-%d'
  ))) %>% 
  arrange(Standort, Anfang_Probenahme, Labor) 

# drop small sites
df <- df %>% 
  group_by(Standort, Labor) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n >= min_obs) 

# drop old site/lab combinations
df <- df %>% 
  group_by(Standort, Labor) %>% 
  nest() %>% 
  group_by(Standort) %>% 
  slice_tail(n = 1) %>% 
  unnest(cols = c(data))

# count NAs per gene and site/lab
df <- df %>%
  group_by(Standort, Labor) %>%
  mutate( across(logvalue,
                 .fns = list(non_nas =~sum(!is.na(.))),
                 .names = "{fn}_{col}" ) ) %>% 
  ungroup()

## compute loess values
df <- loess_fun("logvalue", "non_nas_logvalue")

df_var <- df %>%
  group_by(Standort, Labor) %>%
  mutate(
    # compute mean per site
    m_logvalue = mean(logvalue, na.rm = T),
    # calculate absolute deviations
    dev_abs_15 = abs(loess_logvalue  - logvalue),
    # calculate absolute deviations standardized by mean value
    dev_abs_15_ratio = dev_abs_15/m_logvalue,
  ) %>%
  # calculate "goodness" indicators
  summarise(
    median_abs_error_15 = median(dev_abs_15, na.rm = TRUE),
    median_abs_error_ratio_15 = median(dev_abs_15_ratio, na.rm = TRUE),
    n = n(),
    mean_blod = mean(mean_blod)
  ) %>%
  ungroup()  %>% 
# rename variable
 rename(mae = `median_abs_error_ratio_15`) %>% 
  dplyr::select(-contains("median"))

# save data
saveRDS(df_var,
        here(read_data_here, "mae_measures.RDS"))
