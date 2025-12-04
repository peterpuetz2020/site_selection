# clear environment
rm(list = ls())

# in case not installed, install pacman package
if (!require("pacman"))
  install.packages("pacman")
# install (if not done yet) and load required packages and read in self-written
# functions
pacman::p_load(here)
source(here("Scripts", "setup.R"))

# set relevant variables which should be used for assessing treatment plants
rel_vars <- c("riip_share", "mae", "bloq_share")

# read in data
df <- readRDS(here(read_data_here, "data.RDS"))

# share of values below limit of quantification ---------------------------

# count share of values below loq
df <- df %>% 
  group_by(uwwtd) %>% 
  mutate(bloq_share = mean(below_loq)) %>% 
  ungroup()

# MAE calculation ---------------------------------------------------------
df <- df %>% 
  # use only values below limit of detection
  filter(below_loq == 0) %>% 
  mutate(logvalue = log10(viral_load)) %>% 
  # reformat date variable
  mutate(date = as.Date(format(
    as.POSIXct(date, format = '%Y-%m-%d %H:%M:%S'),
    format = '%Y-%m-%d'
  ))) %>% 
  arrange(site_name, date, lab)

# drop small sites
df <- df %>% 
  group_by(site_name, lab) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n >= min_obs)

# drop old site/lab combinations
df <- df %>% 
  group_by(site_name, lab) %>% 
  nest() %>% 
  group_by(site_name) %>% 
  slice_tail(n = 1) %>% 
  unnest(cols = c(data))

# count NAs per site/lab
df <- df %>%
  group_by(site_name, lab) %>%
  mutate(across(
    logvalue,
    .fns = list(non_nas = ~ sum(!is.na(.))),
    .names = "{fn}_{col}"
  )) %>% 
  ungroup()

## compute loess values
df <- loess_fun("logvalue", "non_nas_logvalue")

df_var <- df %>%
  group_by(site_name, lab) %>%
  mutate(
    # compute mean per site
    m_logvalue = mean(logvalue, na.rm = TRUE),
    # calculate absolute deviations
    dev_abs_15 = abs(loess_logvalue - logvalue),
    # calculate absolute deviations standardized by mean value
    dev_abs_15_ratio = dev_abs_15 / m_logvalue,
  ) %>%
  summarise(
    median_abs_error_15 = median(dev_abs_15, na.rm = TRUE),
    #median_abs_error_ratio_15 = median(dev_abs_15_ratio, na.rm = TRUE),
    n = n()
  ) %>%
  ungroup() %>% 
  # rename variable
  rename(mae = median_abs_error_15) %>% 
  select(-contains("median"))

# merge data
df <- df %>% 
  left_join(df_var)

# riip share calculation --------------------------------------------------
temp <- df %>%
  group_by(site_name, lab) %>%
  mutate(
    date_lag = as.numeric(date - lag(date, 1)),
    R = (viral_load / lag(viral_load))^(1 / date_lag),
    # get R outliers
    R_outlier = identify_outliers_IQR_factor(R),
    # # avoid artificial double R outliers
    # R_outlier = accumulate(R_outlier, 
    #                                     ~ if (.x == 1) 0 else .y),
    # get PV outliers
    PV_outlier = ifelse(
      (R > quantile(R, .75, na.rm = TRUE) & lead(R) < quantile(R, .25, na.rm = TRUE)) |
        (R < quantile(R, .25, na.rm = TRUE) & lead(R) > quantile(R, .75, na.rm = TRUE)),
      1, 0
    ),
    is_outlier = ifelse(R_outlier == 1 | PV_outlier == 1, 1, 0)#,
   # is_outlier = ifelse(is.na(is_outlier), 0, is_outlier)
  ) %>%
  summarise(
    riip_share = mean(is_outlier, na.rm = TRUE)
  )

df <- df %>%
  left_join(temp) %>%
  select(
    site_name, lab, pop_covered, state,
    pop_share, mae, riip_share, bloq_share, kept_sites,
    nominal_load
  ) %>%
  distinct(site_name, .keep_all = TRUE)

# assign ranks
df <- df %>% 
  mutate(across(
    all_of(rel_vars),
    .fns = list(rank = ~ rank(-.)),
    .names = "{col}_{fn}"
  ))

# compute (weighted) scores
df <- df %>%
  rowwise() %>%
  mutate(
    rank_mean = mean(c_across(all_of(
      paste0(rel_vars, "_rank")
    )))
  ) %>%
  ungroup()

df %>% saveRDS(here(read_data_here, "quality_data.RDS"))
