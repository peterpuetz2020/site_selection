# in case not installed, install pacman package
if (!require("pacman"))
  install.packages("pacman")
# install (if not done yet) and load required packages and read in self-written
# functions
pacman::p_load(here, tidyverse)
source(here("Scripts", "setup.R"), encoding = "UTF-8"
)

# read full data
df <- readRDS(here(read_data_here, "all_data.rds")) %>% 
  mutate(Anfang_Probenahme = as.Date(Anfang_Probenahme),
         logvalue = log10(value + 1),
         loq = as.logical(below_lod)) %>% 
  # filter more recent dates as only then enough WWTPS are available
  filter(Anfang_Probenahme >= "2022-11-01",
         Anfang_Probenahme <= "2024-07-15") %>% 
  group_by(uwwtd_code) %>% 
  # add intermediate dates
  pad() %>% 
  fill(Labor, Standort, Bundesland, .direction = "downup")
  
# read in sites that are still in amelag
amelag_sites <-
  import(here(read_data_here, "Standort-Übersicht_29.01.2025.xlsx"),
         sheet = 1,
         setclass = "tibble") %>% 
  select(uwwtd_code = 'UWWTD-Code', amelag = 'AMELAG-finanziert') %>% 
  drop_na()

# select data from sites participating in 2025
df_join <- 
  df %>% 
  left_join(amelag_sites) %>% 
  # indicate 2025 selected sites
  mutate(kept_sites = ifelse(is.na(amelag) | amelag == "nein", "nein", "ja")) %>% 
  select(-amelag)

# create aggregated data for wochenbericht html file
# generate weeks starting on Thursday
thursday_data <-
  df_join %>%
  distinct(Anfang_Probenahme) %>%
  mutate(
    Tag = lubridate::wday(Anfang_Probenahme, week_start = 1),
    is_thursday = ifelse(Tag == 4, 1, 0),
    th_week = cumsum(is_thursday)
  ) %>%
  dplyr::select(Anfang_Probenahme, th_week)

df_agg <- df_join %>%
  left_join(thursday_data) %>% 
  group_by(Standort) %>%
  # complete data
  mutate(loq = ifelse(is.na(logvalue),FALSE, loq)) %>%
  ungroup() 

# add dates with NAs before samples to avoid that 7-days-averages
# drop values if no previous dates are available
new_rows <- df_agg %>%
  group_by(Standort) %>%
  summarise(min_date = min(Anfang_Probenahme, na.rm = TRUE), .groups = 'drop') %>%
  # Create a sequence of dates for the 7 days before the minimum date
  rowwise() %>%
  do(data.frame(Standort = .$Standort, 
                Anfang_Probenahme = seq(.$min_date - days(7), .$min_date - days(1), by = "day"))) %>%  # Set value to NA for new rows
  ungroup()

# combine data
df_agg <- bind_rows(df_agg, new_rows) %>% 
  arrange(Standort, Anfang_Probenahme) %>% 
  mutate(
    KW = lubridate::week(Anfang_Probenahme),
    Jahr = lubridate::year(Anfang_Probenahme)
  ) %>%
  group_by(Standort, th_week) %>%
  # for each site and virus, compute 7-day averages
  mutate_at(vars(contains("value")), ~ mean(., na.rm = TRUE)) %>%
  # remove missings
  filter(!is.na(einwohner), !is.na(logvalue)) %>% 
  # take only one value per week, site
  filter(Anfang_Probenahme == max(Anfang_Probenahme, na.rm = TRUE)) %>%
  ungroup() 

df_agg <- df_agg %>% 
  # calculate unweighted means over the weeks as these unweighted means can be used
  # to calculate differences between site/lab combination from these means - in this
  # way, site/lab specific level regimes can be accounted for
  group_by(th_week) %>% 
  mutate(mean_logvalue = mean(logvalue)) %>% 
  ungroup() %>% 
  # calculate differences from these means
  mutate(logvalue_dev = logvalue-mean_logvalue) %>% 
  # average over these for each virus/site/lab combination
  group_by(Standort, Labor) %>% 
  mutate(logvalue_dev = mean(logvalue_dev)) %>% 
  ungroup() %>% 
  # adjust for deviations
  mutate(logvalue = logvalue - logvalue_dev) %>% 
  # drop variables
  select(-mean_logvalue, -logvalue_dev) %>% 
  # add day
  mutate(Tag = lubridate::wday(Anfang_Probenahme, week_start = 1)
)

# Create an empty list as placeholder
agg_list <- list()

# set seed for replicability
set.seed(22)

# compute (un-)weighted means over all sites for each pathogen
agg_list[[1]] <-
    aggregation(df = df_agg,
                weighting = TRUE) %>% 
    mutate(all_sites = "ja")
  
# same only for amelag sites
  agg_list[[2]] <-
    aggregation(df = df_agg %>% filter(kept_sites == "ja"),
                weighting = TRUE)  %>% 
    mutate(all_sites = "nein")
  
  # combine datasets
df_agg <- map_dfr(agg_list, bind_rows) %>%
  # important
  arrange(all_sites, Anfang_Probenahme)

# calculate loess estimates, see help(loess.as) for details of
# the set options
pred <- df_agg %>%
  group_by(all_sites) %>%
  nest() %>%
  mutate(pred = map(data, ~ predict(
    loess.as(
      .x$obs[!is.na(.x$logvalue)],
      .x$logvalue[!is.na(.x$logvalue)],
      #criterion = "gcv",
      #family = "symmetric",
      degree = 2,
      weights = sqrt(.x$weights[!is.na(.x$logvalue)]),
      #control = loess.control(surface = "direct")
      user.span = .125
    ),
    newdata = data.frame(x = .x$obs),
    se = TRUE))) %>%
  dplyr::select(pred) %>%
  unnest(cols = c(pred))

# store list
pred_list <- pred[, "pred"]$pred

# store number of observations per group
reps <- df_agg %>%
  group_by(all_sites) %>%
  summarise(n = n()) %>%
  pull(n)

df_agg <- df_agg %>%
  # add columns relevant for predictions
  add_column(
    loess_optimized = extract_prediction(lis = pred_list, extract = "fit"),
    loess_optimized_se = extract_prediction(lis = pred_list, extract = "se.fit"),
    loess_optimized_df = extract_prediction(pred_list, "df") %>%
      map2(., reps, ~ rep(.x, .y)) %>%
      unlist()
  ) %>%
  mutate(
    # add pointwise confidence bands
    loess_optimized_pw_lb = loess_optimized - qt(0.975, loess_optimized_df) *
      loess_optimized_se,
    loess_optimized_pw_ub = loess_optimized + qt(0.975, loess_optimized_df) *
      loess_optimized_se
  ) %>%
  add_column(Standort = "Aggregiert") %>%
  dplyr::select(
    Standort,
    all_sites,
    n_non_na,
    Anfang_Probenahme,
    logvalue,
    contains("loess_optimized"),
    anteil_bev,
    -contains("d_df")
  ) %>%
  # back-transform to original scale
  mutate_at(vars(contains("loess")), list(orig = ~ 10 ^ . - 1)) %>%
  mutate(value = 10 ^ logvalue - 1)

df_agg <- df_agg %>%
  mutate_at(vars(value, contains("orig")), ~./1000)

custom_labels <- c("ja" = "All sites", "nein" = "Sites selected for 2025")
Sys.setlocale("LC_TIME", "C")  # Forces English for dates and times
p_trans <- df_agg %>%  
  filter(!is.na(value)) %>% 
  ggplot(aes(x = Anfang_Probenahme, y = value)) +
  geom_ribbon(
    aes(ymin = loess_optimized_pw_lb_orig, ymax = loess_optimized_pw_ub_orig),
    # shadowing cnf intervals
    fill = "lightblue"
  ) +
  geom_point(colour = "grey") +
  geom_line(aes(Anfang_Probenahme, 
                y = loess_optimized_orig),
            linewidth = 1) +
  geom_line(aes(Anfang_Probenahme, y = loess_optimized_pw_lb_orig),linewidth = 1, linetype = 0) +
  geom_line(aes(Anfang_Probenahme, y = loess_optimized_pw_ub_orig),linewidth = 1, linetype = 0) +
  facet_wrap(~all_sites, ncol = 1,
             labeller = labeller(all_sites = custom_labels)) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 15),  # Increase facet title size
    # Increase axis title size
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    # Increase axis tick label size
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
  ) +
  scale_x_date(date_breaks = "2 month", 
               date_labels = "%b %y") +
  labs(y = expression(atop("Viral load in wastewater",atop(paste("in gene copies / liter (in thousand)")))), 
       x = "Date")
p_trans

# compare average confidence band with
ggsave(here(results_here, paste0("curves_compared.svg")), p_trans, width = 12, height = 8)
ggsave(here(results_here, paste0("curves_compared.png")), p_trans, width = 12, height = 8)
