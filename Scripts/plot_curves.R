# in case not installed, install pacman package
if (!require("pacman"))
  install.packages("pacman")
# install (if not done yet) and load required packages and read in self-written
# functions
pacman::p_load(here, tidyverse)
source(here("Scripts", "setup.R"), encoding = "UTF-8")

# read full data
df <- readRDS(here(read_data_here, "data.rds")) %>% 
  mutate(date = as.Date(date),
         log_viral_load = log10(viral_load + 1)) %>% 
  # filter more recent dates as only then enough WWTPS are available
  filter(date >= "2022-10-15") %>% 
  group_by(uwwtd) %>% 
  # add intermediate dates
  pad() %>% 
  fill(lab, site_name, state, .direction = "downup") %>% 
  ungroup() %>% 
  arrange(date)

# create aggregated data for wochenbericht html file
# generate weeks starting on Thursday
thursday_data <-
  df %>%
  distinct(date) %>%
  mutate(
    Tag = lubridate::wday(date, week_start = 1),
    is_thursday = ifelse(Tag == 4, 1, 0),
    th_week = cumsum(is_thursday)
  ) %>%
  dplyr::select(date, th_week)

df_agg <- df %>%
  left_join(thursday_data)

# add dates with NAs before samples to avoid that 7-days-averages
new_rows <- df_agg %>%
  group_by(site_name) %>%
  summarise(min_date = min(date, na.rm = TRUE), .groups = 'drop') %>%
  rowwise() %>%
  do(data.frame(site_name = .$site_name, 
                date = seq(.$min_date - days(7), .$min_date - days(1), by = "day"))) %>%  
  ungroup()

# combine data
df_agg <- bind_rows(df_agg, new_rows) %>% 
  arrange(site_name, date) %>% 
  mutate(
    KW = lubridate::week(date),
    Jahr = lubridate::year(date)
  ) %>%
  group_by(site_name, th_week) %>%
  mutate_at(vars(contains("viral_load")), ~ mean(., na.rm = TRUE)) %>% 
  filter(!is.na(pop_covered), !is.na(log_viral_load)) %>% 
  filter(date == max(date, na.rm = TRUE)) %>%
  ungroup()

df_agg <- df_agg %>% 
  group_by(th_week) %>% 
  mutate(mean_log_viral_load = mean(log_viral_load)) %>% 
  ungroup() %>% 
  mutate(log_viral_load_dev = log_viral_load - mean_log_viral_load) %>% 
  group_by(site_name, lab) %>% 
  mutate(log_viral_load_dev = mean(log_viral_load_dev)) %>% 
  ungroup() %>% 
  mutate(log_viral_load = log_viral_load - log_viral_load_dev) %>% 
  select(-mean_log_viral_load, -log_viral_load_dev) %>% 
  mutate(Tag = lubridate::wday(date, week_start = 1))

# Create an empty list as placeholder
agg_list <- list()

# set seed for replicability
set.seed(22)

# compute (un-)weighted means over all sites
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
  arrange(all_sites, date)

# calculate loess estimates
pred <- df_agg %>%
  group_by(all_sites) %>%
  nest() %>%
  mutate(pred = purrr::map(data, ~ predict(
    loess.as(
      .x$obs[!is.na(.x$log_viral_load)],
      .x$log_viral_load[!is.na(.x$log_viral_load)],
      degree = 2,
      weights = sqrt(.x$weights[!is.na(.x$log_viral_load)]),
      user.span = .13
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
  add_column(
    loess_optimized = extract_prediction(lis = pred_list, extract = "fit"),
    loess_optimized_se = extract_prediction(lis = pred_list, extract = "se.fit"),
    loess_optimized_df = extract_prediction(pred_list, "df") %>%
      map2(., reps, ~ rep(.x, .y)) %>%
      unlist()
  ) %>%
  mutate(
    loess_optimized_pw_lb = loess_optimized - qt(0.975, loess_optimized_df) *
      loess_optimized_se,
    loess_optimized_pw_ub = loess_optimized + qt(0.975, loess_optimized_df) *
      loess_optimized_se
  ) %>%
  add_column(site_name = "Aggregiert") %>%
  dplyr::select(
    site_name,
    all_sites,
    n_non_na,
    date,
    log_viral_load,
    contains("loess_optimized"),
    -contains("d_df")
  ) %>%
  mutate_at(vars(contains("loess")), list(orig = ~ 10 ^ . - 1)) %>%
  mutate(viral_load = 10 ^ log_viral_load - 1)

# compare width of confidence bands
df_agg %>%  
  mutate(width = loess_optimized_pw_ub_orig - loess_optimized_pw_lb_orig) %>% 
  group_by(all_sites) %>% 
  summarise(m = mean(width)) %>% 
  mutate(ratio =  m[2] /  m[1])

df_agg <- df_agg %>%
  mutate_at(vars(viral_load, contains("orig")), ~./1000)

custom_labels <- c("ja" = "All sites", "nein" = "Sites selected for 2025")
Sys.setlocale("LC_TIME", "C")

p_trans <- df_agg %>%  
  ggplot(aes(x = date, y = viral_load)) +
  geom_ribbon(
    aes(ymin = loess_optimized_pw_lb_orig, ymax = loess_optimized_pw_ub_orig),
    fill = "lightblue",
    linemitre = 200
  ) +
  geom_point(colour = "grey") +
  geom_line(aes(date, 
                y = loess_optimized_orig)) +
  facet_wrap(~all_sites, ncol = 1,
             labeller = labeller(all_sites = custom_labels)) +
  scale_x_date(date_breaks = "2 month", 
               date_labels = "%b \n%y", expand = c(0.05,0.05)) +
  labs(y = expression(atop("Viral load in wastewater",
                           atop(paste("in gene copies / liter (in thousand)")))), 
       x = "Date") +
  theme_palatino

p_trans

ggsave(here(results_here, paste0("curves_compared.png")), p_trans, width = 5.25, height = 3.5, dpi = 300)
ggsave(here(results_here, paste0("curves_compared.svg")), p_trans, dpi = 300)
