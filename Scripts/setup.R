# set minimum number of observations required
min_obs = 15

#install packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, broom, fANCOVA, rio, padr, data.table, sp, mapproj, stargazer,
               extrafont, showtext)
# only once
#font_import()
loadfonts(device = "win")  # Use "quartz" for Mac, "win" for Windows

read_data_here <- 
  here(here(), "Data")
results_here <-
  here(here(), "Results")

# function to compute loess values
loess_fun <- function(tar_var = "log_viral_load", span_var = "n_non_na")
{
  temp <- df %>%
    filter(!!sym(span_var)>=min_obs) %>% 
    group_by(site_name, lab) %>%
    mutate(obs = row_number()) %>%
    mutate(across(
      !!sym(tar_var),
      .fns = list(loess = ~ predict(
        loess.as(
          obs[!is.na(.)],
          .[!is.na(.)],
          family = "gaussian",
          degree = 2,
          user.span = 15 / (!!sym(span_var))[1]
        ),
        newdata = data.frame(x = obs)
      )),
      .names = "{fn}_{col}"
    )) %>% 
    ungroup() %>% 
    dplyr::select(site_name, lab, date, paste0("loess_", tar_var))
  
  df <- df %>% left_join(temp)
  return(df)
}

pop = 84000000

# function that computes variance of a weighted mean
var_weighted <- function(x = NULL, wt = NULL) {
  xm = weighted.mean(x, wt)
  return(sum(wt * (x - xm) ^ 2) / (sum(wt) - 1) / sum(wt))
}

# function that aggregates viral loads over all sites
aggregation <- function(df = df_agg,
                        weighting = TRUE) {
  
  # if weighting then use inhabitants as weights, else set weights to 1
  if (weighting)
    df <- df %>%
      mutate(weighting_var = pop_covered) 
  else
    df <- df %>%
      mutate(weighting_var = 1)
  #browser()
  # set seed for replicability
  df <- df %>%
    group_by(Jahr, KW) %>%
    # compute weights for loess curve
    mutate(weights = (1 / var_weighted(x = log_viral_load, wt = weighting_var))) %>%
    # count contributing sites per Wednesday
    mutate(n_non_na = sum(!is.na(viral_load))) %>%
    # compute weighted means
    mutate_at(vars(contains("value")), 
              ~ if (mean(n_non_na) < min_obs) {
                NA
              } else {
                weighted.mean(., (weighting_var), na.rm = TRUE)
              }) %>%
    filter(!is.na(log_viral_load)) %>%
    mutate(
      n_non_na = mean(n_non_na, na.rm = TRUE),
      log_viral_load = mean(log_viral_load, na.rm = TRUE),
      pop_share_country = sum(pop_covered, na.rm = TRUE) / pop,
      weights = mean(weights, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    filter(Tag == 3) %>%
    distinct(th_week, .keep_all = TRUE) %>%
    mutate(weights = weights / mean(weights, na.rm = TRUE)) %>%
    arrange(date) %>%
    dplyr::select(
      date,
      n_non_na,
      log_viral_load,
      pop_share_country,
      weights
    ) %>%
    pad(interval = "day") %>%
    mutate(obs = row_number())
  
  return(df)
}

# function to extract predictions with certain name from list
extract_prediction <- function(lis = NULL, extract = NULL) {
  extracted <- lis[sapply(names(lis), function(x)
    grepl(paste0("^", extract, "$"), x))] %>%
    unlist()
  return(extracted)
}

# fonts for ggplot output
windowsFonts(Palatino = windowsFont("Palatino Linotype"))

# Create a ggplot theme that uses it everywhere
theme_palatino <- theme_minimal(base_size = 10, base_family = "Palatino") +
  theme(
    plot.title      = element_text(family = "Palatino", size = 10),
    plot.subtitle   = element_text(family = "Palatino", size = 10),
    plot.caption    = element_text(family = "Palatino", size = 10),
    axis.title      = element_text(family = "Palatino", size = 10),
    axis.text       = element_text(family = "Palatino", size = 10),
    legend.title    = element_text(family = "Palatino", size = 10),
    legend.text     = element_text(family = "Palatino", size = 10),
    strip.text      = element_text(family = "Palatino", size = 10),
    panel.grid.minor = element_blank()
  )

# function to get R outliers
identify_outliers_IQR_factor <-
  function(x,
           dblfactor = 1.5,
           na.rm = TRUE,
           ...) {
    qnt <- quantile(x, probs = c(.25, .75), na.rm = na.rm, ...)
    H <- dblfactor * IQR(x, na.rm = na.rm)
    y <- rep(0, length(x))
    y[x < (1 / (qnt[2] + H))] <- 1
    y[x > (qnt[2] + H)] <- 1
    y
  }
