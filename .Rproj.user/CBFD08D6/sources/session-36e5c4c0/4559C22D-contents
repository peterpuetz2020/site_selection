# set minimum number of observations required
min_obs = 15

#install packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, broom, fANCOVA, rio, padr)

read_data_here <- 
  here(here(), "Data")
results_here <-
  here(here(), "Results")

# function to compute loess values
loess_fun <- function(tar_var = "logvalue", span_var = "non_nas_logvalue")
{
  temp <- df %>%
    filter(!!sym(span_var)>=min_obs) %>% 
    group_by(Standort, Labor) %>%
    mutate(obs = row_number()) %>%
    mutate( across(!!sym(tar_var),
                   .fns = list(loess = ~predict(
                     loess.as(
                       obs[!is.na(.)],
                       .[!is.na(.)],
                       family = "gaussian",
                       degree = 2,
                       user.span = 15 / (!!sym(span_var))[1]
                     ),
                     newdata = data.frame(x = obs),
                   )),
                   .names = "{fn}_{col}" ) ) %>% 
    ungroup() %>% 
    dplyr::select(Standort, Labor, Anfang_Probenahme, paste0("loess_", tar_var)) 
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
      mutate(weighting_var = einwohner) else
        df <- df %>%
          mutate(weighting_var = 1)
      
      # set seed for replicability
      set.seed(22)
      
      df <- df %>%
        group_by(Jahr, KW) %>%
        # compute weights for loess curve, these are the inverse values of the variance
        # of the (weighted) mean of the observations
        mutate(weights = (1 / var_weighted(x = logvalue, wt = weighting_var))) %>%
        # count contributing sites per Wednesday
        mutate(n_non_na = sum(!is.na(value))) %>%
        # compute weighted means
        mutate_at(vars(contains("value")), # if at least a certain amount of sites provides data
                  ~ if (mean(n_non_na) < min_obs) {
                    NA
                  } else
                  {
                    weighted.mean(., (weighting_var), na.rm = TRUE)
                  }) %>%
        filter(!is.na(logvalue)) %>%
        mutate(
          n_non_na = mean(n_non_na, na.rm = TRUE),
          logvalue = mean(logvalue, na.rm = TRUE),
          anteil_bev = sum(einwohner, na.rm = TRUE) / pop,
          weights = mean(weights, na.rm = TRUE)
        ) %>%
        ungroup() %>%
        filter(Tag == 3) %>%
        distinct(th_week, .keep_all = T) %>%
        # standardize means
        mutate(weights = weights / mean(weights, na.rm = TRUE)) %>%
        arrange(Anfang_Probenahme) %>%
        dplyr::select(Anfang_Probenahme,
                      n_non_na,
                      logvalue,
                      anteil_bev,
                      weights
        ) %>%
        # expand for predictions
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