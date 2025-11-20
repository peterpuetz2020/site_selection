# clear environment
rm(list = ls())

# in case not installed, install pacman package
if (!require("pacman"))
  install.packages("pacman")
# install (if not done yet) and load required packages and read in self-written
# functions
pacman::p_load(here)
source(here("Scripts", "setup.R"))

# options to set ----------------------------------------------------------
# should only amelag treatment plants be selected? Maria and Peter suggest this.
amelag_only = TRUE

# Number of sites to be selected
no_selected <- 161

# Minimum of Nominalbelastung required (once 10,000 was suggested, but
# this is not mandatory)
min_nominalbelastung = 5

# set relevant variables which should be used for assessing treatment plants
rel_vars <- c("riip_share", "mae", "blod_share", "Nominalbelastung")

# set weights for all relevant variables (as ordered above)
weights = c(1, 1, 1, 0)

# set relevant variables that should be reverted (as more is "better")
revert_vars <- c("Nominalbelastung")

# should ranks or categories (1, 2, 3) be used for the criteria which are
# weighted (as mae, riip share and Nominalbelastung)
# use "rank" or "category"
score = "rank"

# data imports and merges -------------------------------------------------
#read prepared data
df_nuts <-
  import(
    here(read_data_here, "20240726_all_plants.csv"),
    encoding = "UTF-8",
    dec=",", sep=";"
  ) |>
  dplyr::select(
    uwwtd_code = pl_eu_key,
    Anlage = pl_name,
    pl_nuts,
    Bundesland = bl_code,
    Nominalbelastung = pl_nload,
    pl_long, 
    pl_lat,
    pl_koord_act
  ) |>
  filter(Bundesland != "EXTRA", Nominalbelastung > min_nominalbelastung)

# amelag sites
existing_sites <-
  import(here(read_data_here, "Standortliste_merged.xlsx"))

df  <- df_nuts  |>
  rename(uwwtd=uwwtd_code) %>% 
  mutate(amelag = ifelse(uwwtd %in% existing_sites$uwwtd, "yes", "no")) %>% 
  as_tibble()

# merge not perfect!
table(df$amelag) #161 sites in Amelag - difference due to 

# inspect non-merging sites
control <- anti_join(existing_sites, df, by="uwwtd")

## adjust uwwtd codes to allow for merging
# remove _KAN-1 and _KAN-2
control$uwwtd <- gsub("_KAN-1|_KAN-2", "", control$uwwtd)
# remove spaces in original dataframe
df$uwwtd <- gsub(" ", "", df$uwwtd)

# join and control after fix
control_df <- inner_join(df,control, by="uwwtd")
control2 <- anti_join(control, control_df, by="uwwtd")

# Introduction of variable uwwtd2, which contains additional 
# information to enable merging at sewage treatment plant level 
# (uwwtd2 only occurs if different from uwwtd)
new_df <- control_df %>%
  left_join(existing_sites, by = "Standort") %>% 
  rename(uwwtd2=uwwtd.y) %>% 
  rename(uwwtd=uwwtd.x) %>% 
  rename(einwohner=einwohner.x) %>% 
  dplyr::select(Standort,uwwtd,uwwtd2, einwohner)

new_df <- new_df %>% 
  inner_join(df, by="uwwtd")

# remove all rows apparent in new_df
df <- df %>%
  anti_join(new_df, by = "uwwtd")
# now add the missing rows
df <- bind_rows(df, new_df)

# add the results
df  <- df  |>
  mutate(amelag = case_when(uwwtd %in% control_df$uwwtd ~ "yes",
                            uwwtd %in% existing_sites$uwwtd ~ "yes",
                            TRUE~"no"))
table(df$amelag)

# select only amelag treatment plants (if set to TRUE above)
if (amelag_only)
  df  <- df  |>
  filter(amelag == "yes")

# fill the rows of df$uwwtd2 if empty:
df$uwwtd2 <- ifelse(is.na(df$uwwtd2), df$uwwtd, df$uwwtd2)

# prepare the dataset for easy merge
df <- df %>%
  rename(uwwtd_code=uwwtd2) 

# merge with median absolute error data
df_mae <- readRDS(here(read_data_here, "mae_measures.RDS"))

# fix df_mae
# remove _KAN-1 and _KAN-2
#df_mae$uwwtd <- gsub("_KAN-1", "", df_mae$uwwtd)
#df_mae$uwwtd<- gsub("_KAN-2", "", df_mae$uwwtd)

# fix wrong keys
df_mae$uwwtd_code[df_mae$uwwtd_code=="DEPT_SN137"] <- "DETP_SN137"
df_mae$uwwtd_code[df_mae$uwwtd_code=="DEPT_SN3016"] <- "DETP_SN3016"
df_mae <- df_mae %>% 
  mutate(uwwtd_code = str_replace_all(uwwtd_code, " ",""))

# consider mismerges
temp <- df %>% 
  anti_join(df_mae, by="uwwtd_code")
temp
# these are sites without data

temp <- df_mae %>% 
  anti_join(df, by="uwwtd_code")
temp
# Augsburg is just missing in df 
# add this manually
df <- df %>% 
  filter(uwwtd_code == "DETP_BYDON-K0001_KAN-1") %>% 
  mutate(uwwtd_code = as.character(ifelse(uwwtd_code == "DETP_BYDON-K0001_KAN-1", "DETP_BYDON-K0001", uwwtd_code))) %>% 
  bind_rows(df)

# merge
df <- df %>% 
  dplyr::select(-uwwtd, -Standort, -einwohner, -pl_nuts, -pl_koord_act, -amelag) %>% 
  distinct(uwwtd_code, .keep_all = T) %>% 
  right_join(df_mae %>% 
               dplyr::select(uwwtd_code,  mae, n, Labor, mean_blod), by="uwwtd_code") %>% 
  as_tibble()

# merge riip share data
df_riip_scores <-
  import(here(read_data_here, "total_scores15.csv"), encoding = "UTF-8") %>% 
  filter(!is.na(riip_share))

# remove spaces in original dataframe
df_riip_scores$uwwtd_code <- gsub(" ", "", df_riip_scores$uwwtd_code)
# change codes manually
df_riip_scores$uwwtd_code[df_riip_scores$uwwtd_code=="DEPT_SN137"] <- "DETP_SN137"
df_riip_scores$uwwtd_code[df_riip_scores$uwwtd_code=="DEPT_SN3016"] <- "DETP_SN3016"

df %>%
  anti_join(df_riip_scores)

df_riip_scores %>%
  anti_join(df)

df <- df %>%
  left_join(df_riip_scores) 

# # merge locations data
loc_data <- import(here(read_data_here, "mapdata.csv")) %>%
  dplyr::select(Standort = Kläranlage, long, lat) %>%
  filter(!is.na(long)) %>%
  as_tibble() %>% 
  add_row(Standort = "Münster", long = 7.6261347, lat = 51.9606649)

df_names <- readRDS(here(read_data_here, "current_amelag_data.rds")) %>%
  rename(uwwtd_code = `uwwtd`)  %>%
  distinct(uwwtd_code, .keep_all = T) %>%
  dplyr::select(uwwtd_code, Standort) %>%
  as_tibble()

loc_data <- loc_data %>%
  left_join(df_names)
# change codes manually
loc_data$uwwtd_code[loc_data$uwwtd_code=="DEPT_SN137"] <- "DETP_SN137"
loc_data$uwwtd_code[loc_data$uwwtd_code=="DEPT_SN3016"] <- "DETP_SN3016"

# merge df with loc_data; Attention: Uses Standorte, that means there s no distinction between: Hamburg Nord and Süd for example
df <- df %>%
  left_join(loc_data, by="uwwtd_code")

# merge not perfect, as Hagen and Bramsche are missing - why???
df %>%
  filter(is.na(long)) %>% 
  pull(Anlage)
# correct manually
df <- df %>% 
  mutate(long = case_when(Anlage == "KA Hagen" ~ 7.4716800,
                          Anlage == "4590142048 KA Bramsche" ~ 8.0015624,
                          .default = long),
         lat = case_when(Anlage == "KA Hagen" ~ 51.3608100,
                         Anlage == "4590142048 KA Bramsche" ~ 52.4077183,
                         .default = lat),
         Standort = case_when(Anlage == "KA Hagen" ~ "Hagen",
                              Anlage == "4590142048 KA Bramsche" ~ "Bramsche",
                              .default = Standort))

# drop sites without values for relevant variables
df <- df %>%
  dplyr::select(Standort, Bundesland, uwwtd_code, Labor, Nominalbelastung,
                long, lat, n, mae, riip_share, blod_share = mean_blod ) %>% 
  filter_at(all_of(rel_vars), all_vars(!is.na(.))) %>%
  filter(!is.na(long)) 

# compute scores ----------------------------------------------------------
# revert variables
df <- df %>%
  mutate(across(all_of(revert_vars),
                ~ . * -1))

# assign categories
df <- df %>% mutate(across(
  all_of(rel_vars),
  .fns = list(category = ~ case_when(
    . < quantile(., 0.25) ~ 3,
    between(., quantile(., 0.25), quantile(., 0.75)) ~ 2,
    . > quantile(., 0.75) ~ 1
  )),
  .names = "{col}_{fn}"
))

# assign ranks
df <- df %>% mutate(across(
  all_of(rel_vars),
  .fns = list(rank = ~ rank(-.)),
  .names = "{col}_{fn}"
))

# compute (weighted) scores
df <- df %>%
  rowwise() %>%
  mutate(
    rank_mean = weighted.mean(c_across(all_of(
      paste0(rel_vars, "_rank")
    )), w = weights),
    category_mean = weighted.mean(c_across(all_of(
      paste0(rel_vars, "_category")
    )), w = weights)
  ) %>%
  ungroup()

# select top sites --------------------------------------------------------
# to select; attention: there is a chance that states are
# missing completely (Saarland is missing)
score_var <- paste0(score, "_mean")

#to ensure that every Bundesland has at least one site:
top_per_state <- df |>
  group_by(Bundesland) |>
  slice_max(!!sym(score_var)) |>
  ungroup()

#antijoin to remove the top selected from each Bundesland
remaining <-  df |>
  anti_join(top_per_state, by = "uwwtd_code")

additional <- remaining |>
  slice_max(!!sym(score_var), n = no_selected - nrow(top_per_state))

# this selection includes at least 1 site per state (to ensure that Saarland is included), and then selects the remaining 64 largest sites (largest by Score)
sites_selected <- bind_rows(top_per_state, additional)

# show results ------------------------------------------------------------
#map data
germany_bl <- readRDS(here(read_data_here, "germany_map_bl.RDS"))

p = ggplot() +
  geom_polygon(
    data = germany_bl,
    aes(x = long, y = lat, group = group),
    colour = "black",
    fill = "white"
  ) + #schwarz/weiß
  # colour = "#4BACDD", fill= "#1273BA")+ #blau
  coord_map() +
  theme_void()

#add Points
pp = p +
  geom_point(data = sites_selected,
             aes(
               x = long,
               y = lat,
               color = "steelblue"
             ),
             size = 3) +
  scale_colour_manual(values = c("steelblue")) +
  theme(
    plot.margin = margin(0, 0, 0, 0),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  ggtitle(
    paste(no_selected, "selected sites"),
    subtitle = paste("Score calculated according to:", score)
  )

pp

# assign name to file
filename_map <-
  paste0("map_", no_selected, "_sites_scoring_by_", score, ".pdf")

# save map
ggsave(
  here(results_here, filename_map),
  width = 30,
  height = 30,
  units = "cm"
)

# create results for bundesländer
results_bula <- sites_selected %>%
  group_by(Bundesland) %>%
  summarise(
    n = n(),
    Mean_score = mean(!!sym(score_var)),
    Mean_Nominalbelastung = -mean(Nominalbelastung),
    Summe_Nominalbelastung = -sum(Nominalbelastung)
  )

# merge with population data for bundesländer
bev_bula <- import(here(read_data_here, "df_nuts.xlsx")) %>%
  distinct(Bundesland, Bev_BL = Bevölkerung_BL)

results_bula <- results_bula %>%
  left_join(bev_bula) %>%
  mutate(Anteil_Bev = Summe_Nominalbelastung / Bev_BL)

# assign name to file
filename_bula <-
  paste0("bula_auswertung_",
         no_selected,
         "_sites_scoring_by_",
         score,
         ".txt")

results_bula %>%
  mutate_if(is.numeric, ~ round(., 2)) %>%
  stargazer(
    summary = FALSE,
    type = "text",
    out = here(results_here, filename_bula)
  )

# save selected sites
# assign name to file
filename_sites <-
  paste0("selected_",
         no_selected,
         "_sites_scoring_by_",
         score,
         ".txt")

sites_selected %>% 
  dplyr::select(Bundesland, Standort, Labor, n, Nominalbelastung, riip_share, blod_share,
                mae, score_var) %>% 
  mutate_if(is.numeric, ~ abs(round(., 4))) %>% 
  arrange(Bundesland, desc(!!sym(paste0(score,"_mean")))) %>% 
  stargazer(
    summary = FALSE,
    type = "text",
    out = here(results_here, filename_sites)
  )

sites_selected %>% 
  dplyr::select(Bundesland, Standort, Labor, Nominalbelastung,n, 
                contains(score)) %>% 
  mutate_if(is.numeric, ~ abs(round(., 4))) %>% 
  arrange(Bundesland, desc(!!sym(paste0(score,"_mean")))) %>% 
  writexl::write_xlsx(path = here(results_here,  paste0("selected_",
                                                        no_selected,
                                                        "_sites_scoring_by_",
                                                        score,
                                                        ".xlsx")))
