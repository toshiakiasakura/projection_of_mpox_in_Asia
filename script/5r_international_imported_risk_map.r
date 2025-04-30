# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .r
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.7
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# ## Visualise the importation probabilities and international travel volumes
# This code will cause an error because of the lack of international flight volume data.

libraries = c("dplyr","magrittr","tidyr","ggplot2","RColorBrewer","zoo","lubridate","tidyverse",
              "ggpattern", "assertthat", "purrr", "igraph", "ggmap", "readxl", "geojsonio", "sf"
             )
for(x in libraries) {library(x, character.only=TRUE, warn.conflicts=FALSE, quietly=TRUE)}

R.Version()$version.string
source("vis_utils.R")

# For sf package.
Sys.setenv("PROJ_LIB"="/opt/conda/share/proj")
Sys.getenv("PROJ_LIB")

shp_df <- read_who_shape_file()

shp <- st_read(file.path("../data/Detailed_Boundary_ADM0_565521753006392799/GLOBAL_ADM0.shp"))
shp <- st_transform(shp, crs=4326)
shp_data_fil <- data.frame(shp) %>% select(WHO_CODE, ISO_3_CODE, CENTER_LON, CENTER_LAT)
df_geo_info <- read_geo_info()

path_flight <- "../data/flight/selected_flight_matrix.csv"
edges_for_plot <- read_clean_flight(path_flight, df_geo_info, thres=10000)

edges_for_plot["log10_values"] %>% summary

# ## Clean importation probabilities

library(ggnewscale)
options(repr.plot.width=15,repr.plot.height=10)

# Function preparations
ret = prepare_scale_bar()
scale_bar_rects = ret$a
scale_bar_labels = ret$b

path <- "../tmp_results/imp_prob_path_sc1_fil.csv"
df_imp <- clean_data(path, df_geo_info)
df_imp_na <- df_imp %>% filter(color=="e_na")

# +
path <- "../tmp_results/imp_prob_path_sc1_fil.csv"
df_imp <- clean_data(path, df_geo_info)
df_imp %<>% mutate(imp_prob_percent = case_when(
    #imp_prob_percent == 100 ~ 70,
    TRUE ~ imp_prob_percent
))

path_map <- "../fig/international_imp_map_sc1_fil.png"
visualise_international_map(path_map, df_imp, edges_for_plot,
    label_set = c("10%", "20%", "30%", "40%"),
    break_set = c(10, 20, 30, 40),
)
# -

df_tmp <- df_imp %>% arrange(desc(imp_prob))
write_csv(df_tmp, "../fig/imp_prob_natsal_4w_IP10.csv")
df_tmp %>% head(10)


