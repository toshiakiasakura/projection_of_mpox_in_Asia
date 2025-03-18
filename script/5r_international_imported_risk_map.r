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
              "ggpattern", "assertthat", "purrr", "igraph", "ggmap", "readxl"
             )
for(x in libraries) {library(x, character.only=TRUE, warn.conflicts=FALSE, quietly=TRUE)}

R.Version()$version.string

# For sf package.
Sys.setenv("PROJ_LIB"="/opt/conda/share/proj")
Sys.getenv("PROJ_LIB")

# Read shape data from the WHO map
library(geojsonio)
library(sf)
shp <- st_read(file.path("../data/Detailed_Boundary_ADM0_565521753006392799/GLOBAL_ADM0.shp"))
shp <- st_transform(shp, crs=4326)
shp_df = sfheaders::sf_to_df(shp)
shp_df %<>% mutate(group_id =
                   paste0(as.character(sfg_id), "_",
                          as.character(multipolygon_id), "_",
                          as.character(polygon_id), "_",
                          as.character(linestring_id)
                         ),
                  )
shp_data <- data.frame(shp)
shp_data_fil <- shp_data %>% select(WHO_CODE, ISO_3_CODE, CENTER_LON, CENTER_LAT)

#### geographical location of each country
read.csv("../data/geo_info.csv") %>% filter(!is.na(lat)) -> df_geo_info
read.csv("../data/iso2_iso3.csv") %>% dplyr::select(alpha.2, alpha.3) %>%
rename(iso=alpha.2, iso_3=alpha.3) -> df_iso
merge(df_geo_info, df_iso, by=c("iso"), all=TRUE) %>% filter(!is.na(country)) -> df_geo_info
df_geo_info %<>%
    mutate(country=case_when(
        country == 'Brunei' ~ "Brunei Darussalam",
        country == "Hong Kong" ~ 'Hong Kong, China',
        country == "Iran" ~ 'Iran, Islamic Republic of',
        country == 'Laos' ~ 'Lao People\'s Democratic Republic',
        country == 'Syria' ~ 'Syrian Arab Republic',
        country == 'South Korea' ~ 'Korea, Republic of',
        country == 'Myanmar [Burma]' ~ 'Myanmar',
        country == 'Vietnam' ~ 'Viet Nam',
        country == 'Taiwan' ~ "Taiwan, Province of China",
        TRUE ~ country,
    ))
df_geo_info <- merge(df_geo_info, shp_data_fil, by.x="iso_3", by.y="ISO_3_CODE", all.x=TRUE)
df_geo_info %<>% mutate(lat = case_when( is.na(CENTER_LAT) ~ lat,
                                        # iso_3 == "MYS" ~ lat,
                                        TRUE ~ CENTER_LAT),
                        lon= case_when( is.na(CENTER_LON) ~ lon,
                                       #iso_3 == "MYS" ~ lon,
                                       TRUE ~ CENTER_LON)
                        )

#### flight volume for international
read_clean_flight <- function(path_flight, df_geo_info, thres=10000){
    flight_matrix <- read.csv(path_flight)
        flight_long <- cbind(flight_matrix[1:2], stack(flight_matrix[3:length(flight_matrix)])) %>%
            rename(target=country, source=ind)
    flight_long$source <- gsub(".", " ", flight_long$source, fixed=TRUE)
    flight_long$source <- gsub("  ", ", ", flight_long$source, fixed=TRUE)
    flight_long$source <- gsub(" s ", "\'s ", flight_long$source, fixed=TRUE)
    flight_long[flight_long == 'Timor Leste'] <- 'Timor-Leste'
    flight_long %>%
        left_join(df_geo_info %>% select(country, lat, lon), by=c("target"="country")) %>%
        rename(ta_y=lat, ta_x=lon) %>%
        #rename(ta_y=CENTER_LAT, ta_x=CENTER_LON) %>%
        left_join(df_geo_info %>% select(country, lat, lon), by=c("source"="country")) %>%
        rename(sr_y=lat, sr_x=lon) %>%
        #rename(sr_y=CENTRE_LAT, sr_x=CENTER_LON) %>%
        mutate(Travel = case_when(
            #target == "Japan" ~ "Japan",
            TRUE ~ "Others",
        )) ->
        edges_for_plot
    edges_for_plot %<>% filter(values > thres)
    edges_for_plot["log10_values"] <- edges_for_plot["values"] %>% log10 - log10(thres)
    edges_for_plot
}

path_flight <- "../data/flight/selected_flight_matrix.csv"
edges_for_plot <- read_clean_flight(path_flight, df_geo_info, thres=10000)

edges_for_plot["log10_values"] %>% summary

# ## Clean importation probabilities

library(ggnewscale)
options(repr.plot.width=15,repr.plot.height=10)

clean_data <- function(path, df_geo_info){
    df_imp <- read.csv(path)
    df_imp["imp_prob_percent"] <- df_imp["imp_prob"]*100
    df_imp %<>% left_join(df_geo_info, by=c("iso_code"="iso_3")) %>% filter(country!="Japan")

    df_obs <- read_excel("../data/mpox_Asia_importation_date.xlsx", sheet="Sheet2")
    df_imp <- left_join(df_imp, df_obs, by ="iso_code")
    df_imp[["obs_imp_flag"]] %<>% replace_na(0)
    df_imp[["asia_imp_flag"]] %<>% replace_na(0)
    df_imp %<>% mutate(color = case_when(
        #iso_code == "KOR" ~ "a_Korea",
        asia_imp_flag == 1 ~ "b_asia_import",
        obs_imp_flag == 1 ~ "c_import",
        TRUE ~ "d_no_import"
    ))
    df_imp
}

# +
visualise_international_map <- function(
    path_map,
    df_imp, edges_for_plot,
    label_set = c("0.0%", "1.0%", "10.0%", "40.0%"),
    break_set = c(0, 1, 10, 40),
    label_set2 = c("10,000", "100,000", "1,000,000", "10,000,000"),
    break_set2 = c(0, 1, 2, 3) # log10 - 4
){
    ggplot() +
        geom_polygon(aes(x = x, y = y, group = group_id), data = shp_df,
                     fill = "#CECECE", color = "#515151", linewidth = 0.15) +

        geom_curve(aes(x = sr_x, y = sr_y, xend = ta_x, yend = ta_y,
                       linewidth = log10_values), #, colour=Travel),
                   data = edges_for_plot, curvature = 0.33, alpha = 0.5,
                   color="#00BFC4",
                  ) +
        #scale_color_manual("Export countries", values = c("#F8766D", "#00BFC4","royalblue3")) +
        scale_linewidth_continuous("Yearly travel volume",
                                   labels=label_set2,
                                   breaks=break_set2,
                                   range = c(0.1, 2.5)) +
        coord_quickmap(xlim = c(30, 140),  ylim = c(-20, 60)) +

        geom_text(aes(x=c(137.9739), y=c(37.53983)),
                  label="★", size=10, family = "HiraKakuPro-W3", color="red") +
        #geom_point(data=df_imp , aes(x=lon, y=lat, size = imp_prob_percent),
        #           shape=21, fill="white", color="black", stroke=0.5) +
        geom_point(data=df_imp ,
                   aes(x=lon, y=lat, size = imp_prob_percent, fill=color),
                   shape=21, color="black", stroke=0.5) +
        scale_fill_manual("Observed importation pattern",
                          values = c("#FF6633", "#FFCC66", "white"),
                           labels=c("≥1 Importation from Asia", "≥1 Importation", "No importation"),
                          ) +
        guides(fill=guide_legend(override.aes = list(size=8))) +
        scale_size_continuous("Simulated importation prob.",
                              labels=label_set, breaks=break_set, range = c(1, 20)) +
        theme(panel.grid = element_blank()) + theme(axis.text = element_blank()) +
        theme(axis.ticks = element_blank()) + theme(axis.title = element_blank()) +
        theme(legend.position = "right", legend.title = element_text(size = 15), legend.text = element_text(size = 15)) +
        theme(panel.grid = element_blank()) +
        theme(panel.background = element_rect(fill = "#596672")) +
        theme(plot.margin = margin(-2, -5, -2, -5, "cm")) +
        guides(color = guide_legend(override.aes = list(linewidth = 5, size=10), order=1),
               size = guide_legend(order = 6),
               linewidth = guide_legend(order = 2))
    ggsave(path_map, width=18.4, height= 10.0, dpi=300)

}

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



