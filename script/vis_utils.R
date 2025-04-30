read_who_shape_file <- function() {
    # Read shape data from the WHO map
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
    shp_df
    #shp_data <- data.frame(shp)
    #shp_data
}

read_geo_info <- function() {
    df_geo_info <- read.csv("../data/geo_info.csv") %>% filter(!is.na(lat))
    df_iso <- read.csv("../data/iso2_iso3.csv") %>% dplyr::select(alpha.2, alpha.3) %>%
        rename(iso=alpha.2, iso_3=alpha.3)
    df_geo_info <- merge(df_geo_info, df_iso, by=c("iso"), all=TRUE) %>% filter(!is.na(country))
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
    df_geo_info
}

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

convert_lon_distance_to_lon <- function(start_lon, start_lat, distance_km) {
  # Approximate kilometers per degree of longitude at the given latitude
  km_per_lon_degree <- 111.320 * cos(abs(start_lat) * pi / 180)

  # Calculate the change in longitude (in degrees)
  delta_lon <- distance_km / km_per_lon_degree

  # Calculate the final longitude
  final_lon <- start_lon + delta_lon

  return(final_lon)
}

prepare_scale_bar <- function(){
    # --- 1. Scale Bar Parameters ---
    total_distance_km <- 3500
    first_half_distance_km <- 1750
    num_segments_first_half <- 4
    scale_bar_location <- data.frame(
      lon_start = 27,
      lon_end = convert_lon_distance_to_lon(35, -18, total_distance_km), # Adjust for desired total length
      lat_bottom = -16,
      lat_top = -15
    )

    # --- 2. Calculate Plot Unit Dimensions (Approximate for Geographic CRS) ---
    ref_lat <- mean(scale_bar_location$lat_bottom, scale_bar_location$lat_top)
    km_per_lon_degree <- 111.320 * cos(ref_lat * pi / 180)
    total_lon_range <- scale_bar_location$lon_end - scale_bar_location$lon_start
    total_plot_distance_km <- total_lon_range * km_per_lon_degree

    first_half_plot_fraction <- first_half_distance_km / total_distance_km
    first_half_lon_range <- total_lon_range * first_half_plot_fraction
    segment_width_lon <- first_half_lon_range / num_segments_first_half

    # --- 3. Create Scale Bar Rectangle Data ---
    scale_bar_rects <- data.frame(
      xmin = scale_bar_location$lon_start + c(0, 1, 2, 3) * segment_width_lon,
      xmax = scale_bar_location$lon_start + c(1, 2, 3, 4) * segment_width_lon,
      ymin = scale_bar_location$lat_bottom,
      ymax = scale_bar_location$lat_top,
      scale_fill = rep(c("black", "white"), 2)
    )

    # Add the black rectangle for the second half
    scale_bar_rects <- rbind(scale_bar_rects, data.frame(
      xmin = scale_bar_location$lon_start + first_half_lon_range,
      xmax = scale_bar_location$lon_end,
      ymin = scale_bar_location$lat_bottom,
      ymax = scale_bar_location$lat_top,
      scale_fill = "black"
    ))

    # --- 4. Create Scale Bar Label Data ---
    scale_bar_labels <- data.frame(
      lon = scale_bar_location$lon_start +
            c(0, 0.5 * first_half_plot_fraction + 0.02, first_half_plot_fraction + 0.1, 1.08) * total_lon_range,
      lat = scale_bar_location$lat_bottom - 0.6,
      label_text = c("0 km", "875 km", "1750 km", "3500 km"),
      hjust = c(0, 0.5, 1, 1)
    )
    list(a=scale_bar_rects, b=scale_bar_labels)
}

clean_data <- function(path, df_geo_info, imp_prob_percent=4.0){
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
    )) %>% # Remove redundant columns
        select(iso_code, country, lon, lat, imp_prob, imp_prob_percent, color)
    # Add data for non-applicable countries for visualisation purpose.
    df_UN_Asia <- read.csv("../data/UN_Asia_list.csv")
    df_imp_na <- df_geo_info %>% filter(iso_3 %in% df_UN_Asia[["Code"]]) %>%
        filter(! iso_3 %in% df_imp[["iso_code"]]) %>%
        filter(iso_3 != "JPN") %>%
        select(iso_3, lat, lon, country) %>%
        rename(iso_code = iso_3) %>%
        mutate(imp_prob_percent = imp_prob_percent,  color="e_na")
    df_imp <- bind_rows(df_imp, df_imp_na)
    df_imp
}

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

        # Plot main points
        geom_curve(aes(x = sr_x, y = sr_y, xend = ta_x, yend = ta_y,
                       linewidth = log10_values), #, colour=Travel),
                   data = edges_for_plot, curvature = 0.33, alpha = 0.5,
                   color="#00BFC4",
                  ) +
        scale_linewidth_continuous("Yearly travel volume",
                                   labels=label_set2,
                                   breaks=break_set2,
                                   range = c(0.1, 2.5)) +
        coord_sf(xlim = c(30, 140),  ylim = c(-20, 60) , crs = 4326) +
        geom_text(aes(x=c(137.9739), y=c(37.53983)),
                  label="★", size=10, family = "HiraKakuPro-W3", color="red") +
        geom_point(data=df_imp ,
                   aes(x=lon, y=lat, size = imp_prob_percent, fill=color),
                   shape=21,
                   color="black", stroke=0.5) +
        scale_fill_manual("Observed importation pattern",
                          values = c("#FF6633", "#FFCC66", "white", "grey"),
                          labels=c("≥1 Importation from Asia", "≥1 Importation", "No importation", "Not included"),
                          ) +
        guides(fill=guide_legend(override.aes = list(size=8))) +
        scale_size_continuous("Simulated importation prob.",
                              labels=label_set, breaks=break_set, range = c(1, 20)) +
        # Plot non-applicable countries
        geom_point(data=df_imp_na, aes(x=lon, y=lat),
                   fill="grey60", stroke=0.3, shape=22, size=3) +
        # Add the scale bar rectangles
        geom_rect(data = scale_bar_rects,
                  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  color = "black", linewidth = 0.2,
                  fill = scale_bar_rects$scale_fill) +
        # Add the scale bar labels
        geom_text(data = scale_bar_labels, aes(x = lon, y = lat, label = label_text, hjust = hjust), size = 3, color="black") +

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