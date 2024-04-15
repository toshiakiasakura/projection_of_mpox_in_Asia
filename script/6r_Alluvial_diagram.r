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

# # Alluvial diagram

libraries = c("dplyr","magrittr","tidyr","ggplot2","RColorBrewer","zoo","lubridate","tidyverse",
              "ggpattern", "assertthat", "purrr", "igraph", "ggmap",
              "htmltools", "htmlwidgets", "RColorBrewer", "viridis", "cmocean",
              "ggalluvial"
             )
for(x in libraries) {library(x, character.only=TRUE, warn.conflicts=FALSE, quietly=TRUE)}

R.Version()$version.string

# +
#packageVersion("networkD3")
# -

path_pop <- "../data/pop_size_edit.csv"
df_pop <- read.csv(path_pop) %>% select(location, iso_code)

# ## Check UN region

dat_UN <- read.csv("../data/UN_Asia_list.csv")
uni_reg <- dat_UN[["Region"]] %>% unique
for (u in uni_reg){
    print("=---")
    dat_UN %>%
        filter(Region == u) %>%
        print
}

df_pop %>% filter(! (iso_code %in% dat_UN[["Code"]]))

dat_UN %>% filter(!(Code %in% df_pop[["iso_code"]]))

# # ggalluvial data cleaning part.

library(ggalluvial)

# +
path = "../tmp_results/exp_imp_gen_pair_path5_IP10_fil.csv"

links <- read.csv(path)
links %<>%
    mutate(cohort= paste0(sim_index,"_",seq_index)) %>%
    filter(import!="Japan") # Exclude Japan from the results
links["freq"] <- 1
# -

links["gen_index"] %>% unique

dat_tmp <- links %>% filter(gen_index == 1)

n_1st <- links %>% filter(gen_index==1) %>% dim
paste0("# of 1st generations, since Japan export multiple times, this number can exceeds 1000: ", n_1st[1])

#- dat_single_gen: One generation data, filtered by sequence id.
#- dat_one: Contains >n_gen + 1 data. Finally appended for the alluvial plot.
one_gen_step <- function(dat_single_gen, dat_one){
    # One generation step
    n_gen <- dat_single_gen[1, "gen_index"]
    index <- 1
    dat_gen_post <- dat_one %>% filter(gen_index==n_gen+1)
    for(i in 1:nrow(dat_single_gen)){
        dat_genM <- dat_gen_post %>% filter(export==dat_single_gen[i, "import"])
        # If the n_gen imported country does not export to any countries (dat_gen_post's export).
        n_row_genM <- nrow(dat_genM)
        if (n_row_genM == 0){
            dat_single_gen[i, "cohort"] <- paste0(dat_single_gen[i, "cohort"], "_", n_gen, "_", index)
            dat_single_gen[i, "freq"] <- 1
            index <- index + 1
            dat_one %<>% bind_rows(dat_single_gen[i, extract_cols])
        } else {
            n_row_genM <- nrow(dat_genM)
            n_uni <- dat_genM[["import"]] %>% unique %>% length
            for (j in 1:n_row_genM){
                dat_single_gen[i, "freq"] <-  dat_genM[j, "freq"]/n_uni
                dat_single_gen[i, "cohort"] <- dat_genM[j, "cohort"]
                dat_one %<>% bind_rows(dat_single_gen[i, extract_cols])
            }
        }
    }
    dat_one
}

links %>% filter(gen_index == 1) %>% nrow

links %>% filter(cohort == "127_3")

# +
# Modifying links dataframe to be compatible with Alluvial plots.
dat_allu <- data.frame()

extract_cols <- c("export", "import", "gen_index", "cohort", "freq")
sim_ids <- links[["sim_index"]] %>% unique
for (l in 1:length(sim_ids)){
    sim <- sim_ids[l]
    dat_sim <- links %>% filter(sim_index == sim)
    seq_ids <- dat_sim[["cohort"]] %>% unique
    n_seq <- length(seq_ids)
    for (k in 1:length(seq_ids)){
        # One seq index
        seq <- seq_ids[k]
        dat <- dat_sim %>% filter(cohort==seq)
        n_gen_max <- dat["gen_index"] %>% max
        dat_one <- dat %>% filter(gen_index == Inf) # Create empty data.frame with columns.
        for( n_gen in seq(n_gen_max, 1, -1)){
            datM <- dat %>% filter(gen_index==n_gen)
            datM["freq"] <- 1 # Create a column, and updated in one_gen_step
            dat_one <- one_gen_step(datM, dat_one)
        }
        # Add Japan's exportation part.
        datM <- dat %>% filter(gen_index==1)
        if (nrow(datM)  > 0){
            datM <- data.frame(
                export=c("Japan"),
                import=c("Japan"),
                gen_index=c(0),
                cohort=c("origin"),
                freq=c(1)
            )
            # Since n_seq corresponds to one exportation from Japan,
            # if we sum over a single seq_id, we obtain 1.0.
            dat_one <- one_gen_step(datM, dat_one)
            # In this part, we have to scale the exportation freq from Japan to be one
            # for each simulation
            dat_one %<>% mutate(freq = case_when(
                gen_index==0 ~ freq/n_seq,
                TRUE ~ freq
            ))
        }
        dat_allu %<>% bind_rows(dat_one)
    }
}
# -

# Re-assign the freq to adjsut the height of the righthand side of links
# to be proportional to actual exportation events.
gen_ind <- 1
dat_allu_new <- data.frame()
for (gen_ind in 1:4){
    gen_export <- dat_allu %>% filter(gen_index == gen_ind)
    gen_import <- dat_allu %>% filter(gen_index == gen_ind + 1)
    # Calculate the proportion of each imported countries stratified by exportation country.
    gen_prop <- gen_import %>%
        group_by(export, import) %>%
        summarise(prob=sum(freq), .groups="keep") %>%
        group_by(export) %>%
        mutate(prop = prob/sum(prob))
    # Calcualte the frequency for the exportation among the export countries.
    gen_import_tmp <- gen_import %>% select(cohort, import) %>% rename(import_gen=import)
    gen_export_tmp <- merge(gen_export, gen_import_tmp, by="cohort", all.x=TRUE) %>%
        mutate(na_freq = case_when(
            !is.na(import_gen) ~ freq,
            is.na(import_gen) ~ 0,
        ),
          na_flag = case_when(
            !is.na(import_gen) ~ 1,
             is.na(import_gen) ~ 0
          ))
    # Exportation probabilities for the export countries.
    gen_export_prob <- gen_export_tmp %>%
        group_by(import) %>%
        summarize(prob_export=sum(na_freq))
    # Number of rows for each exportations.
    gen_export_rows <- gen_export_tmp %>%
        filter(na_flag == 1) %>%
        group_by(import, import_gen) %>%
        summarise(n_row=n(), .groups="keep")
    # Create and assign the individual values.
    gen_master <- merge(gen_prop, gen_export_prob, by.x="export", by.y="import", all.x=TRUE) %>%
        mutate(freq_exp = prob_export*prop)
    gen_master <- merge(gen_master, gen_export_rows,
                        by.x=c("export", "import"),
                        by.y=c("import", "import_gen"),
                        all.x=TRUE
                       ) %>%
        mutate(freq_row = freq_exp/n_row) %>% select(export, import, freq_row)
    dat_rescale <- merge(gen_export_tmp, gen_master,
                         by.x=c("import", "import_gen"),
                         by.y=c("export", "import"),
                         all.x=TRUE,
                        )
    dat_rescale %<>% mutate(freq = case_when(
        na_flag==1 ~ freq_row,
        TRUE ~ freq,
    ))
    dat_allu_new <- rbind(dat_allu_new, dat_rescale)
}
dat_allu_new %<>% select(!c("na_freq", "na_flag", "freq_row", "import_gen"))
# Add Japan and last parts.
dat_allu_jpn <- dat_allu %>% filter(export=="Japan", gen_index==0)
dat_allu_last <- dat_allu %>% filter(gen_index==5)
dat_allu_new <- rbind(dat_allu_new, dat_allu_jpn, dat_allu_last)

# Check the data frame
dat_allu %>% group_by(gen_index) %>% summarise(n=n()) %>% print
dat_allu_new %>% group_by(gen_index) %>% summarise(n=n()) %>% print
dat_allu %>% group_by(gen_index) %>% summarise(n=sum(freq)) %>% print
dat_allu_new %>% group_by(gen_index) %>% summarise(n=sum(freq)) %>% print

dat_allu %>% dim %>% print
dat_allu_new %>% dim %>% print

dat_tmp <- dat_allu %>% filter(gen_index == 1) %>% filter(import == "China")
dat_tmp["freq"] %>% sum %>% print
dat_tmp <- dat_allu_new %>% filter(gen_index == 1) %>% filter(import == "China")
dat_tmp["freq"] %>% sum %>% print

#write.csv(dat_allu_new, file="../tmp_results/dat_allu.csv", row.names=FALSE)
write.csv(dat_allu_new, file="../tmp_results/dat_allu_fil.csv", row.names=FALSE)

# # Alluvial plots

#n_sim <- 5000
#dat_allu <- read.csv("../tmp_results/dat_allu.csv")
n_sim <- 626
dat_allu <- read.csv("../tmp_results/dat_allu_fil.csv")


dat_allu["freq"] <- dat_allu["freq"] /n_sim *100
print("Check the order (maximum is 5?)")
dat_allu[["gen_index"]] %>% unique

is_lodes_form(dat_allu, key=gen_index, id=cohort, value=import)

clean_country <- function(dat_allu){
    # Rename
    dat_allu[dat_allu == "China"] <- "China (mainland)"
    dat_allu[dat_allu == "Hong Kong, China"] <- "Hong Kong"
    dat_allu[dat_allu == "Taiwan, Province of China"] <- "Taiwan"
    dat_allu[dat_allu == "Korea, Republic of"] <- "Republic of Korea"
    dat_allu[dat_allu == "Viet Nam"] <- "Vietnam"
    dat_allu[dat_allu == "Lao People's Democratic Republic"] <- "Laos"
    dat_allu[dat_allu == "Iran, Islamic Republic of"] <- "Iran"
    dat_allu[dat_allu == "Syrian Arab Republic"] <- "Syria"
    dat_allu[dat_allu == "Brunei Darussalam" ] <- "Brunei"
    dat_allu[dat_allu == "Turkey" ] <- "Türkiye"
    dat_allu
}
dat_allu %<>% clean_country

options(repr.plot.width=7,repr.plot.height=5)

path_map <- "../fig/alluvial_fig.png"
ggplot(dat_allu,
       aes(x=gen_index, alluvium=cohort, y=freq, stratum=import)
      ) +
    geom_flow(aes(fill=import), decreasing=TRUE) +
    geom_stratum(aes(fill=import), decreasing=TRUE) +
    stat_stratum(geom = "text", aes(label = import), decreasing = TRUE, min.y=100) +
    #geom_stratum(stat = "alluvium", aes(fill = import), decreasing=TRUE)
    theme(legend.position="none")
#ggsave(path_map, width=18.4, height= 10.0, dpi=150)

library(cmocean)

# +
# Create colour pallate for each region
countries <- c(dat_allu[["export"]], dat_allu[["import"]]) %>% unique %>% sort
dat_cnt <- data.frame(country=countries)
dat_mer <- merge(dat_cnt, dat_UN, by.x = "country", by.y="Name", all.x=TRUE)
#dat_cnt[!(dat_cnt[["country"]] %in% dat_mer[["country"]]), ]
dat_mer %<>% mutate(Region = case_when(
    country == "Türkiye" ~ "Western Asia",
    country == "China (mainland)" ~ "Eastern Asia",
    country == "Hong Kong" ~ "Eastern Asia",
    country == "Republic of Korea" ~ "Eastern Asia",
    country == "Taiwan" ~ "Eastern Asia",
    TRUE ~ Region,
),
                   Code=case_when(
    country == "Türkiye" ~ "TUR",
    country == "China (mainland)" ~ "CHN",
    country == "Hong Kong" ~ "HKG",
    country == "Republic of Korea" ~ "KOR",
    country == "Taiwan" ~ "TWN",
    TRUE ~ Code,
))
#dat_mer <- dat_mer %>% replace(is.na(.), "Eastern Asia")
dat_mer %<>% mutate(region3=case_when(
    (Region == "Eastern Asia") ~ "East",
    (Region == "South-eastern Asia") ~ "South-east" ,
    (Region == "Central Asia") | (Region == "Southern Asia") | (Region == "Western Asia") ~ "Central-west"
))
dat_mer %>% group_by(Region) %>% summarise(n=n()) %>% print
dat_mer %>% group_by(region3) %>% summarise(n=n()) %>% print
# Merge with importation probabilities
path <- "../tmp_results/imp_prob_path5_IP10_fil.csv"
dat_imp <- read.csv(path) %>% select(iso_code, imp_prob)
dat_mer <- merge(dat_mer, dat_imp, by.x="Code", by.y="iso_code", all.x=TRUE)
dat_mer %<>% arrange(region3, desc(imp_prob))

east_col <- cmocean("algae", clip=0.2, direction=1)(6)
south_east_col <- cmocean("ice", clip=0.35, direction=-1)(11)
cent_west_col <- cmocean("amp", clip=0.2, direction=1)(25)
dat_mer["colour"] <- c(cent_west_col, east_col, south_east_col)

dat_mer %<>% arrange(country)
# -

dat_mer %<>% mutate(colour = case_when(
    country == "Japan" ~ east_col[3],
    country == "Taiwan" ~ east_col[4],
    country == "Hong Kong" ~ east_col[5],
    TRUE ~ colour
))

options(repr.plot.width=15,repr.plot.height=10)

# +
path_map <- "../fig/alluvial_fig_side_text.png"
ymin <- 0.04*100
size <- 4
ggplot(dat_allu,
       aes(x=gen_index, alluvium=cohort, y=freq, stratum=import)
      ) +
    geom_flow(aes(fill=import), alpha=0.8, decreasing=TRUE) +
    geom_stratum(aes(fill=import), width=.20, decreasing=TRUE) +
    #stat_stratum(geom = "text", aes(label = import), decreasing = TRUE, min.y=100) +
    ggrepel::geom_text_repel(
        aes(label=ifelse(as.numeric(gen_index) == 0, as.character(import), NA)),
        stat="stratum", decreasing=TRUE, size=size, nudge_x=-0.3, min.y=ymin
    ) +
    ggrepel::geom_text_repel(
        aes(label=ifelse(as.numeric(gen_index) == 1, as.character(import), NA)),
        stat="stratum", decreasing=TRUE, size=size, nudge_x=-0.5, min.y=ymin
    ) +
    ggrepel::geom_text_repel(
        aes(label=ifelse(as.numeric(gen_index) == 2, as.character(import), NA)),
        stat="stratum", decreasing=TRUE, size=size, nudge_x=+0.5, min.y=ymin
    ) +
    ggrepel::geom_text_repel(
        aes(label=ifelse(as.numeric(gen_index) == 3, as.character(import), NA)),
        stat="stratum", decreasing=TRUE, size=size, nudge_x=+0.5, min.y=ymin
    ) +
    ggrepel::geom_text_repel(
        aes(label=ifelse(as.numeric(gen_index) == 4, as.character(import), NA)),
        stat="stratum", decreasing=TRUE, size=size, nudge_x=+0.5, min.y=ymin
    ) +
    ggrepel::geom_text_repel(
        aes(label=ifelse(as.numeric(gen_index) == 5, as.character(import), NA)),
        stat="stratum", decreasing=TRUE, size=size, nudge_x=+0.5, min.y=ymin
    ) +
    scale_fill_manual(values=dat_mer[["colour"]]) +
    #geom_stratum(stat = "alluvium", aes(fill = import), decreasing=TRUE)
    xlab("Generation") + ylab("Simulated probability of international spread events (%, stacked)") +
    theme(legend.position="none",
          #panel.background = element_rect(fill="white"),
          axis.text=element_text(size=16),
          axis.title=element_text(size=20),
         ) +
    scale_x_continuous(breaks=c(0,1,2,3,4,5)) +
    #scale_y_continuous(breaks=c(0, 25, 50, 75, 100, 125, 150)) +
    scale_y_continuous(breaks=c(0, 100, 200, 300, 400)) +
    #scale_y_continuous(breaks=c(0, 10, 20, 30)) + ylim(c(0,30)) +

    annotate("text", label="Eastern Asia (6 countries)",       x=4.5, y=80.0*3, size=size+1) +
    annotate("text", label="South-eastern Asia (11 countries)", x=4.5, y=68.0*3, size=size+1) +
    annotate("text", label="Other Asia (25 countries)",         x=4.5, y=56.0*3, size=size+1)
ggsave(path_map, width=16.4, height= 12.0, dpi=300)

# Colormap was added manually. Heigth 140px for Julia's colormap display, set guide to x=3693px.
# For hatched line, select two areas and merge layer, and use gradient.

# -








