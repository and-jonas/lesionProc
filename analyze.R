
# HEADER ----

rm(list = ls())
Sys.getenv("R_LIBS")
.libPaths("T:/R4Userlibs")
list.of.packages <- c("tidyverse", "data.table", "caret", 
                      "ggsci", "ggpmisc", "ggpubr", "GGally", "quantreg", "progress",
                      "scam")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE, repos = 'https://stat.ethz.ch/CRAN/')

library(data.table)
library(lubridate)
library(caret)
library(ggsci)
library(ggpmisc)
library(ggpubr)
library(GGally)
library(gridExtra)
library(readxl)
library(nls.multstart)
library(quantreg)
library(progress)
library(ggExtra)
library(furrr)
library(tictoc)
library(lme4)
library(nlme)
library(lmerTest)
library(emmeans)
library(tidyverse)

source("Z:/Public/Jonas/004_ESWW007/RScripts/utils.R")

# ============================================================================== -
# 0) Prepare ----
# ============================================================================== -

# # paths
# data_path <- "Z:/Public/Jonas/011_STB_leaf_tracking/output/datasets_test/"
# figure_path <- "Z:/Public/Jonas/011_STB_leaf_tracking/Figures/2/"
# meta_path <- "Z:/Public/Jonas/011_STB_leaf_tracking/output/ts2/"

# paths
in_path <- "Z:/Public/Jonas/011_STB_leaf_tracking/output/"
data_path <- "Z:/Public/Jonas/011_STB_leaf_tracking/data/"
figure_path <- "Z:/Public/Jonas/011_STB_leaf_tracking/Figures/4/"

# base plot
col = pal_jco()(10)[c(1, 8, 9)]
base_plot_batches <- ggplot() +
  scale_color_manual(values =  col, 
                     breaks = c("1", "2", "3"),
                     labels = c("Batch1", "Batch2", "Batch3")) +
  xlab("Time of measurement") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(panel.background = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(),
        legend.title = element_blank(),
        legend.position.inside = c(0.8, 0.1))

# ============================================================================== -
# 1) Get data ----
# ============================================================================== -

# get data sets
lesion_data_esww007 <- readRDS(paste0(in_path, "dataset_ESWW007/data_lesions_aug.rds"))
lesion_data_esww009 <- readRDS(paste0(in_path, "dataset_ESWW009/data_lesions_aug.rds"))
leaf_data_esww007 <- readRDS(paste0(in_path, "dataset_ESWW007/data_leaf_aug.rds"))
leaf_data_esww009 <- readRDS(paste0(in_path, "dataset_ESWW009/data_leaf_aug.rds"))
lesion_data <- bind_rows(lesion_data_esww007, lesion_data_esww009)
leaf_data <- bind_rows(leaf_data_esww007, leaf_data_esww009)

# get meta information
# get information on which leaves in which time span to analyze
meta_esww007 <- read.csv("Z:/Public/Jonas/Data/ESWW007/SingleLeaf/Meta/stb_dominant.csv")
meta_esww009 <- read.csv("Z:/Public/Jonas/Data/ESWW009/SingleLeaf/Meta/stb_dominant_ESWW009.csv")
meta <- bind_rows(meta_esww007, meta_esww009) %>% 
  mutate(date = strsplit(last_image_id, "_") %>% 
           lapply(., "[[", 1) %>% unlist(),
         time = strsplit(last_image_id, "_") %>% 
           lapply(., "[[", 2) %>% unlist(),
         last_valid_timestamp = ymd_hms(paste(date, time, sep = " "))) %>% 
  dplyr::select(leaf_UID, last_valid_timestamp)

# get leaf-level features that are of interest at lesion-level as well
names(leaf_data)
leaf_level_features <- leaf_data %>% 
  dplyr::select(plot_UID, timestamp, leaf_UID, la_healthy_f, placl, rust_density)

# ============================================================================== -
# 2) Filter on Leaves ----
# ============================================================================== -

# select relevant leaves within relevant time span
files <- leaf_data %>% 
  full_join(., meta) %>% 
  dplyr::filter(leaf_UID %in% meta$leaf_UID) %>% 
  dplyr::filter(timestamp <= last_valid_timestamp) %>% 
  # select leaves with at least 90% of the roi present 
  group_by(leaf_UID) %>% 
  mutate(la_tot_init = dplyr::nth(n=2, la_tot)) %>%
  dplyr::filter(la_tot >= 0.90*la_tot_init) %>%
  pull(file_id)
subset <- lesion_data %>% 
  dplyr::filter(file_id %in% files) %>% 
  extract_covars_from_nested(., from = "design", var = c("genotype_name", "exp_UID")) %>% 
  relocate(exp_UID, genotype_name, .after = design)

# add leaf-level data
subset <- left_join(subset, leaf_level_features)

# plot all selected data
p <- list()
for (exp in unique(subset$exp_UID)){
  plotdat <- subset[!is.na(subset$area),]
  plotdat <- plotdat[plotdat$exp_UID == exp,]
  p[[exp]] <- base_plot_batches +
    geom_line(data = plotdat, 
              aes(x = timestamp, y = area, color = as.factor(batch), 
                  group=interaction(plot_UID, leaf_nr, lesion_nr)), 
              alpha = 0.25, linewidth = 0.33) +
    geom_point(data = plotdat[!is.na(plotdat$area),], 
               aes(x = timestamp, y = area, color = as.factor(batch), 
                   group=interaction(plot_UID, leaf_nr, lesion_nr)),  alpha = 0.25) +
    scale_y_continuous(
      limits = c(0, 150), 
      name = bquote("Lesion area (mm" ^2~")")) +
    facet_grid(exp_UID ~ genotype_name) +
    theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))
}
p <- do.call(grid.arrange, c(p, ncol = 1))

png(paste0(figure_path, "lesion_area_time_batches.png"), width = 15, height = 8, units = 'in', res = 400)
plot(p)
dev.off()

# ============================================================================== -
# 3) Add leaf-level features and environmental data ----
# ============================================================================== -

# add heading date
# see Z:/Public/Jonas/004_ESWW007/RScripts/heading.R
hd_esww007 <- read_csv("Z:/Public/Jonas/Data/ESWW007/RefData/BBCH/hd.csv")
hd_esww009 <- read_csv("Z:/Public/Jonas/Data/ESWW009/RefData/BBCH/hd.csv")
hd <- bind_rows(hd_esww007, hd_esww009)%>% 
  dplyr::select(plot_UID, hd)
# estimated heading dates based on previous scorings and images from 06-01 and 06-05
hd[hd$plot_UID == "ESWW0070063",]$hd <- as.Date("2023-06-03")
hd[hd$plot_UID == "ESWW0070183",]$hd <- as.Date("2023-06-03")
hd[hd$plot_UID == "ESWW0090024",]$hd <- as.Date("2024-06-08")
hd[hd$plot_UID == "ESWW0090109",]$hd <- as.Date("2023-06-10")
subset <- subset %>% 
  left_join(., hd, by = "plot_UID")

# get lesion age and "leaf age" (time since leaf was first imaged, time since heading)
subset <- subset %>% 
  group_by(lesion_UID) %>% nest() %>% 
  mutate(first_timestamp_lesion = purrr::map(data, ~dplyr::first(.$timestamp)),
         last_timestamp_lesion = purrr::map(data, ~dplyr::last(.$timestamp))) %>%  
  unnest(c(first_timestamp_lesion, last_timestamp_lesion)) %>% unnest(c(data)) %>% 
  mutate(lesion_age = difftime(timestamp, first_timestamp_lesion, units = "hours"),
         leaf_age = difftime(timestamp, first_timestamp_leaf, units = "hours"),
         leaf_age_hd = difftime(date, hd, units = "days")) %>% 
  relocate(first_timestamp_lesion, last_timestamp_lesion, .after = last_timestamp_leaf)

# get symptomatic and asymptomatic duration
# symptomatic duration is measurement-specific (how much time passed between the measurement and the time point the leaf was first symptomatic)
# asymptomatic duration is leaf-specific (how long did the leaf stay clean since it was first imaged?)
subset <- subset %>% 
  group_by(leaf_UID) %>% 
  # first time a lesion was detected on each leaf
  mutate(first_timestamp_lesion_on_leaf = dplyr::first(timestamp)) %>% 
  relocate(first_timestamp_lesion_on_leaf, .after = last_timestamp_lesion) %>% 
  # convert to duration since event
  mutate(symptomatic_duration = difftime(timestamp, first_timestamp_lesion_on_leaf, units = "hours")) %>% 
  mutate(asymptomatic_duration = difftime(first_timestamp_lesion_on_leaf, first_timestamp_leaf, units = "hours")) %>% 
  mutate(asymptomatic_duration_hd = difftime(as.Date(first_timestamp_lesion_on_leaf), hd, units = "days")) %>% 
  ungroup() %>% 
  dplyr::select(-hd)

# the symptomatic duration and asymptomatic duration must add up to leaf_age for all measurements
ss <- subset %>% dplyr::select(leaf_age, symptomatic_duration, asymptomatic_duration)
ss$totdur <- ss$symptomatic_duration + ss$asymptomatic_duration
plot(ss$leaf_age~ss$totdur)

# add free perimeters as fractions of the total perimeter
subset$frac_x_perimeter_f <- subset$x_perimeter_f / subset$x_perimeter
subset$frac_y_perimeter_f <- subset$y_perimeter_f / subset$y_perimeter
subset$frac_xy_perimeter_f <- subset$xy_perimeter_f / subset$xy_perimeter
subset <- subset %>% 
  relocate(frac_x_perimeter_f, frac_y_perimeter_f, frac_xy_perimeter_f, .after = xy_perimeter_f)

# Convert chronological time to thermal time
# get environmental data
df_covar <- read_csv("Z:/Public/Jonas/Data/ESWW009/MeteoData/covar_data_position_id_/means.csv") %>% 
  dplyr::rename("temp" = "mean_tmp",
                "rh" = "mean_rh") %>% 
  # pre-filter the data set to speed up calculations
  dplyr::filter(timestamp > as.Date("2023-05-25") & timestamp < as.Date("2023-07-01") | 
                  timestamp > as.Date("2024-05-25") & timestamp < as.Date("2024-07-01"))

## Add covar courses
# Generate covar course lookup table
df_covar_lookup <- subset %>% ungroup() %>% 
  dplyr::select(timestamp, first_timestamp_lesion) %>%
  unique() %>% 
  mutate(covar_course_id = seq(n()))

# Merge lookup id to measurements
subset <- left_join(
  subset, df_covar_lookup, by = c("timestamp", "first_timestamp_lesion"))

# Extract temperature data per lookup id
df_covar_lookup <- df_covar_lookup %>%
  group_by(covar_course_id) %>%
  nest() %>%
  mutate(covar_course =
           purrr::map(data, ~extract_covar(df_covar, .$first_timestamp_lesion, .$timestamp))) %>%
  unnest(data) %>%
  dplyr::select(-timestamp, -first_timestamp_lesion)

# Join temperature course id to measurements
subset <- left_join(subset, df_covar_lookup, by = "covar_course_id")
df <- subset %>% ungroup() %>%  
  mutate(lesion_age_gdd = purrr::map_dbl(
    covar_course, 
    ~get_effective_time(., tmin=-4.52, tmax=25.6, topt=23, scale=8.57, dt=1))
  ) %>% 
  relocate(covar_course_id, covar_course, .after = batch)

saveRDS(df, paste0(data_path, "subset_step0.rds"))

# ============================================================================== -

df <- readRDS(paste0(data_path, "subset_step0.rds"))

# Time vs. GDD
p <- base_plot_batches + 
  geom_point(data = df, aes(x = lesion_age, y = lesion_age_gdd, color = as.factor(batch)), alpha = 0.3) +
  ylab("Lesion age (effective)") +
  scale_x_continuous(name = "Lesion age (hours)") +
  facet_wrap(~exp_UID)
png(paste0(figure_path, "thermal_chronological_age.png"), width = 7, height = 4, units = 'in', res = 400)
plot(p)
dev.off()

# plot all selected data
p <- list()
for (exp in unique(df$exp_UID)){
  plotdat <- df[!is.na(df$area),]
  plotdat <- plotdat[plotdat$exp_UID == exp,]
  p[[exp]] <- base_plot_batches +
    geom_line(data = plotdat[!is.na(plotdat$area),], 
              aes(x = lesion_age_gdd, y = area, color = as.factor(batch), 
                  group=interaction(plot_UID, leaf_nr, lesion_nr)), 
              alpha = 0.25, linewidth = 0.33) +
    geom_point(data = plotdat[!is.na(plotdat$area),], 
               aes(x = lesion_age_gdd, y = area, color = as.factor(batch), 
                   group=interaction(plot_UID, leaf_nr, lesion_nr)),  alpha = 0.25) +
    scale_y_continuous(limits = c(0, 100), 
                       name = bquote("Lesion area (mm" ^2~")")) +
    scale_x_continuous(name = "Lesion age (effective)") +
    facet_grid(exp_UID ~ genotype_name)
  }
p <- do.call(grid.arrange, c(p, ncol = 1))

png(paste0(figure_path, "lesion_area_age_batches_gdd.png"), width = 15, height = 8, units = 'in', res = 400)
plot(p)
dev.off()

# ============================================================================== -
# 4) Filter Lesions  ----
# ============================================================================== -

# identify lesions that do not develop sufficient pycnidia for accurate diagnosis
# during the first 5 days of their appearance
rm_lesions_1 <- df %>% ungroup() %>% 
  filter(lesion_age < 120) %>%
  group_by(lesion_UID) %>% 
  dplyr::summarise(mean_n_pycn_early = mean(n_pycn),
                   max_n_pycn_early = max(n_pycn)) %>% 
  filter(mean_n_pycn_early <= 3) %>%  # an average of 3 pycnidia detected in the lesion
  filter(max_n_pycn_early <= 4) %>% # at least 4 pycnidia detected in the lesion in at least one frame
  dplyr::select(lesion_UID) %>% 
  unique() %>% pull(lesion_UID)

# remove those lesions from data set
subset <- df %>% 
  filter(!(lesion_UID %in% rm_lesions_1))

# identify lesions with a very low pycnidia density at advanced stages
rm_lesions_2 <- subset %>% ungroup() %>% 
  filter(lesion_age > 72) %>%
  group_by(lesion_UID) %>% 
  dplyr::summarise(mean_pycn_density = mean(pycn_density_lesion, na.rm = T)) %>% 
  filter(mean_pycn_density <= 0.0005) %>% pull(lesion_UID)

# remove those lesions from data set
subset <- subset %>% 
  filter(!(lesion_UID %in% rm_lesions_2))

# identify lesions that stay very small
# and do not develop at least 2 pycnidia
rm_lesions_3 <- subset %>% 
  filter(lesion_age > 72 & lesion_age < 120) %>%
  group_by(lesion_UID) %>% 
  summarise(mean_n_pycn_early = mean(n_pycn),
            max_n_pycn_early = max(n_pycn),
            max_area_early = max(area)) %>% 
  filter(mean_n_pycn_early <= 3) %>%  # an average of 2 pycnidia detected in the lesion
  filter(max_n_pycn_early <= 4) %>% # at least 4 pycnidia detected in the lesion in at least one frame
  filter(max_area_early < 3) %>% 
  dplyr::select(lesion_UID) %>% 
  unique() %>% pull(lesion_UID)

# remove those lesions from data set
subset <- subset %>% 
  filter(!(lesion_UID %in% rm_lesions_3))

# identify lesions with too high rust density at young age
rm_lesions_4 <- subset %>% 
  group_by(lesion_UID) %>% 
  filter(lesion_age < 96) %>%
  summarise(mean_pycn_density_lesion = mean(pycn_density_lesion),
            mean_rust_density_lesion = mean(rust_density_lesion)) %>% 
  filter(mean_rust_density_lesion > 0.0003) %>% 
  dplyr::select(lesion_UID) %>% 
  unique() %>% pull(lesion_UID)

# remove those lesions from data set
subset <- subset %>% 
  filter(!(lesion_UID %in% rm_lesions_4))
  
# identify lesions with too high rust density at young age
rm_lesions_5 <- subset %>% 
  group_by(lesion_UID) %>% 
  filter(lesion_age < 96) %>%
  summarise(mean_pycn_density_lesion = mean(pycn_density_lesion),
            mean_rust_density_lesion = mean(rust_density_lesion)) %>%   
  filter(3*mean_rust_density_lesion > mean_pycn_density_lesion) %>% 
  dplyr::select(lesion_UID) %>% 
  unique() %>% pull(lesion_UID)

# remove those lesions from data set
subset <- subset %>% 
  filter(!(lesion_UID %in% rm_lesions_5))

saveRDS(subset, paste0(data_path, "subset_step0_for_epiparams.rds"))

# exclude observations of highly "trapped" lesions
# and observations of lesions without an area (e.g., partially outside of roi)
sub <- subset %>% 
  mutate(lag_analyzable_perimeter = dplyr::lag(analyzable_perimeter)) %>%
  # lag is NA for the first detection of a lesion - we keep these,
  # since it should be > 0.5 in most cases
  dplyr::filter(is.na(lag_analyzable_perimeter) | lag_analyzable_perimeter > 0.5) %>%
  dplyr::select(-lag_analyzable_perimeter) %>%  # to avoid conflicts further down
  filter(!is.na(area))

# lesion area vs measurement time point
p <- list()
for (exp in unique(sub$exp_UID)){
  plotdat <- sub[!is.na(sub$area),]
  plotdat <- plotdat[plotdat$exp_UID == exp,]
  p[[exp]] <- base_plot_batches +
    geom_line(data = plotdat, 
              aes(x = timestamp, y = area, color = as.factor(batch), 
                  group=interaction(plot_UID, leaf_nr, lesion_nr)),
              alpha = 0.25, linewidth = 0.33) +
    geom_point(data = plotdat[!is.na(plotdat$area),], 
               aes(x = timestamp, y = area, color = as.factor(batch), 
                   group=interaction(plot_UID, leaf_nr, lesion_nr)),  alpha = 0.25) +
    scale_y_continuous(limits = c(0, 150),
                       name = bquote("Lesion area (mm" ^2~")")) +
    facet_grid(exp_UID ~ genotype_name) +
    theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))
  }
p <- do.call(grid.arrange, c(p, ncol = 1))

png(paste0(figure_path, "01_lesion_development_timestamp.png"), width = 15, height = 6, units = 'in', res = 400)
plot(p)
dev.off()

# lesion area vs lesion age (GDD)
p <- list()
for (exp in unique(subset$exp_UID)){
  plotdat <- subset[!is.na(subset$area),]
  plotdat <- plotdat[plotdat$exp_UID == exp,]
  n_obs_exp <- n_obs[n_obs$exp_UID == exp,] 
  p[[exp]] <- base_plot_batches +
    geom_line(data = plotdat, 
              aes(x = lesion_age_gdd, y = area, color = as.factor(batch), 
                  group=interaction(plot_UID, leaf_nr, lesion_nr)),
              alpha = 0.15, linewidth = 0.33) +
    geom_point(data = plotdat[!is.na(plotdat$area),], 
               aes(x = lesion_age_gdd, y = area, color = as.factor(batch), 
                   group=interaction(plot_UID, leaf_nr, lesion_nr)),  alpha = 0.15) +
    scale_y_continuous(limits = c(0, 100), 
                       name = bquote("Lesion area (mm" ^2~")")) +
    scale_x_continuous(name = "Lesion age (effective)") +
    facet_grid(exp_UID ~ genotype_name) +
    geom_text(data = n_obs_exp, aes(x = 0, y = Inf, label = paste0("leaves: ", n_leaf, "\nlesions: ", n_lesion, "\nintervals: ", n_interval)),
              hjust = 0, vjust = 1.1, inherit.aes = FALSE, size = 2.1, color = "darkblue") +
    theme(axis.text.x=element_text(angle = 0, vjust = 1, hjust=1))
}
p <- do.call(grid.arrange, c(p, ncol = 1))

png(paste0(figure_path, "02_lesion_development_abs_genotypes_gdd.png"), width = 15, height = 5.5, units = 'in', res = 600)
plot(p)
dev.off()

saveRDS(sub, paste0(data_path, "subset_step_0_filtered.rds"))

# ============================================================================== -
# 5) Convert to intervals ----
# ============================================================================== -

sub <- readRDS(paste0(data_path, "subset_step_0_filtered.rds"))

# Get lag variables and interval differences
s <- sub %>% ungroup() %>%
  dplyr::select(-genotype_name) %>%  # conflicting with lag extraction
  group_by(plot_UID, leaf_nr, lesion_nr) %>% nest() %>%
  ungroup() %>%
  mutate(data = purrr::map(data, get_lags_diffs, vars = c(18:65))) %>%
  unnest(c(data))
saveRDS(s, paste0(data_path, "subset_step_0_filtered_withlags.rds"))

# ============================================================================== -
# 6) Filter interval data ----
# ============================================================================== -

# get normalized growth
s <- readRDS(paste0(data_path, "subset_step_0_filtered_withlags.rds"))
s <- s %>% 
  # remove intervals without lag timepoint
  dplyr::filter(!is.na(lag_timestamp)) %>% 
  # remove short intervals
  filter(diff_time > 5) %>% 
  # calculate normalized differences
  mutate(diff_area_norm_chr = diff_area / as.numeric(diff_time),
         diff_area_norm_gdd = diff_area / diff_gdd,
         diff_area_pp_y_norm_chr = diff_area_pp_y/as.numeric(diff_time),
         diff_area_pp_y_norm_gdd = diff_area_pp_y/diff_gdd, 
         diff_area_pp_x_norm_chr = diff_area_pp_x/as.numeric(diff_time),
         diff_area_pp_x_norm_gdd = diff_area_pp_x/diff_gdd,
         diff_area_pp_xy_norm_chr = diff_area_pp_xy/as.numeric(diff_time),
         diff_area_pp_xy_norm_gdd = diff_area_pp_xy/diff_gdd)


# plot cumulative distribution of percentiles for UN-filtered data
x_sorted <- sort(s$diff_area_norm_chr)
percentiles <- ecdf(x_sorted)(x_sorted) * 100
df <- data.frame(value = x_sorted, percentile = percentiles)
A <- ggplot(df, aes(x = value, y = percentile)) +
  geom_line(color = "blue", linewidth = 1.2) +
  labs(x = "Value", y = "Cumulative Percentile") +
  ggtitle("A")

# get quantiles for each growth measure
quants <- s %>% 
  pivot_longer(cols = diff_area_norm_chr:diff_area_pp_xy_norm_gdd) %>% 
  group_by(name) %>% 
  summarise(
    x_min = quantile(value, 0.01),
    x_max = quantile(value, 0.98)
  )
  
# filter observations to be in range for all growth measures
for(colname in quants$name) {
  print(colname)
  lower_limit <- quants[quants$name == colname, ]$x_min
  upper_limit <- quants[quants$name == colname, ]$x_max
  s <- s[s[[colname]] >= lower_limit & s[[colname]] <= upper_limit, ]
}

# get mean environmental variables for intervals
s <- s %>% 
  mutate(mean_interval_temp = purrr::map_dbl(covar_course, ~mean(.$temp, na.rm = T)),
         mean_interval_rh = purrr::map_dbl(covar_course, ~mean(.$rh, na.rm = T)),
         cv_interval_temp = purrr::map_dbl(covar_course, ~sd(.$temp, na.rm = TRUE) / mean(.$temp, na.rm = TRUE)),
         cv_interval_rh = purrr::map_dbl(covar_course, ~sd(.$rh, na.rm = TRUE) / mean(.$rh, na.rm = TRUE)))

# statistics on intervals
interval_stats <- s %>% ungroup() %>% 
  dplyr::summarise(mean_int = mean(diff_time, na.rm = T),
                   sd_int = sd(diff_time, na.rm = T), 
                   min_int = min(diff_time, na.rm = T),
                   max_int = max(diff_time, na.rm = T))

# a problem with one lesion (ESWW0020_9_9); remove
s <- s %>% dplyr::filter(lesion_UID != "ESWW0070020_9_9")

saveRDS(s, paste0(data_path, "subset_step2.rds"))
# data ready for feature selection

# data set size
s <- s %>% extract_covars_from_nested("design", "genotype_name")
n_lesions <- s %>% group_by(lesion_UID) %>% nest()
n_intervals_lesion <- n_lesions %>% ungroup() %>% 
  mutate(n_obs = purrr::map_dbl(data, nrow)) %>% 
  dplyr::summarise(mean_n_obs = mean(n_obs),
                   min_n_obs = min(n_obs),
                   max_n_obs = max(n_obs))
n_leaves <- s %>% group_by(leaf_UID) %>% nest()
n_leafs_geno <- s %>% 
  group_by(exp_UID, genotype_name, leaf_UID) %>% nest() %>% 
  group_by(exp_UID, genotype_name) %>% nest() %>% 
  mutate(n_leaf = purrr::map_dbl(data, nrow)) %>% select(-data)
n_lesions_geno <- s  %>%
  group_by(genotype_name, lesion_UID) %>% nest() %>% 
  group_by(genotype_name) %>% nest() %>% 
  mutate(n_lesion = purrr::map_dbl(data, nrow)) %>% select(-data)
n_interval_geno <- s %>% group_by(genotype_name) %>% nest() %>% 
  mutate(n_interval = purrr::map_dbl(data, nrow)) %>% select(-data)
n_obs <- full_join(n_leafs_geno, n_lesions_geno) %>% 
  full_join(., n_interval_geno)

# Time vs. GDD
p <- base_plot_batches + 
  geom_point(data = s, aes(x = diff_time, y = diff_gdd, color = as.factor(batch)), alpha = 0.2) +
  ylab("Measurement interval (effective)") +
  scale_x_continuous(
    name = "Measurement interval (hours)"
  ) +
  facet_wrap(~exp_UID)
png(paste0(figure_path, "thermal_chronological_interval.png"), width = 7, height = 4, units = 'in', res = 400)
plot(p)
dev.off()

# ============================================================================== -

subset <- readRDS(paste0(data_path, "subset_step2.rds"))

# Sort and compute percentiles
x_sorted <- sort(subset$diff_area_norm_chr)
percentiles <- ecdf(x_sorted)(x_sorted) * 100
df <- tibble(value = x_sorted, percentile = percentiles)
B <- ggplot(df, aes(x = value, y = percentile)) +
  geom_line(color = "blue", linewidth = 1.2) +
  labs(x = "Value", y = "Cumulative Percentile", title = "B") +
  scale_x_continuous(
    limits = c(-0.1, 0.55), 
    breaks = c(0, 0.2, 0.4)
  )
A <- A + 
  geom_vline(xintercept = min(subset$diff_area_norm_chr), color = "red") +
  geom_vline(xintercept = max(subset$diff_area_norm_chr), color = "red")

p <- ggarrange(A, B, nrow = 1)
png(paste0(figure_path, "cum_perc.png"), width = 7, height = 3.5, units = 'in', res = 400)
plot(p)
dev.off()

# plot l_density vs. p_density
p <- ggplot(subset) +
  geom_point(aes(x = mean_l_density, y = mean_p_density)) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlab("Mean lesion pycnidia density") +
  ylab("Mean pycnidiation area pycnidia density")
png(paste0(figure_path, "p_l_pycn_density.png"), width = 7, height = 5, units = 'in', res = 400)
plot(p)
dev.off()

# ============================================================================== -
# 7) Growth ~ Perimeters ----
# ============================================================================== -

# New facet label names for direction variable
dir_labs <- c("X (Major leaf axis)", "Y (Minor leaf axis)", "XY (Major + minor leaf axes)")
names(dir_labs) <- c("x_perimeter", "y_perimeter", "xy_perimeter")

# reshape data for plotting
pdat <- subset %>% 
  dplyr::select(plot_UID, leaf_nr, lesion_nr, 
                lag_xy_perimeter, lag_x_perimeter, lag_y_perimeter,
                lag_xy_perimeter_f, lag_x_perimeter_f, lag_y_perimeter_f,
                diff_area_norm_gdd) %>% 
  pivot_longer(cols = lag_xy_perimeter:lag_y_perimeter_f, names_to = "perimeter") %>% 
  mutate(perimeter_type = ifelse(grepl("_f", perimeter), "Free", "Total")) %>% 
  mutate(perimeter = gsub("lag_", "", perimeter),
         perimeter = gsub("_f", "", perimeter)) %>% 
  mutate(perimeter = as.factor(perimeter),
         perimeter_type = as.factor(perimeter_type))

# can fit loess only if sufficient data is available
# must limit the range of the predictor variable slightly 
# this is means ignoring few extremely large values perimeter values (keep 99.5%)
loess_limits <- pdat %>%
  group_by(perimeter, perimeter_type) %>%
  summarise(
    x_min = quantile(value, 0),
    x_max = quantile(value, 0.99)
  )
mdat <- pdat %>%
  left_join(., loess_limits, by = c("perimeter", "perimeter_type")) %>% 
  filter(value >= x_min & value <= x_max)

# # Fit a median quantile regression model
# qr_models <-  mdat%>% 
#   group_by(perimeter, perimeter_type) %>% nest() %>% ungroup() %>% 
#   mutate(
#     linear_q = purrr::map(data, ~ rq(diff_area_norm_gdd ~ poly(value, 2), data = .)),
#     loess_q = purrr::map(data, ~ lprq(x = .$value,
#                                       y = .$diff_area_norm_gdd,
#                                       h = (range(.$value)[2]-range(.$value)[1])/20,  # 25 works, with m set to default value,
#                                       tau = .5)),
#     new_data = purrr::map(data, ~ data.frame(value = seq(min(.$value), max(.$value), length.out = 100)))
#   )
# 
# # Make predictions for both QR and LPQR models
# qr_models2 <- qr_models %>% 
#   mutate(
#     predicted_qr = purrr::map2(linear_q, new_data, ~ predict(.x, newdata = .y)),
#     predicted_lpqr = purrr::map2(loess_q, new_data, ~ {
#       fit_values <- .x$fv
#       predict_values <- approx(.x$xx, fit_values, xout = .y$value)$y
#       return(predict_values)
#     })
#   )
# 
# # add pseudo R2
# qr_models3 <- qr_models2 %>% 
#   mutate(pseudo_r2 = purrr::map2(data, linear_q, ~ compute_pseudoR2(.x, response = "diff_area_norm_gdd", qr_model = .y)))
# 
# # extract pseudo R2
# qr_models4 <- qr_models3 %>% 
#   mutate(r2 = map_dbl(pseudo_r2, 1))
# # save
# saveRDS(qr_models4, paste0(data_path, "fits_perimeters.rds"))


qr_models4 <- readRDS(paste0(data_path, "fits_perimeters.rds"))

# extract predictions
predictions <- qr_models4 %>% 
  dplyr::select(perimeter, perimeter_type, new_data, predicted_qr, predicted_lpqr) %>% 
  unnest(cols = c(new_data, predicted_qr, predicted_lpqr))

# set axis limits to include 99.9% of the data
# trends not visible otherwise due to extreme outliers
plot_limits <- pdat %>%
  group_by(perimeter) %>%
  filter(perimeter_type == "Total") %>% 
  summarise(
    x_min = quantile(value, 0),
    x_max = quantile(value, 0.99)
  )
pdat <- pdat %>%
  left_join(., plot_limits, by = c("perimeter")) %>% 
  filter(value >= x_min & value <= x_max)

# plot data with model fits
pA <- ggplot() +
  geom_point(data = pdat, aes(x = value, y = diff_area_norm_gdd), alpha = 0.02) +
  # QR loess
  geom_line(data = predictions, aes(y = predicted_lpqr, x = value), color = "green", size = 1) +
  # QR linear
  geom_line(data = predictions, aes(y = predicted_qr, x = value), size = 1, color = "yellow") +
  geom_abline(intercept = 0, slope = 0, lty = 2) +
  scale_y_continuous(
    name = bquote("Lesion growth (" ~ A[t[1]] - A[t[0]] ~ ")"),
    breaks = seq(0, 5, 1)) +
  scale_x_continuous(name = "Perimeter length (pixels)") +
  facet_grid(perimeter_type ~ perimeter, scales = "free",
             labeller = labeller(perimeter = dir_labs)) +
  geom_text(data = qr_models4, aes(x = Inf, y = Inf, label = paste0("R²: ", round(r2, 2))),
            hjust = 1.1, vjust = 1.1, inherit.aes = FALSE)
  # + ggtitle("A")
png(paste0(figure_path, "03_growth_perimeters_models_new_norect.png"), width = 7, height = 5, units = 'in', res = 400)
plot(pA)
dev.off()

# ============================================================================== -
# 8) Growth ~ Init Area ----
# ============================================================================== -

subset <- readRDS(paste0(data_path, "subset_step2.rds"))

init <- subset %>% dplyr::filter(lag_lesion_age == 0)
avg <- subset %>% group_by(lesion_UID) %>% 
  summarise(mean_delta = mean(diff_area_pp_y_norm_gdd),
            mean_age = mean(lesion_age))
both <- left_join(init, avg, by = "lesion_UID")
both <- both %>% extract_covars_from_nested("design", "genotype_name")

# check distribution
# very sparse data at large area
ggplot(both) +
  geom_point(aes(x = area, y=mean_delta))

# can fit loess only if sufficient data is available
# must limit the range of the predictor variable slightly 
# this is means ignoring few extremely large values perimeter values (keep 99.5%)
loess_limits <- both %>%
  summarise(
    x_min = quantile(area, 0),
    x_max = quantile(area, 0.99)
  )
both <- both %>% 
  filter(area >= loess_limits$x_min[[1]] & area <= loess_limits$x_max[[1]])

# Fit a quantile regression model
qr_models <- both %>% nest() %>% 
  mutate(
    linear_q = purrr::map(data, ~ rq(mean_delta ~ area, data = .)),
    loess_q = purrr::map(data, ~ lprq(x = .$area, 
                                      y = .$mean_delta, 
                                      h = (range(.$area)[2]-range(.$area)[1])/20,
                                      tau = .5)),
    new_data = purrr::map(data, ~ data.frame(area = seq(min(.$area), max(.$area), length.out = 1000)))
  )

# Make predictions for both QR and LPQR models
qr_models2 <- qr_models %>% 
  mutate(
    predicted_qr = purrr::map2(linear_q, new_data, ~ predict(.x, newdata = .y)),
    predicted_lpqr = purrr::map2(loess_q, new_data, ~ {
      fit_values <- .x$fv
      predict_values <- approx(.x$xx, fit_values, xout = .y$area)$y
      return(predict_values)
    })
  )

qr_models3 <- qr_models2 %>% 
  mutate(pseudo_r2 = purrr::map2(data, linear_q, ~ compute_pseudoR2(.x, response = "mean_delta", qr_model = .y)))

qr_models4 <- qr_models3 %>% 
  mutate(r2 = map_dbl(pseudo_r2, 1))

predictions <- qr_models2 %>% 
  dplyr::select(new_data, predicted_qr, predicted_lpqr) %>% 
  unnest(cols = c(new_data, predicted_qr, predicted_lpqr))

pB <- ggplot(both) +
  scale_y_continuous(
    name = bquote("Mean Lesion Growth (" ~ A[t[1]] - A[t[0]] ~ ")")) +
  scale_x_continuous(
    name = bquote("Initial Lesion Area (mm" ^2~")")) +
  geom_point(aes(x = area, y=mean_delta), alpha = 0.1) +
  geom_line(data = predictions, aes(y = predicted_lpqr, x = area), color = "green", size = 1) +
  # QR linear
  geom_line(data = predictions, aes(y = predicted_qr, x = area), size = 1, color = "yellow") +
  geom_text(data = qr_models4, aes(x = Inf, y = Inf, label = paste0("R²: ", round(r2, 2))),
            hjust = 1.35, vjust = 2, inherit.aes = FALSE) 
  # + ggtitle("B")

png(paste0(figure_path, "init_growth.png"), width = 3.5, height = 3.5, units = 'in', res = 400)
plot(pB)
dev.off()

# Modify for multi-plot
pB <- pB +
  theme(plot.margin = margin(t = 50, b = 50))  # adjust values as needed)
p <- ggarrange(pA, pB, nrow = 1, widths = c(2.5, 1))
png(paste0(figure_path, "Figure3.png"), width = 9, height = 4.5, units = 'in', res = 400)
plot(p)
dev.off()

# Not related to leaf age
ggplot(both) +
  geom_point(aes(x = leaf_age_hd, y = area)) +
  geom_smooth(aes(x = leaf_age_hd, y = area))

# Not strongly related to lesion age
ggplot(both) +
  geom_point(aes(x = mean_age, y = mean_delta)) +
  geom_smooth(aes(x = mean_age, y = mean_delta))

# interactions?
model <- lm(mean_delta ~ area + mean_age + batch + genotype_name, data = both)
summary(model)
anova(model)

library(car)
Anova(model)
Anova(model, type = "III")

# ============================================================================== -
# 9) Remaining Growth ~ Area ----
# ============================================================================== -

# SHOW FOR ALL DIMENSIONS
mdat0 <- readRDS(paste0(data_path, "subset_step2.rds"))

# reshape data for plotting
pdat <- mdat0 %>% 
  dplyr::select(plot_UID, leaf_nr, lesion_nr, 
                area, diff_area_norm_gdd, 
                diff_area_pp_x_norm_gdd,
                diff_area_pp_y_norm_gdd,
                diff_area_pp_xy_norm_gdd) %>% 
  pivot_longer(cols = 5:8)

# can fit loess only if sufficient data is available
# must limit the range of the predictor variable slightly 
# this is means ignoring few extremely large values perimeter values (keep 99.5%)
loess_limits <- pdat %>%
  group_by(name) %>%
  summarise(
    x_min = quantile(area, 0),
    x_max = quantile(area, 0.99)
  )
pdat <- pdat %>% 
  left_join(loess_limits, by = c("name")) %>%
  filter(area >= x_min & area <= x_max)

# # Fit a quantile regression model
# qr_models <- pdat %>% 
#   group_by(name) %>% nest() %>% ungroup() %>% 
#   mutate(
#     linear_q = purrr::map(data, ~ rq(value ~ area, data = .)),
#     linear_q_sqrt = purrr::map(data, ~ rq(value ~ area + sqrt(area), data = .)),
#     loess_q = purrr::map(data, ~ lprq(x = .$area, 
#                                       y = .$value, 
#                                       h = (range(.$area)[2]-range(.$area)[1])/20, 
#                                       tau = .5)),
#     new_data = purrr::map(data, ~ data.frame(area = seq(min(.$area), max(.$area), length.out = 1000)))
#   )
# 
# # compare linear and sqrt-transformed 
# AIC(qr_models$linear_q[[1]])
# AIC(qr_models$linear_q_sqrt[[1]])
# 
# # Make predictions for both QR and LPQR models
# qr_models2 <- qr_models %>% 
#   mutate(
#     predicted_qr = purrr::map2(linear_q, new_data, ~ predict(.x, newdata = .y)),
#     predicted_qr_sqrt = purrr::map2(linear_q_sqrt, new_data, ~ predict(.x, newdata = .y)),
#     predicted_lpqr = purrr::map2(loess_q, new_data, ~ {
#       fit_values <- .x$fv
#       predict_values <- approx(.x$xx, fit_values, xout = .y$area)$y
#       return(predict_values)
#     })
#   )
# 
# # add pseudo R2
# qr_models3 <- qr_models2 %>% 
#   mutate(pseudo_r2 = purrr::map2(data, linear_q_sqrt, ~ compute_pseudoR2(.x, response = "value", qr_model = .y)))
# 
# # extract pseudo R2
# qr_models4 <- qr_models3 %>% 
#   mutate(r2 = map_dbl(pseudo_r2, 1))
# 
# # save
# saveRDS(qr_models4, paste0(data_path, "fits_area.rds"))
qr_models4 <- readRDS(paste0(data_path, "fits_area.rds"))

# extract predictions
predictions <- qr_models4 %>% 
  dplyr::select(name, new_data, predicted_qr, predicted_qr_sqrt, predicted_lpqr) %>% 
  unnest(cols = c(new_data, predicted_qr, predicted_qr_sqrt, predicted_lpqr))

# New facet label names for direction variable
dir_labs <- c("None", "X (Major leaf axis)", "XY (Major + Minor leaf axes)", "Y (Minor leaf axis")
names(dir_labs) <- c("diff_area_norm_gdd", "diff_area_pp_x_norm_gdd", "diff_area_pp_xy_norm_gdd", "diff_area_pp_y_norm_gdd")
p0 <- ggplot() +
  geom_point(data = pdat, aes(x = area, y = value), alpha = 0.02) +
  # # qr loess
  geom_line(data = predictions, aes(y = predicted_lpqr, x = area), color = "green", size = 1) +
  geom_abline(intercept = 0, slope = 0, lty = 2) +
  scale_y_continuous(
    name = bquote("Lesion growth (" ~ A[t[1]] - A[t[0]] ~ ")")) +
  scale_x_continuous(limits = c(0, 50), 
                     name = bquote("Lesion area at"~t[0]~" (mm" ^2~")")) +
  facet_wrap(~ name, scales = "free", nrow = 1,
             labeller = labeller(name = dir_labs)) +
  geom_text(data = qr_models4, aes(x = Inf, y = Inf, label = paste0("R²: ", round(r2, 2))),
            hjust = 1.1, vjust = 1.1, inherit.aes = FALSE) +
  ggtitle("A")

# sqrt
p <- p0 +
  geom_line(data = predictions, aes(y = predicted_qr_sqrt, x = area), size = 1, color = "yellow") +
png(paste0(figure_path, "remaining_growth_area_sqrt.png"), width = 10, height = 3.5, units = 'in', res = 400)
plot(p)
dev.off()

# linear
p <- p0 +
  geom_line(data = predictions, aes(y = predicted_qr, x = area), size = 1, color = "yellow") +
png(paste0(figure_path, "remaining_growth_area_linear.png"), width = 10, height = 3.5, units = 'in', res = 400)
plot(p)
dev.off()

# ============================================================================== -
# 10) Area ~ Age ----
# ============================================================================== -

data <- readRDS(paste0(data_path, "subset_step2.rds"))

# can fit loess only if sufficient data is available
# must limit the range of the predictor variable slightly 
# this is means ignoring few extremely large values perimeter values (keep 99.5%)
loess_limits <- data %>%
  summarise(
    x_min = quantile(area, 0),
    x_max = quantile(area, 0.99),
    y_min = quantile(lesion_age_gdd, 0),
    y_max = quantile(lesion_age_gdd, 0.99)
  )
mdat <- data %>% 
  filter(area >= loess_limits$x_min & area <= loess_limits$x_max) %>% 
  filter(lesion_age_gdd >= loess_limits$y_min & lesion_age_gdd <= loess_limits$y_max)

# Fit a quantile regression model
qr_models <- mdat %>% 
 nest() %>% ungroup() %>% 
  mutate(
    linear_q = purrr::map(data, ~ rq(area ~ lesion_age_gdd, data = .)),
    loess_q = purrr::map(data, ~ lprq(x = .$lesion_age_gdd, 
                                      y = .$area, 
                                      h = (range(.$lesion_age_gdd)[2]-range(.$lesion_age_gdd)[1])/20, 
                                      tau = .5)),
    new_data = purrr::map(data, ~ data.frame(lesion_age_gdd = seq(min(.$lesion_age_gdd), max(.$lesion_age_gdd), length.out = 1000)))
  )

# Make predictions for both QR and LPQR models
qr_models2 <- qr_models %>% 
  mutate(
    predicted_qr = purrr::map2(linear_q, new_data, ~ predict(.x, newdata = .y)),
    predicted_lpqr = purrr::map2(loess_q, new_data, ~ {
      fit_values <- .x$fv
      predict_values <- approx(.x$xx, fit_values, xout = .y$lesion_age_gdd)$y
      return(predict_values)
    })
  )

# add pseudo R2
qr_models3 <- qr_models2 %>% 
  mutate(pseudo_r2 = purrr::map2(data, linear_q, ~ compute_pseudoR2(.x, response = "lesion_age_gdd", qr_model = .y)))

# extract pseudo R2
qr_models4 <- qr_models3 %>% 
  mutate(r2 = map_dbl(pseudo_r2, 1))

# extract predictions
predictions <- qr_models2 %>% 
  dplyr::select(new_data, predicted_qr, predicted_lpqr) %>% 
  unnest(cols = c(new_data, predicted_qr, predicted_lpqr))

plot <- ggplot(mdat)+
  geom_point(aes(x = lesion_age_gdd, y = area), alpha = 0.02) +
  geom_line(data = predictions, aes(y = predicted_lpqr, x = lesion_age_gdd), color = "green", size = 1) +
  geom_line(data = predictions, aes(y = predicted_qr, x = lesion_age_gdd), size = 1, color = "yellow") + 
  geom_text(data = qr_models4, aes(x = Inf, y = Inf, label = paste0("R²: ", round(r2, 2))),
            hjust = 1.5, vjust = 3, inherit.aes = FALSE) +
  scale_y_continuous(name = bquote("Lesion area (mm" ^2~")")) +
  scale_x_continuous(name = "Lesion age (effective)")
png(paste0(figure_path, "age_area.png"), width = 3.5, height = 3.5, units = 'in', res = 400)
plot(plot)
dev.off()

# ============================================================================== -
# 10) Growth ~ Age ----
# ============================================================================== -

mdat <- readRDS(paste0(data_path, "subset_step2.rds"))

# check distribution
# very sparse data at large area
ggplot(mdat) +
  geom_point(aes(x = lesion_age_gdd, y = diff_area_pp_y_norm_gdd))

# can fit loess only if sufficient data is available
# must limit the range of the predictor variable slightly 
# this is means ignoring few extremely large values perimeter values (keep 99.5%)
loess_limits <- mdat %>%
  summarise(
    x_min = quantile(lesion_age_gdd, 0),
    x_max = quantile(lesion_age_gdd, 0.995)
  )
mdat <- mdat %>% 
  filter(lesion_age_gdd >= loess_limits$x_min[[1]] & lesion_age_gdd <= loess_limits$x_max[[1]])

# df for predictions
new_preds <- mdat %>% ungroup() %>% 
  do(., data.frame(
    lesion_age_gdd = seq(round(min(as.numeric(.$lesion_age_gdd))), max(as.numeric(.$lesion_age_gdd)), by = 0.2),
    stringsAsFactors = FALSE)
  )

# linear quantile regression
rqmodel <- rq(diff_area_pp_y_norm_gdd ~ lesion_age_gdd, tau = .5, data = mdat)

# loess
fit <- lprq(mdat$lesion_age_gdd, mdat$diff_area_pp_y_norm_gdd, 
            h = (range(mdat$lesion_age_gdd)[2]-range(mdat$lesion_age_gdd)[1])/20, 
            tau = .5)
rq_yy <- predict(rqmodel, newdata = new_preds)

## Fit models ================================================================== -

x = "lesion_age_gdd"
y = "diff_area_pp_y_norm_gdd"

# fits <- mdat %>% nest() %>% 
#   mutate(linear_q = purrr::map(.x = data, .f = linear_quantile, x=x, y=y),
#          nls_q_exp = purrr::map(.x = data, .f = nls_quantile_exp, n_samples = 300, x=x, y=y)) %>%  
#   tidyr::pivot_longer(cols = linear_q:nls_q_exp, names_to = "type", values_to = "fit")
# saveRDS(fits, paste0(data_path, "model_fits.rds"))

fits <- readRDS(paste0(data_path, "model_fits.rds"))

# get pseudo-r2
resid_fit <- residuals(fits$fit[[2]])
rho_tau <- function(u, tau) {
  u * (tau - (u < 0))
}
numerator <- sum(rho_tau(resid_fit, tau=.5))
denominator <- sum(rho_tau(mdat[[y]] - median(mdat[[y]]), tau=.5))
pseudo_R2 <- 1 - numerator / denominator

pdat <- fits %>% 
  mutate(preds = purrr::map(fit, broom::augment, newdata = new_preds))
preds <- pdat %>% dplyr::select(type, preds) %>% unnest(preds) %>% 
  dplyr::filter(type ==  "nls_q_exp")
pd <- pdat %>% 
  dplyr::filter(type ==  "nls_q_exp") %>%  ## only one dataset needed
  dplyr::select(type, data) %>% unnest(data)

p1 <- ggplot() +
  geom_point(data = pd, aes(x = lesion_age_gdd, y = diff_area_pp_y_norm_gdd), alpha = 0.025) +
  # add loess fit
  geom_line(aes(x = fit$xx, y = fit$fv), color = "green") +
  # add linear and non-linear fit
  geom_line(data = preds, aes(x = lesion_age_gdd, y = .fitted), color = "yellow") +
  scale_y_continuous(
    name = bquote("Lesion growth (" ~ A[t[1]] - A[t[0]] ~ ")"),
  ) +
  geom_abline(intercept = 0, slope = 0, lty = 2) +
  geom_text(aes(x = Inf, y = Inf, label = paste0("R²: ", round(pseudo_R2, 2))),
            hjust = 1.2, vjust = 2, inherit.aes = FALSE) +
  theme(
    legend.position = "None",
    panel.background = element_rect(fill = "#E5E5E5")
  ) + 
  xlab(bquote("Lesion age at" ~t[1]~"(effective)"* phantom("X")))

png(paste0(figure_path, "05_fits_growth_age_all_density.png"), width = 3.5, height = 3.5, units = 'in', res = 400)
plot(p1)
dev.off()

# ============================================================================== -

# distribution of GDD-nomalized growth is OK
y_dens <- ggplot(pdat) +
  geom_density(aes(y = value)) +
  scale_y_continuous(limits = c(-.002, .01))

# is this not in contradiction with the trend observed for lesion age (see below) ?
pdat_peri <- mdat0 %>% 
  dplyr::select(plot_UID, leaf_nr, lesion_nr, 
                area, lesion_age_gdd, x_perimeter, y_perimeter, xy_perimeter) %>% 
  dplyr::filter(area <= 150)
p <- 5
p1 <- ggpairs(pdat_peri[4:8],
              lower = list(continuous = wrap("points", alpha = 0.02),
                           theme = theme(panel.background = element_blank(), 
                                         panel.grid = element_blank())),
              diag = list(continuous = "densityDiag"),
              upper = list(continuous = wrap("cor", size = 3, color = "black")))

# Generate the correlation plot
p2 <- ggcorr(pdat_peri[4:8], label = TRUE, label_round = 2) +
  scale_fill_distiller(palette = "BrBG", direction = 1, limits = c(-1, 1))
g2 <- ggplotGrob(p2)
colors <- g2$grobs[[6]]$children[[3]]$gp$fill

# Change background color to tiles in the upper triangular matrix of plots 
idx <- 1
for (k1 in 1:(p-1)) {
  for (k2 in (k1+1):p) {
    plt <- getPlot(p1,k1,k2) +
      theme(panel.background = element_rect(fill = colors[idx], color="white"),
            panel.grid.major = element_line(color=colors[idx]))
    p1 <- putPlot(p1,plt,k1,k2)
    idx <- idx+1
  }
}

png(paste0(figure_path, "area_age_peri.png"), width = 6, height = 6, units = 'in', res = 400)
p1
dev.off()

# ============================================================================== -
# 6) Environmental main effects ----
# ============================================================================== -

mdat <- readRDS(paste0(data_path, "subset_step2.rds")) %>% 
  extract_covars_from_nested("design", "genotype_name") %>% 
  mutate(batch_UID = paste(exp_UID, batch, sep = "_")) %>% 
  dplyr::rename(diff = "diff_area_pp_y_norm_chr")

pdat <- mdat %>% 
  dplyr::select(plot_UID, leaf_nr, lesion_nr, 
                diff, 
                cv_interval_rh, 
                mean_interval_rh,
                cv_interval_temp, 
                mean_interval_temp) %>% 
  pivot_longer(cols = 5:8)

lims <- pdat %>%
  group_by(name) %>%
  summarise(
    x_min = quantile(value, 0.005),
    x_max = quantile(value, 0.995),
    y_min = quantile(diff, 0.005),
    y_max = quantile(diff, 0.995)
  )

pdat_sub <- pdat %>% 
  left_join(lims, by = c("name")) %>%
  filter(value >= x_min & value <= x_max & diff >= y_min & diff <= y_max) %>% 
  group_by(name)

# Fit a quantile regression model
qr_models <- pdat_sub %>% 
  group_by(name) %>% nest() %>% ungroup() %>% 
  mutate(linear_q = purrr::map(data, ~ rq(diff ~ value, data = .)), 
         # loess_q = purrr::map(data, ~ lprq(x = .$value, y = .$diff, m = nrow(.)/10, h = 0.1*range(.$value), tau = .5)),
         new_data = purrr::map(data, ~ data.frame(value = seq(min(.$value), max(.$value), length.out = 1000))))

# Make predictions for both QR and LPQR models
qr_models2 <- qr_models %>% 
  mutate(
    predicted_qr = purrr::map2(linear_q, new_data, ~ predict(.x, newdata = .y)),
    # predicted_lpqr = purrr::map2(loess_q, new_data, ~ {
    #   fit_values <- .x$fv
    #   predict_values <- approx(.x$xx, fit_values, xout = .y$value)$y
    #   return(predict_values)
    # })
  )

qr_models3 <- qr_models2 %>% 
  mutate(pseudo_r2 = purrr::map2(data, linear_q, ~ compute_pseudoR2(.x, response = "diff", qr_model = .y)))

qr_models4 <- qr_models3 %>% 
  mutate(r2 = map_dbl(pseudo_r2, 1))

predictions <- qr_models2 %>% 
  dplyr::select(name, new_data, predicted_qr) %>% 
  unnest(cols = c(new_data, predicted_qr))

# New facet label names for direction variable
dir_labs <- c("CV Relative Humidity", "CV Temperature", "Mean Temperature", "Mean Relative Humidity")
names(dir_labs) <- c("cv_interval_rh", "cv_interval_temp", "mean_interval_temp", "mean_interval_rh")

p <- ggplot() +
  geom_point(data = pdat_sub, aes(x = value, y = diff), alpha = 0.02) +
  # # qr loess
  # geom_line(data = predictions, aes(y = predicted_lpqr, x = value), color = "green", size = 1) +
  # QR linear
  geom_line(data = predictions, aes(y = predicted_qr, x = value), size = 1, color = "yellow") +
  geom_abline(intercept = 0, slope = 0, lty = 2) +
  scale_y_continuous(
    name = bquote("Lesion growth (" ~ A[t[1]] - A[t[0]] ~ ")")) +
  scale_x_continuous(name = "Environmental covariate value") +
  facet_wrap(~name,
             labeller = labeller(name = dir_labs),
             scales = "free_x") +
  geom_text(data = qr_models4, aes(x = Inf, y = Inf, label = paste0("Pseudo-R²: ", round(r2, 2))),
            hjust = 1.1, vjust = 1.1, inherit.aes = FALSE) +
  theme(panel.background = element_rect(fill = "#E5E5E5"),
        panel.spacing = unit(1.5, "lines")) 

png(paste0(figure_path, "covars_regr_univariate.png"), width = 6, height = 5, units = 'in', res = 400)
plot(p)
dev.off()

# ============================================================================== -
# Lesion age main effect
# ============================================================================== -

mdat <- readRDS(paste0(data_path, "subset_step2.rds"))

# df for predictions
new_preds <- mdat %>% ungroup() %>% 
  do(., data.frame(
    lesion_age_gdd = seq(round(min(as.numeric(.$lesion_age_gdd))), max(as.numeric(.$lesion_age_gdd)), by = 0.2),
    stringsAsFactors = FALSE)
  )

# linear quantile regression
rqmodel <- rq(diff_area_pp_y_norm_gdd ~ lesion_age_gdd, tau = .5, data = mdat)

# loess
fit <- lprq(mdat$lesion_age_gdd, mdat$diff_area_pp_y_norm_gdd, m = nrow(new_preds), h = 3, tau = .5)
rq_yy <- predict(rqmodel, newdata = new_preds)

## Fit models ================================================================== -

x = "lesion_age_gdd"
y = "diff_area_pp_y_norm_gdd"

fits <- mdat %>% nest() %>% 
  mutate(linear_q = purrr::map(.x = data, .f = linear_quantile, x=x, y=y),
         nls_q_exp = purrr::map(.x = data, .f = nls_quantile_exp, n_samples = 300, x=x, y=y)) %>%  
  tidyr::pivot_longer(cols = linear_q:nls_q_exp, names_to = "type", values_to = "fit")

saveRDS(fits, paste0(data_path, "model_fits.rds"))
fits <- readRDS(paste0(data_path, "model_fits.rds"))

# get pseudo-r2
resid_fit <- residuals(fits$fit[[2]])
rho_tau <- function(u, tau) {
  u * (tau - (u < 0))
}
numerator <- sum(rho_tau(resid_fit, tau=.5))
denominator <- sum(rho_tau(mdat[[y]] - median(mdat[[y]]), tau=.5))
pseudo_R2 <- 1 - numerator / denominator

pdat <- fits %>% 
  mutate(preds = purrr::map(fit, broom::augment, newdata = new_preds))
preds <- pdat %>% dplyr::select(type, preds) %>% unnest(preds) %>% 
  dplyr::filter(type ==  "nls_q_exp")
pd <- pdat %>% 
  dplyr::filter(type ==  "nls_q_exp") %>%  ## only one dataset needed
  dplyr::select(type, data) %>% unnest(data)

lims <- mdat %>%
  summarise(
    y_min = quantile(diff_area_pp_y_norm_gdd, 0.005),
    y_max = quantile(diff_area_pp_y_norm_gdd, 0.995)
  )

p1 <- ggplot() +
  geom_point(data = pd, aes(x = lesion_age_gdd, y = diff_area_pp_y_norm_gdd), alpha = 0.025) +
  # add loess fit
  geom_line(aes(x = fit$xx, y = fit$fv), color = "green") +
  # add linear and non-linear fit
  geom_line(data = preds, aes(x = lesion_age_gdd, y = .fitted), color = "yellow") +
  scale_x_continuous(limits = c(0, 50)) +
  scale_y_continuous(
    name = bquote("Lesion growth (" ~ A[t[1]] - A[t[0]] ~ ")"),
    # limits = c(lims$y_min, lims$y_max),
    limits = c(-0.01, 0.06)
  ) +
  geom_abline(intercept = 0, slope = 0, lty = 2) +
  geom_text(aes(x = Inf, y = Inf, label = paste0("Pseudo-R²: ", round(pseudo_R2, 2))),
            hjust = 1.2, vjust = 2, inherit.aes = FALSE) +
  ggtitle("A") +
  theme(
    legend.position = "None",
    panel.background = element_rect(fill = "#E5E5E5"),
    plot.title = element_text(face = "bold")
  ) + 
  xlab(bquote("Lesion age at" ~t[1]~"(effective)"))

# p <- ggMarginal(p, type = "density", fill = "gray", alpha = 0.5)

png(paste0(figure_path, "05_fits_growth_age_all_density.png"), width = 7, height = 5, units = 'in', res = 400)
plot(p1)
dev.off()

# ============================================================================== -
# Pycn distance
# ============================================================================== -

pd <- readRDS(paste0(data_path, "mdat.rds")) %>% 
  extract_covars_from_nested("design", "genotype_name") %>% 
  mutate(batch_UID = paste(exp_UID, batch, sep = "_"))

ggplot() +
  geom_point(data = pd, aes(x = lag_max_dist, y = diff_area_pp_y_norm_gdd), alpha = 1) +
  geom_smooth(data = pd, aes(x = lag_max_dist, y = diff_area_pp_y_norm_gdd)) +
  geom_smooth(data = pd, method = "lm", aes(x = lag_max_dist, y = diff_area_pp_y_norm_gdd))

# ============================================================================== - 

# Fit a quantile regression model
pdat <- pd %>% 
  mutate(dist = lag_max_dist) %>% 
  filter(!is.na(dist)) 

lims <- pdat %>%
  summarise(
    x_min = quantile(dist, 0.005),
    x_max = quantile(dist, 0.995),
    y_min = quantile(diff_area_pp_y_norm_gdd, 0.005),
    y_max = quantile(diff_area_pp_y_norm_gdd, 0.995)
  )

pdat_sub <- pdat %>% 
  filter(diff_area_pp_y_norm_gdd >= lims$y_min & diff_area_pp_y_norm_gdd <= lims$y_max & dist >= lims$x_min & dist <= lims$x_max)

qr_models <- pdat_sub %>% 
  # sample_n(5000) %>%
  nest() %>% ungroup() %>% 
  mutate(linear_q = purrr::map(data, ~ rq(diff_area_pp_y_norm_gdd ~ dist, data = .)), 
         loess_q = purrr::map(data, ~ lprq(x = .$dist, y = .$diff_area_pp_y_norm_gdd, m = nrow(.)/10, h = 0.2*range(.$dist), tau = .5)),
         new_data = purrr::map(data, ~ data.frame(dist = seq(min(.$dist), max(.$dist), length.out = 1000))))

# Make predictions for both QR and LPQR models
qr_models2 <- qr_models %>% 
  mutate(
    predicted_qr = purrr::map2(linear_q, new_data, ~ predict(.x, newdata = .y)),
    predicted_lpqr = purrr::map2(loess_q, new_data, ~ {
      fit_values <- .x$fv
      predict_values <- approx(.x$xx, fit_values, xout = .y$dist)$y
      return(predict_values)
    })
  )

qr_models3 <- qr_models2 %>% 
  mutate(pseudo_r2 = purrr::map2_dbl(data, linear_q, ~ compute_pseudoR2(.x, response = "diff_area_pp_y_norm_gdd", qr_model = .y)))

predictions <- qr_models2 %>% 
  dplyr::select(new_data, predicted_qr, predicted_lpqr) %>% 
  unnest(cols = c(new_data, predicted_qr, predicted_lpqr))

# New facet label names for direction variable
dir_labs <- c("CV Relative Humidity", "CV Temperature", "Mean Temperature", "Mean Relative Humidity")
names(dir_labs) <- c("cv_interval_rh", "cv_interval_temp", "mean_interval_temp", "mean_interval_rh")

p2 <- ggplot() +
  geom_point(data = pdat_sub, aes(x = dist, y = diff_area_pp_y_norm_gdd), alpha = 0.02) +
  # qr loess
  geom_line(data = predictions, aes(y = predicted_lpqr, x = dist), color = "green") +
  # QR linear
  geom_line(data = predictions, aes(y = predicted_qr, x = dist), color = "yellow") +
  geom_abline(intercept = 0, slope = 0, lty = 2) +
  scale_y_continuous(limits = c(-0.01, 0.06),
                     # limits = c(lims$y_min, lims$y_max),
                     name = bquote("Lesion growth (" ~ A[t[1]] - A[t[0]] ~ ")")) +
  scale_x_continuous(limits = c(lims$x_min, lims$x_max),
                     name = "Maximum distance") +
  geom_text(data = qr_models3, aes(x = Inf, y = Inf, label = paste0("Pseudo-R²: ", round(pseudo_r2, 2))),
            hjust = 1.2, vjust = 2, inherit.aes = FALSE) +
  ggtitle("B") +
  theme(
    panel.background = element_rect(fill = "#E5E5E5"),
    axis.title.y = element_blank(),
    legend.position = "None",
    plot.title = element_text(face = "bold"),
    axis.text.y = element_blank()
  )

p2

png(paste0(figure_path, "dist_diff.png"), width = 5, height = 4, units = 'in', res = 400)
plot(p2)
dev.off()

# SHOW AREA EFFECT ON GROWTH FOR Y ONLY
mdat0 <- readRDS(paste0(data_path, "mdat.rds"))

pdat <- mdat0 %>% 
  dplyr::select(plot_UID, leaf_nr, lesion_nr, 
                area,
                diff_area_pp_y_norm_gdd) %>% 
  pivot_longer(cols = 5:5)

# remove nas and non-finites
pdat <- pdat %>%
  filter(!is.na(value) & !is.na(area) & 
           is.finite(value) & is.finite(area))

facet_limits <- pdat %>%
  group_by(name) %>%
  summarise(
    x_min = quantile(value, 0.005),
    x_max = quantile(value, 0.995),
    y_min = quantile(area, 0.005),
    y_max = quantile(area, 0.995)
  )

pdat_sub <- pdat %>% 
  left_join(facet_limits, by = c("name")) %>%
  filter(value >= x_min & value <= x_max & area >= y_min & area <= y_max)

# Fit a quantile regression model
qr_models <- pdat_sub %>% 
  group_by(name) %>% nest() %>% ungroup() %>% 
  mutate(linear_q = purrr::map(data, ~ rq(value ~ area, data = .)), 
         loess_q = purrr::map(data, ~ lprq(x = .$area, y = .$value, m = nrow(.)/100, h = 10, tau = .5)),
         new_data = purrr::map(data, ~ data.frame(area = seq(min(.$area), max(.$area), length.out = 1000))))

# Make predictions for both QR and LPQR models
qr_models2 <- qr_models %>% 
  mutate(
    predicted_qr = purrr::map2(linear_q, new_data, ~ predict(.x, newdata = .y)),
    predicted_lpqr = purrr::map2(loess_q, new_data, ~ {
      fit_values <- .x$fv
      predict_values <- approx(.x$xx, fit_values, xout = .y$area)$y
      return(predict_values)
    })
  )

qr_models3 <- qr_models2 %>% 
  mutate(pseudo_r2 = purrr::map2_dbl(data, linear_q, ~ compute_pseudoR2(.x, response = "value", qr_model = .y)))

predictions <- qr_models2 %>% 
  dplyr::select(name, new_data, predicted_qr, predicted_lpqr) %>% 
  unnest(cols = c(new_data, predicted_qr, predicted_lpqr))

p0 <- ggplot() +
  geom_point(data = pdat_sub, aes(x = area, y = value), alpha = 0.02) +
  # # qr loess
  geom_line(data = predictions, aes(y = predicted_lpqr, x = area), color = "green") +
  # QR linear
  geom_line(data = predictions, aes(y = predicted_qr, x = area), color = "yellow") +
  geom_abline(intercept = 0, slope = 0, lty = 2) +
  scale_y_continuous(
    limits = c(-0.01, 0.06),
    name = bquote("Lesion growth (" ~ A[t[1]] - A[t[0]] ~ ")")) +
  scale_x_continuous(limits = c(0, 50), 
                     name = bquote("Lesion area (mm" ^2~")")) +
  geom_text(data = qr_models3, aes(x = Inf, y = Inf, label = paste0("Pseudo-R²: ", round(pseudo_r2, 2))),
            hjust = 1.2, vjust = 2, inherit.aes = FALSE) +
  ggtitle("C") +
  theme(panel.background = element_rect(fill = "#E5E5E5"),
        axis.title.y = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.y = element_blank())

# combine lesion age and pycn dist
p <- gridExtra::grid.arrange(p1, p2, p0, nrow = 1, widths = c(1.1, 0.95, 0.95))

png(paste0(figure_path, "age_dist_area_diff.png"), width = 8, height = 4, units = 'in', res = 400)
plot(p)
dev.off()

# ============================================================================== -

ggplot(data = pd) +
  geom_boxplot(aes(x = genotype_name, y = log(lag_max_dist)))

a <- ggplot(data = pd) +
  geom_boxplot(aes(x = batch_UID, y = log(lag_max_dist), fill = as.factor(batch)), alpha = 0.4) +
  scale_fill_manual(
    values = c("1" = col[1], "2" = col[2])) +
  ggtitle("A") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "None", 
        axis.text.x = element_text(angl = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

b <- ggplot(data = pd) +
  geom_boxplot(aes(x = genotype_name, y = log(lag_max_dist), fill = "grey95"), alpha = 0.4) +
  geom_jitter(aes(x = genotype_name, y = log(lag_max_dist), fill = "grey95"), alpha = 0.025) +
  ggtitle("B") +
  scale_fill_manual(
    values = c("1" = col[1], "2" = col[2])) +
  xlab("Host cultivar") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "None",
        axis.text.x = element_text(angl = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

p <- gridExtra::grid.arrange(a, b, nrow = 1, widths = c(1, 2))

png(paste0(figure_path, "pycn_patterns.png"), width = 8, height = 4, units = 'in', res = 400)
plot(p)
dev.off()

cordat <- pd %>% 
  dplyr::select(lag_max_dist, lag_variance_dist, lag_mean_dist, lag_median_dist)
M<-cor(cordat, use = "pairwise.complete.obs")

corrplot(M, addCoef.col = "black", tl.col="black", tl.srt=45,
         method = "color", type="upper", diag=FALSE, 
         number.cex = 1, tl.cex= 1)

mod <- lm(log(lag_max_dist) ~ genotype_name+batch_UID, data = pd)
summary(mod)
hist(residuals(mod), breaks = 30, main = "Histogram of Residuals")

anova(mod)
plot(mod)

# ============================================================================== -
# Asymptomatic
# ============================================================================== -

pd <- readRDS(paste0(data_path, "mdat.rds")) %>% 
  extract_covars_from_nested("design", c("harvest_year", "genotype_name")) %>% 
  mutate(batch_UID = paste(exp_UID, batch, sep = "_")) %>% 
  rename(year = "harvest_year")

ggplot() +
  geom_jitter(data = pd, aes(x = lag_asymptomatic_duration_hd, y = diff_area_pp_y_norm_gdd), alpha = 0.05) +
  geom_smooth(data = pd, aes(x = lag_asymptomatic_duration_hd, y = diff_area_pp_y_norm_gdd)) +
  geom_abline(intercept = 0, slope = 0, lty = 2) +
  scale_y_continuous(limits = c(-0.02, 0.1))

pdat_sub <- pd %>% 
  dplyr::filter(!is.na(lag_asymptomatic_duration_hd)) %>% 
  # sample_n(5000) %>%
  mutate(lag_asymptomatic_duration_hd = as.numeric(lag_asymptomatic_duration_hd))

# Fit a quantile regression model
qr_models <- pdat_sub %>% group_by(year, batch) %>% 
  nest() %>% ungroup() %>% 
  mutate(linear_q = purrr::map(data, ~ rq(diff_area_pp_y_norm_gdd ~ lag_asymptomatic_duration_hd, data = .)), 
         loess_q = purrr::map(data, ~ lprq(x = .$lag_asymptomatic_duration_hd, y = .$diff_area_pp_y_norm_gdd, m = nrow(.)/10, h = 1, tau = .5)),
         new_data = purrr::map(data, ~ data.frame(lag_asymptomatic_duration_hd = seq(min(.$lag_asymptomatic_duration_hd), max(.$lag_asymptomatic_duration_hd), length.out = 1000))))

# Make predictions for both QR and LPQR models
qr_models2 <- qr_models %>% 
  mutate(
    predicted_qr = purrr::map2(linear_q, new_data, ~ predict(.x, newdata = .y)),
    predicted_lpqr = purrr::map2(loess_q, new_data, ~ {
      fit_values <- .x$fv
      predict_values <- approx(.x$xx, fit_values, xout = .y$lag_asymptomatic_duration_hd)$y
      return(predict_values)
    })
  )

qr_models3 <- qr_models2 %>% 
  mutate(pseudo_r2 = purrr::map2_dbl(data, linear_q, ~ compute_pseudoR2(.x, response = "diff_area_pp_y_norm_gdd", qr_model = .y)))

predictions <- qr_models2 %>% 
  dplyr::select(batch, year, new_data, predicted_qr, predicted_lpqr) %>% 
  unnest(cols = c(new_data, predicted_qr, predicted_lpqr))

p <- ggplot() +
  geom_jitter(data = pdat_sub, aes(x = lag_asymptomatic_duration_hd, y = diff_area_pp_y_norm_gdd), alpha = 0.02) +
  # # qr loess
  # geom_line(data = predictions, aes(y = predicted_lpqr, x = lag_asymptomatic_duration_hd), color = "green", size = 1) +
  # QR linear
  geom_line(data = predictions, aes(y = predicted_qr, x = lag_asymptomatic_duration_hd), size = 1, color = "yellow") +
  geom_abline(intercept = 0, slope = 0, lty = 2) +
  scale_y_continuous(
    limits = c(-0.02, 0.05), 
    name = bquote("Lesion growth (" ~ A[t[1]] - A[t[0]] ~ ")")) +
  scale_x_continuous(name = bquote("Lesion area (mm" ^2~")")) +
  facet_grid(year~batch, scales = "free_x") +
  geom_text(data = qr_models3, aes(x = Inf, y = Inf, label = paste0("Pseudo-R²: ", round(pseudo_r2, 2))),
            hjust = 1.1, vjust = 1.1, inherit.aes = FALSE)

png(paste0(figure_path, "leaf_asymptomatic.png"), width = 6, height = 6, units = 'in', res = 400)
plot(p)
dev.off()

# >> to leaves ----

# mdat_leaf <- mdat %>% 
#   group_by(leaf_id) %>% nest() %>% 
#   mutate(n = purrr::map_dbl(data, nrow)) %>% 
#   filter(n > 20)
# 
# num_ticks <- n_groups(mdat_leaf)
# pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) elapsed :elapsed eta :eta",
#                        total = num_ticks)
# 
# # processing serially (slow!)
# fits <- mdat_leaf%>% 
#   mutate(
#     linear_q = purrr::map(.x = data, .f = linear_quantile),
#     nls_q_exp = purrr::map(.x = data, .f = possibly(~ nls_quantile_exp(.x, n_samples = 500), otherwise = NA))) %>%
#   tidyr::pivot_longer(cols = linear_q:nls_q_exp, names_to = "type", values_to = "fit")
# 
# saveRDS(fits, paste0(data_path, "model_fits_leaves.rds"))
# processing in parallel: see 'lesion_level_dynamic_par.R' and 'utils.R'

# ============================================================================== -

# load from processed in parallel
preds <- readRDS(paste0(data_path, "model_preds_leaves.rds"))
pars <- readRDS(paste0(data_path, "model_pars_leaves.rds"))

pd <- readRDS(paste0(data_path, "mdat.rds"))

p <- ggplot() +
  geom_line(data = preds, aes(x = lesion_age_gdd, y = .fitted, group = type,  col = type)) +
  geom_point(data = pd, aes(x = lesion_age_gdd, y = diff_area_norm_gdd_peri_y), alpha = 0.4) +
  scale_y_continuous(
    name = bquote("Normalized growth: " ~ A[t[x]] - A[t[x-1]]),
    limits = c(-0.01, 0.025),
    # breaks = seq(0, 5, 1)
  ) +
  facet_wrap(~leaf_id) +
  theme_bw( base_family = "Helvetica") +
  theme(panel.grid = element_blank()) +
  xlab(bquote("lesion age at" ~t[2]~" (?C days)")) + ylab(bquote("Lesion relative growth: " ~ A[t[x]] / A[t[x-1]])) 
png(paste0(figure_path, "01_fits_growth_age_leaves.png"), width = 25, height = 25, units = 'in', res = 400)
plot(p)
dev.off()

# parameters - leaves
gen <- lesion_data %>% dplyr::select(leaf_id, batch, genotype_name) %>% unique()
pars_aug <- left_join(pars, gen, by = "leaf_id") %>% 
  tidyr::pivot_longer(., cols = y0_fitted:slope_fitted, names_to = "parameter") %>% 
  filter(!is.na(value))
# reorder 
pars_aug$parameter <- factor(pars_aug$parameter, 
                             levels = c("y0_fitted", "yf_fitted", "alpha_fitted",
                                        "intercept_fitted", "slope_fitted"))
p <- ggplot(data = pars_aug) +
  geom_boxplot(aes(x = genotype_name, y = value), outlier.shape = NA) +
  geom_jitter(aes(x = genotype_name, y = value), alpha = 0.2) +
  theme_bw( base_family = "Helvetica") +
  facet_wrap(~parameter, scales = "free")+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))
png(paste0(figure_path, "07_parameters_genotypes.png"), width = 7, height = 5, units = 'in', res = 400)
plot(p)
dev.off()

# parameters - batches
pars_aug$batch <- factor(pars_aug$batch, levels = c("1", "2"))
p <- ggplot(data = pars_aug) +
  geom_boxplot(aes(x = batch, y = value), outlier.shape = NA) +
  geom_jitter(aes(x = batch, y = value), alpha = 0.2) +
  theme_bw( base_family = "Helvetica") +
  facet_wrap(~parameter, scales = "free")+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))
png(paste0(figure_path, "parameters_batches.png"), width = 7, height = 5, units = 'in', res = 400)
plot(p)
dev.off()

# ============================================================================== -

# >> to batches ----

mdat <- readRDS(paste0(data_path, "mdat.rds"))
# extract only data for genotypes in both batches
mdat <- mdat %>% 
  extract_covars_from_nested(., from="design", "genotype_name") %>% 
  filter(genotype_name %in% c("FORNO", "BORNEO"))

mdat_batch <- mdat %>% 
  group_by(batch) %>% nest() %>% 
  mutate(n = purrr::map_dbl(data, nrow))

# processing serially
x = "lesion_age_gdd"
y = "diff_area_norm_gdd_peri_y"
fits <- mdat_batch%>% 
  mutate(
    linear_q = purrr::map(.x = data, .f = linear_quantile, x=x, y=y),
    nls_q_exp = purrr::map(.x = data, .f = possibly(~ nls_quantile_exp(.x, n_samples = 100, x=x, y=y), otherwise = NA))) %>%
  tidyr::pivot_longer(cols = linear_q:nls_q_exp, names_to = "type", values_to = "fit")

# saveRDS(fits, paste0(data_path, "model_fits_batches.rds"))
# fits <- readRDS(paste0(data_path, "model_fits_batches.rds"))

pars <- fits %>% 
  filter(type == "nls_q_exp") %>% 
  mutate(pars = purrr::map(fit, get_qr_pars)) %>% unnest(pars)

pdat <- fits %>%
  filter(type == "nls_q_exp") %>% 
  mutate(preds = purrr::map(fit, broom::augment, newdata = new_preds))

preds <- pdat %>% dplyr::select(batch, type, preds) %>% unnest(preds)
pd <- pdat %>% 
  filter(type == "nls_q_exp") %>%  # only one copy of the data needed
  dplyr::select(batch, data) %>% unnest(data)

col = pal_jco()(2)
p <- ggplot() +
  geom_point(data = pd, aes(x = lesion_age_gdd, y = diff_area_norm_gdd_peri_y, col = as.factor(batch)), alpha = 0.025) +
  geom_line(data = preds, aes(x = lesion_age_gdd, y = .fitted, group = interaction(type, batch),  col = as.factor(batch))) +
  scale_y_continuous(
    name = bquote("GDD-Normalized growth " ~ pp[y]~":" ~ A[t[x]] - A[t[x-1]]),
    limits = c(-0.01, 0.025),
  ) +
  scale_color_manual("batch", values =  col) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw(base_family = "Helvetica") +
  theme(panel.grid = element_blank()) +
  xlab(bquote("lesion age at" ~t[2]~" (?C days)")) + ylab(bquote("Lesion relative growth: " ~ A[t[x]] / A[t[x-1]])) 
png(paste0(figure_path, "06_fits_growth_age_batches.png"), width = 5, height = 3.5, units = 'in', res = 400)
plot(p)
dev.off()

# ============================================================================== -
# 6) heading date ----
# ============================================================================== -

# see Z:/Public/Jonas/004_ESWW007/RScripts/heading.R
hd <- read_csv("Z:/Public/Jonas/Data/ESWW007/RefData/BBCH/hd.csv") %>% 
  mutate(hd_das = yday(hd))
mdat <- readRDS(paste0(data_path, "mdat.rds"))

pd <- left_join(mdat, hd, by = "plot_UID") %>% 
  dplyr::filter(diff_area_norm_gdd_peri_y > -0.2) %>% 
  dplyr::filter(!is.na(hd_das))

fit <- lprq(pd$hd_das, pd$diff_area_norm_gdd_peri_y, m = nrow(new_preds), h = 25, tau = .5)
rq_yy <- predict(rqmodel, newdata = new_preds)

fits_lin <- rq(diff_area_norm_gdd_peri_y ~ hd_das, data = pd)

new_preds_hd <- pd %>% ungroup() %>% 
  do(., data.frame(
    hd_das = seq(round(min(as.numeric(.$hd_das))), max(as.numeric(.$hd_das)), by = 0.5),
    stringsAsFactors = FALSE)
  )
preds <- broom::augment(fits_lin, new_data = new_preds_hd)

p <- base_plot_batches +
  geom_point(data = pd, aes(x = hd_das, y = diff_area_norm_gdd_peri_y), alpha = 0.025) +
  # add QR loess fit
  geom_line(aes(x = fit$xx, y = fit$fv), color = "green", size = 1) +
  # add OLS fit
  geom_smooth(method = "lm", data = pd, aes(x = hd_das, y = diff_area_norm_gdd_peri_y)) +
  # add QR linear fit
  geom_line(data = preds, aes(x = hd_das, y = .fitted), size = 1, lty = 2, col = "red") +
  ggpubr::stat_cor(data = pd, aes(x = hd_das, y = diff_area_norm_gdd_peri_y, label = after_stat(rr.label)), color = "blue", geom = "label") +
  scale_y_continuous(
    name = bquote("Normalized growth: " ~ A[t[x]] - A[t[x-1]]),
    limits = c(-0.01, 0.025),
  ) +
  geom_smooth(method = "lm") +
  xlab("Host cultivar heading time (Days after sowing)")

png(paste0(figure_path, "07_lesion_growth_hd.png"), width = 5, height = 3.5, units = 'in', res = 400)
plot(p)
dev.off()

# ============================================================================== -
# 7) Fit models (area) ---- 
# ============================================================================== -

mdat0 <- readRDS(paste0(data_path, "mdat.rds"))

# df for predictions
new_preds <- mdat0 %>% ungroup() %>% 
  do(., data.frame(
    area = seq(round(min(as.numeric(.$area))), max(as.numeric(.$area)), by = 1),
    stringsAsFactors = FALSE)
  )

pdat <- mdat0 %>% 
  dplyr::select(plot_UID, leaf_nr, lesion_nr, 
                area, diff_area_norm_gdd, 
                diff_area_norm_gdd_peri_x,
                diff_area_norm_gdd_peri_y,
                diff_area_norm_gdd_peri_xy)

# df for predictions
new_preds <- mdat %>% ungroup() %>% 
  do(., data.frame(
    area = seq(round(min(as.numeric(.$area))), max(as.numeric(.$area)), by = 1),
    stringsAsFactors = FALSE)
  )

# linear quantile regression
rqmodel <- rq(diff_area_norm_gdd_peri_xy ~ area, tau = .5, data = mdat)

# loess
fit <- lprq(mdat$area, mdat$diff_area_norm_gdd_peri_xy, m = nrow(new_preds), h = 20, tau = .5)
rq_yy <- predict(rqmodel, newdata = new_preds)

fits <- mdat %>% nest() %>% 
  mutate(linear_q = purrr::map(.x = data, .f = linear_quantile, x = "area", y = "diff_area_norm_gdd_peri_xy")) %>% 
  tidyr::pivot_longer(cols = linear_q:linear_q, names_to = "type", values_to = "fit")

d <- fits %>% 
  mutate(preds = purrr::map(fit, broom::augment, newdata = new_preds))
preds <- d %>% dplyr::select(type, preds) %>% unnest(preds)

p <- base_plot_batches+
  geom_point(data = pdat, aes(x = area, y = diff_area_norm_gdd_peri_xy), alpha = 0.1) +
  geom_smooth(data = pdat, aes(x = area, y = diff_area_norm_gdd_peri_xy), se = F) +
  geom_smooth(method = "lm", data = pdat,  aes(x = area, y = diff_area_norm_gdd_peri_xy), se = F, color = "red") +
  # add loess fit
  geom_line(aes(x = fit$xx, y = fit$fv), color = "green") +
  # add linear fit
  geom_line(data = preds, aes(x = area, y = .fitted, color = type)) +
  ggpubr::stat_cor(data = pdat, aes(x = area, y = diff_area_norm_gdd_peri_xy, label = after_stat(rr.label)), color = "red", geom = "label") +
  xlab("Lesion area") +   
  scale_y_continuous(name = bquote("GDD-normalized growth: " ~ A[t[x]] - A[t[x-1]])) +
  scale_x_continuous(limits = c(0, 150), 
                     name = bquote("Lesion area (mm" ^2~")"))

png(paste0(figure_path, "lesion_growth_area_fits.png"), width = 7, height = 5, units = 'in', res = 400)
plot(p)
dev.off()

# linear seems OK. 

# ============================================================================== -
# Trying random regression ----
# ============================================================================== -

mdat <- readRDS(paste0(data_path, "mdat.rds"))

data <- mdat %>% 
  extract_covars_from_nested(., from="design", "genotype_name")

ggplot(data) +
  geom_histogram(aes(x = diff_area_pp_y_norm_chr))

# # plots
# ggplot(data) +
#   geom_boxplot(aes(x = plot_UID, y = diff_area_pp_y_norm_chr))
# 
# # leaves
# ggplot(data) +
#   geom_boxplot(aes(x = leaf_UID, y = diff_area_pp_y_norm_chr))
# 
# # cultivars
# ggplot(data) +
#   geom_boxplot(aes(x = genotype_name, y = diff_area_pp_y_norm_chr))

moddat <- data %>% 
  dplyr::select(exp_UID, batch, genotype_name, plot_UID, leaf_UID, lesion_UID, diff_area_pp_y_norm_chr, mean_interval_temp, mean_interval_rh, lesion_age) %>% 
  mutate(batch_UID = paste(exp_UID, batch, sep = "_")) %>% 
  dplyr::relocate(batch_UID, .after = exp_UID) %>% 
  dplyr::select(-batch) %>% 
  mutate_at(C(1:6), factor)

ggplot(moddat) +
  geom_boxplot(aes(x = batch_UID, y = diff_area_pp_y_norm_chr))

ggplot(moddat, aes(x = mean_interval_temp, y = diff_area_pp_y_norm_chr)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(genotype_name~batch_UID)

library(asreml)
asreml.options(workspace="8192mb")

ll <- asreml.options()
ll$workspace

model <- asreml(fixed = diff_area_pp_y_norm_chr~mean_interval_temp + genotype_name + mean_interval_temp:genotype_name + exp_UID + batch_UID + lesion_age,
                random = ~ plot_UID/leaf_UID, 
                data = moddat, maxit=20)
plot(model)
wald.asreml(model)

model_linear <- lme(
  diff_area_pp_y_norm_chr ~ genotype_name*mean_interval_rh*mean_interval_temp + batch + lesion_age_gdd, 
  data = data, 
  random = ~ 1 | plot_UID / leaf_nr / lesion_nr
)

res <- residuals(model_linear)
shapiro.test(res[seq(1, length(res), 10)])

summary(model_linear)
anova(model_linear)
plot(model_linear)

hist(residuals(model_linear), breaks = 30, main = "Histogram of Residuals")
qqnorm(residuals(model_linear))  # Q-Q plot
qqline(residuals(model_linear), col = "red")
plot(data$genotype_name, residuals(model_linear), main = "Residuals vs. Genotype")
plot(data$mean_interval_temp, residuals(model_linear), main = "Residuals vs. Temperature")
plot(data$lesion_age_gdd, residuals(model_linear), main = "Residuals vs. Lesion Age")

# transformation
# ensure all values are > 0
data$response <- data$diff_area_pp_y_norm_chr + abs(min(data$diff_area_pp_y_norm_chr)) + 0.0001
min(data$response)

model_linear <- lme(
  log(response) ~ genotype_name*mean_interval_rh*mean_interval_temp + batch + lesion_age_gdd, 
  data = data, 
  random = ~ 1 | plot_UID / leaf_nr / lesion_nr
)

hist(residuals(model_linear), breaks = 30, main = "Histogram of Residuals")
qqnorm(residuals(model_linear))  # Q-Q plot
qqline(residuals(model_linear), col = "red")



plot(model)
summary(model)$varcomp
summary(model, coef=TRUE)$coef.fixed
summary(model, coef=TRUE)$coef.random
wald.asreml(model)

anova(model)


# as in Adhikari et al.
model_linear <- lme(
  diff_area_pp_y_norm_chr ~ genotype_name*mean_interval_rh*mean_interval_temp + batch, 
  # diff_area_pp_y_norm_chr ~ genotype_name*mean_interval_rh*mean_interval_temp + lesion_age_gdd + batch + area, 
  data = data, 
  random = ~ 1 | plot_UID / leaf_nr / lesion_nr
)

VarCorr(model_linear)

plot(fitted(model_linear), residuals(model_linear),
     xlab = "Fitted Values", ylab = "Residuals", 
     main = "Residuals vs. Fitted Values")
abline(h = 0, col = "red")

qqnorm(residuals(model_linear))  # Q-Q plot
qqline(residuals(model_linear), col = "red")

shapiro.test(residuals(model_linear))

hist(residuals(model_linear), breaks = 30, main = "Histogram of Residuals")

plot(data$genotype_name, residuals(model_linear), main = "Residuals vs. Genotype")
plot(data$mean_interval_temp, residuals(model_linear), main = "Residuals vs. Temperature")
plot(data$lesion_age_gdd, residuals(model_linear), main = "Residuals vs. Lesion Age")

summary(model_linear)
anova(model_linear)
summary(model_linear)$tTable
ranef(model_linear)

emmeans(model_linear, ~ lesion_age_gdd | genotype_name)
plot(emmeans(model_linear, ~ genotype_name))

# install.packages("performance")  # If not installed
library(performance)

# Compute R? for mixed model
r2_results <- r2_nakagawa(model_linear)
print(r2_results)




model <- lmer(
  diff_area_ppy ~ genotype_name*lesion_age_gdd + area + 
    (1 | batch) + (1 | plot_UID / leaf_nr / lesion_nr),
  data = data
)

model <- nlmer(
  diff_area_pp_y ~ SSasymp(lesion_age_gdd, yf, alpha, y0) + 
    genotype_name + diff_gdd + area ~ 1 + (1 | batch) + (1 | plot_UID / leaf_nr / lesion_nr),
  data = data,
  start = list(yf = 0.0, alpha = 0.1, y0 = 1) 
)


summary(model)

summary(model)
ranef(model)

model_reduced <- lmer(
  diff_area_pp_y ~ diff_gdd + lesion_age_gdd + area +
    (1 | batch) + (1 | plot_UID / leaf_nr / lesion_nr),
  data = data,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
)

anova(model_reduced, model)

# ============================================================================== -

# get cleaned data (without intervals)
base_plot <- ggplot() +
  xlab("Lesion age (?C days)") +
  scale_color_manual(values =  col) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(panel.background = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(),
        legend.title = element_blank())

data <- readRDS(paste0(data_path, "mdat.rds")) %>% 
  extract_covars_from_nested(., from="design", "genotype_name")

# cross batches for common genotypes
subset <- data %>% 
  filter(genotype_name %in% c("FORNO", "BORNEO"))
col = pal_jco()(2)
base_plot +
  geom_point(data = subset, aes(x = lesion_age_gdd, y = diff_area_norm_gdd_peri_y, color = as.factor(batch)),
             alpha = 0.1) +
  geom_smooth(data = subset, method = "lm", aes(x = lesion_age_gdd, y = diff_area_norm_gdd_peri_y, color = as.factor(batch)), se = F) + 
  geom_smooth(data = subset, aes(x = lesion_age_gdd, y = diff_area_norm_gdd_peri_y, color = as.factor(batch)), se = F)

# as in Adhikari et al.
model_linear <- lme(
  diff_area_ppy ~ genotype_name*diff_gdd*lesion_age_gdd + batch + area, 
  data = data, 
  random = ~ 1 | plot_UID / leaf_UID / lesion_UID
)

summary(model_linear)
anova(model_linear)
summary(model_linear)$tTable
ranef(model_linear)

emmeans(model_linear, ~ lesion_age_gdd | genotype_name)
plot(emmeans(model_linear, ~ genotype_name))

# ============================================================================== -
# data for statslab ----
# ============================================================================== -

# Function to get range for numeric columns
get_range <- function(x) {
  if (is.numeric(x)) {
    paste0(round(min(x, na.rm = TRUE), 2), " - ", round(max(x, na.rm = TRUE), 2))
  } else {
    NA  # No range for non-numeric columns
  }
}

# Function to get mean for numeric columns
get_mean <- function(x) {
  if (is.numeric(x))  {
    round(mean(x, na.rm = TRUE), 2)
  } else {
    NA  # No mean for non-numeric columns
  }
}

data <- readRDS(paste0(data_path, "subset_step2.rds"))

data_statslab <- data %>% 
  dplyr::select(plot_UID:frac_xy_perimeter_f, 
                la_healthy_f:placl, 
                lesion_age:lesion_age_gdd, 
                lag_area:lag_frac_xy_perimeter_f, 
                lag_la_healthy_f:lag_placl,
                lag_lesion_age:lag_symptomatic_duration,
                lag_lesion_age_gdd) %>% 
  extract_covars_from_nested("design", "harvest_year") %>% 
  relocate(harvest_year, .after = batch) %>% 
  dplyr::select(-file_id)

saveRDS(data_statslab, "Z:/Public/Jonas/Share/statslab/data.rds")

# Generate metadata structure
meta_df <- tibble(
  column_name = names(data_statslab),
  description = c("plot unique identifier",
                  "leaf number, randomly assigned to selected leaves",
                  "lesion number, assigned according to order of appearance on leaves",
                  "lesion unique identifier, composed of plot_UID, leaf_nr, and lesion_nr)",
                  "design information for the experimental plot",
                  "experiment uniqe identifier (redundant with harvest_year)",
                  "leaf unique identifier, composed of plot_UID and leaf_nr",
                  "datetime of image acquisition at t",
                  "datetime of image acquisition at t-1",
                  "date of image acquisition at t",
                  "time of image acquisition at t",
                  "datetime of precending image acquisition for the leaf",
                  "datetime of first image acquisition for the leaf",
                  "datetime of last image acquisition for the leaf",
                  "datetime of first measurement of the lesion",
                  "datetime of last measurement of the lesion",
                  "datetime of first lesion detected on the leaf",
                  "measurement batch, 1 = penultimate leaves, 2 = flag leaves",
                  "year in which measurements were made",
                  "unique identifier for a covariate course",
                  "hourly covariate data for the interval",
                  "duration of the interval in chronological time",
                  "duration of the interval in thermal time (growing degree days)",
                  "difference in lesion area measured for the interval, absolute, in mm^2",
                  "difference in lesion area measured for the interval, divided by the length of the xy perimeter",
                  "difference in lesion area measured for the interval, divided by the length of the x perimeter",
                  "difference in lesion area measured for the interval, divided by the length of the y perimeter",
                  "difference in lesion area measured for the interval, divided by the length of the y perimeter, per hour",
                  "difference in lesion area measured for the interval, divided by the length of the y perimeter, per growing degree day",
                  "difference in lesion area measured for the interval, divided by the length of the x perimeter, per hour",
                  "difference in lesion area measured for the interval, divided by the length of the x perimeter, per growing degree day",
                  "difference in lesion area measured for the interval, divided by the length of the xy perimeter, per hour",
                  "difference in lesion area measured for the interval, divided by the length of the xy perimeter, per growing degree day",
                  "mean temperature in the interval",
                  "mean relative humidity in the interval",
                  "coefficient of variation of temperature in the interval",
                  "coefficient of variation of relative humidity in the interval",
                  "differene in lesion area measured for the interval, devided by chronological time",
                  "differene in lesion area measured for the interval, devided by thermal time",
                  "area of the lesion at t",
                  "solidity of the lesion at t",
                  "maximum width of the lesion at t",
                  "maximum height of the lesion at t",
                  "length of the lesion perimeter at t",
                  "length of the free lesion perimeter at t",
                  "length of the occluded lesion perimeter at t",
                  "length of the y perimeter at t",
                  "length of the x perimeter at t",
                  "length of the xy perimeter at t",
                  "length of the occluded y perimter at t",
                  "length of the occluded x perimeter at t",
                  "length of the occluded xy perimeter at t",
                  "length of the free y perimter at t",
                  "length of the free x perimeter at t",
                  "length of the free xy perimeter at t",
                  "fraction of free x perimeter at t",
                  "fraction of free y perimeter at t",
                  "fraction of free xy perimeter at t",
                  "healthy leaf area at t",
                  "percent leaf area covered by lesions at t",
                  "lesion age in chronological time at t",
                  "leaf age in chronological time at t",
                  "leaf age relative to heading date",
                  "duration of the symptomatic phase of the leaf at t",
                  "duration of the asymptomatic phase of the leaf",
                  "duration of the asymptomatic phase of the leaf since heading",
                  "lesion age in thermal time at t",
                  "area of the lesion at t-1",
                  "solidity of the lesion at t-1",
                  "maximum width of the lesion at t-1",
                  "maximum height of the lesion at t-1",
                  "length of the lesion perimeter at t-1",
                  "length of the free lesion perimeter at t-1",
                  "length of the occluded lesion perimeter at t-1",
                  "length of the y perimeter at t-1",
                  "length of the x perimeter at t-1",
                  "length of the xy perimeter at t-1",
                  "length of the occluded y perimter at t-1",
                  "length of the occluded x perimeter at t-1",
                  "length of the occluded xy perimeter at t-1",
                  "length of the free y perimter at t-1",
                  "length of the free x perimeter at t-1",
                  "length of the free xy perimeter at t-1",
                  "fraction of free x perimeter at t-1",
                  "fraction of free y perimeter at t-1",
                  "fraction of free xy perimeter at t-1",
                  "healthy leaf area at t-1",
                  "percent leaf area covered by lesions at t-1",
                  "lesion age in chronological time at t-1",
                  "leaf age in chronological time at t-1",
                  "leaf age relative to heading date at t-1",
                  "duration of the symptomatic phase of the leaf at t-1",
                  "lesion age in thermal time at t-1"),
  data_type = sapply(data_statslab, class),  # Extract data type
  range = map_chr(data_statslab, get_range),  # Extract range
  mean = map_dbl(data_statslab, get_mean)  # Extract mean
)

res <- meta_df %>% 
  mutate(data_type = map_chr(data_type, ~ .[[1]]))

write.csv(res, "Z:/Public/Jonas/Share/statslab/meta.csv", row.names = F)

DDD <- readRDS("Z:/Public/Jonas/Share/statslab/data.rds")

# diff_area_pp_xy = diff_area / lag_xy_perimeter_f,

A <- data$diff_area_pp_xy 
B <- data$diff_area/data$lag_xy_perimeter_f
A == B

A <- data$diff_area_pp_xy_norm_chr 
B <- data$diff_area_pp_xy/as.numeric(data$diff_time)
A == B

ggplot(data_statslab) +
  geom_point(aes(x = diff_area_pp_xy_norm_chr, y=(diff_area_norm_chr/lag_xy_perimeter_f)))
