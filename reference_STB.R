
#=============================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2025  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

#=============================================================================== -

rm(list = ls())
.libPaths("T:/R4Userlibs")

library(SpATS)
library(tidyverse)
library(car)
library(asreml)
library(boot)

data <- read_csv("Z:/Public/Jonas/011_STB_leaf_tracking/data/stb_reference_data.csv")
design <- read_csv("Z:/Public/Jonas/011_STB_leaf_tracking/data/stb_design.csv")

# split data according to assessment time point  
subsets <- split(data, interaction(data$exp_UID, data$time_point))[1:5]

# ============================================================================== -

# >> incidence ----

geno_pred <- list()  # to collect results
for(i in 1:length(subsets)){
  
  # get subset
  dataset <- subsets[[i]] %>% 
    dplyr::select(plot_UID, exp_UID, time_point, leaf, scorer, incidence)
  exp <- unique(dataset$exp_UID)
  des <- design %>% filter(exp_UID == exp)
  
  # join data and design
  # convert all to factors
  # remove check plots
  data <- full_join(dataset, des) %>% 
    mutate(time_point = as.factor(time_point),
           scorer = as.factor(scorer),
           leaf = as.factor(leaf),
           plot_UID = as.factor(plot_UID),
           genotype_name = as.factor(genotype_name),
           block = as.factor(block),
           treatment = as.factor(treatment),
           Xf = as.factor(Xf),
           Yf = as.factor(Yf)) %>% 
    filter(complete.cases(.)) %>% 
    filter(treatment != "Check")
  
  # check which factors have more than one level and should therefore be included as fixed effects
  # genotype, Xf, Yf, and plot_UID are always included 
  df <- data %>%
    dplyr::select(plot_UID, leaf, incidence, Xf, Yf, treatment, genotype_name, leaf, scorer, block, range_lot, row_lot) %>% 
    filter(complete.cases(.)) %>% droplevels()
  check_factors <- function(data) {
    factor_cols <- c("treatment", "leaf", "scorer", "block")
    factor_levels <- sapply(data[, factor_cols, drop = FALSE], function(x) length(levels(x)) > 1)
    return(names(factor_levels)[factor_levels])
  }
  valid_factors <- check_factors(df)
  
  # perform spatial correction
  response <-"incidence"
  data$incidence <- car::logit(data$incidence, adjust=0.01, percents = F)
  random <- "~ Xf + Yf + plot_UID" 
  fixed <- paste("~", paste(valid_factors, collapse = " + "))
  genotype = "genotype_name"
  genotype.as.random <- FALSE  # to get BLUEs
  spats <- SpATS(response = response,
                 random = random, 
                 fixed = fixed,
                 spatial = ~PSANOVA(range_lot, row_lot, nseg = c(5,5), nest.div=c(2,2)),
                 genotype = genotype, 
                 genotype.as.random = genotype.as.random,
                 data = data,
                 control = list(maxit = 100, tolerance = 1e-03, monitoring = 0))
  
  # # heritability
  # print(getHeritability(spats))
  
  # get BLUEs and standard errors
  geno_pred[[i]] <- predict(spats, which = "genotype_name")[c("genotype_name", "predicted.values", "standard.errors")]
  
}

# add time point and year
t <- names(subsets)
for(i in 1:length(geno_pred)){
  geno_pred[[i]]$tp <- strsplit(t[i], "\\.") %>% lapply("[[", 2) %>% unlist()
  geno_pred[[i]]$exp <- strsplit(t[i], "\\.") %>% lapply("[[", 1) %>% unlist()
}
Blues <- bind_rows(geno_pred)

# Second stage: Estimate heritability across time points, and get an overall estimate of the genotype value
Blues$TP <- paste(Blues$exp, Blues$tp, sep = "_")
Blues <- Blues %>%
  mutate(TP = as.factor(TP),
         genotype_name = as.factor(genotype_name))

### With asreml can fit TP:genotype interactions although design matrix is not complete
weights_ <- 1/(Blues$standard.errors^2)
mod1 <- asreml(data=Blues,
       fixed = predicted.values ~ TP,
       random = ~ genotype_name + TP:genotype_name,
       weights = weights_,
       trace=TRUE)
mod1 <- update(mod1)
summary(mod1)$varcomp
vg <- summary(mod1)$varcomp["genotype_name", "component"]
ve <- summary(mod1)$varcomp["units!R", "component"]
heritability <- vg / (vg + ve)

# get genotypic BLUP
#genoblups <- ranef(mod_2stage)$genotype_name
genoblups <- predict(mod1, classify = "genotype_name")$pvals
genoblups$predicted.value <- inv.logit(genoblups$predicted.value)

write.csv(genoblups, "Z:/Public/Jonas/011_STB_leaf_tracking/data/inc_blups.csv")

# ============================================================================== -

# >> conditional severity ----
geno_pred <- list()
for(i in 1:length(subsets)){
  
  # get subset
  dataset <- subsets[[i]] %>% 
    dplyr::select(plot_UID, exp_UID, time_point, leaf, placl_mean, placl_var)
  exp <- unique(dataset$exp_UID)
  des <- design %>% filter(exp_UID == exp)
  
  # join data and design
  # convert all to factors
  # remove check plots
  data <- full_join(dataset, des) %>% 
    mutate(time_point = as.factor(time_point),
           leaf = as.factor(leaf),
           plot_UID = as.factor(plot_UID),
           genotype_name = as.factor(genotype_name),
           treatment = as.factor(treatment),
           block = as.factor(block),
           Xf = as.factor(Xf),
           Yf = as.factor(Yf)) %>% 
    filter(complete.cases(.)) %>% 
    filter(treatment != "Check")

   data$placl_mean <- car::logit(data$placl_mean, adjust=0.01, percents = F)
  
  # check which factors have more than one level and should therefore be included as fixed effects
  # genotype, Xf, Yf, and plot_UID are always included 
  df <- data %>%
    dplyr::select(plot_UID, leaf, placl_mean, placl_var, Xf, Yf, treatment, genotype_name, leaf, block, range_lot, row_lot) %>% 
    filter(complete.cases(.)) %>% droplevels()
  check_factors <- function(data) {
    factor_cols <- c("treatment", "leaf", "block")
    factor_levels <- sapply(data[, factor_cols, drop = FALSE], function(x) length(levels(x)) > 1)
    return(names(factor_levels)[factor_levels])
  }
  valid_factors <- check_factors(df)
  
  # perform spatial correction
  # If multiple levels of leaf ("L1", "L2", "L3), include plot_UID, 
  # since incidence estimates per leaf layer are from the same plot
  response <-"placl_mean"
  if("leaf" %in% valid_factors){
    print("multiple leaf layers")
    random <- "~ Xf + Yf + plot_UID"
  } else {
    random <- "~ Xf + Yf" 
  }
  fixed <- paste("~", paste(valid_factors, collapse = " + ") )
  genotype = "genotype_name"
  genotype.as.random <- F  # to get BLUEs
  weights = 1/data$placl_var  # variance of placl observed in 8 different leaves
  spats <- SpATS(response = response, 
                 random = random, 
                 fixed = fixed,
                 spatial = ~PSANOVA(range_lot, row_lot, nseg = c(5,5), nest.div=c(2,2)),
                 genotype = genotype, 
                 genotype.as.random = genotype.as.random,
                 family = "gaussian",
                 weights = weights,
                 data = data,
                 control = list(maxit = 100, tolerance = 1e-03, monitoring = 0))
  
  # # heritability
  # print(getHeritability(spats))
  
  # get BLUEs and standard errors
  geno_pred[[i]] <- predict(spats, which = "genotype_name")[c("genotype_name", "predicted.values", "standard.errors")]
  
}

t <- names(subsets)
for(i in 1:length(geno_pred)){
  geno_pred[[i]]$tp <- strsplit(t[i], "\\.") %>% lapply("[[", 2) %>% unlist()
  geno_pred[[i]]$exp <- strsplit(t[i], "\\.") %>% lapply("[[", 1) %>% unlist()
}

Blues <- bind_rows(geno_pred) %>% 
  mutate(assess = paste(tp, exp, sep = "_"))

Blues <- Blues %>%
  mutate(assess = as.factor(assess),
         exp = as.factor(exp),
         genotype_name = as.factor(genotype_name))

# Second stage: Estimate heritability across time points, and get an overall estimate of the genotype value
weights_ <- 1/(Blues$standard.errors^2)
mod1 <- asreml(data=Blues,
       fixed = predicted.values ~ exp,
       random = ~ genotype_name + assess:genotype_name,
       weights = weights_,
       trace=TRUE)
mod1 <- update(mod1)
summary(mod1)$varcomp
vg <- summary(mod1)$varcomp["genotype_name", "component"]
ve <- summary(mod1)$varcomp["units!R", "component"]
heritability <- vg / (vg + ve)

genoblups <- predict(mod1, classify = "genotype_name")$pvals
genoblups$predicted.value <- inv.logit(genoblups$predicted.value)
write.csv(genoblups, "Z:/Public/Jonas/011_STB_leaf_tracking/data/sev_blups.csv")

# ============================================================================== -

# >> Visual scoring data ----

# get raw plot visual scoring data
scordat <- read_excel("O:/Evaluation/FIP/2016/WW012/RefTraits/Septoria/2016_FIP_scoringdataFabio.xlsx", 
                      sheet = 1, skip = 3)
# get rid of unwanted rows
scordat <- scordat %>% 
  dplyr::filter(!plot %in% c("TEM", "lot3")) %>%
  dplyr::filter(!is.na(plot)) %>% 
  dplyr::filter(!is.na(nosemis)) %>% 
  mutate(plot = as.numeric(plot))

# still some non-unique plots
# duplicate for 787, but 788 is missing
# one 176 in place of 716
# I assume these are errors from manually preparing the table and correct them
scordat[634, ]$plot <- 716
scordat[706, ]$plot <- 788

# re-format
scordat <- scordat %>% 
  dplyr::select(3, 4, 5, 7, 9) %>% 
  rename("20160520" = 3,
         "20160621" = 4,
         "20160629" = 5)

# calculate AUDPC (what is the reference date ? taken from Fabio)
audpc <- scordat %>% 
  mutate(AUDPC = `20160520`*4/2+((`20160520`+`20160621`)*(44-4)/2)+((`20160621`+`20160629`)*(52-44)/2)) %>% 
  dplyr::select(plot, AUDPC)

# add experimental design
design <- read_xls("O:/Evaluation/FIP/2016/Feldbuch/WW/2015_Design_FPWW012_FPWW016.xls",
                   sheet = 4, skip = 5) %>% 
  dplyr::select(Plot_ID, Plot, Lot, RangeLot, RowLot, Gen_Name, Gen_ID, Checks)
scdat <- left_join(design, audpc, by = c("Plot" = "plot"))
# 720 rows, OK. 

# make grid gap between lots
moddat <- scdat %>% 
  mutate(RangeLot = ifelse(Lot == 3, RangeLot + 18, RangeLot),
         RowLot = ifelse(Lot == 3, RowLot + 24, RowLot))

# convert to factors
moddat <- moddat %>% 
  mutate(RowF = as.factor(RowLot),
         RangeF = as.factor(RangeLot),
         Lot = as.factor(Lot),
         Gen_ID = as.factor(Gen_ID),
         Check = as.factor(Checks)) %>% 
  dplyr::select(-Checks) %>% 
  relocate(AUDPC, .after = Check)

# PARADOR repliace is extreme outlier - remove
moddat[moddat$Plot_ID == "FPWW0120633", ]$AUDPC <- NA

# examine correlation across lots
pdat <- moddat %>% dplyr::select(Lot, Gen_Name, AUDPC)
pdat <- split(pdat, pdat$Lot)
names(pdat[[1]])[3] <- "AUDPC_3"
names(pdat[[2]])[3] <- "AUDPC_4"
pdat <- full_join(pdat[[1]], pdat[[2]], by = "Gen_Name")
ggplot(pdat) +
  geom_point(aes(x = AUDPC_3, y = AUDPC_4))
cor(pdat$AUDPC_3, pdat$AUDPC_4, use = "pairwise.complete.obs")
# No correlation

# perform spatial correction
response <- "AUDPC"
random <- as.formula("~ RowF + RangeF")
fixed <- as.formula("~ Check")
genotype = "Gen_Name"
genotype.as.random <- T  # to get variance estimate
spats <- SpATS(response = response,
               random = random, 
               fixed = fixed,
               spatial = ~PSANOVA(RowLot, RangeLot, nseg = c(12,12), nest.div=c(2,2)),
               genotype = genotype, 
               genotype.as.random = genotype.as.random,
               data = moddat,
               control = list(maxit = 100, tolerance = 1e-03, monitoring = 0))
getHeritability(spats)
plot(spats,  spaTrend = "percentage")

pred <- predict.SpATS(spats, which = "Gen_Name") %>% 
  dplyr::select(Gen_Name, predicted.values) %>% as_tibble() %>% 
  rename("AUDPC" = 2)
write.csv(pred, "Z:/Public/Jonas/011_STB_leaf_tracking/data/AUDPC_blups.csv",
          row.names = F)

# ============================================================================== -

# Correlation with overall QR ----

# correlation between measures of disease severity and lesion growth rate

sev_blups <- read_csv("Z:/Public/Jonas/011_STB_leaf_tracking/data/sev_blups.csv") %>% 
  dplyr::select(genotype_name, predicted.value) %>% 
  rename(pred_sev = predicted.value)
inc_blups <- read_csv("Z:/Public/Jonas/011_STB_leaf_tracking/data/inc_blups.csv") %>% 
  dplyr::select(genotype_name, predicted.value) %>% 
  rename(pred_inc = predicted.value)
audpc_blups <- read_csv("Z:/Public/Jonas/011_STB_leaf_tracking/data/AUDPC_blups.csv") %>% 
  rename(genotype_name = Gen_Name)
growth_blups <- read_csv("Z:/Public/Jonas/011_STB_leaf_tracking/data/growth_blups.csv") %>% 
  dplyr::select(genotype_name, predicted.value) %>% 
  rename(pred_growth = predicted.value)
slope_blups <- read_csv("Z:/Public/Jonas/011_STB_leaf_tracking/data/slopes_blups.csv") %>% 
  mutate(across(where(is.numeric), ~ .x * 1e5))

# merge
pdat <- full_join(sev_blups, inc_blups) %>% 
  full_join(., growth_blups) %>% 
  left_join(., audpc_blups) %>% 
  left_join(., slope_blups)

# outlier 
mod <- lm(pred_growth ~ pred_sev, data = pdat)
summary(mod)
d <- cooks.distance(mod)
m_d <- mean(d[-2])
s_d <- sd(d[-2])

pdat2 <- pdat %>% filter(genotype_name != "AUBUSSON")
plotdata <- list(pdat, pdat2)

# Function to compute R² and p-value and return annotation text
get_lm_label <- function(x, y) {
  model <- lm(y ~ x)
  r2 <- summary(model)$r.squared
  p <- summary(model)$coefficients[2, 4]
  
  # Use scientific format for p-values and bquote for expression rendering
  if (p < 2.2e-16) {
    expr <- bquote(R^2 == .(round(r2, 2)) ~ "\n" ~ italic(P) < 2.2 %*% 10^-16)
  } else if(p < 0.001) {
    # Get base and exponent manually
    sci_p <- formatC(p, format = "e", digits = 2)
    base <- sub("e.*", "", sci_p)
    exp <- sub(".*e", "", sci_p)
    expr <- bquote(R^2 == .(round(r2, 2)) ~ "\n" ~ italic(P) == .(base) %*% 10^.(as.integer(exp)))
  } else {
    expr <- bquote(R^2 == .(round(r2, 2)) ~ "\n" ~ italic(P) == .(round(p, 2)))
  }
  return(expr)
}

data_all <- plotdata[[1]]
data_nooutlier <- plotdata[[2]]

# Plot 1: pred_sev vs pred_inc
# Plot 2: pred_sev vs predicted.value
label1 <- get_lm_label(data_all$pred_sev, data_all$pred_inc)
label1_2 <- get_lm_label(data_nooutlier$pred_sev, data_nooutlier$pred_inc)
p1 <- ggplot(data_all, aes(x = pred_sev, y = pred_inc)) +
  geom_point(aes(color = ifelse(genotype_name == "AUBUSSON", "red", "black"))) + # Set color for AUBUSSON
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_smooth(data = data_nooutlier, method = "lm", se = FALSE, lty = 2, color = "slateblue") +
  annotate("text", x = min(data_all$pred_sev, na.rm = T), y = max(data_all$pred_inc, na.rm = T), 
           label = label1, hjust = 0, vjust = 1, size = 3) +
  annotate("text", x = min(data_all$pred_sev, na.rm = T), y = 0.92*max(data_all$pred_inc, na.rm = T), 
           label = label1_2, hjust = 0, vjust = 1, color = "slateblue", size = 3) +
  labs(x = "Conditional Severity", y = "Incidence") +
  scale_color_identity() +  ggtitle("A") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

# Plot 2: pred_sev vs predicted.value
label2 <- get_lm_label(data_all$pred_sev, data_all$pred_growth)
label2_2 <- get_lm_label(data_nooutlier$pred_sev, data_nooutlier$pred_growth)
p2 <- ggplot(data_all, aes(x = pred_sev, y = pred_growth)) +
  geom_point(aes(color = ifelse(genotype_name == "AUBUSSON", "red", "black"))) + # Set color for AUBUSSON
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_smooth(data = data_nooutlier, method = "lm", se = FALSE, lty = 2, color = "slateblue") +
  annotate("text", x = min(data_all$pred_sev, na.rm = T), y = max(data_all$pred_growth, na.rm = T), 
           label = label2, hjust = 0, vjust = 0, size = 3) +
  annotate("text", x = min(data_all$pred_sev, na.rm = T), y = 0.97*max(data_all$pred_growth, na.rm = T), 
           label = label2_2, hjust = 0, vjust = 0, color = "slateblue", size = 3) +
  labs(x = "Conditional Severity", y = "Lesion Growth") +
  scale_y_continuous(limits = c(0.0006, 0.001)) +
  scale_color_identity() +
  ggtitle("B") +
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"))

# Plot 3: pred_inc vs predicted.value
label3 <- get_lm_label(data_all$pred_inc, data_all$pred_growth)
label3_2 <- get_lm_label(data_nooutlier$pred_inc, data_nooutlier$pred_growth)
p3 <- ggplot(data_all, aes(x = pred_inc, y = pred_growth)) +
  geom_point(aes(color = ifelse(genotype_name == "AUBUSSON", "red", "black"))) + # Set color for AUBUSSON
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_smooth(data = data_nooutlier, method = "lm", se = FALSE, lty = 2, color = "slateblue") +
  annotate("text", x = min(data_all$pred_inc, na.rm = T), y = max(data_all$pred_growth, na.rm = T), 
           label = label3, hjust = 0, vjust = 0, size = 3) +
  annotate("text", x = min(data_all$pred_inc, na.rm = T), y = 0.97*max(data_all$pred_growth, na.rm = T), 
           label = label3_2, hjust = 0, vjust = 0, color = "slateblue", size = 3) +
  labs(x = "Incidence", y = "Lesion Growth") + 
  scale_color_identity() + 
  scale_y_continuous(limits = c(0.0006, 0.001)) +
  ggtitle("C") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(face = "bold"))

# Plot 4: AUDPC vs predicted.value
label4 <- get_lm_label(data_all$AUDPC, data_all$pred_growth)
label4_2 <- get_lm_label(data_nooutlier$AUDPC, data_nooutlier$pred_growth)
p4 <- ggplot(data_all, aes(x = AUDPC, y = pred_growth)) +
  geom_point(aes(color = ifelse(genotype_name == "AUBUSSON", "red", "black"))) + # Set color for AUBUSSON
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_smooth(data = data_nooutlier, method = "lm", se = FALSE, lty = 2, color = "slateblue") +
  annotate("text", x = min(data_all$AUDPC, na.rm = T), y = max(data_all$pred_growth, na.rm = T), 
           label = label4, hjust = 0, vjust = 0, size = 3) +
  annotate("text", x = min(data_all$AUDPC, na.rm = T), y = 0.97*max(data_all$pred_growth, na.rm = T), 
           label = label4_2, hjust = 0, vjust = 0, color = "slateblue", size = 3) +
  labs(x = "AUDPC", y = "Lesion Growth") + 
  scale_color_identity() + 
  scale_y_continuous(limits = c(0.0006, 0.001)) +
  ggtitle("C") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(face = "bold"))

# Plot 5: AUDPC vs Severity
label5 <- get_lm_label(data_all$AUDPC, data_all$pred_sev)
label5_2 <- get_lm_label(data_nooutlier$AUDPC, data_nooutlier$pred_sev)
p5 <- ggplot(data_all, aes(x = AUDPC, y = pred_sev)) +
  geom_point(aes(color = ifelse(genotype_name == "AUBUSSON", "red", "black"))) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_smooth(data = data_nooutlier, method = "lm", se = FALSE, lty = 2, color = "slateblue") +
  annotate("text", x = min(data_all$AUDPC, na.rm = T), y = max(data_all$pred_sev, na.rm = T), 
           label = label5, hjust = 0, vjust = 0, size = 3) +
  annotate("text", x = min(data_all$AUDPC, na.rm = T), y = 0.97*max(data_all$pred_sev, na.rm = T), 
           label = label5_2, hjust = 0, vjust = 0, color = "slateblue", size = 3) +
  labs(x = "AUDPC", y = "Lesion Growth") + 
  scale_color_identity() + 
  ggtitle("C") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(face = "bold"))

all <- ggarrange(p1, p2, p3, nrow = 1, widths = c(1, 1.075, 0.925))
png("Z:/Public/Jonas/011_STB_leaf_tracking/Figures/paper/growth_qr_corr.png", width = 9, height = 3.5, units = 'in', res = 400)
plot(all)
dev.off()

# Plot 6: Slopes vs Severity
pdat_all <- data_all %>% 
  dplyr::select(genotype_name, pred_sev, slope_temp, slope_rh, slope_age) %>% 
  tidyr::pivot_longer(slope_temp:slope_age)
pdat_nooutlier <- data_nooutlier %>% 
  dplyr::select(genotype_name, pred_sev, slope_temp, slope_rh, slope_age) %>% 
  tidyr::pivot_longer(slope_temp:slope_age)
# Create label data frame for all data
label_df_all <- pdat_all %>%
  group_by(name) %>% nest() %>% 
  mutate(label = purrr::map(data, ~get_lm_label(.$value, .$pred_sev)),
         x = purrr::map_dbl(data, ~min(.$value, na.rm = T)),
         y = purrr::map_dbl(data, ~max(.$pred_sev, na.rm = T)-0.02)) %>% 
  dplyr::select(-data)
label_df_nooutlier <- pdat_nooutlier %>%
  group_by(name) %>% nest() %>% 
  mutate(label = purrr::map(data, ~get_lm_label(.$value, .$pred_sev)),
         x = purrr::map_dbl(data, ~min(.$value, na.rm = T)),
         y = purrr::map_dbl(data, ~max(.$pred_sev, na.rm = T)+0.0125)) %>% 
  dplyr::select(-data)
# for facet names
trait_labels <- c(slope_age = "Lesion Age", 
                  slope_rh = "Relative Humdity",
                  slope_temp = "Temperature")
p6 <- ggplot(pdat_all, aes(x = value, y = pred_sev)) +
  geom_point(aes(color = ifelse(genotype_name == "AUBUSSON", "red", "black"))) +
  geom_smooth(method = "lm", se = F, color = "black") +
  geom_smooth(data = pdat_nooutlier, method = "lm", se = F, lty = 2, color = "slateblue") +
  geom_text(data = label_df_all, aes(x = x, y = y, label = label), 
            inherit.aes = T, hjust = 0, vjust = 0, size = 3, parse = T) +
  geom_text(data = label_df_nooutlier, aes(x = x, y = y, label = label), 
            inherit.aes = FALSE, hjust = 0, vjust = 0, size = 3, parse = T, color = "slateblue") +
  labs(x = "Slope (×10⁵)", y = "Conditional Severity") + 
  scale_color_identity() + 
  facet_wrap(~ name, scales = "free_x", labeller = labeller(name = trait_labels)) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))
p6
png("Z:/Public/Jonas/011_STB_leaf_tracking/Figures/paper/slope_qr_corr.png", width = 9, height = 3.5, units = 'in', res = 400)
plot(p6)
dev.off()

# Plot 7: Slopes vs Severity
pdat_all <- data_all %>% 
  dplyr::select(genotype_name, pred_growth, slope_temp, slope_rh, slope_age) %>% 
  tidyr::pivot_longer(slope_temp:slope_age)
pdat_nooutlier <- data_nooutlier %>% 
  dplyr::select(genotype_name, pred_growth, slope_temp, slope_rh, slope_age) %>% 
  tidyr::pivot_longer(slope_temp:slope_age)
# Create label data frame for all data
label_df_all <- pdat_all %>%
  group_by(name) %>% nest() %>% 
  mutate(label = purrr::map(data, ~get_lm_label(.$value, .$pred_growth)),
         x = purrr::map_dbl(data, ~min(.$value, na.rm = T)),
         y = purrr::map_dbl(data, ~max(.$pred_growth, na.rm = T)-0.000035)) %>% 
  dplyr::select(-data)
label_df_nooutlier <- pdat_nooutlier %>%
  group_by(name) %>% nest() %>% 
  mutate(label = purrr::map(data, ~get_lm_label(.$value, .$pred_growth)),
         x = purrr::map_dbl(data, ~min(.$value, na.rm = T)),
         y = purrr::map_dbl(data, ~max(.$pred_growth, na.rm = T)-0.000015)) %>% 
  dplyr::select(-data)
# for facet names
trait_labels <- c(slope_age = "Lesion Age", 
                  slope_rh = "Relative Humdity",
                  slope_temp = "Temperature")
p7 <- ggplot(pdat_all, aes(x = value, y = pred_growth)) +
  geom_point(aes(color = ifelse(genotype_name == "AUBUSSON", "red", "black"))) +
  geom_smooth(method = "lm", se = F, color = "black") +
  geom_smooth(data = pdat_nooutlier, method = "lm", se = F, lty = 2, color = "slateblue") +
  geom_text(data = label_df_all, aes(x = x, y = y, label = label), 
            inherit.aes = T, hjust = 0, vjust = 0, size = 3, parse = T) +
  geom_text(data = label_df_nooutlier, aes(x = x, y = y, label = label), 
            inherit.aes = FALSE, hjust = 0, vjust = 0, size = 3, parse = T, color = "slateblue") +
  labs(x = "Slope (×10⁵)", y = "Genotype Main Effect") + 
  scale_color_identity() + 
  facet_wrap(~ name, scales = "free_x", labeller = labeller(name = trait_labels)) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))
p7
png("Z:/Public/Jonas/011_STB_leaf_tracking/Figures/paper/slope_main_corr.png", width = 9, height = 3.5, units = 'in', res = 400)
plot(p7)
dev.off()

# ============================================================================== -
