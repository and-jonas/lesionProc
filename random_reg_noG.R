#=============================================================================== -

#HEADER ----

# Author: Lukas Roth, ETH Zürich
# Copyright (C) 2025  ETH Zürich, Lukas Roth (lukas.roth@usys.ethz.ch)

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
library(tidyverse)
library(asreml)

# get data and ensure variables are global unique
moddat <- readRDS("Z:/Public/Jonas/011_STB_leaf_tracking/data/growth_data_v2.rds")
moddat <- moddat[complete.cases(moddat), ]
data <- moddat %>%
  mutate(batch = as.factor(paste0(year, "_", as.character(batch))),
         leaf_nr = as.factor(paste0(plot_UID, "_", leaf_nr)),
         lesion_nr = as.factor(paste0(plot_UID, "_", leaf_nr, "_", lesion_nr)),
         log_age = log(age)
  )

# ============================================================================== -

# Genotype agnostic model on a plot level
# the previously selected best model is used

m_noG_rh_temp <- asreml(
  fixed= delta ~ age + mean_temp + mean_rh + area + maxdist,
  random =
    ~ leaf_nr + lesion_nr +
    ~str(~ age:plot_UID + mean_temp:plot_UID + mean_rh:plot_UID, ~diag(3):id(plot_UID)),
  data=data,
  maxiter=20
)
m_noG_rh_temp = update(m_noG_rh_temp)

# Predict plot_UID main effect
df_main <- predict(m_noG_rh_temp, classify = "plot_UID",
                   levels = list(
                     plot = unique(data$plot_UID)
                   )
)$pvals
df_age_main_noG <- df_main %>% select(-std.error) %>% 
  mutate(main = predicted.value)

# Predict plot_UID slopes
## Temperature
df_temp_resp_noG <- predict(m_noG_rh_temp, classify = "plot_UID:mean_temp",
                            levels= list(
                              mean_temp = c(0, 1)
                            )
)$pvals
df_temp_resp_noG <- df_temp_resp_noG %>% select(-std.error) %>% pivot_wider(names_from = mean_temp, values_from = predicted.value) %>%
  mutate(slope = `1` - `0`)

## Relative Humidity
df_rh_resp_noG <- predict(m_noG_rh_temp, classify = "plot_UID:mean_rh",
                          levels= list(
                            mean_rh = c(0, 1)
                          )
)$pvals
df_rh_resp_noG <- df_rh_resp_noG %>% select(-std.error) %>% pivot_wider(names_from = mean_rh, values_from = predicted.value) %>%
  mutate(slope = `1` - `0`)

## Lesion Age
df_age_resp_noG <- predict(m_noG_rh_temp, classify = "plot_UID:age",
                           levels= list(
                             age = c(0, 1)
                           )
)$pvals
df_age_resp_noG <- df_age_resp_noG %>% select(-std.error) %>% pivot_wider(names_from = age, values_from = predicted.value) %>%
  mutate(slope = `1` - `0`)

# add experimental design
df_design <- data %>% select(plot_UID, batch, year, genotype_name) %>% distinct()
df_temp_resp_noG <- inner_join(df_temp_resp_noG, df_design)
df_rh_resp_noG <- inner_join(df_rh_resp_noG, df_design)
df_age_resp_noG <- inner_join(df_age_resp_noG, df_design)

vars <- list(temp = df_temp_resp_noG,
             rh = df_rh_resp_noG,
             age = df_age_resp_noG)

# loop over response variables
for (i in 1:length(vars)) {
  var_name = names(vars)[i]
  df_ = vars[[i]]
  # Model over batch and years
  m_genotype <- asreml(fixed= slope ~ 1,
                       random = ~ batch + year + genotype_name,
                       data=df_,
                       maxiter=20
  )
  df_preds <- predict(m_genotype, classify = "genotype_name", only = "genotype_name")$pvals %>%
    as.data.frame() %>%
    arrange(predicted.value)
  print(df_preds)
  
  # Cullis Heritability
  m_genotype <- update(m_genotype)
  
  # Genetic variance
  vc.g <- summary(m_genotype)$varcomp['genotype_name', 'component']
  # Mean variance of a difference of two genotypic BLUPs
  # obtain squared s.e.d. matrix
  vdBLUP.mat <- predict(m_genotype, classify = "genotype_name", only = "genotype_name", sed = TRUE)$sed^2
  # take mean of upper triangle
  vdBLUP.avg <- mean(vdBLUP.mat[upper.tri(vdBLUP.mat, diag = FALSE)])
  
  H2Cullis <- 1 - (vdBLUP.avg / (vc.g * 2))
  print(paste(var_name, " - H2 Cullis: ", H2Cullis))
}

# main effect
df_main <- inner_join(df_main, df_design)

# main effect
m_genotype <- asreml(fixed= predicted.value ~ batch + year,
                     random = ~ genotype_name,
                     data=df_main,
                     maxiter=20)
df_preds <- predict(m_genotype, classify = "genotype_name", only = "genotype_name")$pvals %>%
  as.data.frame() %>%
  arrange(predicted.value)
print(df_preds)

# Cullis Heritability
m_genotype <- update(m_genotype)
# Genetic variance
vc.g <- summary(m_genotype)$varcomp['genotype_name', 'component']
# Mean variance of a difference of two genotypic BLUPs
# obtain squared s.e.d. matrix
vdBLUP.mat <- predict(m_genotype, classify = "genotype_name", only = "genotype_name", sed = TRUE)$sed^2
# take mean of upper triangle
vdBLUP.avg <- mean(vdBLUP.mat[upper.tri(vdBLUP.mat, diag = FALSE)])

H2Cullis <- 1 - (vdBLUP.avg / (vc.g * 2))
print(paste("genotype main effect", " - H2 Cullis: ", H2Cullis))

# ============================================================================== -
