# prepare data
library(tidyr)
library(data.table)
library(dplyr)
library(purrr)
library(rlist)
library(splitstackshape)

# diversity estimation
library(iNEXT)
library(SpadeR)
library(vegan)
library(paleotree)
library(betapart)

# plotting
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(viridis)

# modeling
library(brms)
library(glmmTMB)
library(MuMIn)
library(sjPlot)
library(webshot)
##
library(effects)
library(ggeffects)
library(extrafont)
library(ggpubr)
library(gridExtra)
library(tidybayes)

## open datasets
alpha = fread("Data\\analysis_and_plots\\table_analysis.csv", header = TRUE)
beta = fread("Data\\analysis_and_plots\\table_analysis_beta.csv", header = TRUE)
##

# number of species and patches in the final dataset analyzed
# data = fread("Data\\fragSAD_predicts_ewers.csv", header = TRUE)
# data_final <- data[data$dataset_label %in% levels(as.factor(alpha$dataset_id)), ]
# levels(as.factor(paste(data_final$frag_id, data_final$dataset_label))) # 1149 patches
# levels(as.factor(data_final$species)) # 4351 species 

# Data preparation 
# generate a column to split the data in a list of tables to calculate slopes within each scenario
alpha$split <- paste(alpha$dataset_id, alpha$simulation_number, alpha$habitat_amount)
alpha_splitted <- split(alpha, alpha$split) 
slopes <- split(alpha, alpha$split)

# calculate for each scenario the slope of the log2(species richness) ~ log10(mean patch area) relationship 
for(i in 1: length(slopes)) {
  slopes[[i]] <- summary(lm(log2(slopes[[i]]$richness+0.1) ~ log10(slopes[[i]]$mean_patch_area)))[[4]][[2]]
}

# we also calculated the slope of scaled evenness ~ scaled mean patch area, but we do not present it in the manuscript as it is not discussed in Chase et al. 2020
slopes_evenness <- split(alpha, alpha$split)
for(i in 1: length(slopes)) {
  slopes_evenness[[i]] <- summary(lm(slopes_evenness[[i]]$evenness_scaled ~ slopes_evenness[[i]]$mean_patch_area_sc))[[4]][[2]]
}

names_slopes <- names(slopes)

slopes <- do.call(rbind, slopes)
slopes_evenness <- do.call(rbind, slopes_evenness)

slopes_beta <- as.data.frame(cbind(slopes, slopes_evenness))
colnames(slopes_beta)[1] <- "slope_richness"
colnames(slopes_beta)[2] <- "slope_evenness"
slopes_beta$split <- names_slopes

#merge slopes and beta diversity
beta$split <- paste(beta$dataset, beta$simulation, beta$habitat_amount)
beta <- merge(beta, slopes_beta, by = "split")

# subset randomly one dataset for the analysis presented in the manuscript; see appendix 1.3
# included larger model with the entire dataset
set.seed(144)
sample(1:100, 1)
alpha_plot <- subset(alpha, simulation_number ==39)
beta_plot <- subset(beta, simulation ==39)

# transform evenness for beta regression following 
# Smithson, M. & Verkuilen, J. A better lemon squeezer? Maximum-likelihood regression with beta-distributed dependent variables. Psychol. Methods 11, 54–71 (2006). DOI: 10.1037/1082-989X.11.1.54
alpha_plot$evenness_t <- (alpha_plot$evenness*(nrow(alpha_plot) - 1) + 0.5)/nrow(alpha_plot)
alpha_plot$log_rich <- log2(alpha_plot$richness)
alpha_plot$mean_patch_area_log <- log10(alpha_plot$mean_patch_area)
alpha_plot$habitat_amount_per_dataset <- paste(alpha_plot$dataset_id, alpha_plot$habitat_amount)

# for the analysis of IUCN declining species, subset only datasets with at least one declining species across all sets of patches
alpha_plot_decl <- split(alpha_plot, alpha_plot$dataset_id)

for(i in 1:length(alpha_plot_decl)){
  if (sum(alpha_plot_decl[[i]]$richness_IUCN) == 0) {
    alpha_plot_decl[[i]] <- NA
  }
}

alpha_plot_decl <- Filter(Negate(anyNA),alpha_plot_decl)
alpha_plot_decl <- do.call(rbind, alpha_plot_decl)
alpha_plot_decl$log_rich_IUCN <- log2(alpha_plot_decl$richness_IUCN+0.1)

########################
## MODELS
########################
# see appendix 1.3 for model results; the tables reported in appendix generated using function tab_model()

# Bayesian mixed effect models for richness
# model 1-6, page 9-14 in appendix 1
model_richness <- brm(formula = log_rich ~  mean_patch_area_log  + (mean_patch_area_log|dataset_id/habitat_amount), 
                      data   = alpha_plot,
                      seed   = 123,
                      chains = 4, cores = 10)

model_richness_decl <- brm(formula = log_rich_IUCN ~  mean_patch_area_log  + (mean_patch_area_log|dataset_id/habitat_amount), 
                           data   = alpha_plot_decl,
                           seed   = 123,
                           chains = 4, cores = 10)

model_richness_taxa <- brm(formula = log_rich ~  mean_patch_area_log * taxa  + (mean_patch_area_log|dataset_id/habitat_amount), 
                           data   = alpha_plot,
                           seed   = 123,
                           chains = 4, cores = 10)

model_richness_time <- brm(formula = log_rich ~  mean_patch_area_log * time.since.fragmentation  + (mean_patch_area_log|dataset_id/habitat_amount), 
                           data   = alpha_plot,
                           seed   = 123,
                           chains = 4, cores = 10)

model_richness_matrix <- brm(formula = log_rich ~  mean_patch_area_log * Matrix.category  + (mean_patch_area_log|dataset_id/habitat_amount), 
                             data   = alpha_plot,
                             seed   = 123,
                             chains = 4, cores = 10)

model_richness_continent <- brm(formula = log_rich ~  mean_patch_area_log * continent  + (mean_patch_area_log|dataset_id/habitat_amount), 
                                data   = alpha_plot,
                                seed   = 123,
                                chains = 4, cores = 10)

# supplementary model for richness based on the full dataset with 995,400 observations
# models 7-11 from pages 15-20
model_richness_sc <- glmmTMB(richness_scaled ~ 0 + mean_patch_area_sc  + (0+mean_patch_area_sc|simulation_number/dataset_id/habitat_amount),# +
                             data = alpha, family = "gaussian", control = glmmTMBControl(parallel = 10))

model_richness_sc_taxa <- glmmTMB(richness_scaled ~ 0 + mean_patch_area_sc * taxa  + (0+mean_patch_area_sc|simulation_number/dataset_id/habitat_amount),# +
                                  data = alpha, family = "gaussian", control = glmmTMBControl(parallel = 10))

model_richness_sc_region <- glmmTMB(richness_scaled ~ 0 + mean_patch_area_sc * continent  + (0+mean_patch_area_sc|simulation_number/dataset_id/habitat_amount),# +
                                    data = alpha, family = "gaussian", control = glmmTMBControl(parallel = 10))

model_richness_sc_matrix <- glmmTMB(richness_scaled ~ 0 + mean_patch_area_sc * Matrix.category  + (0+mean_patch_area_sc|simulation_number/dataset_id/habitat_amount),# +
                                    data = alpha, family = "gaussian", control = glmmTMBControl(parallel = 10))

model_richness_sc_time <- glmmTMB(richness_scaled ~ 0 + mean_patch_area_sc * time.since.fragmentation  + (0+mean_patch_area_sc|simulation_number/dataset_id/habitat_amount),# +
                                  data = alpha, family = "gaussian", control = glmmTMBControl(parallel = 10))

# to evaluate model predictions, use the following line of code
# plot(allEffects(model_richness_sc), type = "response")


# supplementary model for evenness based on the full dataset with 995,400 observations
# models 12-16 from pages 21-26
model_evenness_sc <- glmmTMB(evenness_scaled ~ 0 + mean_patch_area_sc  + (0+mean_patch_area_sc|simulation_number/dataset_id/habitat_amount),# +
                             data = alpha, family = "gaussian", control = glmmTMBControl(parallel = 10))

model_evenness_sc_taxa <- glmmTMB(evenness_scaled ~ 0 + mean_patch_area_sc * taxa  + (0+mean_patch_area_sc|simulation_number/dataset_id/habitat_amount),# +
                                  data = alpha, family = "gaussian", control = glmmTMBControl(parallel = 10))

model_evenness_sc_region <- glmmTMB(evenness_scaled ~ 0 + mean_patch_area_sc * continent  + (0+mean_patch_area_sc|simulation_number/dataset_id/habitat_amount),# +
                                    data = alpha, family = "gaussian", control = glmmTMBControl(parallel = 10))

model_evenness_sc_matrix <- glmmTMB(evenness_scaled ~ 0 + mean_patch_area_sc * Matrix.category  + (0+mean_patch_area_sc|simulation_number/dataset_id/habitat_amount),# +
                                    data = alpha, family = "gaussian", control = glmmTMBControl(parallel = 10))

model_evenness_sc_time <- glmmTMB(evenness_scaled ~ 0 + mean_patch_area_sc * time.since.fragmentation  + (0+mean_patch_area_sc|simulation_number/dataset_id/habitat_amount),# +
                                  data = alpha, family = "gaussian", control = glmmTMBControl(parallel = 10))
# to evaluate model predictions, use the following line of code
# plot(allEffects(model_evenness_sc), type = "response")

# models for relationship between slope of the richness ~ mean patch size relationship and turnover
# model 17, page 28 appendix 1
model_beta <- brm(formula = slope_richness ~  beta.JTU  + (1|dataset),  
                  data   = beta_plot,
                  seed   = 123,
                  chains = 4, 
                  cores = 10)

# model 18, page 29 appendix 1
model_beta_ruz <- brm(formula = slope_richness ~  beta.RUZ.BAL + (1|dataset), 
                                data   = beta_plot,
                                seed   = 123,
                                chains = 4, cores = 10)

# models 19-21, pages 29-31 in appendix 1
alpha_plot$log_rich <- log2(alpha_plot$richness)
alpha_plot$log_area <- log10(alpha_plot$total_patch_area)


model_SAR_null <- glmmTMB(formula = log_rich ~ log_area  + (1|dataset_id),
                           data   = alpha_plot, family = "gaussian", control = glmmTMBControl(parallel = 10))

model_frag <- glmmTMB(formula = log_rich ~ mean_patch_area_sc + (1|dataset_id),
                       data   = alpha_plot, family = "gaussian", control = glmmTMBControl(parallel = 10))

model_SAR <- glmmTMB(formula = log_rich ~ log_area + mean_patch_area_sc + (1|dataset_id),
                      data   = alpha_plot, family = "gaussian", control = glmmTMBControl(parallel = 10))


