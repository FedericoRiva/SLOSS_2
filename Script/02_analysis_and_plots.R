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

# modeling
library(glmmTMB)




## open datasets
alpha = fread("C:\\Users\\feder\\OneDrive\\Desktop\\GITHUB\\SLOSS_2\\Data\\analysis_and_plots\\table_analysis.csv", header = TRUE)
beta = fread("C:\\Users\\feder\\OneDrive\\Desktop\\GITHUB\\SLOSS_2\\Data\\analysis_and_plots\\table_analysis_beta.csv", header = TRUE)
##

alpha_plot <- subset(alpha, simulation_number ==3)
alpha_plot$group <- paste(alpha_plot$dataset_id, alpha_plot$habitat_amount)
# plot 1
ggplot(alpha_plot, aes(mean_patch_area_sc,richness_scaled, col = group)) + # also works without specifying data, x, and y
  geom_point()+
  geom_smooth(method = 'lm', se = F, col = 'black')+
  theme(legend.position="none")

# 9554 sets of patches from 425 scenarios. 

ggplot(alpha_plot, aes(mean_patch_area_sc,richness_scaled, col = group)) + 
  geom_point(size = 0.5, shape = 1) + # change size and colour
  labs(y = "Species richness  (scaled)", x = "Mean patch area (scaled)") + # rename axes
  geom_smooth(method = 'lm', se = F, size = 0.2) + # fit linear regression line 
  geom_smooth(method = 'lm', se = TRUE, aes(group = 1), col = "black", size = 2)+
  theme_bw()+
  theme(legend.position="none")

ggplot(alpha_plot, aes(mean_patch_area_sc,evenness_scaled, col = group)) + 
  geom_point(size = 0.5, shape = 1) + # change size and colour
  labs(y = "Evenness richness  (scaled)", x = "Mean patch area (scaled)") + # rename axes
  geom_smooth(method = 'lm', se = F, size = 0.2) + # fit linear regression line 
  geom_smooth(method = 'lm', se = TRUE, aes(group = 1), col = "black", size = 2)+
  theme_bw()+
  theme(legend.position="none")



model <- glmmTMB(evenness_scaled ~ mean_patch_area_sc + #* habitat_amount + 
                   (1|habitat_amount) + 
                   (1|dataset_id/patch_set_number) + #patch set number is the simulated assemblage, from 1 to 100
                   (1|simulation_number), 
                 data = table_analysis, family = "gaussian")

table_analysis$log_rich <- log(table_analysis$richness + 1)
table_analysis$log_area <- log(table_analysis$total_patch_area)

model <- glmmTMB(log_rich ~ log_area *
                   mean_patch_area_sc + #* habitat_amount + 
                   (1|habitat_amount) + 
                   (1|dataset_id/patch_set_number) + #patch set number is the simulated assemblage, from 1 to 100
                   (1|simulation_number), 
                 data = table_analysis, family = "gaussian")

summary(model)

table_analysis$predicted_richness <- predict(model, type = "response")

plot(allEffects(model), type = "response")





plot_model <- ggpredict(model, 
                        c("log_area [all]", "mean_patch_area_sc [-1, 0, 1]" ),
                        type = "re"
) 

plot(plot_model) + 
  labs(
    x = "log (area)", 
    y = "log (richness + 1)", 
    title = ""  ) + 
  labs(colour = "SD from mean patch size in a scenario")




# check ggplot tutorial here https://aosmith.rbind.io/2018/11/16/plot-fitted-lines/

ggplot(table_analysis, aes(x=log10_mean_patch_area, y=log10ric, color=dataset_id, shape=habitat_amount)) +
  #geom_point() +
  geom_smooth(method = "lm")+
  #geom_line(aes(y = predicted_richness))+
  theme(legend.position="none")

# table_analysis$evenness2 <- (table_analysis$evenness * (nrow(table_analysis)-1) + 0.5) / nrow(table_analysis)
# 
# model <- glmmTMB(evenness2 ~ mean_patch_area * habitat_amount + taxa +
#                    (1|dataset_id/patch_set_number) +
#                    (1|simulation_number), 
#                  data = table_analysis, family=beta_family(link="logit"))
# 
# summary(model)
# AIC(model)


plot_model <- ggpredict(model, c("mean_patch_area", "taxa"), type = "re") 
plot(plot_model)



