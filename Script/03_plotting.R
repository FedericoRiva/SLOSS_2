library(brms)
library(extrafont)
library(ggpubr)
library(gridExtra)
library(tidybayes)

## Fig. 3 is based on models from script 02_models
model_richness <- brm(formula = log_rich ~  mean_patch_area_log  + (mean_patch_area_log|dataset_id/habitat_amount), 
                      data   = alpha_plot,
                      iter = 5000,
                      seed   = 123,
                      chains = 4, cores = 10)

model_richness_decl <- brm(formula = log_rich_IUCN ~  mean_patch_area_log  + (mean_patch_area_log|dataset_id/habitat_amount), 
                           data   = alpha_plot_decl,
                           iter = 5000,
                           seed   = 123,
                           chains = 4, cores = 10)

model_beta <- brm(formula = slope_richness ~  beta.JTU  + (1|dataset),  
                  data   = beta_plot,
                  seed   = 123,
                  chains = 4, 
                  cores = 10)

save.image(file='review_ele.RData')
load('review_ele.RData')

## plot
predict_model_richness <- predict(model_richness)
model_richness_pred <- posterior_samples(model_richness)

predict_model_richness_decl <- predict(model_richness_decl)
model_richness_decl_pred <- posterior_samples(model_richness_decl)

predict_model_beta <- predict(model_beta)
model_beta_pred <- posterior_samples(model_beta)

alpha_plot$group <- paste(alpha_plot$dataset_id, alpha_plot$habitat_amount)
alpha_plot$mean_patch_area_log <- log10(alpha_plot$mean_patch_area)
alpha_plot$log_rich <- log2(alpha_plot$richness)


plot_rich <-  ggplot(data =  alpha_plot, 
       aes(x = mean_patch_area_log, y = log_rich, col = group)) +
  geom_jitter(alpha = 0.15)+
  geom_smooth(method = 'lm', se = F, size = 0.1, linetype = "solid", alpha = 0.1, formula = y~x )+
  geom_abline(intercept = model_richness_pred[1:1500, 1],
              slope     = model_richness_pred[1:1500, 2],
              size = 1/3, alpha = .05, color = "darkgrey") +
  geom_abline(intercept = mean(model_richness_pred[, 1]),
              slope     = mean(model_richness_pred[, 2]),
              size = 1.5, color ="black", linetype = 1) +
  xlab("Mean patch size (hectares)") + ylab("Species richness")+
  theme_few()+
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5), 
                     labels = c("0.1","1","10","100","1.000","10.000","100.000")) +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8), 
                     labels = c("1","2","4","8","16","32","64", "128", "256")) +
  theme(legend.position="none")+
  scale_color_viridis(discrete=TRUE, option="turbo")+
  theme(text=element_text(family="serif",  size=16)) +
  theme(text = element_text(size = 20)) +
  annotate(geom="text", family="serif", size =5,
           x=4, 
           y=8, parse=TRUE, label= "italic(β) == -0.15  (-0.28, -0.02)")
   

plot_rich

plot_slope <-  ggplot(data =  beta_plot, 
       aes(x = beta.JTU, y = slope_richness, col = dataset)) +
  geom_jitter(alpha = 0.25) +
  geom_abline(intercept = model_beta_pred[1:1500, 1],
              slope     = model_beta_pred[1:1500, 2],
              size = 1/3, alpha = .05, color = "darkgrey") +
  geom_abline(intercept = mean(model_beta_pred[, 1]),
              slope     = mean(model_beta_pred[, 2]),
              size = 1.5, color ="black", linetype = 1) +
  xlab("Turnover (Jaccard dissimilarity)") + 
  ylab(expression('Slope of log'[2]*'(species richness) ~ log'[10]*'(mean patch size)'))+
  #ylab("log2(species richness) ~ log10(mean patch size)")+
  ylim(-2,2)+
  theme_few()+
  theme(legend.position="none")+
  scale_color_viridis(discrete=TRUE, option="turbo")+
  theme(text=element_text(family="serif",  size=16)) +
  theme(text = element_text(size = 20)) +
  annotate(geom="text", family="serif", size = 5,
           x=0.78, 
           y=2, parse=TRUE, label= "italic(β) == -0.42 (-0.90, 0.06)")+
  geom_hline(yintercept=0, linetype="dashed", color = "darkred", size = 1) 
  
plot_slope

g <- arrangeGrob(plot_rich, plot_slope, ncol = 2)
g
ggsave("fig3.jpg", path = "Figures", width = 4500, height = 1600, units = "px",  device='jpeg', dpi=300, g)




## fig 4
plot_rich_iucn <-  ggplot(data =  alpha_plot, 
                     aes(x = mean_patch_area_log, y = log_rich, col = group)) +
  geom_smooth(method = 'lm', se = F, size = 0.1, linetype = "blank",  formula = y~x )+
  geom_abline(intercept = model_richness_pred[1:1500, 1],
              slope     = model_richness_pred[1:1500, 2],
              size = 1/3, alpha = .05, color = "darkgrey") +
  geom_abline(intercept = mean(model_richness_pred[, 1]),
              slope     = mean(model_richness_pred[, 2]),
              size = 1.5, color ="black", linetype = 1) +
  
  geom_abline(intercept = model_richness_decl_pred[1:1500, 1],
              slope     = model_richness_decl_pred[1:1500, 2],
              size = 1/3, alpha = .05, color = "red") +
  geom_abline(intercept = mean(model_richness_decl_pred[, 1]),
              slope     = mean(model_richness_decl_pred[, 2]),
              size = 1.5, color ="darkred", linetype = 1) +
  
  xlab("Mean patch size (hectares)") + ylab("Species richness")+
  theme_few()+
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5), 
                     labels = c("0.1","1","10","100","1.000","10.000","100.000")) +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8), 
                     labels = c("1","2","4","8","16","32","64", "128", "256")) +
  theme(legend.position="none")+
  scale_color_viridis(discrete=TRUE, option="turbo")+
  theme(text=element_text(family="serif",  size=16)) +
  theme(text = element_text(size = 20)) +
  annotate(geom="text", family="serif", size =5,
           x=4, 
           y=8, parse=TRUE, label= "italic(β) == -0.15  (-0.28, -0.02)")+
  annotate(geom="text", family="serif", size =5,
           x=4, 
           y=7.2, parse=TRUE, label= "italic(β[IUCN]) == -0.55  (-1.58, 0.47)")


#plot_rich_iucn
ggsave("fig4.jpg", path = "Figures", width = 2250, height = 1600, units = "px",  device='jpeg', dpi=300, plot_rich_iucn)
