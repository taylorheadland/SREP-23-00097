## This script runs the modelling necessary to obtain the results for the trait based analysis, and will reproduce figure 3. 

### Trait based analysis

library(tidyverse)
library(car)
library(broom)
library(PerformanceAnalytics)
library(ggeffects)
library(effects)
library(ggfortify)
library(patchwork)

analysis_data <- readxl::read_xlsx("Data/Trait_data.xlsx")

Tolerance_values <- read_csv("Data/Results/Urban_Tolerance_Index.csv")

#join estimates to trait data 
trait_data <- Tolerance_values %>% 
  select(species, estimate) %>% 
  left_join(analysis_data, by = "species") %>% 
  mutate(feeding_guild = fct_relevel(feeding_guild, "Generalist", "Bird specialist", "Mammal specialist", "Fish specialist"))

lm.mod_ALAN <- lm(scale(estimate) ~ 
                    scale(body_mass) +
                    scale(habitat_breadth) +
                    scale(nest_substrate_breadth) +
                    feeding_guild +
                    migratory_status,
                  data = trait_data)

summary(lm.mod_ALAN)

# Anova 
Anova(lm.mod_ALAN) %>% 
  kableExtra::kable(digits = 3)
confint(lm.mod_ALAN)

# Model evaluation
mod.eval <- autoplot(lm.mod_ALAN) 
#png("Model.eval.png", units = "cm", res = 300, width = 14, height = 12)
#print(mod.eval)
#dev.off()

plot(allEffects(lm.mod_ALAN), rescale.axis = F)

#### correlations ####

cor <- trait_data %>% 
  select(estimate, body_mass, habitat_breadth, nest_substrate_breadth)

chart.Correlation(cor, histogram = TRUE)

# tidying up results

result_ALAN <- tidy(lm.mod_ALAN) %>%
  mutate(lwr_95_confint=confint(lm.mod_ALAN)[,1]) %>%
  mutate(upr_95_confint=confint(lm.mod_ALAN)[,2]) %>%
  mutate(significance=ifelse(p.value <=0.05, "Significant", "Non-significant")) %>%
  mutate(response=ifelse(.$estimate > 0, "positive", "negative")) %>%  
  kableExtra::kable(digits = 3)

#effects plots of models
#body mass
bm_effect <- ggpredict(lm.mod_ALAN, terms = "body_mass")%>% 
  plot(., residuals = TRUE) + 
  labs(title = "",
       x = "body mass",
       y = "Urban tolerance index") +
  theme_bw()

#nest substrate breadth
nsb_effect <- ggpredict(lm.mod_ALAN, terms = "nest_substrate_breadth")%>% 
  plot(., residuals = TRUE) + 
  labs(title = "",
       x = "nest substrate breadth",
       y = "Urban tolerance index") +
  theme_bw()

#habitat breadth
hb_effect <- ggpredict(lm.mod_ALAN, terms = "habitat_breadth")%>% 
  plot(., residuals = TRUE) + 
  labs(title = "",
       x = "habitat breadth",
       y = "Urban tolerance index") +
  theme_bw()

#feeding guild
fg_effect <- ggpredict(lm.mod_ALAN, terms = "feeding_guild")%>% 
  plot(., residuals = TRUE) + 
  labs(title = "",
       x = "feeding guild",
       y = "Urban tolerance index") +
  theme_bw()

#migratory status
ms_effect <- ggpredict(lm.mod_ALAN, terms = "migratory_status") %>% 
  plot(., residuals = TRUE) + 
  labs(title = "",
       x = "migratory status",
       y = "Urban tolerance index") +
  theme_bw()

#create plot

patchwork <- (bm_effect) / (nsb_effect + hb_effect) / (fg_effect + ms_effect)

trait_plots <- patchwork + plot_annotation(tag_levels = 'A') & theme(plot.tag.position = c(0.1, 0.96))

trait_plots

ggsave(trait_plots, filename = "Figures/Figure_3.png", dpi = 300, height = 24, width = 22, units = "cm")