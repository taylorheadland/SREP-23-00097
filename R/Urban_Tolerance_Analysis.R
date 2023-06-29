## This script runs the function to compute the urban tolerance values for each species in the analysis, and reproduces figure 2. 
library(tidyverse)
library(dggridR)
library(mgcv)
library(broom)

# resolve namespace conflicts
select <- dplyr::select

# function to run the model

analysis_function <- function(species_name){
  
  message(paste0("Analyzing ", species_name))
  
  # ebird data
  ebird <- readRDS(paste0("Data/eBird data/eBird_", species_name, "_zf_cropped.RDS")) %>% 
    mutate(protocol_type = factor(protocol_type, 
                                  levels = c("Stationary", "Traveling", "BirdLife Australia 20min-2ha survey", "BirdLife Australia 500m radius search","BirdLife Australia 5 km radius search")))%>%
    # remove observations with no count
    filter(!is.na(observation_count))
  
  # ALAN data
  ALAN <- readr::read_csv(paste0("Data/ALAN data/", species_name, "_ALAN_location_year.csv"), show_col_types = FALSE) %>% 
    mutate(year = as.integer(year)) %>% 
    select(locality_id, year, ALAN_median)
  
  # combine ebird and ALAN data
  ebird_ALAN <- inner_join(ebird, ALAN, by = c("locality_id", "year"))
  
  # generate hexagonal grid with ~ 5 km between cells
  dggs <- dgconstruct(spacing = 5)
  # get hexagonal cell id and week number for each checklist
  checklist_cell <- ebird_ALAN %>% 
    mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum,
           week = week(observation_date))
  # sample one checklist per grid cell per week
  # sample detection/non-detection independently 
  ebird_ss <- checklist_cell %>% 
    group_by(species_observed, year, week, cell) %>% 
    sample_n(size = 1) %>% 
    ungroup() %>% 
    select(-cell, -week)
  
  # gam formula
  k_time <- 7
  gam_formula <- observation_count ~ 
    ALAN_median +
    s(day_of_year, k = 5) + 
    s(duration_minutes, k = 5) + 
    s(effort_distance_km, k = 5) + 
    s(number_observers, k = 5) + 
    s(longitude, latitude) + 
    s(time_observations_started, bs = "cc", k = k_time)
  
  # explicitly specify where the knots should occur for time_observations_started
  # this ensures that the cyclic spline joins the variable at midnight
  # this won't happen by default if there are no data near midnight
  time_knots <- list(time_observations_started = seq(0, 24, length.out = k_time))
  
  # negative binomial GAM
  m_nb <- bam(gam_formula,
              data = ebird_ss, 
              family = "nb",
              knots = time_knots,
              method = "fREML",
              discrete = TRUE,
              nthreads = 10)
  
  #clean data frame
  summary_df <- broom::tidy(m_nb, parametric = TRUE, conf.int = TRUE) %>% 
    mutate(Response=ifelse(estimate >0, "Urban Tolerant", "Urban Avoidant")) %>% 
    mutate(species = species_name) %>% 
    slice(-1)
  
  return(summary_df)
  
}

### #Birds removed due to large confidence intervals relative to the other species
#The majority of these birds have less than 1000 observations post-sub sampling. Even though the Eastern Grass Owl has quite a small confidence interval making the cutoff 1000 will facilitate for easier interpretation of results. 
#Letter-winged Kite
#Grey Falcon
#Sooty Owl
#Red Goshawk
#Lesser Sooty Owl
#Australian Masked Owl
#Rufous Owl
#Eastern Grass Owl
#Tasmanian Boobook
#Australian Masked Owl
#Black-breasted Buzzard

#names of all mainland Australian raptor species within the analysis
list_names <- c(
  "AustralianHobby", 
  "BarkingOwl", 
  "BlackFalcon", 
  "BlackKite", 
  "BlackShoulderedKite", 
  "BrahminyKite", 
  "BrownFalcon", 
  "BrownGoshawk", 
  "CollaredSparrowhawk", 
  "EasternBarnOwl", 
  "GreyGoshawk", 
  "LittleEagle", 
  "NankeenKestrel", 
  "Osprey", 
  "PacificBaza", 
  "PeregrineFalcon", 
  "PowerfulOwl", 
  "SouthernBoobook", 
  "SpottedHarrier", 
  "SquareTailedKite", 
  "SwampHarrier", 
  "WedgeTailedEagle", 
  "WhistlingKite", 
  "WhiteBelliedSeaEagle") 

# create data frame of 100 iterations of the function per species

df <- as.list(1:100) %>% 
  map_df(~ tibble(iteration = .x,
                  species = list_names)) %>% 
  arrange(species)

results_df <- map_df(df$species, analysis_function, .progress = TRUE) %>% 
  bind_cols(df$iteration)

# clean up the names
results_df <- results_df %>% 
  mutate(species = recode(species, "AustralianHobby" = "Australian Hobby",
                          "BarkingOwl" = "Barking Owl",
                          "BlackFalcon" = "Black Falcon", 
                          "BlackKite" = "Black Kite", 
                          "BlackShoulderedKite" = "Black-shouldered Kite", 
                          "BrahminyKite" = "Brahminy Kite", 
                          "BrownFalcon" = "Brown Falcon", 
                          "BrownGoshawk" = "Brown Goshawk", 
                          "CollaredSparrowhawk" = "Collared Sparrowhawk", 
                          "EasternBarnOwl" = "Eastern Barn Owl",
                          "GreyGoshawk" = "Grey Goshawk",
                          "LittleEagle" = "Little Eagle", 
                          "NankeenKestrel" = "Nankeen Kestrel", 
                          "Osprey" = "Eastern Osprey", 
                          "PacificBaza" = "Pacific Baza", 
                          "PeregrineFalcon" = "Peregrine Falcon", 
                          "PowerfulOwl" = "Powerful Owl", 
                          "SouthernBoobook" = "Southern Boobook", 
                          "SpottedHarrier" = "Spotted Harrier", 
                          "SquareTailedKite" = "Square-tailed Kite", 
                          "SwampHarrier" = "Swamp Harrier",
                          "WedgeTailedEagle" = "Wedge-tailed Eagle", 
                          "WhistlingKite" = "Whistling Kite", 
                          "WhiteBelliedSeaEagle" = "White-bellied Sea-Eagle")) %>% 
  rename("iteration" = ...10)

#group by and summarise index
results_clean <-  results_df %>% 
  group_by(species) %>% 
  summarise(estimate = mean(estimate),
            std.error = mean(std.error),
            statistic = mean(statistic),
            p.value = mean(p.value),
            conf.low = mean(conf.low),
            conf.high = mean(conf.high)) %>% 
  mutate(Response=ifelse(estimate >0, "Urban Tolerant", "Urban Avoidant"),
         species = fct_reorder(species, estimate)) %>% 
  arrange(desc(estimate)) %>% 
  ungroup()

# Save urban tolerance index
write_csv(results_clean, "Data/Results/Urban_Tolerance_Index.csv")

# Plot uban tolerance index
Urban_Tolerance <- 
  results_clean %>%
  mutate(species = fct_reorder(species, estimate)) %>% 
  arrange(desc(estimate)) %>%  
  ggplot(aes(x=estimate, y=species)) +
  geom_vline(xintercept = 0, size=0.25, color = "#AAAAAA") +
  geom_hline(aes(yintercept = species), size=0.05, color = "#AAAAAA") + 
  geom_point(aes(colour = Response)) +
  scale_colour_manual(values = c('Urban Tolerant' = "#0072B2", 'Urban Avoidant' = "#D55E00"), limits = c("Urban Tolerant", "Urban Avoidant"))+
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, height = 0.5, colour = Response))+
  ylab("")+
  xlab("Relationship to urbanisation")+
  scale_x_continuous(breaks = scales::breaks_pretty())+
  theme_classic() +
  theme(legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.position = c(0.85,0.2),
        panel.border = element_rect(fill=NA, colour="black"))


ggsave(Urban_Tolerance, filename = "Figures/Figure_2.png", dpi = 300, units = "cm")