## This script runs the code to reproduce figure 1. 

library(tidyverse)
library(sf)
library(scico)
library(ozmaps)

# ebird data
ebird <- readRDS("Data/eBird data/ebird_EasternBarnOwl_zf_cropped.RDS")

# ALAN data
ALAN <- readr::read_csv("Data/ALAN data/EasternBarnOwl_ALAN_location_year.csv")

# combine ebird and ALAN data
ebird_ALAN <- inner_join(ebird, ALAN, by = c("locality_id", "year"))

# load map of Australia
Aus <- st_transform(ozmaps::ozmap_states, 7844)

#plot

map.viirs <- 
  ggplot()+
  # plot Aus shapefile
  geom_sf(data = Aus, fill = "black", colour = "white", size = 0.1) +
  # plot bird data
  geom_point(data = ebird_ALAN, 
             mapping = aes(x=longitude, y=latitude, colour= log10(ALAN_median)),
             size = 0.1, 
             alpha = 0.5) +
  # extra aesthetic options
  coord_sf(xlim = c(110, 155), expand = TRUE) +
  labs(color = expression(log(Median~radiance)~" "(nW~cm^{-2}~sr^{-1}))) + 
  scale_color_scico(palette = "vikO") +
  theme_void()+
  theme(axis.text = element_blank(),
        legend.position = c(0.4,0.1),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = NA, colour = NA),
        text = element_text(size = 9)) +
  guides(color = guide_colourbar(barwidth = 10.3, 
                                 barheight = 0.7,
                                 title.position = "top"))

map.viirs


ggsave(map.viirs, filename = "Figures/Figure_1.png", dpi = 300, units = "cm")