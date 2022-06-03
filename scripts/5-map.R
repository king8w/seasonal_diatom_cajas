#---------------------------------------------------------------------------------------
# Script: 
# Paper: 
# Author: Benito, X.
# e-mail: xavier.benito.granell@gmail.com
#---------------------------------------------------------------------------------------


#Load packages for functions used
library(maps)
library(raster)
library(cowplot)
library(ggrepel)
library(cividis)
library(ggsn)
library(tidyverse)
library(viridis)
library(rnaturalearth)

world_df <- map_data("world")
interest <- c("Ecuador")
ecuador <- world_df %>% filter(str_detect(region, interest))

world <- ne_countries(scale = "medium", returnclass = "sf")

southamerica <- ggplot(data = world) +
  geom_sf() +
  geom_polygon(data=ecuador, aes(x=long, y=lat, group=group), fill="blue") +
  coord_sf(ylim=c(-55,15), xlim=c(-85,-35)) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank(),
        axis.ticks.x=element_blank())
southamerica

#Plot elevation base map (raster library)
#unzip(zipfile = "data/DEM/dem.zip", exdir = "data/DEM")
DEM <- raster("data/DEM/dem2.bil")
ext<-extent(-82,-50,-30,15)
ext_ecuador<-extent(-81,-77,-4,2)
ext_cajas <- extent(-80, -78,-4,-2)
altmod<-crop(DEM,ext)
altmod_ecuador <- crop(DEM,ext_ecuador)
altmod_cajas <- crop(DEM,ext_cajas)

#convert the raster to points for plotting
map.p <- rasterToPoints(altmod)
map.p_ecuador <- rasterToPoints(altmod_ecuador)
map.p_cajas <- rasterToPoints(altmod_cajas)

#Make the points a dataframe for ggplot
df <- data.frame(map.p)
df_ecuador <- data.frame(map.p_ecuador)
df_cajas <- data.frame(map.p_cajas)

#Make appropriate column headings
colnames(df) <- c("long", "lat", "Elevation")
colnames(df_ecuador) <- c("long", "lat", "Elevation")
colnames(df_cajas) <- c("long", "lat", 'Elevation')

## test Cajas map
plot_cajas <- ggplot(data=df_cajas, aes(y=lat, x=long)) +
  geom_raster(aes(fill=Elevation)) +
  theme_bw() +
  coord_equal()
plot_cajas
  
# Prepare dataset to plot
model_LCBD_var_average <- merge(df1, meta, by=0) %>%
  dplyr::select(-Row.names) %>%
  select(-c("Month", "SubCuenca", "rock_type", "Vertiente")) %>%
  select(-c(4,5)) %>%
  group_by(Lake) %>%
  summarise(across(everything(), list(mean)))  
  

## Cajas map
plt_cajas <- ggplot(data=df_cajas, aes(y=lat, x=long))+
  geom_raster(aes(fill=Elevation)) +
  scale_fill_cividis()+
  geom_point(data=model_LCBD_var_average, aes(x=Longitude_1, y=Latitude_1, size=LCBD_1, colour=wetland_1))+
  scale_size(range = c(1,5), breaks = c(0.035,0.015), "LCBD") +
  geom_text_repel(data=model_LCBD_var_average, aes(Longitude_1, Latitude_1, label = Lake),col="white") +
  xlab("Longitude")+
  ylab("Latitude")+
  labs(colour = "Wetland (sqrt %)") +
  labs(fill="Elevation (m)") +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(size=0.1, linetype = "solid", color = "grey"),
        panel.grid.minor = element_line(size=0.1, linetype = "solid", color = "grey"),
        legend.position = "right")+
  north(df_cajas, symbol = 4) + #add north arrow
  scalebar(df_cajas, dist = 20, dist_unit = "km",
           transform = TRUE, model = "WGS84", st.dist = 0.02, st.size = 2)+
  coord_equal()
plt_cajas  

##
composite <- ggdraw() +
  draw_plot(plt_cajas) +
  #draw_plot(LCBDwetland, x = 0.095, y = 0.10, width = .3, height = .3)
  draw_plot(LCBDrich, x=0.15, y=0.1, width = .3, height = .3) +
  draw_plot(southamerica, x=0.64, y=0.75, width = 0.2, height = 0.2)
composite

#save plot
ggsave("outputs/LCBDmap_Wetland_richness2.png", last_plot(), height = 8, width = 10)


