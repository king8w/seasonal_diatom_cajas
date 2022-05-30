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

world <- map_data("world")

southamerica <- ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="lightgrey") +
  theme(legend.position = "right")+
  coord_map("albers", parameters = c(-100, -100),  ylim=c(-30,15), xlim=c(-82,-50)) +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw()
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
  

## Cajas map
plt_cajas <- ggplot(data=df_cajas, aes(y=lat, x=long))+
  geom_raster(aes(fill=Elevation)) +
  scale_fill_cividis()+
  geom_point(data=model_LCBD_var_average, aes(x=Longitude, y=Latitude, size=LCBD, colour=wetland))+
  scale_size(range = c(1,5), breaks = c(0.035,0.015), "LCBD") +
  geom_text_repel(data=model_LCBD_var_average, aes(Longitude, Latitude, label = Lake),col="white") +
  xlab("Longitude")+
  ylab("Latitude")+
  labs(colour = "sqrt wetland (%)") +
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
  draw_plot(LCBDwetland, x = 0.095, y = 0.10, width = .3, height = .3)
  #draw_plot(LCBDrich, x = 0.40, y = 0.10, width = .3, height = .3) 

composite

#save plot
ggsave("outputs/LCBDmap_Wetland_richness.png", last_plot(), height = 8, width = 10)


