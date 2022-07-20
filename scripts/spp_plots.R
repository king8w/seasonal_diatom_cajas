### Create some plots of species distributions across environmental gradients 

# diatoms are counts, species present in more than 20 samples
diat <- read.csv("data/Diatoms_S_2019.csv", row.names = 1, sep=";")
#ecogroups <- read.csv("data/EcoGroups.csv", row.names = 1, sep = ",")
#diat <- diat[,-ncol(diat)] #last column is NAs

# Read in geographical coordinates lakes
spatial_var <- read.csv("data/Spatial_Cajas2019.csv", row.names = 1)
meta <- read.csv("data/metamonth.csv", row.names = 1, sep=";")

# Read most parsimonius CCA variables
model_var <- read.csv("outputs/model_var_v2.csv", row.names=1)

##Select most abundant species across samples
abund <- apply(diat, 2, max)
diat <- diat[, abund>20] # present in >20 samples

# Merge diatom data the most parsimonious variables selected by CCA
diat_env <- merge(diat, model_var, by="row.names")
#diat_env <- merge(ecogroups, model_var, by="row.names")

row.names(diat_env) <- diat_env$Row.names

# Merge diatom data lake metadata
diat_env <- merge(diat_env, meta, by="row.names")
diat_env <- diat_env[,-1]

# Pass row.names 
row.names(diat_env) <- diat_env$Row.names
diat_env <- merge(diat_env, spatial_var, by="row.names")
diat_env <- diat_env[,-1]

# make it long
diatom_env_long <- diat_env %>%
  gather(key = taxa, value = count, -names(model_var), 
         -Row.names, -names(meta), -Latitude, -Longitude) %>%
  filter(str_detect(taxa, "Discostella|Aulacoseira"))

# calculate relative abundance to plot seasonal variation of species by altitude
diat_spp_month_ra <- diatom_env_long %>%
  group_by(taxa, Month, Lake, mix_event,Ca) %>% #here group by variables of interest
  summarise(count = sum(count)) %>%
  #filter(!count == "0" ) %>% #this is to remove empty samples (rows)
  ungroup() %>%
  group_by(Month, Lake) %>%
  mutate(relative_abundance_percent = count / sum(count) * 100) %>%
  ungroup()  
  
# Bar plot
spp.plot <- ggplot(diat_spp_month_ra, aes(fill = taxa, y = relative_abundance_percent, x=fct_reorder2(Lake,relative_abundance_percent,mix_event))) +
  #geom_point(aes(size=Altitude))+  #coord_flip() +
  geom_bar(position = "fill", stat="identity")+  #coord_flip() +
  facet_wrap(Month~., scales = "free") +
  #scale_fill_viridis_d()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Lake (arranged by Altitude (m)")+ylab("Proportion (%)")
spp.plot

# Change species labels
diat_spp_month_ra$taxa <- factor(diat_spp_month_ra$taxa, levels = c("Aulacoseira.alpigena", "Aulacoseira.distans.septentrionalis",
                                                     "Discostella.stelligera"),
                  labels = c("Aulacoseira alpigena", "Aulacoseira distans var. septentrionalis", "Discostella stelligera")
)

# Point-sized plot -- Figure 6 of the manuscript
spp.plot <- ggplot(diat_spp_month_ra, aes(y = taxa, x=fct_reorder2(Lake,relative_abundance_percent,Ca), colour=(10^Ca))) +
  geom_point(aes(size=relative_abundance_percent))+  #coord_flip() +
  facet_wrap(Month~., scales = "free") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 11),
        legend.position = "right")+
  scale_y_discrete(labels=expression(italic("Aulacoseira alpigena"),italic("Aulacoseira distans var. septentrionalis"),
                                     italic("Discostella stelligera")))+
  labs(colour=expression(paste('Ca '^'2+')), size="Relative abundance (%)") +
  xlab("Lakes (arranged by Ca2+)")+ylab("Species") 
  
spp.plot

# Point-sized plot -- Supplementary figure 3
spp.plot <- ggplot(diat_spp_month_ra, aes(y = taxa, x=fct_reorder2(Lake,relative_abundance_percent,mix_event), colour=(10^Ca))) +
  geom_point(aes(size=relative_abundance_percent))+  #coord_flip() +
  facet_wrap(Month~., scales = "free") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 11),
        legend.position = "right")+
  scale_y_discrete(labels=expression(italic("Aulacoseira alpigena"),italic("Aulacoseira distans var. septentrionalis"),
                                     italic("Discostella stelligera")))+
  labs(colour=expression(paste('Ca '^'2+')), size="Relative abundance (%)") +
  xlab("Lakes (arranged by mix events)")+ylab("Species") 

spp.plot

# Save the plot
ggsave("outputs/Discostella_Aulacoseira_Ca_bymixing.png", plot=last_plot(), height=8, width=12,units="in",
        dpi = 400)
