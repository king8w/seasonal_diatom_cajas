#https://stackoverflow.com/questions/24863470/running-multiple-gamm-models-using-for-loop-or-lapply 

d_resp <- d[ c("y", "y1")]
d_pred <- d[, !(colnames(d) %in% c("y", "y1"))]

## create a "matrix" list of dimensions i x j
results_m <- vector("list", length=ncol(d_resp)*ncol(d_pred))
dim(results_m) <- c(ncol(d_resp), ncol(d_pred))

for(i in 1:ncol(d_resp)){
  for(j in 1:ncol(d_pred)){
    results_m[i, j][[1]] <- gamm(d_resp[, i] ~ s(d_pred[, j]))
  }
}

# flatten the "matrix" list
results_l <- do.call("list", results_m)


### Create plots of species distributions across environmental gradients 
# diatoms are counts, species present in more than 20 samples
diat <- read.csv("data/Diatoms_S_2019.csv", row.names = 1, sep=";")
#diat <- diat[,-ncol(diat)] #last column is NAs

# Read in geographical coordinates lakes
spatial_var <- read.csv("data/Spatial_Cajas2019.csv", row.names = 1)
meta <- read.csv("data/metamonth.csv", row.names = 1)

##Select most abundant species across samples
abund <- apply(diat, 2, max)
diat <- diat[, abund>20] # present in >20 samples

# Merge diatom data the most parsimonious variables selected by CCA
diat_env <- merge(diat, full_env, by="row.names")
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
  gather(key = taxa, value = count, -names(model_var), -Row.names, -names(meta), -Latitude, -Longitude) %>%
  filter(str_detect(taxa, "Discostella|Aulacoseira"))

# calculate relative abundance to plot seasonal variation of species by altitude
diat_spp_month_ra <- diatom_env_long %>%
  group_by(taxa, Month, Lake, Ca, Altitude) %>%
  summarise(count = sum(count)) %>%
  #filter(!count == "0" ) %>% #this is to remove empty samples (rows)
  ungroup() %>%
  group_by(Month, Lake) %>%
  mutate(relative_abundance_percent = count / sum(count) * 100) %>%
  ungroup()  
  
# Bar plot
spp.plot <- ggplot(diat_spp_month_ra, aes(fill = taxa, y = relative_abundance_percent, x=fct_reorder2(Lake,relative_abundance_percent,Altitude))) +
  #geom_point(aes(size=Altitude))+  #coord_flip() +
  geom_bar(position = "fill", stat="identity")+  #coord_flip() +
  facet_wrap(Month~., scales = "free") +
  #scale_fill_viridis_d()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Lake (arranged by Altitude (m)")+ylab("Proportion (%)")
spp.plot

# Point-sized plot
spp.plot <- ggplot(diat_spp_month_ra, aes(y = taxa, x=fct_reorder2(Lake,relative_abundance_percent,Altitude), colour=Ca)) +
  geom_point(aes(size=relative_abundance_percent))+  #coord_flip() +
  #geom_bar(position = "fill", stat="identity")+  #coord_flip() +
  facet_wrap(Month~., scales = "free") +
  #scale_fill_viridis_d()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Lake (arranged by Altitude(m))")+ylab("Proportion (%)")
spp.plot
