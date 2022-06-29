#---------------------------------------------------------------------------------------
# Script: Local Contribution to Beta-Diversity analyses
# Paper: Space, not time, drive is driving contemporary planktic diatom species composition in tropical Andean lakes
# Author: Benito, X.
# e-mail: xavier.benito.granell@gmail.com
#---------------------------------------------------------------------------------------

#clear workspace
rm(list=ls(all=TRUE))
dev.off()

#unload all loaded packages
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

#Load libraries for functions used
library(adespatial)
library(tidyverse)
library(mgcv)
library(ggplot2)
library(ggcorrplot)
library(gratia)
library(vegan)
library(psych)

#Load the data for the diatoms
diat <- read.csv("data/Diatoms_S_2019.csv", sep=";",row.names = 1)
head(diat)

#richness calculation
raremax <- min(rowSums(diat)) 
Srare <- rarefy(diat, raremax)

# Read in geographical coordinates lakes
spatial_var <- read.csv("data/Spatial_Cajas2019.csv", row.names = 1)
meta <- read.csv("data/metamonth.csv", row.names = 1, sep=";")

##Transform to relative abundance
total <- apply(diat, 1, sum)
diat <- diat/total*100

##Remove rare species
abund <- apply(diat, 2, max)
n.occur <- apply(diat>0, 2, sum)
diat <- diat[, n.occur>1 & abund>2] #more than 3% of RA and present in >1 sample

#richness calculation
richness<- apply(diat>0,1,sum)

# Perform LCBD
diatomsLCBD <- beta.div(diat, method = "hellinger", nperm = 999, adj = TRUE, sqrt.D=FALSE)

## plot a map of the LCBD indices
LCBD <- data.frame(diatomsLCBD$LCBD)
LCBD$Longitude <- spatial_var$Longitude
LCBD$Latitude <- spatial_var$Latitude
colnames(LCBD)[1] <- c("LCBD")
LCBD$sign <- data.frame(diatomsLCBD$p.LCBD)
colnames(LCBD[,4]) <- "sign"
LCBD$sign.ad <- data.frame(diatomsLCBD$p.adj)

#LCBD - richness relationship
LCBD$richness <- richness
LCBD$Srare <- Srare
plot(LCBD$richness, LCBD$LCBD)
plot(LCBD$Srare, LCBD$LCBD)

cor.test(LCBD$richness, LCBD$LCBD, method = "spearman")
cor.test(LCBD$Srare, LCBD$LCBD, method = "spearman")

LCBD_df <- data.frame(LCBD)

# Plot LCBD-species richness
LCBDrich <- ggplot(data=LCBD_df, aes(x=richness, y=LCBD))+
  geom_point()+
  geom_smooth(method=lm)+
  xlab("Species richness")+
  theme_classic()
LCBDrich

## Prepare data to fit GAM model on LCBD
# Read most parsimonius CCA variables
model_var <- read.csv("outputs/model_var_v2.csv", row.names=1)
#model_var <- read.csv("outputs/model_var.csv", row.names=1) 

#merge by row.names with model var previously subset from CCA
df1 <- merge(LCBD_df, model_var, by=0)
row.names(df1) <- df1$Row.names
df1 <- df1[,!names(df1) %in% c("Row.names")]

model_LCBD_var <- merge(df1, meta, by=0) %>%
  dplyr::select(-Row.names) %>%
  mutate(month=recode(Month,
                      "June"=1,
                      "September"=2,
                      "December"=3,
                      "February"=4)) %>%
  mutate(basin=recode(SubCuenca,"Tomebamba"=1, "Canar"=2, "Balao"=3, "Yanuncay"=4)) %>%
  mutate(geology=recode(rock_type, "andesite"=1, "rhyolite"=2, "dacite"=3)) %>%
  select(-c("Month", "SubCuenca", "rock_type", "Vertiente", "Lake")) %>%
  select(-c(2,3,4,5)) %>%
  select(-c("month", "basin", "geology"))

# Make correlations with environmental variables
# corr <- cor(model_LCBD_var, method = "spearman")
# corr.test(model_LCBD_var, method = "spearman", ci=TRUE)

plot(model_LCBD_var$Fe, model_LCBD_var$LCBD)

corr <- round(cor(model_LCBD_var), 2)
p.mat <- cor_pmat(model_LCBD_var)
head(corr)

ggcorrplot(corr, hc.order = TRUE, p.mat = p.mat, insig = "blank", type = "lower", lab = TRUE,
           ggtheme = ggplot2::theme_classic(), outline.color = "white")

ggsave("outputs/LCBD_correlations_v2.png", plot=last_plot(), height=8, width=10,units="in",
       dpi = 400)

# make some plots
boxplot(model_LCBD_var$LCBD ~ month, data = model_LCBD_var, ylab="LCBD", main="Month")
boxplot(model_LCBD_var$LCBD ~ geology, data = model_LCBD_var, ylab="LCBD", main="Geology")
boxplot(model_LCBD_var$LCBD ~ basin, data = model_LCBD_var, ylab="LCBD", main="Basin")


model_LCBD_var$Lake <- meta$Lake
model_LCBD_var$Month <- meta$Month

LCBD_lake <- ggplot(model_LCBD_var, aes(x = reorder(Lake, LCBD, FUN = median), y = LCBD)) + 
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic(base_size = 14) +
  xlab("Lake") +
  ylab("LCBD")

LCBD_month <- ggplot(model_LCBD_var, aes(x = Month, y = LCBD)) + 
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic(base_size = 14) +
  xlab("Month") +
  ylab("LCBD")


LCBD_anova <- aov(LCBD ~ Lake, data = model_LCBD_var)
LCBD_anova <- aov(LCBD ~ Month, data = model_LCBD_var)

summary(LCBD_anova)

ggsave("outputs/LCBD_month.png", plot=last_plot(), height=8, width=10,units="in",
       dpi = 400)

dev.off()

#Combine plots
plot_composite <- plot_grid(LCBD_month, LCBD_lake, ncol = 1,  
                            rel_heights = c(1, 1), align="v", labels = "AUTO") +
  theme(plot.margin = unit(c(0.5, -1, 0.5, 0.5), "cm"))
plot_composite

ggsave("outputs/LCBD_month_lake.png", plot=last_plot(), height=8, width=10,units="in",
       dpi = 400)
