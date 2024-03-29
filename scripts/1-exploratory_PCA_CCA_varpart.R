#---------------------------------------------------------------------------------------
# Script: Exploratory analysis of environmental variables (PCA) and diatom-environment (CCA)
# Paper: Space, not time, drive is driving contemporary planktic diatom species composition in tropical Andean lakes
# Author: Benito, X.
# e-mail: xavier.benito.granell@gmail.com
#---------------------------------------------------------------------------------------

#clear workspace
rm(list=ls(all=TRUE))
dev.off()

#unload all loaded packages
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# load libraries for functions used
library(psych) #allow to calculate correlation index for dataframes and KMO test
library(vegan)
library(adespatial)
library(usdm)
library(tidyverse)
library(ggcorrplot)

#Read lake metadata for grouping
meta <- read.csv("data/metamonth.csv", row.names = 1, sep=";")
head(meta)
str(meta)

#read in environmental data (original analysis)
env <- read.csv("data/monthlyENV.csv", row.names = 1)
row.names(env)[11] <- "CJ-001-S"
id_env <- row.names(env)
names(env)
head(env)

#read in lake physics data
physics <- read.csv("data/lake_physics_data.csv", sep=";", row.names = 1)
physics$erosion <- as.numeric(gsub(",", ".", gsub("\\.", "", physics$erosion))) #replace commas with dots for decimals
physics$water_bodies <- as.numeric(gsub(",", ".", gsub("\\.", "", physics$water_bodies))) #replace commas with dots for decimals
physics$erosion_prop <- as.numeric(gsub(",", ".", gsub("\\.", "", physics$erosion_prop))) #replace commas with dots for decimals
physics$fetch <- as.numeric(gsub(",", ".", gsub("\\.", "", physics$fetch))) #replace commas with dots for decimals
row.names(physics) <- id_env #this is to fix row.names so all datatables have a common id

# subset environmental variables into chemical and physical
chemical_var <- env[,names(env) %in% c("ChlA_a", "Ca", "Mg", "Na", "K", "Alkalinity", "SO4", "Si", "Fe", "Color", "pH", "Cond", "DO", "waterT")]
chemical_var[,c("secchi_m", "TP", "DOC", "TOC","Fe","waterT")] <- physics[,names(physics) %in% c("secchi_m", "TP", "DOC", "TOC","Fe","waterT")]
row.names(chemical_var) <- id_env #this is to fix row.names so all datatables have a common id

#read in catchment variables data and select those of interest
catchment_var <- read.csv("data/monthlyENV.csv", row.names = 1) %>%
  select(Altitude,Heat,Zmax_m,CA_m2,WRT,Pajonal,PajonalRoca,Roca)

catchment_var$wetland <- physics$wetland
catchment_var$erosion <- physics$erosion
catchment_var$erosion_prop <- physics$erosion_prop
catchment_var$water_bodies <- physics$water_bodies
names(catchment_var)

# Assign variables to the group that belongs to
physics[,c("Heat","Zmax_m", "WRT")] <- catchment_var[,c("Heat","Zmax_m", "WRT")]
physics <- physics[,!names(physics) %in% c("Fe", "TP", "DOC", "TOC", "wetland", "waterT","erosion","water_bodies","erosion_prop","rock_type")]
catchment_var <- catchment_var[,!names(catchment_var) %in% c("Heat","Zmax_m","WRT")] %>%
  mutate(bare_rock=Roca+erosion_prop) #create new variable)
chemical_var <- chemical_var[,!names(chemical_var) %in% c("secchi_m")]

#Read in the phytoplankton richness (non-diatoms)
phyto_richn <- read.csv("data/phyto_richness.csv", row.names=1)
colnames(phyto_richn) <- c("phyto_richness")

# Read in the fDOM variable
fDOM <- read.csv("data/fDOM.csv", row.names=1)

# Combine full_physics and chemical data for data exploration
full_env <- cbind(chemical_var, physics, catchment_var, phyto_richn, fDOM)
names(full_env)
head(full_env)

#Colinearity Panel Function
#panel correlation plots to assess data distribution
panel.hist <- function(x, ...) {     
  usr <- par("usr"); on.exit(par(usr))     
  par(usr = c(usr[1:2], 0, 1.5) )     
  h <- hist(x, plot = FALSE)     
  breaks <- h$breaks; nB <- length(breaks)     
  y <- h$counts; y <- y/max(y)     
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...) 
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {     
  usr <- par("usr"); on.exit(par(usr))     
  par(usr = c(0, 1, 0, 1))     
  r <- abs(cor(x, y, use = "complete"))   
  txt <- format(c(r, 0.123456789), digits = digits)[1]     
  txt <- paste0(prefix, txt)     
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)     
  text(0.5, 0.5, txt, cex = cex.cor * r) }


#Check explanatory variable dataset colinearity
pairs(full_env, diag.panel = panel.hist, upper.panel = panel.smooth, lower.panel = panel.cor, gap = 0, cex.labels = 1, cex=1.5, font.labels = 1)

#transform variables to meet assumptions of homogenity of variances
env_trans <- transform(full_env, Alkalinity=log10(Alkalinity+0.25), Altitude=log10(Altitude+0.25), 
                    CA_m2=log10(CA_m2+0.25), WRT=log10(WRT+0.25), K=log10(K+0.25), 
                    Na=log10(Na+0.25), Mg=log10(Mg+0.25), Ca=log10(Ca+0.25),
                    PajonalRoca=sqrt(PajonalRoca), Pajonal=sqrt(Pajonal), wetland=sqrt(wetland),
                    erosion=log10(erosion+0.25), water_bodies=sqrt(water_bodies),
                    erosion_prop=sqrt(erosion_prop),
                    bare_rock=sqrt(bare_rock),
                    Roca=sqrt(Roca), DO=log10(DO+0.25), Heat=log10(Heat+0.25),
                    Cond=log10(Cond+0.25), Color=log10(Color+0.25), Si=log10(Si+0.25), 
                    SO4=log10(SO4+0.25), ChlA_a=log10(ChlA_a+0.25), TP=log10(TP+0.5),
                    secchi_m=log10(secchi_m+0.25), Fe=log10(Fe+0.25), TOC=log10(TOC+0.25), DOC=log10(DOC+0.25),
                    Zmax_m=log10(Zmax_m+0.25),
                    lake_catch_ratio=log10(lake_catch_ratio+0.25), 
                    catch_volume_ratio=log10(catch_volume_ratio+0.50),
                    mix_event=log10(mix_event+0.25), waterT=log10(waterT+0.25),
                    phyto_richness=sqrt(phyto_richness),
                    fetch=log10(fetch+0.25),
                    fDOM=log10(fDOM+0.25))

#Check explanatory variable dataset colinearity
pairs(env_trans, diag.panel = panel.hist, upper.panel = panel.smooth, lower.panel = panel.cor, gap = 0, cex.labels = 1, cex=1.5, font.labels = 1)

#Check adequacy of PCA ordination
source("scripts/pcor.test.R")

# set PCA df to analyze
PCA_data <- env_trans

# drop  chemical variables that do not have monthly observations (and others)
PCA_data_monthly <- env_trans[,!names(env_trans) %in% c("Heat", "lenght_depth_ratio", "catch_volume_ratio", "WRT",
                                                       "PajonalRoca", "Roca", "Pajonal","erosion","TP","TOC","DOC",
                                                       "phyto_richness", "Color", "fetch", "erosion_prop")]

PCA_data <- PCA_data_monthly #monthly

# Average lake variables by Month factor
PCA_data_yearly$month <- meta$Month
PCA_data_yearly$lake <- meta$Lake

PCA_data_yearly <- PCA_data_yearly %>%
  gather(var, value, -month, -lake) %>%
  group_by(lake, var) %>% 
  summarise(average=mean(value, na.rm=T)) %>%
  spread(var, average) %>%
  ungroup()

# set dataset for PCA yearly
PCA_data <- PCA_data_yearly #yearly
  row.names(PCA_data) <- PCA_data$lake
  PCA_data <- PCA_data[,-1]

# check adequacy of the PCA
Cor.matrix <- NULL
Cor.matrix <- matrix(data = NA, nrow = ncol(PCA_data), ncol = ncol(PCA_data), byrow = FALSE, dimnames = NULL)
colnames(Cor.matrix) <- colnames(PCA_data)

Partial.cor.matrix <- NULL
Partial.cor.matrix <- matrix(data = NA, nrow = ncol(PCA_data), ncol = ncol(PCA_data), byrow = FALSE, dimnames = NULL)
colnames(Partial.cor.matrix) <- colnames(PCA_data)

#Compute Correlation matrix
cor.r <- corr.test(PCA_data, adjust = "none", method = "spearman")$r
cor.p <- corr.test(PCA_data, adjust = "none", method = "spearman")$p
Cor.matrix <- cor.r

for(y1.col in 1:ncol(PCA_data)){
  for(y2.col in 1:ncol(PCA_data)){
    Partial.cor.matrix[y1.col,y2.col] <- pcor.test(PCA_data[,y1.col], PCA_data[,y2.col], PCA_data[,-c(y1.col, y2.col)], use = "mat", method = "spearman", na.rm = F)$estimate
  }
}

#KMO sampling adequacy (KMO > 0.5 means adequate)
KMO(cor.r) 

#Test Barlett (p sign means adequate)
cortest.bartlett(Cor.matrix)

## Run Principal Cmponent Analysis
mod_pca <- rda(PCA_data, scale=TRUE)
summary(mod_pca)

ev <- as.vector(eigenvals(mod_pca, model = "unconstrained")) #extract eigenvalues for then broken stick

evplot <- function(ev) {
  # Broken stick model (MacArthur 1957) Author: Francois Gillet, 25 August 2012
  n <- length(ev)
  bsm <- data.frame(j=seq(1:n), p=0)
  bsm$p[1] <- 1/n
  for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n + 1 - i))
  bsm$p <- 100*bsm$p/n
  # Plot eigenvalues and % of variation for each axis
  op <- par(mfrow=c(2,1))
  barplot(ev, main="Eigenvalues", col="bisque", las=2)
  abline(h=mean(ev), col="red")
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, 
          main="% variation", col=c("bisque",2), las=2)
  legend("topright", c("% eigenvalue", "Broken stick model"), 
         pch=15, col=c("bisque",2), bty="n")
  par(op)
}

br <- evplot(ev) #the first two axes are the only ones eplaining more variability than expected by chance

# How much variability each CCA axis explain
axis.expl <- function(mod, axes = 1:2) {
  if(is.null(mod$CCA)) {
    sapply(axes, function(i) {
      100*mod$CA$eig[i]/mod$tot.chi
    })
  } else {
    sapply(axes, function(i) {
      100*mod$CCA$eig[i]/mod$tot.chi
    })
  }
}

(labs <- axis.expl(mod_pca))

# general plot
plot(mod_pca, scaling=3)

#Factor scores (samples)
#Create data frame with factor scores, month and lake groupings
PCA.scores <- data.frame(PCA1=scores(mod_pca, display = "sites")[,1], 
                         PCA2=scores(mod_pca, display = "sites")[,2])

PCA.scores <- merge(PCA.scores, meta, by="row.names")
row.names(PCA.scores) <- PCA.scores$Row.names
colnames(PCA.scores) <- c("Row.names", "PCA1","PCA2","month","lake","vertiente","basin","geology")
PCA.scores <- PCA.scores[,-1]

#custom plot
png("outputs/PCA_monthly_Cajas_diatoms_v6_geology.png", width=10, height=8, units="in", res=300)
win.metafile("outputs/PCA_monthly_Cajas_diatoms_v7_fDOM.wmf", width=10, height=8, res=300)

par(mfrow=c(1,2))
par(mar=c(5,4,3,3)) #sets the bottom, left, top and right margins respectively of the plot region in number of lines of text. 

plot(PCA.scores$PCA1, PCA.scores$PCA2, type = "n", xlab=paste("PCA1","(",round(labs[1],2),"%",")"), ylab=paste("PCA2","(",round(labs[2],2),"%",")"))
title("Sites")
abline(h=0, col="grey")
abline(v=0, col="grey")

# manual coding for geological units
points(PCA.scores[(PCA.scores$geology=="andesite" & PCA.scores$month=="June"), 1:2], col="forestgreen", pch=17)
points(PCA.scores[(PCA.scores$geology=="andesite" & PCA.scores$month=="September"), 1:2], col="forestgreen", pch=15)
points(PCA.scores[(PCA.scores$geology=="andesite" & PCA.scores$month=="December"), 1:2],  col="forestgreen", pch=18)
points(PCA.scores[(PCA.scores$geology=="andesite" & PCA.scores$month=="February"), 1:2],  col="forestgreen", pch=19)

points(PCA.scores[(PCA.scores$geology=="rhyolite" & PCA.scores$month=="June"), 1:2], col="blue", pch=17)
points(PCA.scores[(PCA.scores$geology=="rhyolite" & PCA.scores$month=="September"), 1:2], col="blue", pch=15)
points(PCA.scores[(PCA.scores$geology=="rhyolite" & PCA.scores$month=="December"), 1:2],  col="blue", pch=18)
points(PCA.scores[(PCA.scores$geology=="rhyolite" & PCA.scores$month=="February"), 1:2],  col="blue", pch=19)

points(PCA.scores[(PCA.scores$geology=="dacite" & PCA.scores$month=="June"), 1:2], col="orange", pch=17)
points(PCA.scores[(PCA.scores$geology=="dacite" & PCA.scores$month=="September"), 1:2], col="orange", pch=15)
points(PCA.scores[(PCA.scores$geology=="dacite" & PCA.scores$month=="December"), 1:2],  col="orange", pch=18)
points(PCA.scores[(PCA.scores$geology=="dacite" & PCA.scores$month=="February"), 1:2],  col="orange", pch=19)

#legend for monhtly dataset
legend("topleft",c("June", "September", "December", "February",
                      "Andesite", "Rhyolite", "Dacite")
       , cex=.8, pch=c(17,15,18,19,NA,NA,NA,NA),
       text.col=c("black", "black", "black", "black",
                  "forestgreen", "blue", "orange"), 
       ncol = 2, xpd = TRUE)

# Add lake name lakes
text(PCA.scores[,1:2], labels=PCA.scores$lake, pos = 1, cex = 0.6, offset = 0.3)

#Variables
comp1 <- as.numeric(scores(mod_pca, display = "species")[,1])
comp2 <- as.numeric(scores(mod_pca, display = "species")[,2])

#Labels
labels <- colnames(PCA_data)
plot(comp1, comp2, pch=16, col="black", xlab=paste("PCA1","(",round(labs[1],2),"%",")"), ylab="", scaling=3)
title("Variables")
abline(h=0, col="grey")
abline(v=0, col="grey")
text(comp1, comp2, labels = labels, pos = 1, cex = 0.5, offset = 0.2)

#save the plot
dev.off()

#Correlations between PCA components and environmental variables
PCA_result <- cbind(PCA.scores, PCA_data)

# Check if a new variable "month" has significnat correlations with PCA axes
PCA_result <- PCA_result %>% mutate(month=recode(month,
                                                 "June"=1,
                                                 "September"=2,
                                                 "December"=3,
                                                 "February"=4)) %>%
  mutate(basin=recode(basin,"Tomebamba"=1, "Canar"=2, "Balao"=3, "Yanuncay"=4)) %>%
  mutate(geology=recode(geology, "andesite"=1, "rhyolite"=2, "dacite"=3))


cor(PCA_result[,1:2], PCA_result[,c(3,6,7,8:ncol(PCA_result))], use = "complete")
cor.r <- corr.test(PCA_result[,1:2], PCA_result[,c(3,6,7,8:ncol(PCA_result))], method = "spearman")$r
cor.p <- corr.test(PCA_result[,1:2], PCA_result[,c(3,6,7,8:ncol(PCA_result))], method = "spearman")$p

cor_PCA_results <- rbind(cor.r, cor.p)
row.names(cor_PCA_results) <- c("PCA1r", "PCA2r", "PCA1p", "PCA2p")
write.table(cor_PCA_results, "outputs/cor_PCA_results.txt")

heatmap(cor.r, Colv = NA, Rowv = NA, scale="column")

### Diatom analysis
#Load the data for the diatoms (response variable)
diat <- read.csv("data/Diatoms_S_2019.csv", row.names = 1, sep=";")
row.names(diat)

##Transform to relative abundance
total <- apply(diat, 1, sum)
diat <- diat/total*100

## species richness
rich <- apply(diat>0, 1, sum)

##Remove rare species
abund <- apply(diat, 2, max)
n.occur <- apply(diat>0, 2, sum)
diat <- diat[, n.occur>1 & abund>2] #more than 2% of RA and present in >1 sample

##  Constrained multivariate ordination
### Hellinger transformation
diat <- decostand(diat, method="hellinger")

par(mfrow=c(5,4))
ccaResult <- list()
for (i in 1:length(PCA_data)) {
  mod <- cca(diat~PCA_data[,i], na=na.omit, scale=TRUE)
  ccaResult$mod[[i]] <- mod
  #plot(ccaResult$mod[[i]], main=colnames(PCA_data[i]))
  ccaResult$anova[[i]] <- anova(ccaResult$mod[[i]])
  ccaResult$ratio[[i]] <- mod[["CCA"]][["eig"]][["CCA1"]]/mod[["CA"]][["eig"]][["CA1"]] #ratio from the CCAconstrained  axis (k1) to the first eigenvalue from the unconstrained axis (k2)
  #print(anova(ccaResult$mod[[i]]))
}

names(ccaResult$mod) <- colnames(PCA_data)
names(ccaResult$anova) <- colnames(PCA_data)
names(ccaResult$ratio) <- colnames(PCA_data)

#extract anova results
anova.cca <- do.call(rbind, ccaResult$anov) 
anova.cca <- filter(anova.cca, !grepl("Residual",row.names(anova.cca)))
vec <- anova.cca[,4]
names(vec) <- row.names(anova.cca)
barplot(vec, col = "grey60", cex.names=0.5, las=1, horiz=T)
abline(v=0.05)

dev.off()
win.metafile("outputs/CCA1_CA1ratio.wmf", width=10, height=8, res=300)
png("outputs/CCA1_CA1ratio.png", width = 12, height = 8)
vec <- t(t(sort(data.frame(ccaResult$ratio), decreasing = FALSE)))
barplot(vec, col = "grey60", cex.names=0.7, las=1, xlab="λ1/λ2 ratio", horiz=T)
dev.off()

##multivariate cca with variables having individual statistical significant effects on diatom data (monthly model)
#exclude variables
var <- c("DO","pH","SO4","Na","ChlA_a","Si", "erosion_prop") 
model_var <- PCA_data[,!names(PCA_data) %in% var]

vifstep(model_var) #check out for multicollinearity

#subset variables and save most parsimonius model variables selected by individual CCAs for later
model_var <- model_var[,!names(model_var) %in% c("CA_m2", "Cond" ,"water_bodies", "fDOM")]

write.csv(model_var, "outputs/model_var_v2.csv")  

# Perform CCA
mod_cca <- cca(diat~., data=model_var, scale=TRUE)
plot(mod_cca, scaling = 3)

#Plot eigenvalues and percentages of variation of an ordination object
anova(mod_cca, by="axis")
ev <- as.vector(eigenvals(mod_cca, model = "constrained")) #extract eigenvalues for then broken stick

evplot <- function(ev) {
  # Broken stick model (MacArthur 1957) Author: Francois Gillet, 25 August 2012
  n <- length(ev)
  bsm <- data.frame(j=seq(1:n), p=0)
  bsm$p[1] <- 1/n
  for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n + 1 - i))
  bsm$p <- 100*bsm$p/n
  # Plot eigenvalues and % of variation for each axis
  op <- par(mfrow=c(2,1))
  barplot(ev, main="Eigenvalues", col="bisque", las=2)
  abline(h=mean(ev), col="red")
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, 
          main="% variation", col=c("bisque",2), las=2)
  legend("topright", c("% eigenvalue", "Broken stick model"), 
         pch=15, col=c("bisque",2), bty="n")
  par(op)
}

br <- evplot(ev) #the first two axes are the only ones eplaining more variability than expected by chance

# How much variability each CCA axis explain
axis.expl <- function(mod, axes = 1:2) {
  if(is.null(mod$CCA)) {
    sapply(axes, function(i) {
      100*mod$CA$eig[i]/mod$tot.chi
    })
  } else {
    sapply(axes, function(i) {
      100*mod$CCA$eig[i]/mod$tot.chi
    })
  }
}

(labs <- axis.expl(mod_cca))

#extract site scores for plotting
scrs <- mod_cca$CCA$wa[,1:2]

#merge by row.names
CCA.scores <- merge(scrs, meta, by=0)

#plot CCA
plot(mod_cca, type="n", xlab=paste("CCA1","(",round(labs[1],2),"%",")"), ylab=paste("CCA2","(",round(labs[2],2),"%",")"))

#manual coding for months
points(CCA.scores[(CCA.scores$Month=="June"), 2:3], col="black", pch=24, bg="grey")
points(CCA.scores[(CCA.scores$Month=="September"), 2:3], col="black", pch=22, bg="black")
points(CCA.scores[(CCA.scores$Month=="December"), 2:3],  col="black", pch=23, bg="grey")
points(CCA.scores[(CCA.scores$Month=="February"), 2:3],  col="black", pch=8)

text(mod_cca, dis="bp", col="grey")
text(CCA.scores[,2:3], labels = CCA.scores$Lake, pos = 2, cex = 0.8, offset = 0.3)

# add species scores
points(mod_cca$CCA$v[,1:2])
text(mod_cca$CCA$v[,1:2])

#Correlations between CCA components and environmental variables
var_env_CCA <- merge(model_var, scrs, by=0)
row.names(var_env_CCA) <- var_env_CCA$Row.names
CCA_result <- merge(meta, var_env_CCA,by=0)

# Check if a new variable "month" has significnat correlations with PCA axes
CCA_result <- CCA_result %>% mutate(month=recode(Month,
                                                 "June"=1,
                                                 "September"=2,
                                                 "December"=3,
                                                 "February"=4)) %>%
  mutate(geology=recode(rock_type,"andesite"=1,"rhyolite"=2,"dacite"=3))

cor(CCA_result[,c("CCA1", "CCA2")], CCA_result[,c(8:21,23,24)], use = "complete")
cor.p <- corr.test(CCA_result[,c("CCA1", "CCA2")], CCA_result[,c(8:21,23,24)], method = "spearman")$p
cor.r <- corr.test(CCA_result[,c("CCA1", "CCA2")], CCA_result[,c(8:21,23,24)], method = "spearman")$r

cor_CCA_results <- rbind(cor.p, cor.r)
row.names(cor_CCA_results) <- c("PCA1p", "PCA2p", "PCA1r", "PCA2r")

write.table(cor_CCA_results, file = "outputs/CCA_correlations.txt")


# Triplot
dev.off()
png("outputs/CCA_cajas_diatoms_new_v6_geology.png", width=12, height=8, units="in", res=300)
win.metafile("outputs/CCA_cajas_diatoms_v6_geology.wmf", width=10, height=8, res=300)

op <- par(no.readonly = TRUE)
par(fig = c(0, 0.5, 0, 0.95))
plot(mod_cca, type="n", xlab=paste("CCA1","(",round(labs[1],2),"%",")"), ylab=paste("CCA2","(",round(labs[2],2),"%",")"))
abline(h=0, col="grey")
abline(v=0, col="grey")
title("Sites")

## Plot by month and geology
points(CCA.scores[(CCA.scores$rock_type=="andesite" & CCA.scores$Month=="June"), 2:3], col="forestgreen", pch=17)
points(CCA.scores[(CCA.scores$rock_type=="andesite" & CCA.scores$Month=="September"), 2:3], col="forestgreen", pch=15)
points(CCA.scores[(CCA.scores$rock_type=="andesite" & CCA.scores$Month=="December"), 2:3], col="forestgreen", pch=18)
points(CCA.scores[(CCA.scores$rock_type=="andesite" & CCA.scores$Month=="February"), 2:3], col="forestgreen", pch=19)

points(CCA.scores[(CCA.scores$rock_type=="rhyolite" & CCA.scores$Month=="June"), 2:3], col="blue", pch=17)
points(CCA.scores[(CCA.scores$rock_type=="rhyolite" & CCA.scores$Month=="September"), 2:3], col="blue", pch=15)
points(CCA.scores[(CCA.scores$rock_type=="rhyolite" & CCA.scores$Month=="December"), 2:3], col="blue", pch=18)
points(CCA.scores[(CCA.scores$rock_type=="rhyolite" & CCA.scores$Month=="February"), 2:3], col="blue", pch=19)

points(CCA.scores[(CCA.scores$rock_type=="dacite" & CCA.scores$Month=="June"), 2:3], col="orange", pch=17)
points(CCA.scores[(CCA.scores$rock_type=="dacite" & CCA.scores$Month=="September"), 2:3], col="orange", pch=15)
points(CCA.scores[(CCA.scores$rock_type=="dacite" & CCA.scores$Month=="December"), 2:3], col="orange", pch=18)
points(CCA.scores[(CCA.scores$rock_type=="dacite" & CCA.scores$Month=="February"), 2:3], col="orange", pch=19)

legend("bottomleft",c("June", "September", "December", "February",
                      "Andesite", "Rhyolite", "Dacite")
       , cex=.8, pch=c(17,15,18,19,NA,NA,NA,NA),
       text.col=c("black", "black", "black", "black",
                  "forestgreen", "blue", "orange"), 
       ncol = 2, xpd = TRUE)

## add lake names
text(CCA.scores[,2:3], labels = CCA.scores$Lake, pos = 1, cex = 0.5, offset = 0.3)

par(fig = c(0.5, 1, 0, 0.95), new = TRUE)
plot(mod_cca, type="n", xlab=paste("CCA1","(",round(labs[1],2),"%",")"), ylab="")
abline(h=0, col="grey")
abline(v=0, col="grey")
title("Species")
#points(mod_cca$CCA$v[,1:2])

#customize species scores by ecological groups
#from TAXA
eco.groups <- c("benthic", "benthic", "planktic", "planktic", "planktic", "benthic", "benthic", "benthic", "planktic", "benthic",
                "benthic", "benthic", "benthic", "tychoplanktic", "tychoplanktic", "benthic", "benthic", "benthic", "benthic", "benthic",
                "benthic", "benthic", "benthic", "benthic", "benthic", "tychoplanktic", "tychoplanktic", "tychoplanktic", "tychoplanktic",
                "tychoplanktic", "tychoplanktic", "tychoplanktic", "tychoplanktic", "tychoplanktic")

text(scores(mod_cca$CCA$v[,1]), scores(mod_cca$CCA$v[,2]), labels = row.names(mod_cca$CCA$v[,1:2]), pos = 1, cex = 0.5, offset = 0.2)
points(scores(mod_cca$CCA$v[,1]), scores(mod_cca$CCA$v[,2]), pch=as.numeric(as.factor(eco.groups)), cex = 0.8)

legend("bottomleft",c("Benthic", "Planktic", "Tychoplanktic"),
       cex=.7, pch=c(1,2,3),
       ncol = 1, xpd = TRUE)

par(fig = c(0.51, 0.72, 0.7, 0.99), new = TRUE)
par(mar = c(0, 0, 0, 0))
bp <- scores(mod_cca, display = "bp", scaling = 3)
plot(bp, type = "n", axes = FALSE, asp = 1)
u <- par("usr")
rect(u[1], u[3], u[2], u[4], col = "white")
arrows(0, 0, bp[, 1] * 0.8, bp[, 2] * 0.8, length = 0.1)
text(bp[, 1] * 0.9, bp[, 2] * 0.9, cex = 0.8, labels = rownames(bp),xpd = NA)
box()
par(op)
#save plot
dev.off()

### Visually explore species-environment relationships
op <- par(no.readonly = TRUE)
par(mfrow = c(6, 6))
par(mar = c(2.5, 3.5, 1, 0.5))
par(mgp = c(1.5, 0.5, 0))
par(oma = c(0, 0, 3, 0))

var <- c("waterT")
for (i in 1:length(diat)) {
  plot(diat[, i] ~ env_trans[, var], ann = F, cex = 0.6, cex.axis = 0.8,
       tcl = -0.4, las = 1, pch = 19, )
  mtext(colnames(diat)[i], side = 3, line = 0.2, cex = 0.6)
}
mtext("Diatom distribution vs. log 10 Secchi disk (m)", outer = TRUE,
      side = 3, cex = 1.2, line = 1)
par(op)

### Model spp-environment with GAMs
source("scripts/functions/model.gams.R")

# Select variable to model with GAM
var <- c("Fe")
env_var <- model_var[,var]

# Run the function
models.gam(diat,env_var, title = "Diatom distribution vs log10 Iron")


### CCA forward-step model selection
# Environmental variables
mod_caa.R2a <- RsquareAdj(mod_cca)$adj.r.squared
env <- PCA_data[,names(PCA_data) %in% c("Ca", "Mg", "K", "Alkalinity", "waterT", "Fe")]

mod_cca.fw.env <- forward.sel(diat, env, adjR2thresh = mod_caa.R2a)
mod_cca.fw.env$R2a.fullmodel <- mod_caa.R2a

env.sel <- env[, which(names(env) %in% mod_cca.fw.env$variables)] 

# Spatial variables
# read mems from 2-MoranEigenvectorMaps.R
spatial <- read.csv("outputs/mems.csv", row.names = 1)

mod_cca.fw.spatial <- forward.sel(diat, spatial, adjR2thresh = mod_caa.R2a)
mod_cca.fw.spatial$R2a.fullmodel <- mod_caa.R2a

spatial.sel <- spatial[, which(names(spatial) %in% mod_cca.fw.spatial$variables)] 

# Physical variables
physical <- PCA_data[,names(PCA_data) %in% c("lake_catch_ratio", "mix_event", "secchi_m", "Zmax_m","Altitude","wetland","bare_rock")]
mod_cca.fw.physical <- forward.sel(diat, physical, adjR2thresh = mod_caa.R2a)
mod_cca.fw.physical$R2a.fullmodel <- mod_caa.R2a

physical.sel <- physical[, which(names(physical) %in% mod_cca.fw.physical$variables)] 

# Save forward-selected variables
mod_cca_fw_full <- rbind(mod_cca.fw.env, mod_cca.fw.spatial, mod_cca.fw.physical)
write.csv(mod_cca_fw_full, "outputs/mod_cca_fw_full.csv")

#Varpart
varpart <- varpart(diat, env.sel, spatial.sel, physical.sel)

plot(varpart, Xnames=c("Environment","Spatial", "Physical"))
#text(locator(), c("Environmental","Spatial", "Physical"), cex=1)

#Test pure effects
anova.cca(rda(diat, env.sel, cbind(spatial.sel, physical.sel)), perm.max = 999) ## test pure environmental (signif 0.05)
anova.cca(rda(diat, spatial.sel, cbind(env.sel, physical.sel)), perm.max = 999) ## test pure spatial (signif 0.005)  
anova.cca(rda(diat, physical.sel, cbind(env.sel, spatial.sel)), perm.max = 999) ## test pure physical effect (signif 0.005)  

