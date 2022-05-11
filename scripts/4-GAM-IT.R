## GAM Models of LCBD or species richness against chemical, physical and spatial factors
# contact: xavier.benito.granell@gmail.com

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
richness<- apply(diat>0,1,sum)
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
model_var <- read.csv("outputs/model_var.csv", row.names=1) 

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
corr <- cor(model_LCBD_var, method = "spearman")
corr.test(model_LCBD_var, method = "spearman", ci=TRUE)

plot(model_LCBD_var$Fe, model_LCBD_var$LCBD)

corr <- round(cor(model_LCBD_var), 2)
p.mat <- cor_pmat(model_LCBD_var)
head(corr)

ggcorrplot(corr, hc.order = TRUE, p.mat = p.mat, insig = "blank", type = "lower")

ggsave("outputs/LCBD_correlations.png", plot=last_plot(), height=8, width=10,units="in",
       dpi = 400)

# make some plots
boxplot(model_LCBD_var$LCBD ~ month, data = model_LCBD_var, ylab="LCBD", main="Month")
boxplot(model_LCBD_var$LCBD ~ geology, data = model_LCBD_var, ylab="LCBD", main="Geology")
boxplot(model_LCBD_var$LCBD ~ basin, data = model_LCBD_var, ylab="LCBD", main="Basin")


model_LCBD_var$Lake <- meta$Lake

ggplot(model_LCBD_var, aes(x = reorder(Lake, LCBD, FUN = median), y = LCBD)) + 
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic(base_size = 14) +
  xlab("Lake") +
  ylab("LCBD")

LCBD_anova <- aov(LCBD ~ Lake, data = model_LCBD_var)
summary(LCBD_anova)

ggsave("outputs/LCBD_lakes.png", plot=last_plot(), height=8, width=10,units="in",
       dpi = 400)


dev.off()

# Linear model + spatial smooths of LCBD
set.seed(10) #set a seed so this is repeatable
mod1 <- gam(LCBD ~ s(Latitude, Longitude, k=5) + Ca + Mg + Alkalinity + SO4 + Si + Altitude + secchi_m + Fe + erosion_prop + lake_catch_ratio + mix_event +
              waterT,
            method="REML", data=model_LCBD_var, select=TRUE, 
            family=gaussian(link = "log"))

visreg(mod1,'lake_catch_ratio', gg=T)+
  ylab("n partial residuals")+
  theme_classic()

summary(mod1)
appraise(mod1)
gam.check(mod1)

# Linear model + spatial smooths of species richness
mod2 <- gam(richness ~ s(Latitude, Longitude, k=5) +Ca + Mg + Alkalinity + SO4 + Si + Altitude + secchi_m + Fe + erosion_prop + lake_catch_ratio + mix_event +
              waterT,
            method="REML", data=model_LCBD_var, select=TRUE, 
            family=gaussian(link = "log"))

summary(mod2)
appraise(mod2)
gam.check(mod2)


#### Summary of models
## Pierre's idea!!
summary(mod1)$p.table

modPred_res <-as.data.frame(rbind(
  round(summary(mod2)$p.table, digits = 4)[-1,],
  round(summary(mod1)$p.table, digits=4)[-1,]))

modPred_res$predictor <- rep(c(names(model_var))) #doesnt work
modPred_res$expl.var <- rep(c("LCBD", "richness"), each=12)
colnames(modPred_res)<- c("coefficient", "SE", "z_value", "P_value","predictor","expl.var")

modPred_res$coefficient<- as.numeric(as.character(modPred_res$coefficient))
modPred_res$SE<- as.numeric(as.character(modPred_res$SE))
modPred_res$sig <- "ns"
modPred_res$sig[modPred_res$P_value < 0.1] <- "P<0.1"
modPred_res$sig[modPred_res$P_value < 0.05] <- "P<0.05"
modPred_res$sig[modPred_res$P_value < 0.01] <- "P<0.01"

str(modPred_res)

# Plot
pltLCBD <- ggplot(modPred_res, aes(x=predictor, y=coefficient, col=sig))+
  geom_pointrange(aes(ymin=coefficient-(1.96*SE), ymax=coefficient+(1.96*SE)))+
  geom_hline(yintercept = 0, linetype=2)+
  facet_wrap(~expl.var,nrow=1, scales='free_y')+
  scale_color_manual(values=c('red', 'green', 'orange', 'blue'))+
  #scale_color_brewer(palette='')+
  theme_bw()+
  theme(axis.text.x=element_text(angle=60, hjust=1), axis.title.x=element_blank())
pltLCBD

# ggsave("figures/Fig2_gamLCBD_modplot.png", plot=pltLCBD, height=8, width=10,units="in",
#        dpi = 400)

############################################
#### GAM-IT-single model inference. It doesn't average 
## Fisher et al. 2020 (https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.4134)
library(MuMIn)
library(FSSgam)

# Prepare dataset 
#Load the data for the diatoms (response variable)
diat <- read.csv("data/Diatoms_S_2019.csv", row.names = 1, sep = ";")
#diat <- diat[,-ncol(diat)] #last column is NAs

# Read in geographical coordinates lakes
spatial_var <- read.csv("data/Spatial_Cajas2019.csv", row.names = 1)
meta <- read.csv("data/metamonth.csv", row.names = 1)

##Select most abundant species across samples
abund <- apply(diat, 2, max)
diat <- diat[, abund>20] # present in >20 samples

# Read most parsimonius CCA variables
#model_var <- read.csv("outputs/model_var_interm.csv", row.names=1)
model_var <- read.csv("outputs/model_var_v2.csv", row.names=1)

#cat.preds <- c("SubCuenca", "Lake")
cont.preds <- c("Ca", "Mg", "SO4", "Alkalinity", "Si", "Altitude", "erosion_prop",
                "Fe", "secchi_m", "lake_catch_ratio", "mix_event", "waterT")
cont.preds <- names(model_var)
null.cars <- c("Longitude", "Latitude")
month <- rep(c(1,2,3,4), 1:nrow(model_LCBD_var), each = 1, len = nrow(model_LCBD_var))

# Run base model
set.seed(10) #set a seed so this is repeatable
mod1 <- gam(LCBD ~ s(month, bs="cc", k=4),
            method="REML", data=model_LCBD_var, select=TRUE, 
            family=gaussian)

model.set <- generate.model.set(use.dat=model_LCBD_var,
                                test.fit=mod1,
                                cyclic.vars=month,
                                null.terms="s(Longitude, Latitude, k=5)",
                                pred.vars.cont=cont.preds)

out.list <- fit.model.set(model.set)

# examine the output
names(out.list)
out.list$failed.models
length(out.list$success.models)
mod.table=out.list$mod.data.out
mod.table=mod.table[order(mod.table$AICc),]
head(mod.table)

# examine the list of failed models
length(out.list$failed.models)
length(out.list$success.models)

# look at the model selection table
mod.table=out.list$mod.data.out
mod.table=mod.table[order(mod.table$AICc),]
head(mod.table)
#write.csv(mod.table[,-2],"modfits.csv")

barplot(out.list$variable.importance$bic$variable.weights.raw,las=2,
        ylab="Relative variable importance")


# extract the best model
mod.table=mod.table[order(mod.table$AIC),]
head(mod.table)

best.model=out.list$success.models[[as.character(mod.table$modname[1])]]
plot(best.model$gam,all.terms=T,pages=1)

gam.check(best.model$gam)
summary(best.model)

# check the predictor correlation matrix
model.set$predictor.correlations

# now run the same thing using the non.linear correlation matrix
model.set=generate.model.set(use.dat=model_LCBD_var,
                             test.fit=mod1,
                             pred.vars.cont=cont.preds,
                             #null.terms="s(Longitude, Latitude, k=5)",
                             non.linear.correlations=TRUE)

model.set$predictor.correlations
out.list=fit.model.set(model.set)
mod.table=out.list$mod.data.out
mod.table=mod.table[order(mod.table$AICc),]
head(mod.table)

#extract coefficients from the model
pred <- "Ca"
# https://github.com/samclifford/mgcv.helper/tree/master/R 
#library(mgcv.helper)
library(dplyr)
bond <- lapply(out.list$success.models, function(z) c(z$coeff, confint(z, parm=pred)))



# Plot all the best models
var.imp <- out.list$variable.importance$aic$variable.weights.raw
all.less.2AICc <- mod.table[which(mod.table$delta.AICc<4),]
top.all <- all.less.2AICc

# plot the all best models
par(mfrow=c(5,4))
par(mar=c(1,2,2,1))
for(r in 1:nrow(all.less.2AICc)){
  best.model.name <- as.character(all.less.2AICc$modname[r])
  best.model <- out.list$success.models[[best.model.name]]
  if(best.model.name!="null"){
    plot(best.model, all.terms=T, residuals=T,pch=16)
    mtext(side=3,text="resp.vars[i]",outer=T)}
}


### Now run the GAM-IT model across species (individual responses)
#### GAM-IT-single model inference. It doesn't average 
## Fisher et al. 2020 (https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.4134)
library(MuMIn)
library(FSSgam)

# Prepare data
# diatoms are counts, species present in more than 20 samples
diat <- read.csv("data/Diatoms_S_2019.csv", row.names = 1, sep = ";")

# Read most parsimonius CCA variables
#model_var <- read.csv("outputs/model_var_interm.csv", row.names=1)
model_var <- read.csv("outputs/model_var_v2.csv", row.names=1)

# Read in geographical coordinates lakes
spatial_var <- read.csv("data/Spatial_Cajas2019.csv", row.names = 1)
meta <- read.csv("data/metamonth.csv", row.names = 1, sep = ";")

##Select most abundant species across samples
abund <- apply(diat, 2, max)
diat <- diat[, abund>20] # present in >20 samples

# Merge diatom data the most parsimonious variables selected by CCA
diat <- diat[order(row.names(diat)),]
model_var <- model_var[order(row.names(model_var)),]
meta <- meta[order(row.names(meta)),]
spatial_var <- spatial_var[order(row.names(spatial_var)),]

diat_env <- cbind(diat, model_var, meta, spatial_var)

# make it long
diatom_env_long <- diat_env %>%
    gather(key = taxa, value = count, -names(model_var), -names(meta), -Latitude, -Longitude)

# calculate relative abundance to plot seasonal variation of species
diat_spp_month_ra <- diatom_env_long %>%
  group_by(taxa, Month, Lake) %>%
  summarise(count = sum(count)) %>%
  #filter(!count == "0" ) %>% #this is to remove empty samples (rows)
  ungroup() %>%
  group_by(Month, Lake) %>%
  mutate(relative_abundance_percent = count / sum(count) * 100) %>%
  ungroup()

# Plot
spp.plot <- ggplot(diat_spp_month_ra, aes(fill = taxa, y = relative_abundance_percent, x=Month)) +
  geom_bar(position = "fill", stat="identity")+  #coord_flip() +
  facet_wrap(Lake~., scales = "free") +
  #scale_fill_viridis_d()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
spp.plot

# ggsave("outputs/spp_months_plot.png", 
#        plot=spp.plot, height=8, width=10,units="in",
#        dpi = 300)

# set predictors
# cont.preds <- c("Ca", "Mg", "SO4", "Alkalinity", "Si", "Altitude", "Fe",
#                 "erosion_prop", "secchi_m", "lake_catch_ratio", "mix_event", "waterT")
cont.preds <- names(model_var)
resp.vars <- names(diat_env[,names(diat_env) %in% names(diat)])
null.vars <- c("Longitude", "Latitude")

# month <- rep(c("December","February","January","September"), 1:nrow(diat_env), each = 1, len = nrow(diat_env))
# month <- rep(c(1,2,3,4), 1:nrow(diat_env), each = 1, len = nrow(diat_env))

diat_env <- cbind(diat, model_var, meta, spatial_var)


#setwd("C:/Users/xbenito/Documents/R/Cajas/outputs") #Set wd for storing the results
setwd("C:/Users/xbenito/Documents/R/seasonal_diatom_cajas/outputs") 

# get rid of NA's and unused columns
use.dat <- na.omit(diat_env[,c(null.vars,cont.preds,resp.vars)])

out.all=list()
var.imp=list()
fss.all=list()
top.all=list()
best.model=list()

i=1
pdf(file="mod_fits_all.pdf",onefile=T)
for(i in 1:length(resp.vars)){
  
  use.dat$response=use.dat[,resp.vars[i]]

  #test.fit model for the particular diatom spp i
  # model1 <- uGamm(use.dat$response ~ s(month, bs="cc", k=4),
  #                 family=gaussian(),
  #                 data=use.dat, lme4 = TRUE)
  
  model1 <- uGamm(use.dat$response ~ Ca + Mg + K + Alkalinity + Fe + waterT + lake_catch_ratio + bare_rock +
                  + mix_event + secchi_m + Zmax_m + Altitude + wetland + s(Latitude, Longitude, k=5),
                  family=gaussian(),
                  data=use.dat, lme4 = TRUE)

  
  model.set=generate.model.set(use.dat=use.dat,
                               test.fit=model1,
                               pred.vars.cont=cont.preds,
                               #cyclic.vars=month,k=4,
                               null.terms = "s(Latitude, Longitude, k=5)")
  
  out.list=fit.model.set(model.set)
  fss.all=c(fss.all,list(out.list))
  mod.table=out.list$mod.data.out
  
  mod.table=mod.table[order(mod.table$AICc),]
  out.i=mod.table
  out.all=c(out.all,list(out.i))
  var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw))
  all.less.2AICc=mod.table[which(mod.table$delta.AICc<2),]
  #out.list.top=success.spp$all.less.2AICc$modname #take less2AICc model names
  top.all=c(top.all,list(all.less.2AICc))
  
  #plot all best models
  par(oma=c(1,1,4,1))
  for(r in 1:nrow(all.less.2AICc)){
  best.model.name <- as.character(all.less.2AICc$modname[r])
  best.model[r] <- fss.all[[resp.vars[i]]]$success.models[[best.model.name]]
  plot.best.model <- out.list$success.models[[best.model.name]]
    if(plot.best.model!="null"){
      plot.gam(plot.best.model$gam, all.terms=T, pages=1,residuals=T,pch=16)
      mtext(side=3,text=resp.vars[i],outer=T)}
  }
}
dev.off()

names(out.all)=resp.vars
names(var.imp)=resp.vars
names(top.all)=resp.vars
names(fss.all)=resp.vars
names(best.model)=resp.vars

all.mod.fits=do.call("rbind",out.all)
all.var.imp=do.call("rbind",var.imp)
top.mod.fits=do.call("rbind",top.all)

# examine the list of failed models
length(out.list$failed.models)
length(out.list$success.models)

# look at the model selection table
mod.table=out.list$mod.data.out
mod.table=mod.table[order(mod.table$AICc),]
head(mod.table)
#write.csv(mod.table[,-2],"modfits.csv")

barplot(out.list$variable.importance$bic$variable.weights.raw,las=2,
        ylab="Relative variable importance")

# extract the best model
mod.table=mod.table[order(mod.table$AIC),]
head(mod.table)

best.model=out.list$success.models[[as.character(mod.table$modname[1])]]
plot(best.model$gam,all.terms=T,pages=1)

gam.check(best.model$gam)
summary(best.model)


## Extract model averaged coefficients and point-wise confidence intervals
# ##All
# Result <- list()
# for (i in seq_along(fss.all)) {
#   Result[[i]] <- lapply(fss.all[[i]]$success.models,
#                         function(z) c(z$gam$coeff, confint(z$gam, parm=z$terms)))
# }
# names(Result) <- names(fss.all)

# Subset the top models for each species
# top.mod.fits.achaffine <- top.mod.fits %>%
#   filter(str_detect(row.names(top.mod.fits), "affine"))
# 
# top.mod.fits.admi <- top.mod.fits %>%
#   filter(str_detect(row.names(top.mod.fits), "minutissimum"))
# 
# top.mod.fits.alpigena <- top.mod.fits %>%
#   filter(str_detect(row.names(top.mod.fits), "alpigena"))
# 
# # Extract GAM coefficients and point-wise CI of the list
# top.mod.fits.achaffine <- lapply(fss.all[["Achnanthidium.affine"]][["success.models"]][top.mod.fits.achaffine$modname], 
#                                  function(z) c(z$gam$coeff, confint(z$gam, parm=z$terms, n=100)))
# 
# top.mod.fits.admi <- lapply(fss.all[["Achnanthidium.minutissimum"]][["success.models"]][top.mod.fits.admi$modname], 
#                                  function(z) c(z$gam$coeff, confint(z$gam, parm=z$terms, n=100)))
# 
# 
# top.mod.fits.alpigena <- lapply(fss.all$Aulacoseira.alpigena[["success.models"]][top.mod.fits.alpigena$modname], 
#                             function(z) c(z$coeff, confint(z, parm=z$terms, n=100)))



## 
require(car)
require(doBy)
require(gplots)
require(RColorBrewer)

pdf(file="var_importance_v2.pdf",height=5,width=7,pointsize=10)
heatmap.2(all.var.imp,notecex=0.4,  dendrogram ="none",
          col=colorRampPalette(c("yellow","orange","red"))(30),
          trace="none",key.title = "",keysize=2,
          notecol="black",key=T,
          sepcolor = "black",margins=c(12,14), lhei=c(3,10),lwid=c(3,10),
          Rowv=FALSE,Colv=FALSE)
dev.off()

write.csv(all.mod.fits,"all_model_fits_v2.csv")
write.csv(top.mod.fits,"top_model_fits_v2.csv")
write.csv(all.var.imp,"all.var.imp_v2.csv")

# save results
saveRDS(fss.all, "spp_GAM_IT_v2.rds")

# Make a nicer plot of variance importance scores
dat.taxa <-read.csv("outputs/all.var.imp_v2.csv", row.names = 1)
all.var.imp <- read.csv("outputs/all_model_fits_v2.csv", row.names = 1)

all.var.imp <- data.frame(dat.taxa)
all.var.imp$taxa <- row.names(all.var.imp)

data.plt <- all.var.imp %>% gather(key=predictor, value=importance, -taxa)
str(data.plt)

# colour ramps-
#re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)
re <- colorRampPalette(brewer.pal(8, "Blues"))(25)

# Labels
legend_title<-"Variable Importance"

# Annotations of the top model for each species
## V3 (bare rock model)
# add up \U2191 and down arrows \U2193; \U2193 straight
dat.taxa.label<-data.plt %>%
  mutate(label=NA) %>%
  mutate(label=ifelse(predictor=="Fe"&taxa=="Achnanthidium.affine","X \U2191",
                      ifelse(predictor=="Zmax_m"&taxa=="Achnanthidium.affine","X \U2193",label)))%>%
  mutate(label=ifelse(predictor=="Fe"&taxa=="Achnanthidium.minutissimum","X \U2191",
                      ifelse(predictor=="wetland"&taxa=="Achnanthidium.minutissimum","X \U2193",
                             ifelse(predictor=="Zmax_m"&taxa=="Achnanthidium.minutissimum","X \U2193",label))))%>%
  mutate(label=ifelse(predictor=="Ca"&taxa=="Aulacoseira.alpigena","X \U2193",
                      ifelse(predictor=="Fe"&taxa=="Aulacoseira.alpigena","X \U2191",
                             ifelse(predictor=="waterT"&taxa=="Aulacoseira.alpigena","X \U2193",label)))) %>%
  mutate(label=ifelse(predictor=="bare_rock"&taxa=="Aulacoseira.distans.septentrionalis","X \U2191",
                      ifelse(predictor=="Zmax_m"&taxa=="Aulacoseira.distans.septentrionalis","X \U2191",
                             ifelse(predictor=="Mg"&taxa=="Aulacoseira.distans.septentrionalis","X \U2193",label)))) %>%
  mutate(label=ifelse(predictor=="Fe"&taxa=="Diatoma.tenuis","X \U2193",
                      ifelse(predictor=="K"&taxa=="Diatoma.tenuis", "X \U2191",
                             ifelse(predictor=="secchi_m"&taxa=="Diatoma.tenuis", "X \U2191",label)))) %>%
  mutate(label=ifelse(predictor=="mix_event"&taxa=="Discostella.stelligera","X \U2191",
                      ifelse(predictor=="Fe"&taxa=="Discostella.stelligera","X \U2191",
                             ifelse(predictor=="Zmax_m"&taxa=="Discostella.stelligera","X \U2193 \U2191",label)))) %>%
  mutate(label=ifelse(predictor=="Mg"&taxa=="Fragilaria.tenera","X \U2193",
                             ifelse(predictor=="wetland"&taxa=="Fragilaria.tenera","X \U2191",
                                    ifelse(predictor=="Zmax_m"&taxa=="Fragilaria.tenera","X \U2193 \U2191",label)))) %>%
  mutate(label=ifelse(predictor=="K"&taxa=="Navicula.notha","X \U2191 \U2193",
                      ifelse(predictor=="mix_event"&taxa=="Navicula.notha","X \U2193",label))) %>%
  mutate(label=ifelse(predictor=="bare_rock"&taxa=="Navicula.radiosa","X \U2191 \U2193",
                      ifelse(predictor=="Zmax_m"&taxa=="Navicula.radiosa","X \U2191",label))) %>%
  mutate(label=ifelse(predictor=="lake_catch_ratio"&taxa=="Nitzschia.Atucyacu.1","X \U2193",label)) %>%
  mutate(label=ifelse(predictor=="Ca"&taxa=="Nitzschia.cf.oberheimiana","X \U2191 \U2193",
                      ifelse(predictor=="mix_event"&taxa=="Nitzschia.cf.oberheimiana","X \U2193 \U2191",label))) %>%
  mutate(label=ifelse(predictor=="Zmax_m"&taxa=="Pseudostaurosira.laucensis","X \U2191",
                      ifelse(predictor=="waterT"&taxa=="Pseudostaurosira.laucensis","X \U2193",
                             ifelse(predictor=="Fe"&taxa=="Pseudostaurosira.laucensis","X \U2193",label)))) %>%
  mutate(label=ifelse(predictor=="Mg"&taxa=="Pseudostaurosira.santaremensis","X \U2193",
                      ifelse(predictor=="Zmax_m"&taxa=="Pseudostaurosira.santaremensis","X \U2193",
                             ifelse(predictor=="mix_event"&taxa=="Pseudostaurosira.santaremensis","X \U2191",label)))) %>%
  mutate(label=ifelse(predictor=="Ca"&taxa=="Staurosirella.pinnata","X \U2193",
                      ifelse(predictor=="waterT"&taxa=="Staurosirella.pinnata","X \U2193",label))) %>%
  mutate(label=ifelse(predictor=="mix_event"&taxa=="Staurosirella.sp.2", "X \U2193",
                      ifelse(predictor=="Zmax_m"&taxa=="Staurosirella.sp.2","X \U2193",label))) %>%
  mutate(label=ifelse(predictor=="bare_rock"&taxa=="Tabellaria.fenestrata","X \U2191 \U2193",
                      ifelse(predictor=="Ca"&taxa=="Tabellaria.fenestrata","X \U2191",label))) %>%
  mutate(label=ifelse(predictor=="mix_event"&taxa=="Tabellaria.flocculosa","X \U2193",
                      ifelse(predictor=="Zmax_m"&taxa=="Tabellaria.flocculosa","X \U2193",label))) %>%
  mutate(label=ifelse(predictor=="Altitude"&taxa=="Ulnaria.delicatissima","X \U2191",
                      ifelse(predictor=="Zmax_m"&taxa=="Ulnaria.delicatissima","X \U2193 \U2191",label))) 

# Plotting theme
Theme1 <-
  theme( 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill="white"),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.text = element_text(size=8),
    legend.title = element_text(size=8, face="bold"),
    legend.position = "top",
    legend.direction="horizontal",
    text=element_text(size=10),
    #strip.text.y = element_text(size = 10,angle = 0),
    #axis.title.x=element_text(vjust=0.3, size=10),
    #axis.title.y=element_text(vjust=0.6, angle=90, size=10),
    axis.text.x=element_text(size=10,angle=45, hjust=1,vjust=1),
    axis.text.y=element_text(size=10,face="italic"),
    axis.line.x=element_line(colour="black", size=0.5,linetype='solid'),
    axis.line.y=element_line(colour="black", size=0.5,linetype='solid'),
    strip.background = element_blank())

# Plot gg.importance.score
gg.importance.scores <- ggplot(dat.taxa.label, aes(x=predictor,y=taxa,fill=importance))+
  geom_tile(show.legend=T) +
  scale_fill_gradientn(legend_title,colours=c("white", re), na.value = "grey98",
                       limits = c(0, max(dat.taxa.label$importance)))+
  scale_x_discrete(limits=c("Ca",
                            "Mg",
                            "K",
                            "Fe",
                            "waterT",
                            "lake_catch_ratio",
                            "mix_event",
                            "secchi_m",
                            "Zmax_m",
                            "Altitude",
                            "wetland",
                            "bare_rock"),
                   labels=c(
                     "Calcium (µeq /L)",
                     "Magnesium (µeq /L)",
                     "Potassium (µeq /L)",
                     "Iron (µeq /L)",
                     "Water temperature (º)",
                     "Lake/ Catchment area ratio",
                     "Mix event (degrees)",
                     "Secchi depth (m)",
                     "Lake depth max (m)",
                     "Altitude (m)",
                     "Wetland (%)",
                     "Bare rock (%)"))+
  # scale_y_discrete(limits = c(""),
  #                  labels=c(""))+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  Theme1+
  geom_text(aes(label=label))
gg.importance.scores

# Save the plot
ggsave("outputs/GAM_IT_varimportance_v3.png", plot=gg.importance.scores, height=8, width=10,units="in",
       dpi = 400)


# # Annotations of the top model for each species
# ## V2
# # add up \U2191 and down arrows \U2193; \U2193 straight
# dat.taxa.label<-data.plt %>%
#   mutate(label=NA) %>%
#   mutate(label=ifelse(predictor=="erosion_prop"&taxa=="Achnanthidium.affine","X \U2193",
#                       ifelse(predictor=="SO4"&taxa=="Achnanthidium.affine","X \U2191",label)))%>%
#   mutate(label=ifelse(predictor=="Alkalinity"&taxa=="Achnanthidium.minutissimum","X \U2193",
#                       ifelse(predictor=="SO4"&taxa=="Achnanthidium.minutissimum","X \U2193",ifelse(predictor=="erosion_prop"&taxa=="Achnanthidium.minutissimum","X \U2192",label))))%>%
#   mutate(label=ifelse(predictor=="Ca"&taxa=="Aulacoseira.alpigena","X \U2193",
#                       ifelse(predictor=="Fe"&taxa=="Aulacoseira.alpigena","X \U2191 \U2193",
#                              ifelse(predictor=="secchi_m"&taxa=="Aulacoseira.alpigena","X \U2193",label)))) %>%
#   mutate(label=ifelse(predictor=="erosion_prop"&taxa=="Aulacoseira.distans.septentrionalis","X \U2191",
#                       ifelse(predictor=="lake_catch_ratio"&taxa=="Aulacoseira.distans.septentrionalis","X \U2193",
#                              ifelse(predictor=="Fe"&taxa=="Aulacoseira.distans.septentrionalis","X \U2191",label)))) %>%
#   mutate(label=ifelse(predictor=="erosion_prop"&taxa=="Diatoma.tenuis","X \U2191",label)) %>%
#   mutate(label=ifelse(predictor=="mix_event"&taxa=="Discostella.stelligera","X \U2191",
#                       ifelse(predictor=="Fe"&taxa=="Discostella.stelligera","X \U2191",
#                              ifelse(predictor=="erosion_prop"&taxa=="Discostella.stelligera","X \U2191",label)))) %>%
#   mutate(label=ifelse(predictor=="erosion_prop"&taxa=="Fragilaria.tenera","X \U2191",
#                       ifelse(predictor=="Fe"&taxa=="Fragilaria.tenera","X \U2193",
#                              ifelse(predictor=="Mg"&taxa=="Fragilaria.tenera","X \U2193",label)))) %>%
#   mutate(label=ifelse(predictor=="Ca"&taxa=="Navicula.notha","X \U2193",
#                       ifelse(predictor=="Fe"&taxa=="Navicula.notha","X \U2193",
#                              ifelse(predictor=="waterT"&taxa=="Navicula.notha","X \U2193",label)))) %>%
#   mutate(label=ifelse(predictor=="mix_event"&taxa=="Navicula.radiosa","X \U2193",
#                       ifelse(predictor=="Fe"&taxa=="Navicula.radiosa","X \U2193",label))) %>%
#   mutate(label=ifelse(predictor=="lake_catch_ratio"&taxa=="Nitzschia.Atucyacu.1","X \U2193",
#                       ifelse(predictor=="Si"&taxa=="Nitzschia.Atucyacu.1","X \U2191",label))) %>%
#   mutate(label=ifelse(predictor=="erosion_prop"&taxa=="Nitzschia.cf.oberheimiana","X \U2193",
#                       ifelse(predictor=="Fe"&taxa=="Nitzschia.cf.oberheimiana","X \U2193 \U2191",
#                              ifelse(predictor=="mix_event"&taxa=="Nitzschia.cf.oberheimiana","X \U2193 \U2191",label)))) %>%
#   mutate(label=ifelse(predictor=="erosion_prop"&taxa=="Pseudostaurosira.laucensis","X \U2193 \U2191",
#                       ifelse(predictor=="waterT"&taxa=="Pseudostaurosira.laucensis","X \U2193",
#                              ifelse(predictor=="Fe"&taxa=="Pseudostaurosira.laucensis","X \U2191",label)))) %>%
#   mutate(label=ifelse(predictor=="Mg"&taxa=="Pseudostaurosira.santaremensis","X \U2193",
#                       ifelse(predictor=="Si"&taxa=="Pseudostaurosira.santaremensis","X \U2191",
#                              ifelse(predictor=="mix_event"&taxa=="Pseudostaurosira.santaremensis","X \U2191",label)))) %>%
#   mutate(label=ifelse(predictor=="Ca"&taxa=="Staurosirella.pinnata","X \U2193",
#                       ifelse(predictor=="waterT"&taxa=="Staurosirella.pinnata","X \U2193",label))) %>%
#   mutate(label=ifelse(predictor=="mix_event"&taxa=="Staurosirella.sp.2", "X \U2193",label)) %>%
#   mutate(label=ifelse(predictor=="Si"&taxa=="Tabellaria.fenestrata","X \U2191",label)) %>%
#   mutate(label=ifelse(predictor=="erosion_prop"&taxa=="Tabellaria.flocculosa","X \U2191",
#                       ifelse(predictor=="secchi_m"&taxa=="Tabellaria.flocculosa","X \U2193",label))) %>%
#   mutate(label=ifelse(predictor=="Altitude"&taxa=="Ulnaria.delicatissima","X \U2191",
#                       ifelse(predictor=="erosion_prop"&taxa=="Ulnaria.delicatissima","X \U2191",label))) 
# 
# # Plotting theme
# Theme1 <-
#   theme( 
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.background = element_rect(fill="white"),
#     legend.key = element_blank(), # switch off the rectangle around symbols in the legend
#     legend.text = element_text(size=8),
#     legend.title = element_text(size=8, face="bold"),
#     legend.position = "top",
#     legend.direction="horizontal",
#     text=element_text(size=10),
#     #strip.text.y = element_text(size = 10,angle = 0),
#     #axis.title.x=element_text(vjust=0.3, size=10),
#     #axis.title.y=element_text(vjust=0.6, angle=90, size=10),
#     axis.text.x=element_text(size=10,angle=45, hjust=1,vjust=1),
#     axis.text.y=element_text(size=10,face="italic"),
#     axis.line.x=element_line(colour="black", size=0.5,linetype='solid'),
#     axis.line.y=element_line(colour="black", size=0.5,linetype='solid'),
#     strip.background = element_blank())
# 
# # Plot gg.importance.score
# gg.importance.scores <- ggplot(dat.taxa.label, aes(x=predictor,y=taxa,fill=importance))+
#   geom_tile(show.legend=T) +
#   scale_fill_gradientn(legend_title,colours=c("white", re), na.value = "grey98",
#                        limits = c(0, max(dat.taxa.label$importance)))+
#   scale_x_discrete(limits=c("Ca",
#                             "Mg",
#                             "Fe",
#                             "Alkalinity",
#                             "SO4",
#                             "Si",
#                             "waterT",
#                             "Altitude",
#                             "secchi_m",
#                             "erosion_prop",
#                             "lake_catch_ratio",
#                             "mix_event"),
#                    labels=c(
#                      "Calcium (µeq /L)",
#                      "Magnesium (µeq /L)",
#                      "Iron (µeq /L)",
#                      "Alkalinity (µmol /L)",
#                      "Sulfate (µeq /L)",
#                      "Silica (µeq /L)",
#                      "Water temperature (º)",
#                      "Altitude (m)",
#                      "Secchi depth (m)",
#                      "Erosion (%)",
#                      "Lake/ Catchment area ratio",
#                      "Mix event (degrees)"))+
#   # scale_y_discrete(limits = c(""),
#   #                  labels=c(""))+
#   xlab(NULL)+
#   ylab(NULL)+
#   theme_classic()+
#   Theme1+
#   geom_text(aes(label=label))
# gg.importance.scores
# 
# # Save the plot
# ggsave("outputs/GAM_IT_varimportance_v2.png", plot=gg.importance.scores, height=8, width=10,units="in",
#        dpi = 400)
# 

# 
# # Annotations of the top model for each species
# ## V1
# dat.taxa.label<-data.plt %>%
#   mutate(label=NA) %>%
#   mutate(label=ifelse(predictor=="SO4"&taxa=="Achnanthidium.affine", "X", label)) %>%
#   mutate(label=ifelse(predictor=="Alkalinity"&taxa=="Achnanthidium.minutissimum","X",
#                       ifelse(predictor=="SO4"&taxa=="Achnanthidium.minutissimum","X",ifelse(predictor=="Zmax_m"&taxa=="Achnanthidium.minutissimum","X",label))))%>%
#   mutate(label=ifelse(predictor=="Ca"&taxa=="Aulacoseira.alpigena","X",
#                       ifelse(predictor=="Si"&taxa=="Aulacoseira.alpigena","X",
#                       ifelse(predictor=="waterT"&taxa=="Aulacoseira.alpigena","X",label)))) %>%
#   mutate(label=ifelse(predictor=="mix_event"&taxa=="Aulacoseira.distans.septentrionalis","X",
#                       ifelse(predictor=="Zmax_m"&taxa=="Aulacoseira.distans.septentrionalis","X",label))) %>%
#   mutate(label=ifelse(predictor=="lake_catch_ratio"&taxa=="Diatoma.tenuis","X",
#                       ifelse(predictor=="Si"&taxa=="Diatoma.tenuis","X",label))) %>%
#   mutate(label=ifelse(predictor=="mix_event"&taxa=="Discostella.stelligera","X",
#                       ifelse(predictor=="Si"&taxa=="Discostella.stelligera","X",label))) %>%
#   mutate(label=ifelse(predictor=="wetland"&taxa=="Fragilaria.tenera","X",
#                       ifelse(predictor=="Zmax_m"&taxa=="Fragilaria.tenera","X",label))) %>%
#   mutate(label=ifelse(predictor=="K"&taxa=="Navicula.notha","X",
#                       ifelse(predictor=="lake_catch_ratio"&taxa=="Navicula.notha","X",
#                       ifelse(predictor=="Si"&taxa=="Navicula.notha","X",label)))) %>%
#   mutate(label=ifelse(predictor=="Zmax_m"&taxa=="Navicula.radiosa","X",
#                       ifelse(predictor=="wetland"&taxa=="Navicula.radiosa","X",label))) %>%
#   mutate(label=ifelse(predictor=="lake_catch_ratio"&taxa=="Nitzschia.Atucyacu.1","X",
#                       ifelse(predictor=="Si"&taxa=="Nitzschia.Atucyacu.1","X",label))) %>%
#   mutate(label=ifelse(predictor=="Alkalinity"&taxa=="Nitzschia.cf.oberheimiana","X",
#                       ifelse(predictor=="SO4"&taxa=="Nitzschia.cf.oberheimiana","X",
#                       ifelse(predictor=="Zmax_m"&taxa=="Nitzschia.cf.oberheimiana","X",label)))) %>%
#   mutate(label=ifelse(predictor=="Alkalinity"&taxa=="Pseudostaurosira.laucensis","X",
#                       ifelse(predictor=="Zmax_m"&taxa=="Pseudostaurosira.laucensis","X",
#                       ifelse(predictor=="waterT"&taxa=="Pseudostaurosira.laucensis","X",label)))) %>%
#   mutate(label=ifelse(predictor=="Ca"&taxa=="Pseudostaurosira.santamarensis","X",
#                       ifelse(predictor=="waterT"&taxa=="Pseudostaurosira.santamarensis","X",
#                       ifelse(predictor=="Si"&taxa=="Pseudostaurosira.laucensis","X",label)))) %>%
#   mutate(label=ifelse(predictor=="Ca"&taxa=="Staurosirella_pinnata","X",
#                       ifelse(predictor=="waterT"&taxa=="Staurosirella_pinnata","X",label))) %>%
#   mutate(label=ifelse(predictor=="mix_event"&taxa=="Staurosirella.sp.2", "X",
#                       ifelse(predictor=="Zmax_m"&taxa=="Staurosirella.sp.2","X",label))) %>%
#   mutate(label=ifelse(predictor=="Si"&taxa=="Tabellaria.fenestrata","X", label)) %>%
#   mutate(label=ifelse(predictor=="mix_event"&taxa=="Tabellaria.flocculosa","X",
#                       ifelse(predictor=="Zmax_m"&taxa=="Tabellaria.flocculosa","X",label))) %>%
#   mutate(label=ifelse(predictor=="Altitude"&taxa=="Ulnaria.delicatissima","X",
#                       ifelse(predictor=="Zmax_m"&taxa=="Ulnaria.delicatissima","X",label))) 
#   
# # Plotting theme
# Theme1 <-
#   theme( 
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.background = element_rect(fill="white"),
#     legend.key = element_blank(), # switch off the rectangle around symbols in the legend
#     legend.text = element_text(size=8),
#     legend.title = element_text(size=8, face="bold"),
#     legend.position = "top",
#     legend.direction="horizontal",
#     text=element_text(size=10),
#     #strip.text.y = element_text(size = 10,angle = 0),
#     #axis.title.x=element_text(vjust=0.3, size=10),
#     #axis.title.y=element_text(vjust=0.6, angle=90, size=10),
#     axis.text.x=element_text(size=10,angle=45, hjust=1,vjust=1),
#     axis.text.y=element_text(size=10,face="italic"),
#     axis.line.x=element_line(colour="black", size=0.5,linetype='solid'),
#     axis.line.y=element_line(colour="black", size=0.5,linetype='solid'),
#     strip.background = element_blank())
# 
# # Plot gg.importance.score
# gg.importance.scores <- ggplot(dat.taxa.label, aes(x=predictor,y=taxa,fill=importance))+
#   geom_tile(show.legend=T) +
#    scale_fill_gradientn(legend_title,colours=c("white", re), na.value = "grey98",
#                         limits = c(0, max(dat.taxa.label$importance)))+
#   scale_x_discrete(limits=c("Ca",
#                             "SO4",
#                             "Alkalinity",
#                             "Si",
#                             "Altitude",
#                             "Zmax_m",
#                             "wetland",
#                             "secchi_m",
#                             "lake_catch_ratio",
#                             "mix_event",
#                             "waterT"),
#                    labels=c(
#                      "Calcium (µeq /L)",
#                      "Sulfate (µeq /L)",
#                      "Alkalinity (µeq /L)",
#                      "Silica (µmol /L)",
#                      "Altitude (m)",
#                      "Z max (m)",
#                      "Wetland (%)",
#                      "Secchi disk (m)",
#                      "Lake/ Catchment area ratio",
#                      "Mix event (degrees)",
#                      "water temperature (º)"))+
#   # scale_y_discrete(limits = c(""),
#   #                  labels=c(""))+
#   xlab(NULL)+
#   ylab(NULL)+
#   theme_classic()+
#   Theme1+
#   geom_text(aes(label=label))
# gg.importance.scores
# 
# # Save the plot
# ggsave("GAM_IT_varimportance.png", plot=gg.importance.scores, height=8, width=10,units="in",
#         dpi = 400)
# 


##### Manually make the most parsimonious GAM models for each taxa ----
spp_GAM_IT_v2 <- readRDS("~/R/seasonal_diatom_cajas/outputs/spp_GAM_IT_v2.rds")

### now  make a nice plot of the most interesting models
library(gridExtra)
library(grid)

# Theme-
Theme1 <-
  theme( # use theme_get() to see available options
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # legend.background = element_rect(fill="white"),
    legend.background = element_blank(),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.text = element_text(size=15),
    legend.title = element_blank(),
    legend.position = c(0.2, 0.8),
    text=element_text(size=15),
    strip.text.y = element_text(size = 15,angle = 0),
    axis.title.x=element_text(vjust=0.3, size=15),
    axis.title.y=element_text(vjust=0.6, angle=90, size=15),
    axis.text.x=element_text(size=15),
    axis.text.y=element_text(size=15),
    axis.line.x=element_line(colour="black", size=0.5,linetype='solid'),
    axis.line.y=element_line(colour="black", size=0.5,linetype='solid'),
    strip.background = element_blank())



# Model Achnanthidium.affine Fe + Zmax
dat.af <- diatom_env_long %>% filter(taxa=="Achnanthidium.affine")
mod.ds <- gam(count ~ s(Fe,k=5, bs="cr") + s(Zmax_m,k=5,  bs="cr"), 
           data=dat.af, method="REML")

# Make partial residual plots and save
ggsave("outputs/Discostella_Aulacoseira_Ca_bymixing.png", plot=draw(mod.ds, residuals = TRUE), height=8, width=10,units="in",
       dpi = 400)

# Predict Fe from model of Achnanthidium affine
mod <- mod.ds
testdata <- expand.grid(Fe=seq(min(mod$model$Fe),max(mod$model$Fe)),
                        Zmax_m=seq(min(mod$model$Zmax_m),max(mod$model$Zmax_m), length.out = length(mod$model$Fe))) 
  
fits <- predict.gam(mod, newdata=testdata, type='response', se.fit=T)
head(fits,2)
predicts.aff <- testdata %>% data.frame(fits)
  

# PLOTS Achnanthidium affine Fe + Zmax ----
ggmod.aff.Fe <- ggplot() +
  ylab(" ")+
  #geom_point(data=dat.af,aes(x=Zmax_m,y=count),  alpha=0.75, size=2,show.legend=FALSE)+
  geom_line(data=predicts.aff, aes(x=Zmax_m,y=fit), alpha=0.5)+
  geom_line(data=predicts.aff, aes(x=Zmax_m,y=fit - se.fit),linetype="dashed",alpha=0.5)+
  geom_line(data=predicts.aff, aes(x=Zmax_m,y=fit + se.fit),linetype="dashed",alpha=0.5)+
  theme_classic()+
  #Theme1+
  annotate("text", x = -Inf, y=Inf, label = "(e)",vjust = 1, hjust = -.1,size=5)
  #annotate("text", x = -Inf, y=Inf, label = "  Pagurus novizelandiae",vjust = 1, hjust = -.1,size=5,fontface="italic")
  #geom_blank(data=dat.af,aes(x=lobster,y=response*1.05))#to nudge data off annotations
ggmod.aff.Fe



