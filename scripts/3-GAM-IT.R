#---------------------------------------------------------------------------------------
# Script: Diatom species level models against chemical, physical and spatial variables modeled with GAM-IT 
# Paper: Space, not time, drive is driving contemporary planktic diatom species composition in tropical Andean lakes
# Author: Benito, X.
# e-mail: xavier.benito.granell@gmail.com
#---------------------------------------------------------------------------------------

#Load libraries
library(MuMIn)
library(FSSgam) #Fisher et al. 2020 (https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.4134)

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



