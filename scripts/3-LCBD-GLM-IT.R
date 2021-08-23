
############################################
### GLM-IT
#Load functions (from Carles Alcaraz)
source("scripts/functions/IC.model.inference.R") 

#Load GLM variables and select from database
model_var <- c("LCBD", "Mg", "K", "Alkalinity", "Cond", "Si", "Altitude", "Zmax_m", "Pajonal", 
               "lake_catch_ratio", "mix_event")  
modelGLM_var <- model_LCBD_var[, (names(model_LCBD_var) %in% model_var)]

#Perform GLM_IT analysis
#Define GLM Full model
Global.model <- lm(LCBD ~ . , data = modelGLM_var) 

#Perform normality tests
Normality.results <- NULL
Normality.results <- matrix(data = NA, nrow = 5, ncol = 4, byrow = FALSE, dimnames = NULL)
colnames(Normality.results) <- c("Method", "Statistic", "Statistic.value", "P.value")
Normality.results[1, 1] <- ad.test(Global.model$residuals)$method
Normality.results[1, 2] <- names(ad.test(Global.model$residuals)[[1]])
Normality.results[1, 3] <- ad.test(Global.model$residuals)$statistic
Normality.results[1, 4] <- ad.test(Global.model$residuals)$p.value
Normality.results[2, 1] <- cvm.test(Global.model$residuals)$method
Normality.results[2, 2] <- names(cvm.test(Global.model$residuals)[[1]])
Normality.results[2, 3] <- cvm.test(Global.model$residuals)$statistic
Normality.results[2, 4] <- cvm.test(Global.model$residuals)$p.value
Normality.results[3, 1] <- lillie.test(Global.model$residuals)$method
Normality.results[3, 2] <- names(lillie.test(Global.model$residuals)[[1]])
Normality.results[3, 3] <- lillie.test(Global.model$residuals)$statistic
Normality.results[3, 4] <- lillie.test(Global.model$residuals)$p.value
Normality.results[4, 1] <- pearson.test(Global.model$residuals)$method
Normality.results[4, 2] <- names(pearson.test(Global.model$residuals)[[1]])
Normality.results[4, 3] <- pearson.test(Global.model$residuals)$statistic
Normality.results[4, 4] <- pearson.test(Global.model$residuals)$p.value
Normality.results[5, 1] <- sf.test(Global.model$residuals)$method
Normality.results[5, 2] <- names(sf.test(Global.model$residuals)[[1]])
Normality.results[5, 3] <- sf.test(Global.model$residuals)$statistic
Normality.results[5, 4] <- sf.test(Global.model$residuals)$p.value
Normality.results
#write.table(Normality.results, file = "Normality.results.txt", sep = "\t", col.names = TRUE, row.names = FALSE, na = "")


#Load model parameters
IC.rank <- c("AICc")  
IC.VIF <- TRUE
C.value <- 0.95
D.value <- 7

#Define all possible combinations
Mod.reg <- IC.selection(Global.model, IC.sel = IC.rank, C.hat = 1000, IC.max.var = NA, IC.fix.var = NULL, beta = FALSE, IC.save = "IC.model.results.txt", IC.coef = FALSE)

#Define all possible combinations
Mod.reg <- IC.selection(Global.model, IC.sel = IC.rank, C.hat = 1000, IC.max.var = NA, IC.fix.var = NULL, beta = FALSE, IC.save = "IC.model.results.txt", IC.coef = FALSE)

#Perform Likelihood ratio test to evaluate models that perform better than the null model
Mod.reg <- Likelihood.test(Mod.reg, Global.model, LH.deletion = TRUE, LH.alpha = 0.05, LH.save = "Likelihood.selected.models.txt", LH.analysis = "Likelihood.results.txt") 

#Perform VIF test to evaluate collinearity effects
Mod.reg <- VIF.evaluation(Mod.reg, Global.model, beta = FALSE, VIF.deletion = TRUE, VIF.value = 5, VIF.save = "VIF.selected.models.txt", VIF.analysis = "VIF.results.txt")

#Select models
#Selection based on delta AIC
Mod.reg <- Subset.models(Mod.reg, subset = Delta < D.value, sub.save = "IC.models.subset.txt")
#Selection based on AIC weight
#Mod.reg <- Subset.models(Mod.reg, selection = cumsum(Weight) <= C.value, sub.save = "IC.models.subset.txt")		

#Extract selected models
Models <- Extract.models(Mod.reg)

#Model averaging
Mod.reg <- Average.models(Models, IC.sel = IC.rank, C.hat = 1000, beta = FALSE, na.values = c("NA"), alpha = 0.05, Avm.save = "Average.model.results.txt")
print(Mod.reg)

#Compute model evaluation data: residuals and fitted values
Model.eval <- Models.evaluation(Global.model, Models, Mod.reg, IC.sel = IC.rank, C.hat = 1000, na.values = c("NA"), Eval.save = "Models.evaluation.results.txt")

#Compute parameter estimation bias
Bias.models(Mod.reg, Global.model, Bias.save = "Models.parameters.bias.txt")

#Perform model evaluation graphs
Graph.data <- model_LCBD_var
#Evaluation.plot(X.var = "Order_M" , Y.var = "LCBD", Graph.data, Model.eval, Graphs = c(1, 2, 3, 4, 5), Output = "Chloro Vs Time.png", Title = "Temporal Vs Fitted models", X.Title = "Temporal variation", Y.Title = "Log(Chlorophyll (?g / L))", Legend = 2, Xv = 2, Yv = 0.2)
Evaluation.plot(X.var = "LCBD", Y.var = "LCBD", Graph.data, Model.eval, Graphs = c(1, 2, 3, 4, 5), Output = "Chloro Vs Models.png", Title = "Chlorophyll Vs Fitted models", X.Title = "Log(Chlorophyll (?g / L))", Y.Title = "Log(Chlorophyll (?g / L))", Legend = 1, Xv = 0.2, Yv = 0.2)


############################################
### GLM-IT
## Approach 2 from Feld et al. 2016 STOTEN multistressor
model_LCBD_var_2 <- transform(model_LCBD_var, SubCuenca=factor(SubCuenca))

mod1 <- lmer(LCBD ~ Mg + K + Alkalinity + Cond + Si + Altitude + Zmax_m + Pajonal + lake_catch_ratio +
                mix_event + (1|SubCuenca),  data=model_LCBD_var_2) 

# Calculate all possible models with different combinations
options(na.action="na.fail") #necessary to run dredge()
gmod1 <- dredge(mod1, rank="AICc", extra = "adjR^2")

# Select a set of "top models"
g_model1_set1c <- get.models (gmod1, subset=delta<=2) # subset with delta AIC<2

# AVerage models
g_model1_av <- model.avg(g_model1_set1c, revised.var=TRUE)
summary(g_model1_av)

#Write model outputs
#MA.est.table<-round(g_model1_av$coefficients[,c(1,2:14),10])

#Model validation
model1.res <- resid(mod1, type="pearson")
hist(model1.res)

scatter.smooth(fitted(mod1), model1.res)

#Check spatial autocorrelation of model residuals
#Extract Latitude and Longitude
xy.coord <- data.frame(spatial_var$Latitude, spatial_var$Longitude)
my.coord <- coordinates (xy.coord)

#no duplicate rows
my.coord.nd <- my.coord[!duplicated(my.coord), ]

duplicated(my.coord)
sel <- c(113)

model1.res2 <- as.matrix(model1.res)
model1.res2 <- model1.res[-sel,]

coord.nb <- tri2nb(my.coord.nd, row.names = NULL)
my.coord.listw <- nb2listw (coord.nb, glist=NULL, style="W",
                            zero.policy=FALSE)

moran.test (as.vector (model1.res2), my.coord.listw)




