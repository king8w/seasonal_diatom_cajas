######################################
#### MORAN'S Eigenvector Maps  #######
######################################

#### Code  from Dray et al. 2012 #####
library(spacemakeR)

# Read in geographical coordinates lakes
spatial_var <- read.csv("data/Spatial_Cajas2019.csv", row.names = 1)

#Extract Latitude and Longitude
xy.coord <- data.frame(spatial_var$Latitude, spatial_var$Longitude)

## create the Gabriel graph
xy.coord <- as.matrix(xy.coord)

nb1 <- graph2nb(gabrielneigh(xy.coord, nnmult = 5), sym = T)

lw1 <- nb2listw(nb1)

# create the MEMs
U <- scores.listw(lw1)

# test the significance of Moran's I for MEMs
mctests <- test.scores(U, lw1, 999)
colnames(U$vectors) <- paste("MEM", 1:ncol(U$vectors), sep="")
U.sel <- U
U.sel$vectors <- U.sel$vectors[, mctests[,2]<0.05] ## keep only MEMS with signif values

#Extract MEMS
mems <- U.sel$vectors
spatial <- data.frame(mems)

# save MEMs
write.csv(spatial, "outputs/mems.csv")
