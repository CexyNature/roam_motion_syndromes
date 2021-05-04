library(adehabitatHR)
library(adehabitatLT)
library(lubridate)
library(ggfortify)

rm(list=ls())

# Some of the code below is coming from: 
# https://github.com/dpseidel/MovEco-R-Workshop/blob/b55ef39dcc474c94589dfe6b23fdf9756e9f0bc8/Materials/Day6/Syndromes.md

suffix <- "_CH1"

# Load simulated data
files <- list.files(path = "data/simulated_individuals/", pattern= paste(suffix, ".Rds"), full.names = TRUE)
syndromes <- lapply(files, readRDS)

central <- as.ltraj(xy = syndromes[[1]][,c("Longitude", "Latitude")],
                    id = "central",
                    date = syndromes[[1]]$Datetime,
                    infolocs = syndromes[[1]][, c("Dist_from_Center", "Angle_from Center")])

migrant <- as.ltraj(xy = syndromes[[2]][,c("Longitude", "Latitude")],
                    id = "migrant",
                    date = syndromes[[2]]$Datetime,
                    infolocs = syndromes[[2]][, c("Dist_from_Center", "Angle_from Center")])

nomad <- as.ltraj(xy = syndromes[[3]][,c("Longitude", "Latitude")],
                    id = "nomad",
                    date = syndromes[[3]]$Datetime,
                    infolocs = syndromes[[3]][, c("Dist_from_Center", "Angle_from Center")])

territorial <- as.ltraj(xy = syndromes[[4]][,c("Longitude", "Latitude")],
                    id = "territorial",
                    date = syndromes[[4]]$Datetime,
                    infolocs = syndromes[[4]][, c("Dist_from_Center", "Angle_from Center")])

# Mean turn angle correlation (TAC)
central.tac <- matrix(ncol=1, nrow=length(central))
for (i in 1:length(central)){
  SA <- adehabitatLT::acfang.ltraj(central[i], which = "relative", plot = FALSE, lag = 1)
  central.tac[i,] <- 1/(mean(SA[[1]][1,], na.rm = TRUE))
}

migrant.tac <- matrix(ncol=1, nrow=length(migrant))
for (i in 1:length(migrant)){
  SA <- adehabitatLT::acfang.ltraj(migrant[i], which = "relative", plot = FALSE, lag = 1)
  migrant.tac[i,] <- 1/(mean(SA[[1]][1,], na.rm = TRUE))
}

nomad.tac <- matrix(ncol=1, nrow=length(nomad))
for (i in 1:length(nomad)){
  SA <- adehabitatLT::acfang.ltraj(nomad[i], which = "relative", plot = FALSE, lag = 1)
  nomad.tac[i,] <- 1/(mean(SA[[1]][1,], na.rm = TRUE))
}

territorial.tac <- matrix(ncol=1, nrow=length(territorial))
for (i in 1:length(territorial)){
  SA <- adehabitatLT::acfang.ltraj(territorial[i], which = "relative", plot = FALSE, lag = 1)
  territorial.tac[i,] <- 1/(mean(SA[[1]][1,], na.rm = TRUE))
}

# # Cesar's Residency time (RT)
# central.rt <- residenceTime(central, 
#               radius =  mean(central[[1]]$dist, na.rm = TRUE), 
#               maxt = 12, 
#               addinfo = FALSE,
#               units = "hours")
# central.rt <-mean(central.rt[[1]]$RT.634.0561, na.rm=T)
# 
# 
# migrant.rt <- residenceTime(migrant, 
#                             radius =  mean(migrant[[1]]$dist, na.rm = TRUE), 
#                             maxt = 12, 
#                             addinfo = FALSE,
#                             units = "hours")
# migrant.rt <-mean(migrant.rt[[1]]$RT.265.2344, na.rm=T)
# 
# nomad.rt <- residenceTime(nomad, 
#                             radius =  mean(nomad[[1]]$dist, na.rm = TRUE), 
#                             maxt = 12, 
#                             addinfo = FALSE,
#                             units = "hours")
# nomad.rt <-mean(nomad.rt[[1]]$RT.254.0528, na.rm=T)
# 
# territorial.rt <- residenceTime(territorial, 
#                           radius =  mean(territorial[[1]]$dist, na.rm = TRUE), 
#                           maxt = 12, 
#                           addinfo = FALSE,
#                           units = "hours")
# territorial.rt <-mean(territorial.rt[[1]]$RT.699.6798, na.rm=T)

# Residency Time (RT) and Time to Return (TtoR or T2R)
RTandT2R <- function(x, radius, maxt, units="hour", addinfo = F){
  fR <- function(x, dframe, radius, maxt, units=units){
    tmp <- dframe[c(x:nrow(dframe)),]
    dists <- sqrt((tmp$x - tmp$x[1])^2 + (tmp$y - tmp$y[1])^2)
    dists <- as.numeric(dists<=radius)
    ext <- which(dists[-length(dists)] > dists[-1])+1
    entr <-  which(dists[-length(dists)] < dists[-1])+1
    bts <- difftime(tmp$date[entr], tmp$date[ext[c(1:length(entr))]], units=units)    
    tmp1 <- as.numeric(difftime(tmp$date[ext[(as.numeric(bts)>maxt)][1]], tmp$date[1], units=units)) #first exit
    if (is.na(tmp1) & length(ext)>0) tmp1 <- as.numeric(difftime(tmp$date[ext[length(ext)]], tmp$date[1], units=units))  
    tmp2 <- as.numeric(difftime(tmp$date[entr[(as.numeric(bts)>maxt)][1]], tmp$date[1], units=units)) #first re-entry
    return(c(tmp1, tmp2))
  } 
  res <- data.frame(do.call(rbind, lapply(c(1:nrow(x)), fR, dframe=x, radius=radius, maxt=maxt, units=units)))
  names(res) <- c(paste("RT", radius, maxt, sep="_"), paste("T2R", radius, maxt, sep="_"))
  
  if (!addinfo) return(res)
  if (addinfo) {
    attributes(x)$infolocs <- cbind(attributes(x)$infolocs, res)
    return(x) 
  }
}


central.rt_t2r <- RTandT2R(central[[1]], 
                        radius = mean(central[[1]]$dist, na.rm = TRUE), 
                        maxt= 12, 
                        units="hour", 
                        addinfo = F)

central.rt <- mean(central.rt_t2r[,1], na.rm = TRUE)
central.t2r <- mean(central.rt_t2r[,2], na.rm = TRUE)


migrant.rt_t2r <- RTandT2R(migrant[[1]], 
                           radius = mean(migrant[[1]]$dist, na.rm = TRUE), 
                           maxt= 12, 
                           units="hour", 
                           addinfo = F)

migrant.rt <- mean(migrant.rt_t2r[,1], na.rm = TRUE)
migrant.t2r <- mean(migrant.rt_t2r[,2], na.rm = TRUE)

nomad.rt_t2r <- RTandT2R(nomad[[1]], 
                           radius = mean(nomad[[1]]$dist, na.rm = TRUE), 
                           maxt= 12, 
                           units="hour", 
                           addinfo = F)

nomad.rt <- mean(nomad.rt_t2r[,1], na.rm = TRUE)
nomad.t2r <- mean(nomad.rt_t2r[,2], na.rm = TRUE)

territorial.rt_t2r <- RTandT2R(territorial[[1]], 
                           radius = mean(territorial[[1]]$dist, na.rm = TRUE), 
                           maxt= 12, 
                           units="hour", 
                           addinfo = F)

territorial.rt <- mean(territorial.rt_t2r[,1], na.rm = TRUE)
territorial.t2r <- mean(territorial.rt_t2r[,2], na.rm = TRUE)


# Volume of Intersection (VI)
central.vi <- matrix(ncol=1, nrow=length(central)) # create an empty matrix to populate with a for-loop
for (i in 1:length(central)){
  central[[i]] <- central[[i]][complete.cases(central[[i]][,c("x","y")]),] # remove NAs from coordinates
  central[[i]]$month <- month(central[[i]]$date)
  kudoverlap_monthly <- adehabitatHR::kerneloverlap(SpatialPointsDataFrame(central[[i]][,c("x","y")], 
                                                                           data=data.frame(id=central[[i]]$month)), 
                                                    grid=200, method="VI")
  mw <- matrix(nr = nrow(kudoverlap_monthly), nc = nrow(kudoverlap_monthly))
  mw<-(row(mw) == col(mw) - 1) + 0 
  monthval<-kudoverlap_monthly * mw
  avg_kudoverlap_monthly<-sum(monthval)/(nrow(monthval)-1) # average month-to-month volume of intersection 
  central.vi[i,] <- c(avg_kudoverlap_monthly) 
}

central.vi<-data.frame(central.vi)
central.vi

migrant.vi <- matrix(ncol=1, nrow=length(migrant)) # create an empty matrix to populate with a for-loop
for (i in 1:length(migrant)){
  migrant[[i]] <- migrant[[i]][complete.cases(migrant[[i]][,c("x","y")]),] # remove NAs from coordinates
  migrant[[i]]$month <- month(migrant[[i]]$date)
  kudoverlap_monthly <- adehabitatHR::kerneloverlap(SpatialPointsDataFrame(migrant[[i]][,c("x","y")], 
                                                                           data=data.frame(id=migrant[[i]]$month)), 
                                                    grid=200, method="VI")
  mw <- matrix(nr = nrow(kudoverlap_monthly), nc = nrow(kudoverlap_monthly))
  mw<-(row(mw) == col(mw) - 1) + 0 
  monthval<-kudoverlap_monthly * mw
  avg_kudoverlap_monthly<-sum(monthval)/(nrow(monthval)-1) # average month-to-month volume of intersection 
  migrant.vi[i,] <- c(avg_kudoverlap_monthly) 
}

migrant.vi<-data.frame(migrant.vi)
migrant.vi


nomad.vi <- matrix(ncol=1, nrow=length(nomad)) # create an empty matrix to populate with a for-loop
for (i in 1:length(nomad)){
  nomad[[i]] <- nomad[[i]][complete.cases(nomad[[i]][,c("x","y")]),] # remove NAs from coordinates
  nomad[[i]]$month <- month(nomad[[i]]$date)
  kudoverlap_monthly <- adehabitatHR::kerneloverlap(SpatialPointsDataFrame(nomad[[i]][,c("x","y")], 
                                                                           data=data.frame(id=nomad[[i]]$month)), 
                                                    grid=200, method="VI")
  mw <- matrix(nr = nrow(kudoverlap_monthly), nc = nrow(kudoverlap_monthly))
  mw<-(row(mw) == col(mw) - 1) + 0 
  monthval<-kudoverlap_monthly * mw
  avg_kudoverlap_monthly<-sum(monthval)/(nrow(monthval)-1) # average month-to-month volume of intersection 
  nomad.vi[i,] <- c(avg_kudoverlap_monthly) 
}

nomad.vi<-data.frame(nomad.vi)
nomad.vi

territorial.vi <- matrix(ncol=1, nrow=length(territorial)) # create an empty matrix to populate with a for-loop
for (i in 1:length(territorial)){
  territorial[[i]] <- territorial[[i]][complete.cases(territorial[[i]][,c("x","y")]),] # remove NAs from coordinates
  territorial[[i]]$month <- month(territorial[[i]]$date)
  kudoverlap_monthly <- adehabitatHR::kerneloverlap(SpatialPointsDataFrame(territorial[[i]][,c("x","y")], 
                                                                           data=data.frame(id=territorial[[i]]$month)), 
                                                    grid=200, method="VI")
  mw <- matrix(nr = nrow(kudoverlap_monthly), nc = nrow(kudoverlap_monthly))
  mw<-(row(mw) == col(mw) - 1) + 0 
  monthval<-kudoverlap_monthly * mw
  avg_kudoverlap_monthly<-sum(monthval)/(nrow(monthval)-1) # average month-to-month volume of intersection 
  territorial.vi[i,] <- c(avg_kudoverlap_monthly) 
}

territorial.vi<-data.frame(territorial.vi)
territorial.vi

# Maximum net squared displacement (MNSD)
central.mnsd <- max(central[[1]]$R2n, na.rm = TRUE) / min(central[[1]]$R2n[central[[1]]$R2n > 0], na.rm = TRUE)
migrant.mnsd <- max(migrant[[1]]$R2n, na.rm = TRUE) /  min(migrant[[1]]$R2n[migrant[[1]]$R2n > 0], na.rm = TRUE)
nomad.mnsd <- max(nomad[[1]]$R2n, na.rm = TRUE) /  min(nomad[[1]]$R2n[nomad[[1]]$R2n > 0], na.rm = TRUE)
territorial.mnsd <- max(territorial[[1]]$R2n, na.rm = TRUE) /  min(territorial[[1]]$R2n[territorial[[1]]$R2n > 0], na.rm = TRUE)


# Create final df
syndrome_metrics <- data.frame(syndrome = c("central", "migrant", "nomad", "territorial"),
                               group = suffix,
                               tac = c(central.tac[1,1], migrant.tac[1,1], nomad.tac[1,1], territorial.tac[1,1]),
                               rt = c(central.rt, migrant.rt, nomad.rt, territorial.rt),
                               t2r = c(central.t2r, migrant.t2r, nomad.t2r, territorial.t2r),
                               vi = c(central.vi[1,1], migrant.vi[1,1], nomad.vi[1,1], territorial.vi[1,1]),
                               mnsd = c(central.mnsd, migrant.mnsd, nomad.mnsd, territorial.mnsd))

filename_group_df <- paste("data/simulated_group/group", suffix, ".Rds")
saveRDS(syndrome_metrics, filename_group_df)

# PCA
filename_pca <- paste("figures/simulated_group_pca/pca_group", suffix, ".png")

pca <- prcomp(syndrome_metrics[,-1:-2],
              center = TRUE,
              scale = TRUE)

autoplot(pca, data = syndrome_metrics, colour = "syndrome", size = 3) +
  labs(title = paste("PCA Group: ", suffix),
       color = "Syndrome") +
  theme(panel.grid.major = element_line(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "#ffffff",color = NA),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 0, color = "#4e4d47"),
        legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size=5)))

dev.copy(png, filename_pca, width=6, height=6, units="in", res=500)
dev.off()
