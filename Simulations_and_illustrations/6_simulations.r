library(doParallel)
library(raster)
library(spatstat)
library(mvtnorm)
library(truncnorm)
library(mobsim)
library(plyr)
library(maptools)
library(rgeos)
library(rgdal)
library(randomForest)
library(gbm)
library(dplyr)
library(tidyverse)
library(viridis)
library(sampling)


# ------------------------------------------------------------------------------

# determine how many simulations have already been done 
out.file <- "data_simulated_10Jun2021.csv"
out.data <- read.csv(out.file)

if(out.file %in% list.files())
{
  N.done <- max(out.data$sim.id)
} else(N.done = 0)

# ------------------------------------------------------------------------------

source("functions_analysis.r")
source("functions_simulation.r")
source("functions_sampling.r")
source("functions_graphics.r")

# -----------------------------------------------------------------------------


# register number of cores for parallel processing
registerDoParallel(cores = 4)

SIM.PARAMS <- expand.grid(n.spec = c(13, 26, 52), # random forest can't handle more than 53 sp.
                          n.ind.clust = c(100, 500),
                          n.clust.tot = c(100, 500),
                          immigr.rate = c(0, 0.2, 0.7),
                          SAD.trend = c(-0.7, -0.3, 0, 0.3, 0.7),
                          RSD.trend = c(-0.7, -0.3, 0, 0.3, 0.7),
                          prop1 = c(0.1, 0.3, 0.5),
                          prop2 = c(0.1, 0.3, 0.5),
                          overlap = c(TRUE, FALSE))

#SIM.PARAMS <- expand.grid(n.spec = c(26, 52), # random forest can't handle more than 53 sp.
#                          n.ind.clust = c(100, 500),
#                          n.clust.tot = c(100, 500),
#                          immigr.rate = c(0.2, 0.7),
#                          SAD.trend = c(-0.7,  0, 0.7),
#                          RSD.trend = c(-0.7,  0, 0.7),
#                          prop1 = c( 0.3, 0.5),
#                          prop2 = c(0.3, 0.5),
#                          overlap = c(TRUE, FALSE))

examined.grains <-  c(4, 8, 16, 32)
N.tree = 500

res <- list()

# initialize the file in case we start from scratch
if(N.done == 0)
{
  res.header <- c("sim.id", "Area", "S1", "S2", "Gain", "Loss", "Delta", "est.Delta", "Formula",
                  "Method", "n.spec", "n.ind.clust", "n.clust.tot", "immigr.rate","SAD.trend",
                  "RSD.trend", "prop1", "prop2", "overlap")
  resf <- data.frame(matrix(nrow=0, ncol=length(res.header)))
  names(resf) <- res.header
  write.table(resf, file=out.file, sep=",")
}

# RUN THE SIMULATIONS
for(i in (N.done+1):nrow(SIM.PARAMS))
{
  set.seed(12345)
  
  message(paste("Simulation", i, "out of", nrow(SIM.PARAMS)))

  message("1 Simulating community")
  print(SIM.PARAMS[i,])
  
  prms <- set.simulation(n.spec = SIM.PARAMS$n.spec[i], 
                         n.ind.clust = SIM.PARAMS$n.ind.clust[i],
                         n.clust.tot = SIM.PARAMS$n.clust.tot[i],
                         immigr.rate= SIM.PARAMS$immigr.rate[i],
                         SAD.trend = SIM.PARAMS$SAD.trend[i],
                         RSD.trend = SIM.PARAMS$RSD.trend[i])
  
  # do the actual simulations
  comm <- simulate.comms(prms)

  # create polygon layers
  message("2 Generating polygons")
  PLS <- make.tess.polygons(prop1 =  SIM.PARAMS$prop1[i], 
                            prop2 =  SIM.PARAMS$prop2[i], 
                            overlap =  SIM.PARAMS$overlap[i])

  
  # convert ppp patterns to sp objects
  message("3 Converting ppp to sp objects")
  smp1 <- count.S.in.poly(comm = pplist.to.sp(comm$PPlist1), poly = PLS$PLS1)
  smp2 <- count.S.in.poly(comm = pplist.to.sp(comm$PPlist2), poly = PLS$PLS2)
  smp12 <- rbind(data.frame(smp1, Time = 1), data.frame(smp2, Time = 2))

  message("4 True SEGARs (Species Extinction- or Gain-Area Relationships)")
  SEGARS <- SEGAR(grains = examined.grains,
                  comm.list = comm, 
                  area.span = range(smp12$Area),
                  return.rasters = FALSE)
  
  
  # returned object in case there is an error
  error.fun <- function(cond)
  {
    res.i <-  data.frame(Area = NA, 
                         S1 = NA, 
                         S2 = NA,
                         Gain = NA,
                         Loss = NA,
                         Delta = NA,
                         est.Delta = NA,
                         Formula = NA)
    assign("res.i", res.i, envir = .GlobalEnv )
    #return(res.i)
  }
  
  # APPLY SAR-based ML algorithms
  tryCatch(
  {
    message("5a SAR-based RandomForest models")
    
    SAR.based.rf <- fit.variants.rf(data = smp12,
                                    true.SEGAR = SEGARS,
                                    grains = examined.grains, 
                                    area.span = range(smp12$Area),
                                    n.trees = N.tree) 
    SAR.based.rf <- data.frame(SAR.based.rf, Method = "Random forest" )
    
    message("5b SAR-based gbm models")
    SAR.based.gbm <- fit.variants.gbm(data = smp12,
                                      true.SEGAR = SEGARS,
                                      grains = examined.grains, 
                                      area.span = range(smp12$Area),
                                      n.trees = N.tree) 
    SAR.based.gbm <- data.frame(SAR.based.gbm, Method = "Boosted trees" )
    
    
    # APPLY OAR-based ML algorithms
    message("6a OAR-based RandomForest model")
    X <- count.spec.by.spec(comm = comm, poly = PLS)

    OAR.model.rf <- randomForest(as.factor(occ) ~ x + y + Area + Time + spec, 
                                 data = X, 
                                 mtry = 3,
                                 ntree = N.tree)

    OAR.preds.rf <- predict.OAR.mean(OAR.model.rf, 
                                area.span = range(smp12$Area))
    OAR.preds.rf <- select(OAR.preds.rf$wide, Area, est.Delta = Delta.pred)
    OAR.preds.rf <- left_join(SEGARS, OAR.preds.rf, by = "Area")
    OAR.preds.rf <- data.frame(OAR.preds.rf, 
                               Formula = "OAR-based", Method = "Random forest")
   
    message("6b OAR-based gbm model")
    OAR.model.gbm <- gbm(occ ~ x + y + Area + Time + as.factor(spec), 
                         data = X, 
                         interaction.depth = 5, 
                         distribution = "bernoulli",
                         n.minobsinnode = 3,
                         n.trees = N.tree)
    OAR.preds.gbm <- predict.OAR.mean(OAR.model.gbm, 
                                     area.span = range(smp12$Area))
    OAR.preds.gbm <- select(OAR.preds.gbm$wide, Area, est.Delta = Delta.pred)
    OAR.preds.gbm <- left_join(SEGARS, OAR.preds.gbm, by = "Area")
    OAR.preds.gbm <- data.frame(OAR.preds.gbm, 
                               Formula = "OAR-based", Method = "Boosted trees")
    

    res.i <- rbind(SAR.based.rf, OAR.preds.rf, SAR.based.gbm, OAR.preds.gbm)
  },
  error = error.fun)
  
  message("8 Stepping outside the tryCatch fun")
  res.i <- data.frame(sim.id = i, res.i, SIM.PARAMS[i,])
  
  write.table(res.i, file= out.file,
              sep=",", 
              append=TRUE,
              row.names = FALSE, col.names = FALSE)
  
  
  #res[[i]] <- res.i
  print(res.i)
}




