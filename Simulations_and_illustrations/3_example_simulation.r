source("functions_analysis.r")
source("functions_simulation.r")
source("functions_sampling.r")
source("functions_graphics.r")


# create species-abundance and range-size distributions
prms <- set.simulation(n.spec = 20, 
                       n.ind.clust = 50,
                       n.clust.tot = 100,
                       immigr.rate = 0.9,
                       SAD.trend =  -0.5,
                       RSD.trend =  -0.9)

# do the actual simulations
comm <- simulate.comms(prms)

# plot the simulating species
comm.plot(comm)



# OVERLAYING THE POLYGONS ------------------------------------------------------

# create polygon layers
PLS <- make.tess.polygons(prop1 = 0.5, prop2 = 0.5, overlap = TRUE)

# convert ppp patterns to sp objects
smp1 <- count.S.in.poly(comm = pplist.to.sp(comm$PPlist1), poly = PLS$PLS1)
smp2 <- count.S.in.poly(comm = pplist.to.sp(comm$PPlist2), poly = PLS$PLS2)
smp12 <- rbind(data.frame(smp1, Time = 1), data.frame(smp2, Time = 2))


# TRUTH ------------------------------------------------------------------------
# calculate "true" scaling relationships

registerDoParallel(cores = 4)

#system.time(SAR(comm.list = comm, area.span = range(smp12$Area)))
#system.time(SEGAR(comm.list = comm, area.span = range(smp12$Area)))

SEGARS <- SEGAR(comm.list = comm, 
                area.span = range(smp12$Area),
                return.rasters = FALSE)

# FITTING THE SAR MODEL --------------------------------------------------------

n.trees = 500

SAR.model <- gbm(S ~ Area + Time + x + y , 
                 data = smp12, 
                 interaction.depth = 4, 
                 distribution = "poisson",
                 n.minobsinnode = 3,
                 n.trees = n.trees)

fit.variants.gbm(data=smp12, true.SEGAR=SEGARS, area.span=range(smp12$Area))

# PREDICTING FROM THE SAR MODEL ------------------------------------------------

SARS.pred <- predict.SAR.mean(SAR.model, 
                              grains = c(32, 16, 8, 4), 
                              area.span = range(smp12$Area))

obs.pred.SAR <- left_join(x = SEGARS, y = SARS.pred$wide, by = "Area")

# OBSERVED VS PREDICTED --------------------------------------------------------

oi.sar <- ggplot(data = obs.pred.SAR, aes(x = Delta, y = Delta.pred)) +
  geom_point( size = 3) +
  geom_abline(intercept=0, slope=1) + 
  geom_hline(yintercept=0, linetype = 2) +
  geom_vline(xintercept=0, linetype = 2) +
  coord_fixed() +
  labs(title = SARS.pred$modelclass,
       subtitle = "SAR method", x="True Delta S", y="Estimated Delta S") +
  theme_bw() 


# FITTING THE INDIVIDUAL MODEL -------------------------------------------------

### BOOSTED TREES ###

# species-by-site detailed polygon data
X <- count.spec.by.spec(comm = comm, poly = PLS)

MODEL.I <- gbm(occ ~ x + y + Area + Time + spec, 
               data = X, 
               interaction.depth = 5, 
               distribution = "bernoulli",
               n.minobsinnode = 3,
               n.trees = n.trees)

I.preds <- predict.I.mean(MODEL.I, area.span = range(smp12$Area))

obs.pred.I <- left_join(x = SEGARS, y = I.preds$wide, by = "Area")

oi.I <- ggplot(data = obs.pred.I, aes(x = Delta, y = Delta.pred)) +
  geom_point( size = 3) +
  geom_abline(intercept=0, slope=1) + 
  geom_hline(yintercept=0, linetype = 2) +
  geom_vline(xintercept=0, linetype = 2) +
  coord_fixed() +
  labs(subtitle = "OAR method", x="True Delta S", y="Estimated Delta S") +
  theme_bw() 

grid.arrange(oi.sar, oi.I, ncol=2)




# ------------------------------------------------------------------------------
### RANDOM FORESTS ###

SAR.model <- randomForest(S ~ Area + Time + x + y , 
                          data = smp12, 
                          mtry = 3,
                          ntree = n.trees)

fit.variants.rf(data=smp12, true.SEGAR=SEGARS, area.span=range(smp12$Area))


SARS.pred <- predict.SAR.mean(SAR.model, 
                              grains = c(32, 16, 8, 4), 
                              area.span = range(smp12$Area))
obs.pred.SAR <- left_join(x = SEGARS, y = SARS.pred$wide, by = "Area")

oi.sar <- ggplot(data = obs.pred.SAR, aes(x = Delta, y = Delta.pred)) +
  geom_point( size = 3) +
  geom_abline(intercept=0, slope=1) + 
  geom_hline(yintercept=0, linetype = 2) +
  geom_vline(xintercept=0, linetype = 2) +
  coord_fixed() +
  labs(title = SARS.pred$modelclass,
       subtitle = "SAR method", x="True Delta S", y="Estimated Delta S") +
  theme_bw() 


# species-by-site detailed polygon data
X <- count.spec.by.spec(comm = comm, poly = PLS)

MODEL.I <- randomForest(as.factor(occ) ~ x + y + Area + Time + spec, 
               data = X, 
               mtry = 3, 
               ntree = n.trees)

I.preds <- predict.I.mean(MODEL.I, area.span = range(smp12$Area))

oi.I <- ggplot(data = obs.pred.I, aes(x = Delta, y = Delta.pred)) +
  geom_point( size = 3) +
  geom_abline(intercept=0, slope=1) + 
  geom_hline(yintercept=0, linetype = 2) +
  geom_vline(xintercept=0, linetype = 2) +
  coord_fixed() +
  labs(subtitle = "OAR method", x="True Delta S", y="Estimated Delta S") +
  theme_bw() 


grid.arrange(oi.sar, oi.I, ncol=2)

