
source("functions_analysis.r")
source("functions_simulation.r")
source("functions_sampling.r")
source("functions_graphics.r")


# -----------------------------------------------------------------------------

set.seed(12345)

# create the polygon layers
PLS <- make.tess.polygons(prop1 = 0.3, prop2 = 0.3, overlap = FALSE)

pdf("../graphics/SAR_polygon_examples.pdf", width =  6, height = 3)
polygon.plot(PLS)
dev.off()

# set simulation paramters
# COMPLEX SCENARIO
prms1 <- set.simulation(n.spec = 60, 
                        n.ind.clust = 80,
                        n.clust.tot = 400,
                        immigr.rate= 0.9,
                        SAD.trend =  -0.5,
                        RSD.trend =  0.9)

# LOSS SCENARIO
prms2 <- set.simulation(n.spec = 60, 
                        n.ind.clust = 80,
                        n.clust.tot = 400,
                        immigr.rate= 0,
                        SAD.trend =  -0.5,
                        RSD.trend =  -0.5)

# GAIN SCENARIO
prms3 <- set.simulation(n.spec = 60, 
                        n.ind.clust = 80,
                        n.clust.tot = 400,
                        immigr.rate= 0.5,
                        SAD.trend =  0.5,
                        RSD.trend =  0.5)

prms.list <- list(prms1, prms2, prms3)


pdf("../graphics/SAR_example_scenarios.pdf", width =  8, height = 4)
for(i in 1:3)
{
  print(i)  
  
  prms <- prms.list[[i]]
  
  # run the main simulation
  comm <- simulate.comms(prms)
  
  # calculate "true" scaling relationships
  SARS <- SEGAR(comm.list = comm)
  
  list2env(comm, envir=.GlobalEnv); names(comm)
  list2env(prms, envir=.GlobalEnv); names(prms)
  list2env(SARS, envir=.GlobalEnv); names(SARS)
  
  
  # convert ppp patterns to sp objects
  comm1 <- pplist.to.sp(PPlist1)
  comm2 <- pplist.to.sp(PPlist2)
  
  # count species richness in the polygons
  smp1 <- count.S.in.poly(comm = comm1, poly = PLS$PLS1)
  smp2 <- count.S.in.poly(comm = comm2, poly = PLS$PLS2)
  smp12 <- rbind(data.frame(smp1, Time = 1), data.frame(smp2, Time = 2))
  
  # species-by-site detailed polygon data
  X <- count.spec.by.spec(comm = comm, poly = PLS)
  
  # the truth in all of the polygons
  truth.poly <- count.SEG.in.poly(comm = comm, poly = PLS$PLS.complete) 
  
  # ------------------------------------------------------------------------------
  SAR.model <- gbm(S ~ Area + Time + x + y , 
                   data = smp12, 
                   interaction.depth = 4, 
                   distribution = "poisson",
                   n.minobsinnode = 3,
                   n.trees = 500)
  
  # ------------------------------------------------------------------------------
  # predictions in all polygons
  S1.pred <- predict(SAR.model, newdata = data.frame(truth.poly, Time = 1), 
                     n.trees = 500, type = "response")
  
  S2.pred <- predict(SAR.model, newdata = data.frame(truth.poly, Time = 2), 
                     n.trees = 500, type = "response")
  Delta.pred <- S2.pred - S1.pred
  poly.all <- data.frame(truth.poly, S1.pred, S2.pred, Delta.pred)
  
  # RAW PRED VS OBS
  #ggplot(data = poly.all, aes(x = Delta, y = Delta.pred)) +
  #  geom_point() +  
  #  geom_abline(intercept=0, slope=1) + 
  #  geom_hline(yintercept=0, linetype = 2) +
  #  geom_vline(xintercept=0, linetype = 2) +
  #  coord_fixed() 
  
  # ------------------------------------------------------------------------------
  # predicted mean values
  # --------------------
  # adjust the data to be usable in the predictions and to not extrapolate
  Truth <- data.frame(Area = SARS$Delta$Area,
                      Delta = SARS$Delta$Delta, 
                      S.mean1 = SARS$SAR1$S.mean, 
                      S.mean2 = SARS$SAR2$S.mean)
  Delta <- Delta[Delta$Area >= min(smp12$Area) ,]
  Delta <- Delta[Delta$Area <= max(smp12$Area) ,]
  
  
  GRID <- prepare.grid(area.span = c(min(smp12$Area), max(smp12$Area))) 
  
  
  S.pred1 <- predict(SAR.model, 
                     newdata = GRID1, 
                     n.trees = 500,
                     type = "response")
  S.pred2 <- predict(SAR.model, 
                     newdata = GRID2, 
                     n.trees = 500,
                     type = "response")
  
  GRID <- data.frame(GRID1, S.pred1, S.pred2, Delta.pred = S.pred2 - S.pred1)
  summary(GRID)
  
  SAR12.pred <- plyr::ddply(.data = GRID,
                            .variables=c("Area", "Time"),
                            .fun = summarise,
                            S.mean.pred1 = mean(S.pred1),
                            S.mean.pred2 = mean(S.pred2),
                            Delta.mean.pred = mean(Delta.pred))
  SAR12.pred <- select(SAR12.pred, -Time)
  
  SAR <- left_join(x = SAR12.pred, y = Truth, by = "Area")
  SAR <- na.omit(SAR)
  
  # ------------------------------------------------------------------------------
  # PLOTTING
  # 1. The true SARs
  for.plot.true <- rbind(data.frame(Area = SARS$Delta$Area, S = SARS$SAR1$S.mean),
                         data.frame(Area = SARS$Delta$Area, S = SARS$SAR2$S.mean)) 
  for.plot.true <- data.frame(for.plot.true, 
                              Time = rep(1:2, each=nrow(for.plot.true)/2),
                              Type = "Truth")
  for.plot.est <- rbind(data.frame(Area = SAR$Area, S = SAR$S.mean.pred1),
                        data.frame(Area = SAR$Area, S = SAR$S.mean.pred2))
  for.plot.est <- data.frame(for.plot.est,
                             Time = rep(1:2, each=nrow(for.plot.est)/2),
                             Type = "Estimate")
  for.plot <- rbind(for.plot.true, for.plot.est)
  
  
  SAR.plot <- ggplot(data = for.plot, aes(x = Area, y = S)) +
    geom_vline(data = SAR, xintercept = SAR$Area, colour = "darkgrey") + 
    geom_line(aes(colour = as.factor(Time), linetype = Type), size = 1) +
    scale_x_log10() + scale_y_log10() +
    labs(title = "Species-Area Relationships", colour = "Time") + 
    theme_bw() +
    theme(legend.position = c(0.8,0.3)) +
    scale_colour_brewer(palette="Set1") 
  
  
  # 2. Estimated deltas
  Delta.plot <- ggplot(data = SAR, aes(x = Delta, y = Delta.mean.pred)) + 
  geom_point(aes(), size = 3) +
  scale_colour_viridis(trans= "log10") +
  geom_abline(intercept=0, slope=1) + 
  geom_hline(yintercept=0, linetype = 2) +
  geom_vline(xintercept=0, linetype = 2) +
  coord_fixed() +
  labs(title = "Observed vs Predicted Delta S", x="True Delta S", y="Estimated Delta S") +
  theme_bw() +
              theme(legend.position = c(0.8,0.3))


print(grid.arrange(SAR.plot, Delta.plot, ncol=2 ))

}
dev.off()
dev.off()
dev.off()
dev.off()

