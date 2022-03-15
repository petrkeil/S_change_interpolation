
source("simulation_functions.r")
source("ML_functions.r")

# -----------------------------------------------------------------------------

# set simulation paramters
prms <- set.simulation2(n.spec = 30, 
                        n.ind.clust = 40,
                        n.clust.tot = 200,
                        immigr.rate= 0,
                        SAD.trend =  0.5,
                        RSD.trend =  0.9)

# run the main simulation
comm <- simulate.comms(prms)

# calculate "true" scaling relationships
SARS <- SEGAR(comm.list = comm)

Truth <- data.frame(Area = SARS$Delta$Area,
                    S1 = SARS$SAR1$S.mean, 
                    S2 = SARS$SAR2$S.mean,
                    Gain = SARS$Gain$Gain.mean,
                    Loss = SARS$Loss$Loss.mean,
                    Delta = SARS$Delta$Delta)
Delta <- Delta[Delta$Area >= min(smp12$Area) ,]
Delta <- Delta[Delta$Area <= max(smp12$Area) ,]


list2env(comm, envir=.GlobalEnv); names(comm)
list2env(prms, envir=.GlobalEnv); names(prms)
list2env(SARS, envir=.GlobalEnv); names(SARS)

# create the polygon layers
PLS <- make.tess.polygons(prop1 = 0.4, prop2 = 0.4, overlap= FALSE)

# convert ppp patterns to sp objects
comm1 <- pplist.to.sp(PPlist1)
comm2 <- pplist.to.sp(PPlist2)

# count species richness in the polygons
smp1 <- count.S.in.poly(comm = comm1, poly = PLS$PLS1)
smp2 <- count.S.in.poly(comm2, PLS$PLS2)
smp12 <- rbind(data.frame(smp1, Time = 1), data.frame(smp2, Time = 2))

# species-by-site detailed polygon data
X <- count.spec.by.spec(comm = comm, poly = PLS)


# ------------------------------------------------------------------------------
# BOOSTED TREES - LOSS AND GAIN

MODEL1 <- gbm(occ ~ x + y + Area + Time + spec, 
              data = X, 
              interaction.depth = 5, 
              distribution = "bernoulli",
              n.minobsinnode = 3,
              n.trees = 500)


# ------------------------------------------------------------------------------

GRID <- prepare.grid() 
GRID <- data.frame(GRID, Cell = 1:nrow(GRID))
GRID <- GRID[GRID$Area >= min(smp12$Area) ,]
GRID <- GRID[GRID$Area <= max(smp12$Area) ,]
GRID <- data.frame(GRID, spec=rep(unique(X$spec), each = nrow(GRID)))

GRID1 <- data.frame(GRID, Time = rep(1, each = nrow(GRID)))
GRID2 <- data.frame(GRID, Time = rep(2, each = nrow(GRID)))


# ------------------------------------------------------------------------------

P1 <- predict(MODEL1, newdata = GRID1, n.trees= 500, type = "response")
P2 <- predict(MODEL1, newdata = GRID2, n.trees= 500, type = "response")

P.ext = P1 * (1-P2) # probability of extinction
P.gain = P2 * (1-P1) # probability of gain

T.all <- data.frame(P.ext, P.gain, P1, P2, GRID)

T.sum <- ddply(.data = T.all, 
               .variables=c("Cell"),
               .fun=summarize,
               est.S1 = sum(P1),
               est.S2 = sum(P2),
               est.Gain = sum(P.gain),
               est.Loss = sum(P.ext),
               est.Delta = est.Gain - est.Loss, 
               Area = mean(Area))

T.mean <- ddply(.data = T.sum,
                .variables = "Area",
                .fun = summarize,
                est.S1 = mean(est.S1),
                est.S2 = mean(est.S2),
                est.Gain = mean(est.Gain),
                est.Loss = mean(est.Loss),
                est.Delta = mean(est.Delta),
                Area = mean(Area))

all <- na.omit(left_join(Truth, T.mean, by= "Area"))

par(mfrow=c(1,4))
plot(S1 ~ Area, data = all, type = "l", col = "#e41a1c",
     ylim = c(min(c(all$S1, all$S2)), max(c(all$S1, all$S2))),
     ylab = "S", lwd=2)
lines(all$Area, all$S2, type = "l", col = "#377eb8", lwd = 2)
legend("bottomright", legend=c("Time 1", "Time 2"), lwd = c(2,2),
       col = c("#e41a1c","#377eb8"))

plot(Gain ~ Area, data = all, type = "b", col = "#4daf4a", 
     ylim = c(min(c(all$Gain, -all$Loss)), max(c(all$Gain, -all$Loss))))
lines(est.Gain ~ Area, data = all, col = "#4daf4a", lty = 2, type = "b")
abline(a = 0, b = 0)

plot(-Loss ~ Area, data = all, col = "red", type = "b",
     ylim = c(min(c(all$Gain, -all$Loss)), max(c(all$Gain, -all$Loss))))
lines(-est.Loss ~ Area, data = all,  col = "red", lty = 2, type = "b")
abline(a = 0, b = 0)

plot(Delta ~ Area, data = all, col = "grey", type = "b",
     ylim = c(min(c(all$Gain, -all$Loss)), max(c(all$Gain, -all$Loss))))
lines(est.Delta ~ Area, data = all, col = "grey", lty = 2, type = "b")
abline(a = 0, b = 0)


