
# Prepare regular prediction grids
prepare.grid <- function(grains = c(64, 32, 16, 8, 4, 2), 
                         area.span = NULL)
{
  res.multigrain <- list()
  for(grain in grains)
  {
    rast <- raster(as.im(square(), dimyx=c(grain, grain)))
    Area <- res(rast)[1]^2
    coords <- data.frame(coordinates(rast))
    res <- data.frame(coords, Area)
    res.multigrain[[grain]] <- res
  }  
  
  GRID <- do.call(what=rbind, res.multigrain)
  GRID1 <- data.frame(GRID, Time = rep(1, each = nrow(GRID)))
  GRID2 <- data.frame(GRID, Time = rep(2, each = nrow(GRID)))
  GRID <- rbind(GRID1, GRID2)
  
  if(is.null(area.span) == FALSE){
    GRID <- GRID[GRID$Area >= area.span[1] ,]
    GRID <- GRID[GRID$Area <= area.span[2] ,]
  }
  
  return(GRID)
}
# prepare.grid()

# ------------------------------------------------------------------------------

predict.SAR.mean <- function(SAR.model, 
                        grains = c(32, 16, 8, 4), 
                        area.span = NULL,
                        return.rasters = FALSE)
{
  GRID <- prepare.grid(grains, area.span)
  
  if("gbm" %in% class(SAR.model))
  {
    S.pred <- predict(SAR.model, 
                      newdata = GRID, 
                      n.trees = SAR.model$n.trees,
                      type = "response")
    modelclass <- "gbm"
  }
  if("randomForest" %in% class(SAR.model))
  {
    S.pred <- predict(SAR.model, 
                      newdata = GRID)
    modelclass <- "randomForest"
  }
  
  
  pred.GRID <- data.frame(GRID, S.pred)
  
  SAR12.long <- plyr::ddply(.data = pred.GRID,
                            .variables=c("Area", "Time"),
                            .fun = summarise,
                            Area = max(Area),
                            Time = max(Time),
                            S.pred = mean(S.pred))
  
  p1 <- SAR12.long[SAR12.long$Time == 1, ]
  p2 <- SAR12.long[SAR12.long$Time == 2, ]
  p12 <- data.frame(Area = p1$Area, 
                    S.pred1 = p1$S.pred,
                    S.pred2 = p2$S.pred,
                    Delta.pred = p2$S.pred - p1$S.pred)
  res <- list(modelclass = modelclass,
              wide = p12,
              long = SAR12.long)
  return(res)
}

# ------------------------------------------------------------------------------

predict.OAR.mean <- function(model.OAR, 
                             grains = c(32, 16, 8, 4), 
                             area.span = NULL)
{
  GRID <- prepare.grid(grains, area.span)
  GRID <- data.frame(GRID, Cell = 1:nrow(GRID))
  
  GRID.I <- data.frame(GRID, spec=rep(unique(X$spec), each = nrow(GRID)))
  G1 <- GRID.I[GRID.I$Time == 1,]
  G2 <- GRID.I[GRID.I$Time == 2,]
  
  if("gbm" %in% class(model.OAR))
  {
    P1 <- predict(model.OAR, 
                  newdata = G1, 
                  n.trees= model.OAR$n.trees,
                  type = "response")
    P2 <- predict(model.OAR, 
                  newdata = G2, 
                  n.trees= model.OAR$n.trees, 
                  type = "response")
    modelclass <- "gbm"
  }
  
  if("randomForest" %in% class(model.OAR))
  {
    P1 <- predict(model.OAR, 
                  newdata = G1,
                  type = "prob")[,2]
    P2 <- predict(model.OAR, 
                  newdata = G2,
                  type = "prob")[,2]
    modelclass <- "randomForest"
  }

  
  P.ext = P1 * (1-P2) # probability of extinction
  P.gain = P2 * (1-P1)

  T.all <- data.frame(P.ext, P.gain, P1, P2, G1) %>% select(-Time)
  
  T.sum <- ddply(.data = T.all, 
                 .variables=c("Cell"),
                 .fun=summarize,
                 est.S1 = sum(P1),
                 est.S2 = sum(P2),
                 est.Gain = sum(P.gain),
                 est.Loss = sum(P.ext),
                 est.Delta = est.Gain - est.Loss, 
                 Area = mean(Area))
  
  T.wide <- ddply(.data = T.sum,
                  .variables = "Area",
                  .fun = summarize,
                  Area = mean(Area),
                  Spred.1 = mean(est.S1),
                  Spred.2 = mean(est.S2),
                  Gain.pred = mean(est.Gain),
                  Loss.pred = mean(est.Loss),
                  Delta.pred = mean(est.Delta))
  
  T.long <- select(T.wide, -Gain.pred, -Loss.pred, -Delta.pred) %>%
    gather(key = Time, "Spred.2", -Area); T.long
  T.long$Time <- rep(c(1,2), each = nrow(T.wide))
  T.long <- dplyr::rename(T.long, Spred = Spred.2)
  
  res <- list(modelclass = modelclass,
              wide = T.wide,
              long = T.long)
  return(res)
}




# ------------------------------------------------------------------------------

fit.variants.gbm <- function(data, 
                             true.SEGAR, 
                             grains = c(32, 16, 8, 4), 
                             n.trees = 200,
                             area.span = NULL) 
{
  if(is.null(area.span) == FALSE)
  {
    areas <- (1/grains)^2
    grains <- grains[areas >= area.span[1]]
    grains <- grains[areas <= area.span[2]]
  }
  
  # minimum number of observations per node
  NM <- 3
  
  newdat <- prepare.grid(grains, area.span = area.span)
  # -----------------------
  message("- model 1")
  MODEL1 <- gbm(S ~ Area + Time + x + y , 
                data = data, 
                interaction.depth = 4, 
                distribution = "poisson",
                n.minobsinnode = NM,
                n.trees = n.trees)

  est.S1 <- predict(MODEL1, 
                    newdata = newdat[newdat$Time == 1, ],
                    n.trees = n.trees,
                    type = "response")
  
  est.S2 <- predict(MODEL1, 
                    newdata = newdat[newdat$Time == 2, ],
                    n.trees = n.trees,
                    type = "response")
  
  est.Delta <- est.S2 - est.S1
  new <- data.frame(newdat, est.Delta = est.Delta)
  new <- ddply(.data = new, 
               .variables="Area",
               .fun=summarize,
               est.Delta = mean(est.Delta))
  
  PRD1 <- left_join(true.SEGAR, new, by = "Area")
  PRD1 <- data.frame(PRD1, Formula = "Area + Time + x + y")
  
  # ---------------------------------------------------------
  message("- model 2")
  MODEL2 <- gbm(S ~ Area + Time, 
                data = data, 
                interaction.depth = 4, 
                distribution = "poisson",
                n.minobsinnode = NM,
                n.trees = n.trees)

  est.S1 <- predict(MODEL2, 
                    newdata = newdat[newdat$Time == 1, ],
                    n.trees = n.trees,
                    type = "response")
  est.S2 <- predict(MODEL2, 
                    newdata = newdat[newdat$Time == 2, ],
                    n.trees = n.trees,
                    type = "response")
  est.Delta <- est.S2 - est.S1
  new <- data.frame(newdat, est.Delta = est.Delta)
  new <- ddply(.data = new, 
               .variables="Area",
               .fun=summarize,
               est.Delta = mean(est.Delta))
  
  PRD2 <- left_join(true.SEGAR, new, by = "Area")
  PRD2 <- data.frame(PRD2, Formula = "Area + Time ")
  
  # ---------------------------------------------------------
  message("- model 3")
  MODEL3t1 <- gbm(S ~ Area, 
                  data = data[data$Time == 1,], 
                  interaction.depth = 4, 
                  distribution = "poisson",
                  n.minobsinnode = NM,
                  n.trees = n.trees)
  
  MODEL3t2 <- gbm(S ~ Area, 
                  data = data[data$Time == 2,], 
                  interaction.depth = 4, 
                  distribution = "poisson",
                  n.minobsinnode = NM,
                  n.trees = n.trees)

  est.S1 <- predict(MODEL3t1, 
                    newdata = newdat[newdat$Time == 1, ],
                    n.trees = n.trees,
                    type = "response")
  est.S2 <- predict(MODEL3t2, 
                    newdata = newdat[newdat$Time == 2, ],
                    n.trees = n.trees,
                    type = "response")
  est.Delta <- est.S2 - est.S1
  new <- data.frame(newdat, est.Delta = est.Delta)
  new <- ddply(.data = new, 
               .variables="Area",
               .fun=summarize,
               est.Delta = mean(est.Delta))
  
  PRD3 <- left_join(true.SEGAR, new, by = "Area")
  PRD3 <- data.frame(PRD3, Formula = "Area, Time separate")
  
  # ---------------------------------------------------------
  message("- model 4")
  MODEL4t1 <- gbm(S ~ Area + x + y, 
                  data = data[data$Time == 1,], 
                  interaction.depth = 4, 
                  distribution = "poisson",
                  n.minobsinnode = NM,
                  n.trees = n.trees)
  MODEL4t2 <- gbm(S ~ Area + x + y, 
                  data = data[data$Time == 2,], 
                  interaction.depth = 4, 
                  distribution = "poisson",
                  n.minobsinnode = NM,
                  n.trees = n.trees)

  
  est.S1 <- predict(MODEL4t1, 
                    newdata = newdat[newdat$Time == 1, ],
                    n.trees = n.trees,
                    type = "response")
  est.S2 <- predict(MODEL4t2, 
                    newdata = newdat[newdat$Time == 1, ],
                    n.trees = n.trees,
                    type = "response")
  est.Delta <- est.S2 - est.S1
  new <- data.frame(newdat, est.Delta = est.Delta)
  new <- ddply(.data = new, 
               .variables="Area",
               .fun=summarize,
               est.Delta = mean(est.Delta))
  
  PRD4 <- left_join(true.SEGAR, new, by = "Area")
  PRD4 <- data.frame(PRD4, Formula = "Area + x + y, Time separate")
  
  # ---------------------------------------------------------
  PRD.all <- rbind(PRD1, PRD2, PRD3, PRD4)
  PRD.all <- PRD.all[PRD.all$Area <= max(data$Area),]
  return(PRD.all)
}



# ------------------------------------------------------------------------------

fit.variants.rf <- function(data, 
                            true.SEGAR, 
                            grains = c(32, 16, 8, 4), 
                            n.trees = 500,
                            area.span = NULL) 
{
  if(is.null(area.span) == FALSE)
  {
    areas <- (1/grains)^2
    grains <- grains[areas >= area.span[1]]
    grains <- grains[areas <= area.span[2]]
  }
  
  # number of variables to omit in each rf resampling
  vars.omit <- 1
  
  newdat <- prepare.grid(grains, area.span = area.span)
  # -----------------------
  message("- model 1")
  MODEL1 <- randomForest(S ~ Area + Time + x + y , 
                data = data, 
                mtry = 4 - vars.omit, 
                ntree = n.trees)
  
  est.S1 <- predict(MODEL1, 
                    newdata = newdat[newdat$Time == 1, ])
  
  est.S2 <- predict(MODEL1, 
                    newdata = newdat[newdat$Time == 2, ])
  
  est.Delta <- est.S2 - est.S1
  new <- data.frame(newdat, est.Delta = est.Delta)
  new <- ddply(.data = new, 
               .variables="Area",
               .fun=summarize,
               est.Delta = mean(est.Delta))
  
  PRD1 <- left_join(true.SEGAR, new, by = "Area")
  PRD1 <- data.frame(PRD1, Formula = "Area + Time + x + y")
  
  # ---------------------------------------------------------
  message("- model 2")
  MODEL2 <-  randomForest(S ~ Area + Time, 
                                    data = data, 
                                    mtry = 2 - vars.omit, 
                                    ntree = n.trees)
  
  est.S1 <- predict(MODEL2, 
                    newdata = newdat[newdat$Time == 1, ])
  est.S2 <- predict(MODEL2, 
                    newdata = newdat[newdat$Time == 2, ])
  est.Delta <- est.S2 - est.S1
  new <- data.frame(newdat, est.Delta = est.Delta)
  new <- ddply(.data = new, 
               .variables="Area",
               .fun=summarize,
               est.Delta = mean(est.Delta))
  
  PRD2 <- left_join(true.SEGAR, new, by = "Area")
  PRD2 <- data.frame(PRD2, Formula = "Area + Time ")
  
  # ---------------------------------------------------------
  message("- model 3")
  MODEL3t1 <- randomForest(S ~ Area, 
                  data = data[data$Time == 1,],
                  mtry = 1, 
                  ntree = n.trees)
  
  MODEL3t2 <- randomForest(S ~ Area, 
                  data = data[data$Time == 2,],
                  mtry = 1, 
                  ntree = n.trees)
  
  est.S1 <- predict(MODEL3t1, 
                    newdata = newdat[newdat$Time == 1, ])
  est.S2 <- predict(MODEL3t2, 
                    newdata = newdat[newdat$Time == 2, ])
  
  est.Delta <- est.S2 - est.S1
  new <- data.frame(newdat, est.Delta = est.Delta)
  new <- ddply(.data = new, 
               .variables="Area",
               .fun=summarize,
               est.Delta = mean(est.Delta))
  
  PRD3 <- left_join(true.SEGAR, new, by = "Area")
  PRD3 <- data.frame(PRD3, Formula = "Area, Time separate")
  
  # ---------------------------------------------------------
  message("- model 4")
  MODEL4t1 <- randomForest(S ~ Area + x + y,
                  data = data[data$Time == 1,],
                  mtry = 3 - vars.omit, 
                  ntree = n.trees)
  MODEL4t2 <- randomForest(S ~ Area + x + y, 
                  data = data[data$Time == 2,],
                  mtry = 3 - vars.omit, 
                  ntree = n.trees)
  
  
  est.S1 <- predict(MODEL4t1, 
                    newdata = newdat[newdat$Time == 1, ])
  est.S2 <- predict(MODEL4t2, 
                    newdata = newdat[newdat$Time == 1, ])
  est.Delta <- est.S2 - est.S1
  new <- data.frame(newdat, est.Delta = est.Delta)
  new <- ddply(.data = new, 
               .variables="Area",
               .fun=summarize,
               est.Delta = mean(est.Delta))
  
  PRD4 <- left_join(true.SEGAR, new, by = "Area")
  PRD4 <- data.frame(PRD4, Formula = "Area + x + y, Time separate")
  
  # ---------------------------------------------------------
  PRD.all <- rbind(PRD1, PRD2, PRD3, PRD4)
  PRD.all <- PRD.all[PRD.all$Area <= max(data$Area),]
  return(PRD.all)
}


# ------------------------------------------------------------------------------

predict.lossgain.raster <- function(brt.model, rast, specnames, n.trees = 200)
{
  Area <- (res(rast)[1])^2
  xy <- coordinates(rast)
  Col <- colFromX(rast, xy[,'x'])
  Row <- rowFromY(rast, xy[,'y'])
  Cell <- 1:ncell(rast)
  xyA <- data.frame(Area, xy, Col, Row, Cell)
  Time.spec1 <- expand.grid(Time = 1, spec = specnames)
  Time.spec2 <- expand.grid(Time = 2, spec = specnames)
  new.1 <- merge(xyA, Time.spec1, by = NULL)
  new.2 <- merge(xyA, Time.spec2, by = NULL)
  
  # prediction
  prds1 <- data.frame(yhat = predict(brt.model, 
                                     newdata = new.1, 
                                     type = "response", 
                                     n.trees = n.trees), new.1)
  
  prds2 <- data.frame(yhat = predict(brt.model, 
                                     newdata = new.2, 
                                     type = "response", 
                                     n.trees = n.trees), new.2)

  P1 <- prds1$yhat
  P2 <- prds2$yhat
  
  P.ext = P1 * (1-P2) # probability of extinction
  P.gain = P2 * (1-P1) # probability of gain

  T.all <- data.frame(P.ext, P.gain, P1, P2, new.1)
  
  T.sum <- ddply(.data = T.all, 
                 .variables=c("Cell"),
                 .fun=summarize,
                 S1 = sum(P1),
                 S2 = sum(P2),
                 est.Gain = sum(P.gain),
                 est.Loss = (-1)*sum(P.ext),
                 est.Delta = est.Gain + est.Loss)
  
  TT <- left_join(xyA, T.sum, by = "Cell")
  
  S1 <- S2 <- Loss <- Gain <- Delta <- rast
  S1[as.matrix(cbind(TT$Row, TT$Col))] <- TT$S1
  S2[as.matrix(cbind(TT$Row, TT$Col))] <- TT$S2
  Loss[as.matrix(cbind(TT$Row, TT$Col))] <- TT$est.Loss
  Gain[as.matrix(cbind(TT$Row, TT$Col))] <- TT$est.Gain
  Delta[as.matrix(cbind(TT$Row, TT$Col))] <- TT$est.Delta
  
  res <- stack(list(est.S1 = S1,
                    est.S2 = S2,
                    est.Loss = Loss,
                    est.Gain = Gain,
                    est.Delta = Delta))
  return(res)

}

##rast <- raster(as.im(square(), dimyx=c(8, 8)))
#xx.pred <- predict.lossgain.raster(brt.model = MODEL, 
##                        rast=rast, 
#                        specnames = as.character(unique(X$spec)))

#plot(xx.pred)

# ------------------------------------------------------------------------------
