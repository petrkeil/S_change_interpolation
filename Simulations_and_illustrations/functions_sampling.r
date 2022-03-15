
# ------------------------------------------------------------------------------
#' Calculates the "true" SAR and extinction- and gain-area relationships, 
#' based on continuous regular grids.

#' @param grains Numeric integer vector giving numbers of grid cells along a side
#' of a square domain. 
#' @param comm.list An object returned by the simulate.comms function. A list of lists
#' of density images and point patterns.


require(foreach) # use library for parallel computing
require(doParallel)


# --------------
GLS <- function(i, PPlist, G)
{
  rst <- raster(as.im(PPlist[[i]], dimyx = c(G, G)) )
  rst[rst > 1] <- 1
  return(rst)
}
# --------------

SEGAR <- function(grains = c(32, 16, 8, 4, 2), 
                  comm.list,   
                  area.span = NULL,
                  return.rasters = FALSE)
{
  if(is.null(area.span) == FALSE)
  {
    areas <- (1/grains)^2
    good1 <- areas >= area.span[1]
    good2 <- areas <= area.span[2]
    good <- good1*good2 == 1
    grains <- grains[good]
  }

  PPlist1 <- comm.list$PPlist1
  PPlist2 <- comm.list$PPlist2
  
  S1 <- S2 <- Loss <- Gain <- D <- list() # empty containers
  
  for(G in grains) # for each grain
  {
    Glist1 <- foreach(i = 1:length(PPlist1), .combine = stack) %dopar% 
    {
      GLS(i, PPlist1, G)    
    }
    
    Glist2 <- foreach(i = 1:length(PPlist2), .combine = stack) %dopar% 
    {
      GLS(i, PPlist2, G)    
    }

    gains <- losses <- Glist1 - Glist2
    gains[gains > -1] <- 0
    gains[gains == -1] <- 1
    losses[losses < 1] <- 0
    losses[losses == 1] <- 1
    
    # sum up the individual rasters at a given grain
    S1[[as.character(G)]] <- sum(Glist1)
    S2[[as.character(G)]] <- sum(Glist2)
    Gain[[as.character(G)]] <- sum(gains)
    Loss[[as.character(G)]] <- sum(losses)
    D[[as.character(G)]] <-  S2[[as.character(G)]] - S1[[as.character(G)]]
  }
  
  # calculate the *mean* scaling curves
  myfun <- function(RAST) mean(RAST[])
  SAR1.mean <- ldply(lapply(X=S1, FUN=myfun), .id = "Area")
  SAR2.mean <- ldply(lapply(X=S2, FUN=myfun), .id = "Area")
  Gain.mean <- ldply(lapply(X=Gain, FUN=myfun), .id = "Area")
  Loss.mean <- ldply(lapply(X=Loss, FUN=myfun), .id = "Area")
  Delta.mean <- data.frame(Area = Gain.mean$Area, Delta = Gain.mean[,2] - Loss.mean[,2])

  # convert number of cells to area
  SAR1.mean$Area <- SAR2.mean$Area <- Gain.mean$Area <- 
  Loss.mean$Area <- Delta.mean$Area <- 
  (1/as.numeric(as.character(SAR1.mean$Area)) )^2
  
  names(SAR1.mean)[2] <- "S1"
  names(SAR2.mean)[2] <- "S2"
  names(Gain.mean)[2] <- "Gain"
  names(Loss.mean)[2] <- "Loss"
  names(Delta.mean)[2] <- "Delta"
  
  SEGARS <- list(SAR1.mean, 
                 SAR2.mean, 
                 Gain.mean, 
                 Loss.mean,
                 Delta.mean)
  
  SEGARS <- Reduce(function(x, y) merge(x, y, all=TRUE, by = "Area"), SEGARS)
  
  if(return.rasters)
  {
    SEGARS <- list(Delta.rast = D, 
                   Gain.rast = Gain, 
                   Loss.rast = Loss,
                   means = SEGARS)  
  }
  
  return(SEGARS)
}


# ------------------------------------------------------------------------------
#' Count number of species in each polygon in a polygon layer

#' @param comm List of spatstat's point pattern 'ppp' objects, 
#' each reprsenting one species. 
#' @param poly sp's SpatialPolygons-class object with all the polygons. 

count.S.in.poly <- function(comm, poly)
{
  inventories <- over(x = poly, y = comm, returnList = TRUE)
  
  countfun <- function(x) nrow(unique(x))
  
  S <- lapply(X=inventories, FUN = countfun)
  S <- data.frame(S = unlist(S))
  
  Area <- data.frame(Area = gArea(poly, byid = TRUE))
  coords <- data.frame(coordinates(poly))
  names(coords) <- c("x","y")
  .id = rownames(S)
  
  return(data.frame(.id, coords, Area, S))
}


# ------------------------------------------------------------------------------

poly.over.sp <- function(sp, poly, sp.names)
{
  # if there are 0 points in the pattern
  if(length(sp) == 0)
  {
    ALL <- data.frame(.id = numeric(0),
                      spec = character(0),
                      occ = numeric(0),
                      x = numeric(0),
                      y = numeric(0),
                      Area = numeric(0))
    return(ALL)
  }
  
  sp@data <- data.frame(sp@data, occ = 1)
  
  X <- over(x = poly, y = sp, returnList = TRUE)
  for(poly.id in names(poly))
  {
    X[[poly.id]] <- unique(merge(X[[poly.id]], sp.names, by="spec", all.y = TRUE))
  }
  X <- ldply(X)
  X$occ[is.na(X$occ)] <- 0
  
  Area <- data.frame(Area = gArea(poly, byid = TRUE))
  Area <- data.frame(.id = rownames(Area), Area)
  coords <- data.frame(round(coordinates(poly),4))
  names(coords) <- c("x","y")
  coords <- data.frame(.id = rownames(coords), coords)
  ALL <- merge(X, coords, by=".id", all.x = TRUE)
  ALL <- merge(ALL, Area, by =".id", all.x = TRUE)  
  return(ALL)
}

#poly.over.sp(sp2, PLS2, sp.names)

# ------------------------------------------------------------------------------

count.spec.by.spec <- function(comm, poly)
{
  PPlist1 <- comm$PPlist1
  PPlist2 <- comm$PPlist2
  PLS1 <- poly$PLS1
  PLS2 <- poly$PLS2
  sp.names <- data.frame(spec = paste("sp", 1:length(PPlist1), sep=""))
  
  sp1 <- pplist.to.sp(PPlist1)
  sp2 <- pplist.to.sp(PPlist2)
  
  T1 <- poly.over.sp(sp1, PLS1, sp.names)
  T1 <- data.frame(T1, Time = rep(1, times = nrow(T1)))
  
  T2 <- poly.over.sp(sp2, PLS2, sp.names)
  T2 <- data.frame(T2, Time = rep(2, times = nrow(T2)))
  
  return(rbind(T1, T2))
}
#X <- count.spec.by.spec(comm = comm, poly = PLS)
# ------------------------------------------------------------------------------
