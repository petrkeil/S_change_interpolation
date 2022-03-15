# Functions for simulation and machine learning-based analysis of biodiversity
# change in spatially misaligned polygons.

# Author: Petr Keil 
# pkeil@seznam.cz

# ------------------------------------------------------------------------------

# LIBRARIES
 
library(mvtnorm)
library(truncnorm)
library(mobsim)
# ------------------------------------------------------------------------------
library(raster)
library(spatstat)
library(sp)
library(maptools) # for converstion of tesseletions to SpatialPolygons
library(rgeos)
# ------------------------------------------------------------------------------
library(randomForest)
library(gbm) # gradient boosting - boosted regression trees
library(pdp) # partial dependence plots for machine learning
# ------------------------------------------------------------------------------
library(ggplot2)
library(tidyverse)
library(plyr)




# ------------------------------------------------------------------------------
#' Bivariate probablity density function in a 2D square

#' @param x x Coordinates of points whose PD should be evaluated.
#' @param y y Coordinates of points whose PD should be evaluated.
#' @param var Variance of the PDF, which is the same in x and y direction. Note: covariance is 0.
#' @param x.centr x coordinate of the mean in the 2D space
#' @param y.centr y coordinate of the mean in the 2D space


dpoint.MVN <- function(x, y, var, x.centr, y.centr)
{
  sigma = matrix(c(var, 0, 0, var), nrow=2, ncol=2)
  dmvnorm(x=cbind(x,y), mean = c(x.centr, y.centr), sigma)
}

#dpoint.MVN (0.1, 0.1, var=0.1, x.centr=0.1, y.centr=0.1)

# ------------------------------------------------------------------------------
#' Function that generates N clusters of n points in a 2D domain, each
#' cluster with a given variance around its mother point.

#' @param var bivariate variance of a single cluster
#' @param N number of clusters
#' @param n number of points per cluster
#' @export

image.MVN <- function(var, N, centers = NULL)
{
  if(is.null(centers)) # are centers of clusters provided?
  { 
    # centers of the clusters (= "mother points")
    x <- runif(N);  y <- runif(N)
    centers <- data.frame(x, y)
  }
  
  if(N == 0) # if we want a zero density image (with no clusters)
  {
    D <- as.im(dpoint.MVN, 
               var = var, 
               x.centr = 0.5, 
               y.centr = 0.5, 
               W = square())
    D[] <- 0
  }
  
  if (N >= 1) # initialize with the first cluster
  {
    D <- as.im(dpoint.MVN, 
               var = var, 
               x.centr = centers$x[1], 
               y.centr = centers$y[1], 
               W = square())
  }
  
  if (N > 1) # and continue adding clusters, if appropriate
  {
    # generate each cluster, one by one
    for(i in 2:N)
    {
      # generate the point cluster with density given by dpoint.MVN, and with
      # given coordinates
      D <- D + as.im(dpoint.MVN, 
                     var = var, 
                     x.centr = centers$x[i], 
                     y.centr = centers$y[i], 
                     W = square())
    }
  }
  
  # save parameters of the simulation
  D[["centers"]] <- centers
  D[["var"]] <- var
  D[["N"]] <- N
  
  return(D)
}

# D <- image.MVN(var=0.001,  N=2)
# plot(D)

# ------------------------------------------------------------------------------
#' Add or remove clusters from a density D image object

#' @param D An object returned by image.MVN function. This is a 
#' spatstat 'im' object, enhanced by information on clusters. 
#' @param N Number of clusters to remove or add. If N < 0, clusters are removed,
#' if N > 0, clusters are added.

change.clusters <- function(D, N)
{
  if(N == 0){ return(D) } # if no cluster is added
  
  if(N > 0) # add clusters
  {
    # choose cluster centers to add
    D.to.add <- image.MVN(var = D$var, N = N)
    
    D2 <- D +  D.to.add
    D2[["centers"]] <- rbind(D$centers, D.to.add$centers)
    D2[["var"]] <- D$var
    D2[["N"]] <- D$N + N
  }
  
  if(N < 0) # remove clusters
  {
    N <- N*(-1)
    if(N == D$N) # if complete extinction
    { 
      D2 <- D; D2[] <- 0 
      D2[["centers"]] <- NULL
      D2[["var"]] <- NULL
      D2[["N"]] <- 0
    } 
    
    else
    {
      # choose cluster centers to remove
      id.to.remove <- sample(x = 1:nrow(D$centers), size = N)
      centers.to.remove <- D$centers[id.to.remove,]
      id.to.keep   <- which( 1:nrow(D$centers) %in% id.to.remove == FALSE )
      sub.D <- image.MVN(var = D$var, N = N, centers = centers.to.remove)
      
      # remove the clusters
      D2 <- D - sub.D
      
      D2[["centers"]] <- rbind(D$centers[id.to.keep,])
      D2[["var"]] <- D$var
      D2[["N"]] <- D$N - N
    }
  }
  
  #get rid of potential tiny accidental negative probs 
  D2[D2 < 0] <- 0
  
  return(D2)
}

# par(mfrow = c(1,2)); plot(D); plot(change.clusters(D, 200))


# ------------------------------------------------------------------------------
# SIMULATE SPATIALLY IMPLICIT DISTRIBUTIONS OF THE SIMULATIONS

set.simulation <- function(n.spec = 100, 
                           n.ind.clust = 10, 
                           n.clust.tot = 50,
                           immigr.rate = 0.5,
                           SAD.trend = 0,
                           RSD.trend = 0)
{ 
  # STEP 1 - SAD1, RSD1, VARD - species-abundance, 
  #                             range-size distributions,
  #                             and var distributions
  RSD1 <- as.numeric(mobsim::sim_sad(s_pool = n.spec, 
                                     n_sim = n.clust.tot, 
                                     sad_type = "lnorm"))
  RSD1 <- c(RSD1, rep(0, times = n.spec - length(RSD1)))
  SAD1 <- rpois(length(RSD1), lambda = RSD1*n.ind.clust) # abundances are proportional to range sizes
  
  # sizes of clusters
  VARD <- rexp(n = n.spec, rate = 40)
  
  # immigration
  IMMIGR.candidates <- 1* ( RSD1 == 0 )
  IMMIGR <- rbinom(n = length(IMMIGR.candidates), prob=immigr.rate, size= 1) * IMMIGR.candidates
  
  # delta coefficients and total losses
  deltaRSD <- rtruncnorm(n = n.spec, a = -1, b = 1, mean = RSD.trend, sd = 0.1)
  diffRSD <- round(RSD1 * deltaRSD) + IMMIGR
  
  deltaSAD <- rtruncnorm(n = n.spec, a = -1, b = 1, mean = SAD.trend, sd = 0.1)
  diffSAD <- round(SAD1 * deltaSAD + rpois(length(RSD1), lambda = IMMIGR * n.ind.clust)) 
  
  # STEP 2 - SAD2 and RSD2 - species-abundance and range-size distributions
  SAD2 <- SAD1 + diffSAD
  RSD2 <- RSD1 + diffRSD + IMMIGR
  
  params <- list(#PARAMS,
    SAD1 = SAD1, SAD2 = SAD2, 
    RSD1 = RSD1, RSD2 = RSD2, 
    VARD = VARD, 
    diffRSD = diffRSD)
  
  return(params)
}
# prms <- set.simulation2(n.spec = 100)



# ------------------------------------------------------------------------------
#' Simulate multi-species point pattern communities in 2 time periods

#' @param params Parameters returned by the set.simulation functions. 

simulate.comms <- function(params, spat.gradient = 1)
{
  IMGlist1 <- IMGlist2 <- list()
  PPlist1  <- PPlist2 <- list()
  
  # extract the specific parameters 
  VARD <- params$VARD
  RSD1 <- params$RSD1
  RSD2 <- params$RSD2
  SAD1 <- params$SAD1
  SAD2 <- params$SAD2
  diffRSD <- params$diffRSD
  
  for(i in 1:length(RSD1)) # for each species
  {
    
    IMG1 <- image.MVN(var = VARD[i], N = RSD1[i], centers = NULL)
    
    if(RSD1[i] > 0) # if there are at least some clusters
    {
      PP1   <- rpoint(n = SAD1[i], f = IMG1)
    }
    if(RSD1[i] == 0) # if there are no clusters 
    {
      PP1   <- rpoint(n = 0, f = 0, win = owin(xrange = IMG1$xrange, IMG1$yrange))
    }
    
    IMG2 <- change.clusters(IMG1, N = diffRSD[i])
    
    if(RSD2[i] > 0) # if there are at least some clusters
    {
      PP2 <-  rpoint(n = SAD2[i], f = IMG2)
    }
    if(RSD2[i] == 0) # if there are no clusters 
    {
      PP2   <- rpoint(n = 0, f = 0, win = owin(xrange = IMG1$xrange, IMG1$yrange))
    }
    
    IMGlist1[[i]] <- IMG1
    IMGlist2[[i]] <- IMG2
    PPlist1[[i]] <- PP1
    PPlist2[[i]] <- PP2
  }
  
  comm.list <- list(PPlist1 = PPlist1,   
                    PPlist2 = PPlist2,
                    IMGlist1 = IMGlist1, 
                    IMGlist2 = IMGlist2,
                    params = params)
  
  return(comm.list)
}

# comm <- simulate.comms(prms)

# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
#' Convert spatstat's point patterns to sp's SpatialPoints

#' @param pplist List of spatstat's ppp objects. 

pplist.to.sp <- function(pplist)
{
  sp.list <- list()
  for(i in 1:length(pplist))
  {
    #print(i)
    if(pplist[[i]]$n > 0)
    {
      PP.sp <- as(object=pplist[[i]], "SpatialPoints")
      PP.sp <- SpatialPointsDataFrame(PP.sp, 
                                      data = data.frame(spec = rep(paste("sp", i, sep = ""), 
                                                                 times = length(PP.sp))))
    }
    else
    {
      PP.sp <- SpatialPoints(data.frame(x = 0, y = 0))[-1,]
      PP.sp <- SpatialPointsDataFrame(PP.sp,
                                      data = data.frame(spec = rep(paste("sp", i, sep = ""), 
                                                                   times = length(PP.sp))))
    }
    
    sp.list[[i]] <- PP.sp
  }  
  if(sum(unlist(lapply(sp.list, FUN=length))) > 0)
  {
    sp.list <- do.call(rbind, sp.list)
  }
  else
  {
    sp.list <- SpatialPoints(data.frame(x = 0, y = 0))[-1,]
    sp.list <- SpatialPointsDataFrame(sp.list,
                                      data = data.frame(spec = rep(paste("sp", i, sep = ""), 
                                                                 times = length(PP.sp))))  
  }
  
  return(sp.list)
}


# pplist.to.sp(pplist = comm$PPlist2)

# ------------------------------------------------------------------------------
#' Create set of N1 and N2 polygons in time 1 and time 2 respectively

#' @param prop1 Value between 0 and 1. A fraction (proportion) of polygons at time 1.
#' @param prop2 Value between 0 and 1. A fraction of polygons at time 2.
#' @param overlap Logical. Are the polygons sampled with replacement? 

make.tess.polygons <- function(prop1, prop2, overlap = FALSE)
{
  repeat # try to generate the polygons until success - this is necessary since at low
    # proportions the algorithm may fail 
  { 
    success <- TRUE
    tryCatch(
      {
        X <- rThomas(kappa = 5,  scale = 0.08, mu = 50)
        
        X <- ppp(x=X$x, y =X$y, window=square())
        PLS <- dirichlet(X)# ; plot(PLS)
        
        N1 <- round(PLS$n * prop1)
        N2 <- round(PLS$n * prop2)
        
        # if N1 or N2 is zero, sample 1
        if(N1 == 0) N1 <- 1
        if(N2 == 0) N2 <- 1
        
        # sample tiles
        all.id <- 1:length(PLS$tiles)
        id1 <- sample(all.id, size=N1) # the first sample in time 1
        
        if(overlap) # allow for overlap
        {
          id2 <- sample(all.id, size=N2)
        }
        if(overlap == FALSE) # spatially independent samples
        {
          id2 <- sample(which(all.id %in% id1 == FALSE), size=N2)
        }
        
        PLS1 <- PLS[id1]
        PLS2 <- PLS[id2]
        
        # convert tesselation to polygon object
        PLS       <- as(PLS, "SpatialPolygons")
        PLS1.poly <- as(PLS1, "SpatialPolygons")
        PLS2.poly <- as(PLS2, "SpatialPolygons")
        
        res.poly <- list(PLS1 = PLS1.poly, PLS2 = PLS2.poly, PLS.complete = PLS)
        
      }, error = function(cond) assign("success", FALSE, envir = .GlobalEnv ))
    if(success) break
  }
  
  return(res.poly)
}
 x <- make.tess.polygons(prop1 = 0.5, prop2 = 0.1, overlap = TRUE)
# par(mfrow = c(1,2)); plot(x$PLS1); plot(x$PLS2)













