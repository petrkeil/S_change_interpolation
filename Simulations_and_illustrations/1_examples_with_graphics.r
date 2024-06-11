# MAKING ILLUSTRATIVE MAPS OF THE SIMULATION OF A SPECIES COMMUNITY IN TIME

# Author: Petr Keil 

# Description: Code that simulates a point pattern changing between at two time 
# periods. Each pattern represents a community of n.spec species. 
# The two communities are then plotted and exported as a .png file. Finally, 
# there is a little Linux script that creates a .gif animation out of the 
# two .png images.

# ------------------------------------------------------------------------------
# Loading the simulation functions from other files
source("functions_simulation.r")
source("functions_sampling.r")

# ------------------------------------------------------------------------------
# Setting up a list that contains parameters for the community dynamics simulation.
# This includes: Species Abundance Distributions (SADs) at two time periods,
# Range Size Distributions (RSDs) at two time periods.
prms <- set.simulation(n.spec = 24, 
                       n.ind.clust = 10,
                       n.clust.tot = 100,
                       immigr.rate= 1,
                       SAD.trend = 0,
                       RSD.trend = -0.80)

# Simulate the two communities
comm <- simulate.comms(prms)

# Calculate species-area relationships at both time periods, as well as the
# Loss, Gain, and Delta S
SARS <- SEGAR(comm.list = comm)

list2env(comm, envir=.GlobalEnv)
list2env(prms, envir=.GlobalEnv)
list2env(SARS, envir=.GlobalEnv)

# ------------------------------------------------------------------------------
# Plotting a spatially explicit image of each species, in both time periods.

png("../Graphics/time1.png", width=1400, height = 1400, res = 100)
  par(mfrow = c(5,5), mai=c(0.1, 0.1, 0.1, 1))
  for(i in 1:length(RSD1))
  {
    if(IMGlist1[[i]]$N == 0)
    {
      plot(IMGlist1[[i]], ribbon = FALSE, main=NULL, col = hcl.colors(100)[1])
    }
    else{ plot(IMGlist1[[i]], ribbon = FALSE, main=NULL, col = hcl.colors(100)) }
    plot(PPlist1[[i]], add=TRUE, col= "white", pch=19)
  }
  
  text(0.5, 0.5, label="Time1", col="white", cex = 2)
  par(mai=c(0.8,0.5,0.5,0.5))
  plot(SARS$Area, SARS$S1, type = "b", col = "black", lwd=2, ylim = c(0.01, 24))
  lines(SARS$Area, SARS$S2, type = "b", col="grey")
dev.off()


png("../Graphics/time2.png", width=1400, height = 1400, res = 100)
  par(mfrow = c(5,5), mai=c(0.1, 0.1, 0.1, 1))
  for(i in 1:length(RSD1))
  {
    if(IMGlist2[[i]]$N == 0)
    {
      plot(IMGlist2[[i]], ribbon = FALSE, main=NULL, col = hcl.colors(100)[1])
    }
    else{ plot(IMGlist2[[i]], ribbon = FALSE, main=NULL, col = hcl.colors(100)) }
    plot(PPlist2[[i]], add=TRUE, col= "white", pch=19)
  }
  text(0.5, 0.5, label="Time2", col="white", cex = 2)
  par(mai=c(0.8,0.5,0.5,0.5))
  plot(SARS$Area, SARS$S1, type = "b", col = "black", lwd=2, ylim = c(0.01, 24))
  lines(SARS$Area, SARS$S2, type = "b", col="grey")
dev.off()

# ------------------------------------------------------------------------------
# This is a script that converts the two .png images to an animation. It uses
# Ubuntu Linux's "convert" tool to do the animation.

system("
  cd ../graphics
  convert time1.png time2.png time1.png -delay 0 -morph 10 anim.gif
")

