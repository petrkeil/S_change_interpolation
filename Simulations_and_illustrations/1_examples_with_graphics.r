source("simulation_functions.r")

prms <- set.simulation2(n.spec = 24, 
                        n.ind.clust = 10,
                        n.clust.tot = 100,
                        immigr.rate= 1,
                        SAD.trend = 0,
                        RSD.trend = -0.80)

comm <- simulate.comms(prms)
SARS <- SEGAR(comm.list = comm)

list2env(comm, envir=.GlobalEnv)
list2env(prms, envir=.GlobalEnv)
list2env(SARS, envir=.GlobalEnv)

# ------------------------------------------------------------------------------
# extract things for figures

# PLOTTING THE SIMULATION
png("../graphics/time1.png", width=1400, height = 1400, res = 100)
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
  plot(SAR1, type = "b", col = "black", lwd=2, ylim = c(0.01, 24))
  lines(SAR2, type = "b", col="grey")
dev.off()


png("../graphics/time2.png", width=1400, height = 1400, res = 100)
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
  plot(SAR1, type = "b", col = "grey", ylim = c(0.01, 24))
  lines(SAR2, type = "b", col="black", lwd=2)
dev.off()

# ------------------------------------------------------------------------------

system("
  cd ../graphics
  convert time1.png time2.png time1.png -delay 0 -morph 10 anim.gif
")

