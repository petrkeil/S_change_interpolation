comm.plot <- function(comm)
{
  extr <- function(X) data.frame(x = X$x, y = X$y)
  names(comm$PPlist1) <- 1:length(comm$PPlist1)
  names(comm$PPlist2) <- 1:length(comm$PPlist2)
  x1 <- ldply(.data = comm$PPlist1,
              .fun = extr,
              .id = "Species")
  x2 <- ldply(.data = comm$PPlist2,
              .fun = extr,
              .id = "Species")
  x1 <- data.frame(x1, Time = 1)
  x2 <- data.frame(x2, Time = 2)
  x <- rbind(x1, x2)
  x$Species <- paste("Species", x$Species)
  
  ggplot(data = x, aes(x = x, y = y)) +
    geom_point(aes(colour = as.factor(Time)), shape = 1) + 
    facet_wrap(.~Species) + 
    theme_minimal() +
    coord_fixed() + 
    theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
          axis.text.x=element_blank(),  axis.text.y=element_blank(),
          axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
    labs(colour = "Time")
  
}

# ------------------------------------------------------------------------------

polygon.ggplot <- function(PLS)
{
  
  pol1 <- ggplot(data =  fortify(PLS$PLS1), aes(x = long, y = lat, group = group)) + 
    geom_polygon(fill = "lightgrey", colour = "black") + theme_minimal() +
    coord_fixed() +  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                           axis.text.x=element_blank(),  axis.text.y=element_blank(),
                           axis.ticks.x=element_blank(), axis.ticks.y=element_blank())+
    labs(title = "Time 1")
  
  pol2 <- ggplot(data =  fortify(PLS$PLS2), aes(x = long, y = lat, group = group)) + 
    geom_polygon(fill = "lightgrey", colour = "black") + theme_minimal() +
    coord_fixed() +  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                           axis.text.x=element_blank(),  axis.text.y=element_blank(),
                           axis.ticks.x=element_blank(), axis.ticks.y=element_blank())+
    labs(title = "Time 2")
  
  grid.arrange(pol1, pol2, ncol=2)
}

# ------------------------------------------------------------------------------