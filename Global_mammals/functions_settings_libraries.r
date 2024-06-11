# Functions and libraries that will be used in the analysis of empirical 
# mammal data.

# Author: Petr Keil

rm(list = ls()) # clean the environment  

# GIS packages
library(raster)
library(sp)
library(rgdal)
library(rgeos)

# machine learning
library(gbm)
library(randomForest)

# other  
library(ggplot2)
library(tidyverse)
library(reshape)
library(plyr)
library(gridExtra)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(latex2exp)

# useful projections
BEHRMANN <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
WGS84 <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# blank theme for ggplot2 (for maps without figure borders etc.)
blank.theme <- theme(axis.line=element_blank(),axis.text.x=element_blank(),
                     axis.text.y=element_blank(),axis.ticks=element_blank(),
                     axis.title.x=element_blank(),
                     axis.title.y=element_blank(),
                     legend.position=c(0.68, 0),
                     legend.direction = "horizontal",
                     #legend.title = element_blank(),
                     plot.title = element_text(face=quote(bold)),
                     #legend.title.align = 0,
                     panel.background=element_blank(),
                     panel.border=element_blank(),panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),plot.background=element_blank())
