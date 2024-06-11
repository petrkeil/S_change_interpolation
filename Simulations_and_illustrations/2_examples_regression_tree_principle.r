# CODE FOR FIGURE 1 OF THE MAIN TEXT

# Author: Petr Keil 

# Description: This is a code that produces graphics for all the panels of 
# figure 1 of the main text.

# ------------------------------------------------------------------------------
# necessary libraries

require(rpart) # regression trees
require(rpart.plot) # plotting regression trees
require(randomForest) 
require(gbm)
require(viridis) # colour palettes
require(RColorBrewer) # more colour palettes
require(gridExtra)
require(ggplot2)
require(pdp) # partial dependence plots from machine learning 

set.seed(5067) # set random seed for reproducibility

# ------------------------------------------------------------------------------
# setting up the variables and parameters to be used in the figure

Time <- c(1990, 2010) # two time periods (years)
Area <-  10^seq(-2, 1, by = 0.05) # sequence of areas
dat <- expand.grid(Time = Time, Area = Area) # time-area variable

C = 20 # intercept (in log-log) of the species-area relationship (SAR)
Z = -9.8 + dat$Time * 0.005 # time-dependent slope of the SAR
S <- C*(dat$Area^Z) # mean number of species
S.obs = rpois(n = nrow(dat), lambda = S) # observed number of species

# merging everything in a data.frame
dat <- data.frame(dat, 
                  S,
                  S.obs, 
                  lambda = S)

# ------------------------------------------------------------------------------
# FITTING A REGRESSION TREE AND A GLM OBJECT

# fitting the regression tree 
tr <- rpart(S ~ Area + Time, 
            data = dat)

# fitting the GLM
glm.object <- glm(S ~ log(Area)*Time, 
                  data = dat, 
                  family = poisson)

# fitting the random forest
rf <- randomForest(S ~ Area + Time, 
                   data = dat, 
                   mtry = 2, 
                   ntree = 500) 

# fitting the boosted tree
boost <- gbm(S ~ Area + Time, 
             data = dat, 
             n.trees = 500)

# extracting the predictions and merging them to a data.frame
dat <- data.frame(dat, 
                  S.pred.tr  = predict(tr), 
                  S.pred.glm = predict(glm.object, type = "response"),
                  S.pred.rf  = predict(rf),
                  S.pred.gbm = predict(boost,
                                       newdata = dat,
                                       type = "response", n.trees = n.trees))

# ------------------------------------------------------------------------------
# Panel (a) 2D CONTOUR PLOT

pdf("../Graphics/regtree_marginal.pdf", width = 5, height= 4)
  part.tr <- pdp::partial(tr, pred.var=c("Area", "Time"))
  plotPartial(part.tr, levelplot = TRUE, zlab = "yhat", drape = TRUE, 
              colorkey = TRUE, screen = list(z = 40, x = -60))
dev.off()

# ------------------------------------------------------------------------------
# Panel (b) REGRESSION TREE

# plot the regression tree
pdf("../Graphics/regtree_plot.pdf", width = 5, height= 4)
  rpart.plot(tr, type = 5, box.palette = viridis(20))
dev.off()


# ------------------------------------------------------------------------------
# panel (c) GLM FITTED CURVES

p.glm <- ggplot(dat, aes(x= Area, y = S)) +
  geom_point(aes(x=Area, y = S.obs, shape= as.factor(Time), colour = S.obs), size = 2) +
  scale_shape_manual(values=c(19, 1))+
  geom_line(aes(x = Area, y = S.pred.glm, linetype = as.factor(Time)), size = 1) +
  scale_colour_viridis()+
  labs(title="GLM", subtitle="glm(S ~ log(Area)*Time, family = poisson)")  +
  theme(legend.position = "none"); p.glm

# ------------------------------------------------------------------------------
# panel (d) REGRESSION TREE FITTED CURVES

p.tr <- ggplot(dat, aes(x= Area, y = S)) +
  geom_point(aes(x=Area, y = S.obs, shape= as.factor(Time), colour = S), size = 2) +
  scale_shape_manual(values=c(19, 1))+
  scale_colour_viridis()+
  geom_line(aes(x = Area, y = S.pred.tr, linetype = as.factor(Time)), size = 1) +
  theme(legend.position = "none") +
  labs(title="Regression tree", subtitle="rpart(S ~ Area + Time)"); p.tr

# ------------------------------------------------------------------------------
# panel (e) RANDOM FOREST FITTED CURVES

p.rf <- ggplot(dat, aes(x= Area, y = S)) +
  geom_point(aes(x=Area, y = S.obs, shape= as.factor(Time), colour = S), size = 2) +
  scale_shape_manual(values=c(19, 1))+
  scale_colour_viridis()+
  geom_line(aes(x = Area, y = S.pred.rf, linetype = as.factor(Time)), size = 1) +
  theme(legend.position = "none") +
  labs(title="Random forest", subtitle="randomForest(S ~ Area + Time)"); p.rf

# ------------------------------------------------------------------------------
# panel (f) BOOSTED TREES FITTED CURVES

p.gbm <- ggplot(dat, aes(x= Area, y = S)) +
  geom_point(aes(x=Area, y = S.obs, shape= as.factor(Time), colour = S), size = 2) +
  scale_shape_manual(values=c(19, 1))+
  scale_colour_viridis()+
  geom_line(aes(x = Area, y = S.pred.gbm, linetype = as.factor(Time)), size = 1) +
  theme(legend.position = "none") +
  labs(title="Boosted trees", subtitle="gbm(S ~ Area + Time)"); p.gbm

# ------------------------------------------------------------------------------
# Arranging the panels in a grid and exporting them to a file.

grid.arrange(p.glm, p.tr, p.rf, p.gbm, ncol=2, nrow=2)

pdf("../Graphics/ml_panels.pdf", width = 7, height = 7)
  grid.arrange(p.glm, p.tr, p.rf, p.gbm, ncol=2, nrow=2)
dev.off()















