require(gridExtra)
require(RColorBrewer)
require(ggplot2)
require(pdp) # partial dependence plots

# ------------------------------------------------------------------------------

Time <- c(1990, 2010)
Area <- c(0.01, 0.1, 1, 10)
Area <-  10^seq(-2, 1, by = 0.05)
dat <- expand.grid(Time = Time, Area = Area)

set.seed(5067)
C = 20
Z = -9.8 + dat$Time * 0.005
S <- C*(dat$Area^Z)
S.obs = rpois(n = nrow(dat), lambda = S)

dat <- data.frame(dat, 
                  S,
                  S.obs, lambda = S)

glm(S.obs ~ Area, family = "poisson", data=dat)

p1 <- ggplot(dat, aes(x= Area, y = S)) +
      geom_point(aes(x=Area, y = S.obs, colour = as.factor(Time)), 
                 shape = 1) +
      geom_line(aes(colour = as.factor(Time)), size = 1) +
      ggtitle("Simulated truth") +
      theme_bw() +
     # scale_x_log10() + scale_y_log10() +
      theme(legend.position = "none") +
      scale_colour_brewer(palette="Set1") 
p1.log <- ggplot(dat, aes(x= Area, y = S)) +
      geom_point(aes(x=Area, y = S.obs, colour = as.factor(Time)), 
                 shape = 1) +
      geom_line(aes(colour = as.factor(Time)), size = 1) +    
      ggtitle("Simulated truth") +
      theme_bw() +
      scale_x_log10() + scale_y_log10() +
      theme(legend.position = "none") +
      scale_colour_brewer(palette="Set1") 

# ------------------------------------------------------------------------------
# FITTING GLM

glm.object <- glm(S ~ log(Area)*Time, data = dat, family = poisson)


part.glm <- partial(glm.object, pred.var=c("Area", "Time"))
plotPartial(part.glm, levelplot = TRUE, zlab = "yhat", drape = TRUE, 
            colorkey = TRUE, screen = list(z = 40, x = -60))


# ------------------------------------------------------------------------------
# FITTING REGRESSION TREE
require(rpart)
require(rpart.plot)
require(viridis)

tr <- rpart(S ~ Area + Time, data = dat)
plot(tr)
text(tr)
labels(tr)


# 2D contour partial plot
part.tr <- partial(tr, pred.var=c("Area", "Time"))

pdf("../Graphics/regtree_marginal.pdf", width = 5, height= 4)
plotPartial(part.tr, levelplot = TRUE, zlab = "yhat", drape = TRUE, 
            colorkey = TRUE, screen = list(z = 40, x = -60))
dev.off()


# plot the regression tree
pdf("../Graphics/regtree_plot.pdf", width = 5, height= 4)
  rpart.plot(tr,type = 5, box.palette = viridis(20))
dev.off()

predict(tr)
dat <- data.frame(dat, S.pred.tr = predict(tr), S.pred.glm = predict(glm.object, type = "response"))

p.tr <- ggplot(dat, aes(x= Area, y = S)) +
  geom_point(aes(x=Area, y = S.obs, shape= as.factor(Time), colour = S), size = 2) +
  scale_shape_manual(values=c(19, 1))+
  scale_colour_viridis()+
    geom_line(aes(x = Area, y = S.pred.tr, linetype = as.factor(Time)), size = 1) +
  theme(legend.position = "none") +
  labs(title="Regression tree", subtitle="rpart(S ~ Area + Time)"); p.tr

p.tr.log <- ggplot(dat, aes(x= Area, y = S)) +
  geom_point(aes(x=Area, y = S.obs, shape= as.factor(Time), colour = S), size = 2) +
  scale_shape_manual(values=c(19, 1))+
  scale_colour_viridis()+
  geom_line(aes(x = Area, y = S.pred.tr, linetype = as.factor(Time)), size = 1) +
  labs(title="rpart(S ~ Area + Time)", subtitle=expression(paste(Log[10], " axes")))  +
  scale_x_log10() + scale_y_log10() +
  theme(legend.position = "none"); p.tr.log

p.glm.log <- ggplot(dat, aes(x= Area, y = S)) +
  geom_point(aes(x=Area, y = S.obs, shape= as.factor(Time), colour = S), size = 2) +
  scale_shape_manual(values=c(19, 1))+
  scale_colour_viridis()+
  geom_line(aes(x = Area, y = S.pred.glm, linetype = as.factor(Time)), size = 1) +
  labs(title="glm(S ~ log(Area)*Time, family = poisson)", subtitle=expression(paste(Log[10], " axes")))  +
  scale_x_log10() + scale_y_log10() +
  theme(legend.position = "none"); p.glm.log

p.glm <- ggplot(dat, aes(x= Area, y = S)) +
  geom_point(aes(x=Area, y = S.obs, shape= as.factor(Time), colour = S.obs), size = 2) +
  scale_shape_manual(values=c(19, 1))+
  geom_line(aes(x = Area, y = S.pred.glm, linetype = as.factor(Time)), size = 1) +
  scale_colour_viridis()+
  labs(title="GLM", subtitle="glm(S ~ log(Area)*Time, family = poisson)")  +
  theme(legend.position = "none"); p.glm


# ------------------------------------------------------------------------------
# RANDOM FOREST
require(randomForest)

rf <- randomForest(S ~ Area + Time, 
                   data = dat, 
                   mtry = 2, # I only have 2 variables, so I am getting smoother fit with this
                   ntree = 500) 
importance(rf)
varImpPlot(rf, type = 2)
plot(rf)

dat <- data.frame(dat, S.pred.rf = predict(rf))

p.rf <- ggplot(dat, aes(x= Area, y = S)) +
  geom_point(aes(x=Area, y = S.obs, shape= as.factor(Time), colour = S), size = 2) +
  scale_shape_manual(values=c(19, 1))+
  scale_colour_viridis()+
  geom_line(aes(x = Area, y = S.pred.rf, linetype = as.factor(Time)), size = 1) +
  theme(legend.position = "none") +
  labs(title="Random forest", subtitle="randomForest(S ~ Area + Time)"); p.rf

p.rf.log <- ggplot(dat, aes(x= Area, y = S)) +
  geom_point(aes(x=Area, y = S.obs, shape= as.factor(Time), colour = S), size = 2) +
  scale_shape_manual(values=c(19, 1))+
  scale_colour_viridis()+
  geom_line(aes(x = Area, y = S.pred.rf, linetype = as.factor(Time)), size = 1) +
  labs(title="randomForest(S ~ Area + Time)", subtitle=expression(paste(Log[10], " axes")))  +
  scale_x_log10() + scale_y_log10() +
  theme(legend.position = "none"); p.rf.log


# ------------------------------------------------------------------------------
# BOOSTED TREE
require(gbm)

n.trees = 500
boost <- gbm(S ~ Area + Time, 
             data = dat, 
             n.trees = n.trees)

dat <- data.frame(dat, S.pred.gbm = predict(boost,
                                            newdata = dat,
                                            type = "response", n.trees = n.trees))

p.gbm <- ggplot(dat, aes(x= Area, y = S)) +
  geom_point(aes(x=Area, y = S.obs, shape= as.factor(Time), colour = S), size = 2) +
  scale_shape_manual(values=c(19, 1))+
  scale_colour_viridis()+
  geom_line(aes(x = Area, y = S.pred.gbm, linetype = as.factor(Time)), size = 1) +
  theme(legend.position = "none") +
  labs(title="Boosted trees", subtitle="gbm(S ~ Area + Time)"); p.gbm

p.gbm.log <- ggplot(dat, aes(x= Area, y = S)) +
  geom_point(aes(x=Area, y = S.obs, shape= as.factor(Time), colour = S), size = 2) +
  scale_shape_manual(values=c(19, 1))+
  scale_colour_viridis()+
  geom_line(aes(x = Area, y = S.pred.gbm, linetype = as.factor(Time)), size = 1) +
  labs(title="gbm(S ~ Area + Time)", subtitle=expression(paste(Log[10], " axes")))  +
  scale_x_log10() + scale_y_log10() +
  theme(legend.position = "none"); p.gbm.log


grid.arrange(p1, p.tr, p.rf, p.gbm, 
             p1.log, p.tr.log, p.rf.log, p.gbm.log,
             ncol = 4, nrow = 2)


# ------------------------------------------------------------------------------
pdf("../Graphics/ml_panels.pdf", width = 9, height = 4.8)
grid.arrange(p1, p.tr, p.rf, p.gbm, 
             p1.log, p.tr.log, p.rf.log, p.gbm.log,
             ncol = 4, nrow = 2)
dev.off()

pdf("../Graphics/ml_panels2.pdf", width = 7, height = 7)
grid.arrange(p.glm, p.tr, p.rf, p.gbm, ncol=2, nrow=2)
dev.off()

# ------------------------------------------------------------------------------
# CONFIDENCE INTERVALS, STD ERRORS, UNCERTAINTY
#
# https://github.com/swager/randomForestCI

library(devtools); install_github("swager/randomForestCI")
library(randomForestCI)
library(ranger)

rf <- randomForest(S ~ Area + Time, 
                   data = dat, 
                   mtry = 2, # I only have 2 variables, so I am getting smoother fit with this
                   ntree = 500,
                   keep.inbag = TRUE) 
ij = randomForestInfJack(rf, newdata = dat, calibrate = TRUE)

importance(rf)
varImpPlot(rf, type = 2)
plot(rf)

rf <- ranger(S ~ Area + Time, 
             data = dat, 
             mtry = 2, # I only have 2 variables, so I am getting smoother fit with this
             num.trees = 500,
             keep.inbag = TRUE,
             quantreg = TRUE) 
yhat <- predict(rf, data = dat, se.method = "infjack", 
                type = "se")
yhat <- predict(rf, 
                data = dat,
                type = "quantiles", 
                quantiles = c(0.025, 0.5, 0.975))

# -------------------------------------------------------------------------------------
# FIGURE FOR THE BOX

set.seed(891234)
n = 6
#x <- runif(n, 0, 5)
x <- 1:n
y <- 2*x + 3 + rnorm(n)
plot(x, y, pch = 19)
abline(a = 3, b = 2)
n

pdf("../Graphics/GLM_vs_trees.pdf", width = 4, height = 4)
  plot(x, y, pch = 19, ylim = c(4,16))
  abline(a = 3, b = 2, lwd = 2, col="blue")
  abline(v = 2.5)
  abline(v = 4.5)
dev.off()















