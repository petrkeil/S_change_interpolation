library(doParallel)
library(raster)
library(spatstat)
library(mvtnorm)
library(truncnorm)
library(mobsim)
library(plyr)
library(maptools)
library(rgeos)
library(rgdal)
library(randomForest)
library(gbm)
library(dplyr)
library(tidyverse)
library(viridis)
library(sampling)
library(latex2exp)
library(gridExtra)
library(hexbin)

# -----------------------------------------------------------------------------

raw.res <- read.csv("data_simulated_10Jun2021.csv") # this needs to be unzipped first

raw.res <- na.omit(raw.res)
raw.res$Area <- round(raw.res$Area, 4)
raw.res$Grain_fact <- as.factor(paste("Grain = ", raw.res$Area))
raw.res$Grain <- raw.res$Area
raw.res$Formula <- as.character(raw.res$Formula)
raw.res$prop2 <- as.numeric(raw.res$prop2)

raw.res$Error <- abs(raw.res$Delta - raw.res$est.Delta)
raw.res$Bias <- (raw.res$Delta - raw.res$est.Delta)

raw.res$temp.bias <-  raw.res$prop1 - raw.res$prop2
raw.res$min.effort <- pmin(raw.res$prop1, raw.res$prop2)


raw.res$max.effort <- pmax(raw.res$prop1, raw.res$prop2)
raw.res$mean.effort <- (raw.res$prop1 + raw.res$prop2)/2

# remove the super extreme outliers 
raw.res <- raw.res[raw.res$est.Delta > -25,]
raw.res <- raw.res[raw.res$est.Delta < 20,]
hist(raw.res$est.Delta)

# change formula names
Formula2 <- raw.res$Formula
Formula2[Formula2 == "Area + Time "] <- "(a) SAR, A+t"
Formula2[Formula2 == "Area + Time + x + y"] <- "(b) SAR, A+t+X+Y"
Formula2[Formula2 == "Area + x + y, Time separate"] <- "(c) SAR, A+X+Y"
Formula2[Formula2 == "Area, Time separate"] <- "(d) SAR, A"
Formula2[Formula2 == "OAR-based"] <- "(e) OAR, A+t+X+Y"
raw.res <- data.frame(raw.res, Formula2)

# only use the non-overlapping polygons
#raw.res <- raw.res[raw.res$overlap == FALSE, ]

cors <- raw.res %>% 
  group_by(Method, Formula2) %>% 
  dplyr::summarise(r = round(cor(Delta, est.Delta), 2)) %>%
  data.frame()
cors$r <- paste("r =", cors$r)


fig.obspred <- ggplot(data = raw.res, aes(x = Delta, y = est.Delta)) + 
  # geom_point(pch = 1, alpha = 0.3) + 
  # geom_density_2d_filled(aes(fill = stat((count)))) +
  geom_hex(aes(fill = stat(log10(count)))) +
   scale_fill_continuous(name  = TeX(r'($\log_{10} \; \count$)')) +
  #stat_density_2d(geom = "raster", aes(fill = after_stat(density))) + 
  # scale_fill_gradient(trans = "log10") +
  #(trans = "log10", name = "Count") +
  geom_vline(xintercept = 0, linetype = 2, colour = "darkgrey") + 
  geom_hline(yintercept = 0, linetype = 2, 
             colour = "darkgrey") +
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  coord_fixed(ratio = 1) +
  # xlab(expression(True ~ mean ~ Delta ~ S)) +
  # ylab(expression(Predicted ~ mean ~ Delta ~ S)) +
  xlab(TeX(r'($\True \; \hat{S}_{\Delta}$)')) +
  ylab(TeX(r'($\Predicted \; \hat{S}_{\Delta}$)')) +
  facet_grid(Method ~ Formula2) + theme_bw() +
  geom_text(data = cors, mapping = aes(x = -10, y = 15, label = r)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
fig.obspred


pdf("../Graphics/obs_vs_pred.pdf", width = 10, height = 6)
fig.obspred
dev.off()


# ------------------------------------------------------------------------------
sub.res <- raw.res[raw.res$Method == "Boosted trees", ]
#sub.res <- raw.res[raw.res$Formula == "OAR-based",]
sub.res <- sub.res[sub.res$Formula == "Area + Time + x + y",]

# stratified sample by Error
dD <- density(sub.res$Error)
prob <- 1/approx(dD$x, dD$y, xout = sub.res$Error)$y
strat.sampl <- sort(sample(x = 1:nrow(sub.res), replace = FALSE, prob = prob, size = 4000))
sub.ERROR <- sub.res[strat.sampl,]

# stratified sample by Bias
dD <- density(sub.res$Bias)
prob <- 1/approx(dD$x, dD$y, xout = sub.res$Bias)$y
strat.sampl <- sort(sample(x = 1:nrow(sub.res), replace = FALSE, prob = prob, size = 4000))
sub.BIAS <- sub.res[strat.sampl,]


# ------------------------------------------------------------------------------

fig.grain_error <- ggplot(data = sub.ERROR, aes(x = as.factor(Grain), y = Error)) + 
  geom_point(pch=1)  + 
  geom_violin(draw_quantiles=c(0.5), fill = "grey") +
  scale_y_log10() + theme_bw() +
  xlab("Grain") + 
  ylab("Prediction error")  +
  labs(title = "(a)")
fig.grain_error

fig.grain_bias <- ggplot(data = sub.BIAS, aes(x = as.factor(Area), y = Bias)) + 
  geom_point(pch=1)  + 
  geom_hline(yintercept= 0) +
  geom_violin(draw_quantiles=c(0.5), fill = "grey") +
  xlab("Grain") +
  ylab("Prediction bias") +
  labs(title = "(b)") +
  theme_bw()
fig.grain_bias

pdf("../Graphics/grain_perf.pdf", width = 8, height = 4)
grid.arrange(fig.grain_error, fig.grain_bias, ncol= 2)
dev.off()

# ------------------------------------------------------------------------------

fig.tbias_error <- ggplot(data = sub.ERROR, aes(x = as.factor(temp.bias), y = Error)) +
  geom_hline(yintercept= 0) +
  geom_violin(draw_quantiles=c(0.5), fill = "grey") +
  scale_y_log10() + 
  #facet_grid(.~Grain_fact) + 
  xlab("Temporal sampling bias = P1 - P2") +
  ylab("Prediction error")+
  theme_bw() + labs(title="(a)")

fig.tbias_bias <- ggplot(data = sub.BIAS, aes(x = as.factor(temp.bias), y = Bias)) +
  geom_hline(yintercept= 0) +
  geom_violin(draw_quantiles=c(0.5), fill = "grey") +
  # scale_y_log10() + 
  # facet_grid(.~Grain_fact) + 
  xlab("Temporal sampling bias = P1 - P2") +
  ylab("Prediction bias") + labs(title="(b)") +
  theme_bw()

pdf("../Graphics/tbias_perf.pdf", width = 8, height = 4)
grid.arrange(fig.tbias_error, fig.tbias_bias, ncol=2)
dev.off()


# ------------------------------------------------------------------------------

fig.mineff_error <- ggplot(data = sub.ERROR, aes(x = as.factor(mean.effort), y = Error)) + 
  geom_violin(draw_quantiles=c(0.5), fill = "grey") + 
  scale_y_log10()  + 
  #facet_grid(.~Grain_fact ) +
  ylab("Prediction error") + 
  xlab("Mean sampling effort") +
  labs(title="(a)") +
  theme_bw()

fig.mineff_bias <- ggplot(data = sub.BIAS, aes(x = as.factor(mean.effort), y = Bias)) +
  geom_hline(yintercept= 0) +
  geom_violin(draw_quantiles=c(0.5), fill = "grey") + 
  # scale_y_log10()  + 
  #facet_grid(.~Grain_fact ) +
  ylab("Prediction bias") + 
  xlab("Mean sampling effort") +
  labs(title="(b)") +
  theme_bw()

pdf("../Graphics/mineff_perf.pdf", width = 8, height = 4)
grid.arrange(fig.mineff_error, fig.mineff_bias, nrow=1, ncol=2)
dev.off()


# OVERLAP YES OR NO

ggplot(data = sub.ERROR, aes(x = overlap, y = Error)) + 
  geom_violin(draw_quantiles=c(0.5), fill = "grey") + 
  scale_y_log10()  +
  ylab("Prediction error") + 
  xlab("Spatial overlap") +
  labs(title="(a)") +
  theme_bw()

ggplot(data = sub.BIAS, aes(x = overlap, y = Bias)) + 
  geom_violin(draw_quantiles=c(0.5), fill = "grey") + 
  scale_y_log10()  +
  ylab("Prediction bias") + 
  xlab("Spatial overlap") +
  labs(title="(b)") +
  theme_bw()

# Delta



