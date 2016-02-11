setwd("~/Box Sync/My Education/OxWaSP/Modules/Module 6/ABC/Report")
library(abctools)
library(gk)
library(ggplot2)
library(reshape2)
source("palette.R")

# g-and-k distribution, Large dataset
load("gkdataL.rda")
# theta <- gkdataL[,1:4]

# function to generate N-tiles summary order statistics from n simulations
ord_stats_obs <- function(N, n) {
  os <- (n/N) * (1:(N-1))
  obs <- rgk.orderstats(n, os, theta = c(3, 1, 2, 0.5))
  obs
}

# generate the observations (corresponding to gkdataL)
set.seed(1)
obs <- ord_stats_obs(100, 10000)

# use cross-validation to find good tolerance level for different ABC methods
# cv for
cv.res.rej <- cv4abc(param = gkdataL[,1:4], sumstat = gkdataL[,-(1:4)], nval=10, tols=c(0.001, 0.005, 0.01, 0.05, 0.1), method="rejection")
summary(cv.res.rej) # 0.001 best
# since lowest tested is best, try lower
cv.res.rej <- cv4abc(param = gkdataL[,1:4], sumstat = gkdataL[,-(1:4)], nval=10, tols=c(0.0001, 0.0005, 0.001), method="rejection")
summary(cv.res.rej) # 0.001 best for g,k, 0.0001 best for A,B
plot(cv.res.rej)

# cv for regression
cv.res.reg <- cv4abc(param = gkdataL[,1:4], sumstat = gkdataL[,-(1:4)], nval=10, tols=c(0.001, 0.005, 0.01, 0.05, 0.1), method="loclinear")
summary(cv.res.reg) # 0.005 best
plot(cv.res.reg)

# function for plotting output of ABC run
gk_param_plot <- function(post_samples) {
  df <- as.data.frame(post_samples)
  melted <- melt(df)
  colnames(melted) <- c("parameter", "value")
  
  p <- qplot(value, data = melted, geom = "histogram") 
  p <- p + facet_wrap( ~ parameter, scales = "free")
  vline.data <- data.frame(z=c(3, 1, 2, 0.5), parameter = c("A", "B", "g", "k"))
  p + geom_vline(aes(xintercept = z), vline.data, colour="red") + theme_bw()
}


# ABC with tol levels from cross-validation
# rejection
abc.gk.rej.L <- abc(target = obs,
                  param = gkdataL[,1:4],
                  sumstat = gkdataL[,-(1:4)],
                  tol = 0.001,
                  method = "rejection")
pdf("GK_REJ_L_HIST.pdf", height=4)
gk_param_plot(abc.gk.rej.L$unadj.values)
dev.off()

# regression
abc.gk.reg.L <- abc(target = obs,
                    param = gkdataL[,1:4],
                    sumstat = gkdataL[,-(1:4)],
                    tol = 0.005,
                    method = "loclinear")
pdf("GK_REG_L_HIST.pdf", height=4)
gk_param_plot(abc.gk.reg.L$adj.values)
dev.off()

# repeat for small data set
load("gkdata.rda")
set.seed(1)
obs_S <- ord_stats_obs(8, 1000)

# use cross-validation to find good tolerance level for different ABC methods
# cv for rejection
cv.res.rej.S <- cv4abc(param = gkdata[,1:4], sumstat = gkdata[,-(1:4)], nval=10, tols=c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1), method="rejection")
summary(cv.res.rej.S) # 0.005 best for B, 0.0001 best for A, g, k
plot(cv.res.rej.S)

# cv for regression
cv.res.reg.S <- cv4abc(param = gkdata[,1:4], sumstat = gkdata[,-(1:4)], nval=10, tols=c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1), method="loclinear")
summary(cv.res.reg.S) # 0.01 or 0.005 best
plot(cv.res.reg.S)

# ABC with tol levels from cross-validation
# rejection
abc.gk.rej.S <- abc(target = obs_S,
                    param = gkdata[,1:4],
                    sumstat = gkdata[,-(1:4)],
                    tol = 0.005,
                    method = "rejection")
pdf("GK_REJ_S_HIST.pdf", height=4)
gk_param_plot(abc.gk.rej.S$unadj.values)
dev.off()

# regression
abc.gk.reg.S <- abc(target = obs_S,
                    param = gkdata[,1:4],
                    sumstat = gkdata[,-(1:4)],
                    tol = 0.01,
                    method = "loclinear")
pdf("GK_REG_S_HIST.pdf", height=4)
gk_param_plot(abc.gk.reg.S$adj.values)
dev.off()

## Doing much better with small dataset

# Now try semi-auto ABC
mytf <- list(function(x){cbind(x, x ^ 2, x ^ 3, x ^ 4)})

# large dataset
saabc.gk <- semiauto.abc(obs = obs,
                         param = gkdataL[,1:4],
                         sumstats = gkdataL[,-(1:4)],
                         satr = mytf,
                         tol = 0.01, overlap = T, saprop = 1,
                         abcprop = 1, method = "loclinear",
                         final.dens = T
)
pdf("GK_REG_L_HIST_SA.pdf", height=4)
gk_param_plot(saabc.gk$post.sample)
dev.off()

# Generate 2d density of the accepted parameters
dens_gk_ab <- kde2d(saabc.gk$post.sample[, 1],
                    saabc.gk$post.sample[, 2])

# Contour plot of A vs B
pdf("GK_REG_L_CONT_AB.pdf", height=4)
filled.contour(dens_gk_ab, xlab = "A", ylab = "B",plot.axes = {points(3, 1); axis(1); axis(2)},
               color.palette = function(x){palette("#FF0000", "#FFFF00", x)})
dev.off()

dens_gk_gk <- kde2d(saabc.gk$post.sample[, 3],
                    saabc.gk$post.sample[, 4])
pdf("GK_REG_L_CONT_GK.pdf", height=4)
filled.contour(dens_gk_gk, xlab = "g", ylab = "k", plot.axes = { points(2, 0.5); axis(1); axis(2) }, 
               color.palette = function(x){palette("#FF0000", "#FFFF00", x)})
dev.off()

##############################################################

# cross-validation
cv.res.rej <- cv4abc(param = gkdataL[,1:4], sumstat = gkdataL[,-(1:4)], nval=10, tols=c(.005,.01, 0.05, 0.1), method="rejection")
summary(cv.res.rej)
plot(cv.res.rej)

osS <- (1000/8) * (1:7)

obsS <- rgk.orderstats(1000, osS, theta = c(3, 1, 2, 0.5))

abc.gk.rej <- abc(target = obsS,
              param = gkdata[,1:4],
              sumstat = gkdata[,-(1:4)],
              tol = 0.005,
              method = "rejection")

summary(abc.gk.rej)
par(mfrow=c(2,2))
hist(abc.gk.rej)

cv.res.reg <- cv4abc(param = gkdataL[,1:4], sumstat = gkdataL[,-(1:4)], nval=10, tols=c(0.001, .005,.01, 0.05, 0.1), method="loclinear")
summary(cv.res.reg)
plot(cv.res.reg)

abc.gk.reg <- abc(target = obsS,
                  param = gkdata[,1:4],
                  sumstat = gkdata[,-(1:4)],
                  tol = 0.005,
                  method = "loclinear")

summary(abc.gk.reg)
par(mfrow=c(2,2))
hist(abc.gk.reg)


set.seed(1)
entchoice.gk <- mincrit(obs = obs_stats,
                         param = theta,
                         sumstats = t(sumstats),
                         crit = nn.ent, limit = 5, do.only = 1:30,
                         tol = 1e-3, overlap = T, saprop = 1,
                         abcprop = 1, method = "rejection",
                         final.dens = T
)



# g-k dataset example
load("gkdata.rda")
dim(gkdata)
head(gkdata)
params <- gkdata[, 1:4]
octiles <- gkdata[, 5:11]
dim(octiles)

os <- (1000/8)*(1:7)

obs_stats <- rgk.orderstats(1000, os, theta = c(3, 1, 2, 0.5))

tfs <- list(function(x){cbind(x, x^2, x^3, x^4)})

set.seed(1)
saabc.gk <- semiauto.abc(obs = obs_stats, param = params, sumstats = octiles,
                      satr = tfs, overlap = TRUE, saprop = 1,
                      abcprop = 1, tol = 0.005, method = "loclinear",
                      final.dens = TRUE, plot=TRUE)

# Generate 2d density of the accepted parameters
dens_gk_ab <- kde2d(saabc.gk$post.sample[, 1],
                    saabc.gk$post.sample[, 2])

# Contour plot of A vs B
filled.contour(dens_gk_ab, xlab = "A", ylab = "B",plot.axes = {points(3, 1); axis(1); axis(2)},
               color.palette = function(x){palette("#FF0000", "#FFFF00", x)})


dens_gk_gk <- kde2d(saabc.gk$post.sample[, 3],
                    saabc.gk$post.sample[, 4])
filled.contour(dens_gk_gk, xlab = "g", ylab = "k", plot.axes = { points(2, 0.5); axis(1); axis(2) }, 
               color.palette = function(x){palette("#FF0000", "#FFFF00", x)})

dim(saabc.gk$post.sample)
head(saabc.gk$post.sample)


# 50 data sets
obs_50 <- t(replicate(50, rgk.orderstats(1000, os, theta = c(3, 1, 2, 0.5))))
head(obs_50)


dim(saabc.gk$post.sample)

# MSE
A <- matrix(1:20, ncol = 4)
t(t(A) - c(3,1,2,0.5))
A
post <- saabc.gk$post.sample
colMeans(post)

colMeans(t(t(post) - c(3,1,2,0.5))^2)

post.abc <- abc.gk.reg$adj.values

colMeans(t(t(post.abc) - c(3,1,2,0.5))^2)
post.abc.df <- as.data.frame(post.abc)
dim(post.abc.df)

melted <- melt(post.abc.df)
melted$parameter <- mf_labeller('variable',melted$variable)
colnames(melted) <- c("parameter", "value")
head(melted)


p <- qplot(value, data = melted, geom = "histogram") 
p <- p + facet_wrap( ~ parameter, scales = "free")
p
vline.data <- data.frame(z=c(3, 1, 2, 0.5), parameter = c("A", "B", "g", "k"))
p + geom_vline(aes(xintercept = z), vline.data, colour="red")

head(melted)

ggplot(melt(post.abc.df)) + geom_histogram()

gk_param_plot <- function(post_samples) {
  df <- as.data.frame(post_samples)
  melted <- melt(df)
  colnames(melted) <- c("parameter", "value")
  
  p <- qplot(value, data = melted, geom = "histogram") 
  p <- p + facet_wrap( ~ parameter, scales = "free")
  vline.data <- data.frame(z=c(3, 1, 2, 0.5), parameter = c("A", "B", "g", "k"))
  p + geom_vline(aes(xintercept = z), vline.data, colour="red")
}

gk_param_plot(post.abc)
