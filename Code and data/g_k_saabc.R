# Recreating examples from Fearnhead and Prangle (2011)
  
# Inference for g and k distribution
library(abctools)
library(gk)

source("palette.R")

n <- 10e4
set.seed(1)
A <- runif(n, 0, 10)
B <- runif(n, 0, 10)
g <- runif(n, 0, 10)
k <- runif(n, 0, 10)
  
theta <- cbind(A, B, g, k)
  
# Simulate data with 1000 iid draws of each
#X <- sapply(1:n, function(i) { rgk(1000, theta = theta[i,])})
  
os <- (1000/100) * (1:99)
  

  
# Generate the 'observed' summary statistics
obs_stats <- rgk.orderstats(10e4, os, theta = c(3, 1, 2, 0.5))

# Generate the simulated summary statistics
sumstats <- sapply(1:n, function(i) { rgk.orderstats(10e4, orderstats = os,
                                                     theta = theta[i,])})


# The functions for the data transformations  
mytf <- list(function(x){cbind(x, x ^ 2, x ^ 3, x ^ 4)})

# Run the semi-auto ABC
saabc.gk <- semiauto.abc(obs = obs_stats,
                        param = theta,
                        sumstats = t(sumstats),
                        satr = mytf,
                        tol = 1e-3, overlap = T, saprop = 1,
                        abcprop = 1, method = "rejection",
                        final.dens = T
                        )

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

