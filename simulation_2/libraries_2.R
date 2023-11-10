## fdaPDE ----------------------------------------------------------------------
library(fdaPDE)       # Library for Spatio-Temporal Density Estimation with PDE Regularization (STDE-PDE)
# version needed: 1.1-16 or later
# most updated library version always available at https://github.com/fdaPDE/fdaPDE
# periodically released on CRAN: https://cran.r-project.org/web/packages/fdaPDE/index.html

## OTHER LIBRARIES -------------------------------------------------------------
library(mvtnorm)      # Library for multivariate Gaussian distributions
library(stam)         # Library for Spatio-Temporal Kernel Density Estimation (STKDE)
library(gss)          # Library for Spatio-Temporal SPLINE Density Estimation (SPLINE)
library(spatstat)     # Library for ppp (planar point pattern) objects
library(ks)           # Library for Kernel Density Estimation (KDE)
library(TDA)          # Library for k-Nearest Neighbors Density Estimation (KNNDE)
library(pracma)       # Library for two-dimensional data interpolation (for visualization purposes)
library(philentropy)  # Library for Kullback-Leibler distance
library(RColorBrewer) # Library for color palettes
library(rgeos)        # Library for Interface to Geometry Engine
library(deldir)       # Library for Delaunay Triangulation and Voronoi Tessellation
library(fields)       # Library for Visualization of Spatial Data
library(Directional)  # Library for Kent distributions
library(retistruct)   # Library for Kernel Density Estimation (KDE) over manifolds
library(plotrix)      # Library for gradient.rect()
library(latex2exp)    # Library for LateX Expressions

rm(list=ls())
graphics.off()

source("helper_functions_2.R")