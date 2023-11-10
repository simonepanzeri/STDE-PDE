## fdaPDE ----------------------------------------------------------------------
library(fdaPDE)       # Library for Spatio-Temporal Density Estimation with PDE Regularization (STDE-PDE)
# version needed: 1.1-16 or later
# most updated library version always available at https://github.com/fdaPDE/fdaPDE
# periodically released on CRAN: https://cran.r-project.org/web/packages/fdaPDE/index.html

## OTHER LIBRARIES -------------------------------------------------------------
library(mvtnorm)      # Library for multivariate Gaussian distributions
library(stam)         # Library for Spatio-Temporal Kernel Density Estimation (STKDE)
library(gss)          # Library for Spatio-Temporal SPLINE Density Estimation (SPLINE)
library(lgcp)         # Library for Log-Gaussian Cox Processes (LGCP)
library(spatstat)     # Library for ppp (planar point pattern) objects
library(ks)           # Library for Kernel Density Estimation (KDE)
library(TDA)          # Library for k-Nearest Neighbors Density Estimation (KNNDE)
library(pracma)       # Library for two-dimensional data interpolation (for visualization purposes)
library(philentropy)  # Library for Kullback-Leibler distance
library(RColorBrewer) # Library for color palettes
library(sparr)        # Library for kernel estimation of spatio-temporal relative risk
library(rgeos)        # Library for Interface to Geometry Engine
library(deldir)       # Library for Delaunay Triangulation and Voronoi Tessellation
library(fields)       # Library for Visualization of Spatial Data
library(INLA)         # Library for Integrated Nested Laplace Approximation
library(inlabru)      # Library for Integrated Nested Laplace Approximation
library(plotly)       # Library for Plots
library(dplyr)        # Library for "%>%"
library(processx)     # Library for Exporting Plots as Static Images
library(latex2exp)    # Library for LateX Expressions
library(plotrix)      # Library for Various Plotting Functions
library(htmlwidgets)  # Library for Exporting in .html
library(webshot)      # Library for Exporting in .png (screenshot from .html)
library(cubature)     # Library for Adaptive Multivariate Integration over Hypercubes

rm(list=ls())
graphics.off()

options(warn = -1)

source("helper_functions_1.R")
