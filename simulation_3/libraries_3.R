## fdaPDE ----------------------------------------------------------------------
library(fdaPDE)       # Library for Spatio-Temporal Density Estimation with PDE Regularization (STDE-PDE)
# version needed: 1.1-16 or later
# most updated library versin always available at https://github.com/fdaPDE/fdaPDE
# periodically released on CRAN: https://cran.r-project.org/web/packages/fdaPDE/index.html

## OTHER LIBRARIES -------------------------------------------------------------
library(Directional)  # Library for Kent distributions
library(rgl)          # Library for titles in 3D plots
library(retistruct)   # Library for spatial kde on sphere
library(TDA)          # Library for k-Nearest Neighbors Desnity Estimation (KNNDE)
library(philentropy)  # Library for Leibler-Kullback distance
library(RColorBrewer) # Library for color palettes
library(latex2exp)    # Library for LateX Espressions
library(plotrix)      # Library for gradient.rect()

rm(list=ls())
graphics.off()

load("mesh/simulation3.fullPoints.proj.RData")
source("helper_functions_3.R")
