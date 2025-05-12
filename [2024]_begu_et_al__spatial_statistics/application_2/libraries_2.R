## fdaPDE ----------------------------------------------------------------------
library(fdaPDE)       # Library for Spatio-Temporal Density Estimation with PDE Regularization (STDE-PDE)
# version needed: 1.1-16 or later
# most updated library version always available at https://github.com/fdaPDE/fdaPDE
# periodically released on CRAN: https://cran.r-project.org/web/packages/fdaPDE/index.html

## OTHER LIBRARIES -------------------------------------------------------------
library(Directional)  # Library for Kent distribution
library(rgl)          # Library for titles in 3D plots
library(retistruct)   # Library for spatial kde on sphere
library(TDA)          # Library for knnDE
library(philentropy)  # Library for Leibler-Kullback distance
library(lubridate)    # Library for ymd_hms (date-time conversions)
library(stam)         # Library for stkde (not CRAN anymore)
library(maptools)     # Library for coastlines visualization
library(dichromat)    # Library for colorRampPalette
library(stars)        # Library for st_transform
library(tidyverse)    # Library for ggplot
library(plotrix)      # Library for gradient.rect
library(latex2exp)    # Library for LateX Expressions

rm(list=ls())
graphics.off()

source("helper_functions_2.R")