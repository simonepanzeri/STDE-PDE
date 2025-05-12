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
library(sp)           # Library for Spatial Coordinates and Spatial Objects
library(leaflet)      # Library for Highly Customizable Mapping
library(RgoogleMaps)  # Library for Maps
library(marmap)       # Library for Adding City Locations and Names to Maps
library(maps)         # Library for Maps
library(pracma)       # Library for Two-Dimensional Data Interpolation
library(area)         # Library for Size of Cells

# library(PBSmapping)
# library(stpp)
# library(ggmap)         # for mapping points on maps
# library(gplots)        # for col2hex() function
# library(RColorBrewer)  # for color palettes
# library(sf)            # for working with spatial data
# library(openintro)     # for the abbr2state() function
# library(tidyverse)     # for data cleaning and plotting
# library(ggthemes)      # for more themes (including theme_map())
# library(magick)        # for cropping images
# library(imager)        # for loading images
# library(splancs)       # Library for Kernel Smoothing and Spatial/Spatio-Temporal Point Patterns
# library(loa)
# library(latticeExtra)
theme_set(theme_minimal())

rm(list=ls())
graphics.off()

source("helper_functions_1.R")

options(warn = -1)
