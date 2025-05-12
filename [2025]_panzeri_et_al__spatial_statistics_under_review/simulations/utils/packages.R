if(!require(pacman, quietly = TRUE)) install.packages("pacman")
pacman::p_load("devtools")

# removed from CRAN
if(!require(maptools, quietly = TRUE)){
  devtools::install_version("maptools", version = "1.1-8",
                            repos = "https://cran.stat.unipd.it/")
}

# removed from CRAN
if(!require(maptools, quietly = TRUE)){
  devtools::install_version("shp2graph", version = "0-5",
                            repos = "https://cran.stat.unipd.it/")
}

pacman::p_load("fdaPDE", "spatstat", "shp2graph", "sf", "sfnetworks",
               "tidygraph", "stlnpp", "gridExtra", "maptools", "latex2exp",
               "RColorBrewer")

# install.packages("fdaPDE", lib = "C:/Users/simon/Documents")
# library("fdaPDE", lib.loc = "C:/Users/simon/Documents/R-lib-local/fdaPDE_intensity",
#         quietly = TRUE)
