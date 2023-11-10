#'############################################################################'#
#'##################### SPATIO-TEMPORAL DENSITY ESTIMATION ###################'# 
#'######### APPLICATION 1 - GASTROINTESTINAL INFECTIONS IN SOUTHAMPTON #######'#
#'############################################################################'#

## LIBRARIES AND FUNCTIONS -----------------------------------------------------
source("libraries_1.R")

## IMPORT DATA -----------------------------------------------------------------
# Load the Data
load("data/data.RData")

# # Rescale the Spatial Domain in the Square [0,1]x[0,1]
# hx <- xyt$window$xrange[1]
# hy <- xyt$window$yrange[1]
# 
# ScaleFactor <- 1/diff(xyt$window$xrange)
# 
# xyt$x <- xyt$x - hx
# xyt$x <- xyt$x * ScaleFactor
# 
# xyt$y <- xyt$y - hy
# xyt$y <- xyt$y * ScaleFactor
# 
# xyt$window$xrange <- xyt$window$xrange - hx
# xyt$window$xrange <- xyt$window$xrange * ScaleFactor
# 
# xyt$window$yrange <- xyt$window$yrange - hy
# xyt$window$yrange <- xyt$window$yrange * ScaleFactor
# 
# xyt$window$bdry[[1]]$x <- xyt$window$bdry[[1]]$x - hx
# xyt$window$bdry[[1]]$x <- xyt$window$bdry[[1]]$x * ScaleFactor
# 
# xyt$window$bdry[[1]]$y <- xyt$window$bdry[[1]]$y - hy
# xyt$window$bdry[[1]]$y <- xyt$window$bdry[[1]]$y * ScaleFactor
# 
# xyt$window$bdry[[1]]$area <- xyt$window$bdry[[1]]$area * ScaleFactor^2
# 
# save(xyt, hx, hy, ScaleFactor, file = "data/data_rescaled.RData")
# 
# rm(list = ls())

## SPATIO-TEMPORAL DISCRETIZATION ----------------------------------------------
# Spatial 2D Mesh over the Southampton Region for Estimation
boundary_nodes <- cbind(as.data.frame(xyt$window$bdry[[1]]$x),as.data.frame(xyt$window$bdry[[1]]$y))
boundary_segments <- cbind(1:dim(boundary_nodes)[1], c(2:dim(boundary_nodes)[1], 1))
mesh <- create.mesh.2D(nodes = boundary_nodes, segments = boundary_segments)
mesh <- refine.mesh.2D(mesh, maximum_area = 4.25, minimum_angle = 25)
FEMbasis <- create.FEM.basis(mesh)

x11()
plot(mesh)
dev.off()

# Boundary Nodes
boundary_nodes <- as.matrix(boundary_nodes)

# Fine Spatial 2D Mesh over the Southampton Region for Evaluation
# mesh.eval <- refine.mesh.2D(mesh, maximum_area = 0.5, minimum_angle = 25)
# FEMbasis.eval <- create.FEM.basis(mesh.eval)
# 
# x11()
# plot(mesh.eval)
# dev.off()

# Grid Size for Estimation
n <- 500

# Fine Grid for Evaluation
X <- seq(min(boundary_nodes[,1]), max(boundary_nodes[,1]), length.out = n)
Y <- seq(min(boundary_nodes[,2]), max(boundary_nodes[,2]), length.out = n)
grid <- expand.grid(X, Y)
mesh.eval <- create.mesh.2D(grid)
FEMbasis.eval <- create.FEM.basis(mesh.eval)

# Clean Grid on the Southampton Region
grid_clean <- NULL
idx_remove <- NULL
for(p in 1:nrow(mesh.eval$nodes)){
  if(point.in.polygon(mesh.eval$nodes[p,1], mesh.eval$nodes[p,2], boundary_nodes[,1], boundary_nodes[,2])>0){
    grid_clean <- rbind(grid_clean, mesh.eval$nodes[p,])
  }
  else{
    idx_remove <- c(idx_remove, p)
  }
}

x11()
plot(rbind(boundary_nodes, boundary_nodes[1,]), type = "l",
     xlab = "Eastings", ylab = "Northings")
points(grid_clean)
dev.off()

# Average Triangle Area
# tr_area <- 0
# for(tr in 1:nrow(mesh.eval$triangles)){
#   tr_area <- tr_area + triangle_area(rbind(mesh.eval$nodes[mesh.eval$triangles[tr,1],], mesh.eval$nodes[mesh.eval$triangles[tr,2],], mesh.eval$nodes[mesh.eval$triangles[tr,3],]))
# 
# }
# avg_tr_area <- tr_area / nrow(mesh.eval$triangles)

# Temporal 1D Mesh over [0,10]
mesh_time <- c(0, seq(from = 1.5, to = 9.5, by = 2), 10)

# Locations
locations <- cbind(xyt$x, xyt$y)
  
# Times
times <- xyt$t
discrete_times <- times %/% 1

# Sample Size
N <- nrow(locations)

## ESTIMATION PROCEDURE --------------------------------------------------------

### STDE-PDE: Spatio-Temporal Density Estimation with PDE Regulariz. -----------
# Smoothing Parameters
lambda <- 0.1
# alternative: lambda <- 10^seq(from = -2, to = 0, by = 1)

lambda_time <- 0.001
# alternative: lambda_time <- 10^seq(from = -4, to = -2, by = 1)

# Solution
# [If lambda and/or lambda_time are vectors, to select the best proposals by
# 10-folds CV, please specify: preprocess_method = "RightCV"]
solution_STDEPDE <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                mesh_time = mesh_time, lambda = lambda, tol1 = 1e-7, 
                                lambda_time = lambda_time, fvec = NULL, heatStep = 0.1,
                                heatIter = 10, print = T, nfolds = 10, nsimulations = 10000,
                                step_method = "Wolfe_Method", direction_method = "L-BFGS5",
                                preprocess_method = "NoCrossValidation")

FEMfunction_STDEPDE <- FEM.time(solution_STDEPDE$g, mesh_time, FEMbasis, FLAG_PARABOLIC = FALSE)

### LGCP: Log-Gaussian Cox Process ---------------------------------------------
# Spatio-Temporal Planar Point Process
win <- xyt$window
tlim <- c(0, 10)
xy_LGCP <- ppp(x = locations[,1], y = locations[,2], window = win)
xyt_LGCP <- stppp(xy_LGCP, times, tlim = tlim)
xyt_LGCP[["n"]] <- N
xyt_LGCP <- integerise(xyt_LGCP)

# Kernel Smoothed Intensity Function
den <- density.ppp(xyt_LGCP, sigma = 1)
sar <- spatialAtRisk(den)
mut <- muEst(xyt_LGCP)

exceed <- exceedProbs(c(1.5,2,3))

# Temporary Folder
dir.create(file.path(getwd(), "LGCP"), showWarnings = FALSE)
tmpdr <- paste0(getwd(), "/LGCP")

# Grid Size
n_LGCP <- 64

# Solution
solution_LGCP <- lgcpPredict(xyt = xyt_LGCP, T = as.integer(9), laglength = as.integer(9),
                             model.parameters = lgcppars(sigma = 1.6, phi = 1.9, theta = 1.4),
                             gridsize = n_LGCP, spatial.intensity = sar,
                             temporal.intensity = mut,
                             mcmc.control = mcmcpars(mala.length = 120000, burnin = 20000, retain = 100,
                                                   adaptivescheme = andrieuthomsh(inith = 1, alpha = 0.5, C = 1, targetacceptance = 0.574)),
                             output.control = setoutput(gridfunction = dump2dir(dirname = tmpdr, forceSave = TRUE),
                                                      gridmeans = MonteCarloAverage("exceed"))) # to be run locally
intensity_LGCP <- intens(solution_LGCP)
grid_LGCP <- expand.grid(intensity_LGCP$xvals, intensity_LGCP$yvals)
FEMbasis_LGCP <- create.FEM.basis(create.mesh.2D(grid_LGCP))

### INLA: SPDE Spatio-Temporal Model -------------------------------------------
# Temporal Discretization
k <- length(mesh_time)
mesh_time_INLA <- INLA::inla.mesh.1d(mesh_time)

# Spatial Discretization
min_bdy_x <- min(boundary_nodes[,1])
max_bdy_x <- max(boundary_nodes[,1])
min_bdy_y <- min(boundary_nodes[,2])
max_bdy_y <- max(boundary_nodes[,2])

domain_INLA <- SpatialPolygons(list(Polygons(list(Polygon(boundary_nodes)),"0")))

mesh_INLA <- INLA::inla.mesh.2d(boundary = INLA::inla.sp2segment(domain_INLA),
                                max.edge = c(6,6), cutoff = 0)
# alternative: pass locations as additional argument to build the mesh through R-INLA

# Plot of the Spatial Mesh
# x11()
# plot(mesh_INLA)
# points(boundary_nodes, col = "red", pch = 19)

# SPDE Model
spde <- INLA::inla.spde2.pcmatern(mesh = mesh_INLA, prior.range = c(5, 0.01),
                                  prior.sigma = c(1, 0.01))
m <- spde$n.spde

# Spatio-Temporal Projection Matrix: Kronecker Product between the Spatial
# Projector Matrix and the Group Projector one, i.e., the Temporal Dimension
Ast <- INLA::inla.spde.make.A(mesh = mesh_INLA, loc = as.matrix(locations),
                              n.group = length(mesh_time_INLA$n),
                              group = times, group.mesh = mesh_time_INLA)

# Space-Time Index Set
idx_INLA <- INLA::inla.spde.make.index("s", spde$n.spde, n.group = mesh_time_INLA$n)

# Dual Mesh
dmesh <- book.mesh.dual(mesh_INLA)

# Intersection with each Polygon from the Dual Mesh
w <- sapply(1:length(dmesh), function(i) {
  if (gIntersects(dmesh[i,], domain_INLA))
    return(gArea(gIntersection(dmesh[i,], domain_INLA)))
  else return(0)
})

# Spatio-Temporal Volume
st.vol <- rep(w, k) * rep(diag(INLA::inla.mesh.fem(mesh_time_INLA)$c0), m)

# Data Stack
y <- rep(0:1, c(mesh_time_INLA$n * m, N))
expected <- c(st.vol, rep(0, N))
stk <- INLA::inla.stack(
  data = list(y = y, expect = expected), 
  A = list(rbind(Diagonal(n = k * m), Ast), 1), 
  effects = list(idx_INLA, list(a0 = rep(1, k * m + N))))

# Model Fitting with the Gaussian Approximation
pcrho <- list(prior = "pccor1", param = c(0.7, 0.7))
form <- y ~ 0 + a0 + f(s, model = spde, group = s.group, 
                       control.group = list(model = "ar1",
                                            hyper = list(theta = pcrho)))

res <- INLA::inla(form, family = "poisson", 
                  data = INLA::inla.stack.data(stk), E = expect,
                  control.predictor = list(A = INLA::inla.stack.A(stk)),
                  control.inla = list(strategy = "adaptive"))

# NOTE: the exponential of the intercept plus the random effect at each
#       space-time integration point is the relative risk at each of these
#       points. This relative risk times the space-time volume will give the
#       expected number of points (E(n)) at each one of these space-time
#       locations. Summing over them will give a value that approaches the
#       number of observations

eta.at.integration.points <- res$summary.fix[1,1] + res$summary.ran$s$mean
#c(n = n, "E(n)" = sum(st.vol * exp(eta.at.integration.points)))

# Grid Size
n_INLA <- n

# Projection over a Grid for each Time Knot
r0 <- diff(range(boundary_nodes[, 1])) / diff(range(boundary_nodes[, 2]))
prj <- INLA::inla.mesh.projector(mesh_INLA, xlim = range(boundary_nodes[, 1]),
                                 ylim = range(boundary_nodes[, 2]), dims = c(n_INLA, n_INLA)) 
ov <- over(SpatialPoints(prj$lattice$loc), domain_INLA)
m.prj <- lapply(1:k, function(j) {
  r <- INLA::inla.mesh.project(prj, res$summary.ran$s$mean[1:m + (j - 1) * m])
  # alternative: r <- INLA::inla.mesh.project(prj, exp(res_INLA$summary.fix$mean + res_INLA$summary.ran$s$mean[1:m + (j - 1) * m]) * st.vol[1:m + (j - 1) * m] / N)
  r[is.na(ov)] <- NA
  return(r) 
})

# Plot of the Fitted Latent Field at each Time Knot
# igr <- apply(abs(outer(times, mesh_time_INLA$loc, "-")), 1, which.min)
# zlm <- range(unlist(m.prj), na.rm = TRUE)
# x11()
# par(mfrow = c(2, 3), mar = c(0, 0, 0.9, 0))
# for (j in 1:k) {
# image2D(x = prj$x, y = prj$y, z = m.prj[[j]], asp = 1,
#         xlab = "", zlim = zlm, axes = FALSE, col = heat.colors(100),
#         main = paste0("Time: ", j))
# #points(locations[igr == j, 1:2], pch = 19)
# }

# Solution
x_INLA <- prj$x
y_INLA <- prj$y
grid_INLA <- expand.grid(x_INLA, y_INLA)
solution_INLA <- m.prj
FEMbasis_INLA <- create.FEM.basis(create.mesh.2D(grid_INLA))

### INLA-aggregated: SPDE Spatio-Temporal Aggregated Model (Large Dataset) -----
# Temporal Discretization
k <- length(mesh_time)
mesh_time_INLA <- INLA::inla.mesh.1d(mesh_time)

# Spatial Discretization
domain_INLA <- SpatialPolygons(list(Polygons(list(Polygon(boundary_nodes)),"0")))
mesh_INLA <- INLA::inla.mesh.2d(boundary = INLA::inla.sp2segment(domain_INLA),
                                max.edge = c(6,6), cutoff = 0)
# alternative: pass locations as additional argument to build the mesh through R-INLA

# Spatio-Temporal Planar Point Pattern
xyt_INLA <- stppp(xy_LGCP, times, tlim = range(times))

# Space-Time Aggregation
dd <- deldir(mesh_INLA$loc[, 1], mesh_INLA$loc[, 2])
tiles <- tile.list(dd)
polys <- SpatialPolygons(lapply(1:length(tiles), function(i) {
  p <- cbind(tiles[[i]]$x, tiles[[i]]$y)
  n <- nrow(p)
  Polygons(list(Polygon(p[c(1:n, 1), ])), i)
}))
area <- factor(over(SpatialPoints(cbind(xyt_INLA$x, xyt_INLA$y)), polys),
               levels = 1:length(polys))

t.breaks <- sort(c(mesh_time_INLA$loc[c(1, k)],
                   mesh_time_INLA$loc[2:k - 1] / 2 + mesh_time_INLA$loc[2:k] / 2))
time <- factor(findInterval(xyt_INLA$t, t.breaks), levels = 1:(length(t.breaks) - 1))

# Distribution of Data Points on the Time Knots
# table(time)

agg.dat <- as.data.frame(table(area, time))
for(j in 1:2){
  agg.dat[[j]] <- as.integer(as.character(agg.dat[[j]])) 
}

# Areas of the Polygons
w.areas <- sapply(1:length(tiles), function(i) {
  p <- cbind(tiles[[i]]$x, tiles[[i]]$y)
  n <- nrow(p)
  pl <- SpatialPolygons(
    list(Polygons(list(Polygon(p[c(1:n, 1),])), i)))
  if (gIntersects(pl, domain_INLA))
    return(gArea(gIntersection(pl, domain_INLA)))
  else return(0)
})

# Time Width
w.t <- diag(INLA::inla.mesh.fem(mesh_time_INLA)$c0)

# Space-Time Volume (Area Unit per Time Unit) at each Polygon and Time Knot
e0 <- w.areas[agg.dat$area] * (w.t[agg.dat$time])

# Projector Matrix
A.st <- INLA::inla.spde.make.A(mesh_INLA, mesh_INLA$loc[agg.dat$area, ],
                               group = agg.dat$time, mesh.group = mesh_time_INLA)

# SPDE Model Object
spde <- INLA::inla.spde2.pcmatern(mesh_INLA, prior.sigma = c(1,0.01),
                                  prior.range = c(0.05,0.01))

# Space-Time Index Set
idx_INLA <- INLA::inla.spde.make.index("s", spde$n.spde, n.group = k)

# Data Stack
stk <- INLA::inla.stack(data = list(y = agg.dat$Freq, exposure = e0),
                        A = list(A.st, 1),
                        effects = list(idx_INLA, list(b0 = rep(1, nrow(agg.dat)))))

# PC Prior on Correlation
pcrho <- list(theta = list(prior = "pccor1", param = c(0.7, 0.7)))

# Model Formula: Intercept + Spatial Effect + Temporal Effect
formula <- y ~ 0 + b0 + f(s, model = spde, group = s.group,
                          control.group = list(model = "ar1", hyper = pcrho))

# Model Fitting
res <- INLA::inla(formula, family = "poisson",
                  data = INLA::inla.stack.data(stk), E = exposure, 
                  control.predictor = list(A = INLA::inla.stack.A(stk)),
                  control.inla = list(strategy ="adaptive"))

# Value of mu and Summary for the Intercept
# cbind(True = mu, res$summary.fixed[, 1:6])

# Summary for the Hyperparameters
# cbind(True = c(range, sigma, rho), res$summary.hyperpar[, c(1, 2, 3, 5)])

# Grid Size
n_INLA_aggregated <- n

# Projection over a Grid for each Time Knot
r0 <- diff(range(boundary_nodes[, 1])) / diff(range(boundary_nodes[, 2]))
prj <- INLA::inla.mesh.projector(mesh_INLA, xlim = c(min_bdy_x, max_bdy_x), 
                                 ylim = c(min_bdy_y, max_bdy_y),
                                 dims = c(n_INLA_aggregated, n_INLA_aggregated))
g.no.in <- is.na(over(SpatialPoints(prj$lattice$loc), domain_INLA))
t.mean <- lapply(1:k, function(j) {
  z.j <- res$summary.ran$s$mean[idx_INLA$s.group == j]
  # alternative: z.j <- exp(res_INLA_aggregated$summary.fix$mean + res_INLA_aggregated$summary.ran$s$mean[idx_INLA$s.group == j]) * w.areas / N
  z <- INLA::inla.mesh.project(prj, z.j)
  z[g.no.in] <- NA
  return(z)
})

# Plot of the Spatial Surface at each Time Knot
# zlims <- range(unlist(t.mean), na.rm = TRUE)
# x11()
# par(mfrow = c(2, 3), mar = c(0.1, 0.1, 1, 0.1))
# for (j in 1:k) {
#   image(prj$x, prj$y, t.mean[[j]], axes = FALSE, zlim = zlims,
#         col = heat.colors(100), main = paste0("Time knot: ", j))
#   #points(xyt$x[time == j], xyt$y[time == j], cex = 0.1, cex.main = 0.95)
# }

# Solution
x_INLA_aggregated <- prj$x
y_INLA_aggregated <- prj$y
grid_INLA_aggregated <- expand.grid(x_INLA_aggregated, y_INLA_aggregated)
solution_INLA_aggregated <- t.mean
FEMbasis_INLA_aggregated <- create.FEM.basis(create.mesh.2D(grid_INLA_aggregated))

### STKDE: Spatio-Temporal Kernel Density Estimation ---------------------------
# Spatio-Temporal Planar Point Pattern
xy_STKDE <- ppp(x = locations[,1], y = locations[,2], window = win)

# Grid Size
n_STKDE <- n

# Time Refinement
h <- 1/6
t <- seq(from = h, to = 10-h, by = 2*h)

# Solution
solution_STKDE <- spattemp.density(pp = xy_STKDE, h = NULL, tt = times,
                                   lambda = NULL, tlim = c(min(t)-h, max(t)+h),
                                   sedge = "none", tedge = "none", 
                                   sres = n_STKDE, tres = length(t))

x_STKDE <- solution_STKDE$spatial.z$xcol
y_STKDE <- solution_STKDE$spatial.z$yrow
grid_STKDE <- expand.grid(x_STKDE, y_STKDE)
FEMbasis_STKDE <- create.FEM.basis(create.mesh.2D(grid_STKDE))

### STKDE-discrete: Spatio-Temporal Kernel Density Estimation ------------------
# Grid Size
n_STKDE_discrete <- n

# Solution
x11(width = 25, height = 14)
solution_STKDE_discrete <- stkde(xlong = locations[,1], ylat = locations[,2], ztime = discrete_times,
                                 xgrids = n_STKDE_discrete, ygrids = n_STKDE_discrete, breaks = 0.05, alpha = 0.05, nrowspar = 5, bwmethod = "normal-reference") # alternative: bwmethod = "cv.ml"
x_STKDE_discrete <- seq(min(locations[,1]), max(locations[,1]), length.out = n_STKDE_discrete)
y_STKDE_discrete <- seq(min(locations[,2]), max(locations[,2]), length.out = n_STKDE_discrete)
grid_STKDE_discrete <- expand.grid(x_STKDE_discrete, y_STKDE_discrete)
FEMbasis_STKDE_discrete <- create.FEM.basis(create.mesh.2D(grid_STKDE_discrete))
dev.off()

### STKDE-separable: STKDE for 1st-Order Separable Spatio-Temporal Point Processes ------
# Spatial Component
f_space_STKDE_separable <- TDA::kde(locations, mesh.eval$nodes, 0.25)
idx <- which(f_space_STKDE_separable == 0)
f_space_STKDE_separable[idx] <- min(f_space_STKDE_separable[-idx])

# Temporal Component
f_time_STKDE_separable <- ks::kde(x = as.matrix(times), gridsize = 50)

### KNNSTDE-separable: k-Nearest Neighbors Spatio-Temporal Density Estimation -----------
###  for 1st-Order Separable Spatio-Temporal Point Processes
# Spatial Component
f_space_KNNSTDE_separable <- knnDE(as.matrix(locations), as.matrix(mesh.eval$nodes), 250)

# Temporal Component
f_time_KNNSTDE_separable <- knnDE(as.matrix(times), as.matrix(seq(0-0.25, 10+0.25, length.out = 25)), 50)

## ESTIMATES -------------------------------------------------------------------
# Time half width
h <- 1/6

# Time instants at which the solution is evaluated (midpoints of each month)
t <- seq(from = h, to = 10-h, by = 2*h)

# Discrete time instants at which the solution is evaluated (midpoints of each 3 month-period)
t_discrete <- seq(0.5, 9.5, by = 1)

mean_sol_STDEPDE <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
mean_sol_LGCP <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t_discrete))
mean_sol_INLA <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(mesh_time))
mean_sol_INLA_aggregated <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(mesh_time))
mean_sol_STKDE <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
mean_sol_STKDE_discrete <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t_discrete))
mean_sol_STKDE_separable <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
mean_sol_KNNSTDE_separable <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))

for (time_index in 1:length(t)) {
  
  # STDE-PDE
  evaluation_STDEPDE <- eval.FEM.time(FEM.time = FEMfunction_STDEPDE, locations = mesh.eval$nodes, time.instants = t[time_index])
  evaluation_STDEPDE <- exp(evaluation_STDEPDE)
  evaluation_STDEPDE[idx_remove] <- NA
  evaluation_STDEPDE <- evaluation_STDEPDE/sum(evaluation_STDEPDE, na.rm = TRUE)
  
  if(sum(abs(t[time_index] - mesh_time) < 1e-06)){
    
    time_index_INLA <- which(abs(t[time_index] - mesh_time) < 1e-06)
  
    # INLA
    coeff <- c(exp(solution_INLA[[time_index_INLA]]))
    FEM_INLA <- FEM(coeff, FEMbasis_INLA)
    evaluation_INLA <- eval.FEM(FEM_INLA, mesh.eval$nodes)
    evaluation_INLA[idx_remove] <- NA
    evaluation_INLA <- evaluation_INLA/sum(evaluation_INLA, na.rm = TRUE)
    
    # INLA-aggregated
    coeff <- c(exp(solution_INLA_aggregated[[time_index_INLA]]))
    FEM_INLA_aggregated <- FEM(coeff, FEMbasis_INLA_aggregated)
    evaluation_INLA_aggregated <- eval.FEM(FEM_INLA_aggregated, mesh.eval$nodes)
    evaluation_INLA_aggregated[idx_remove] <- NA
    evaluation_INLA_aggregated <- evaluation_INLA_aggregated/sum(evaluation_INLA_aggregated, na.rm = TRUE)
    
  }
  
  # STKDE
  coeff <- c(t(solution_STKDE$z[[time_index]]$v))
  FEM_STKDE <- FEM(coeff, FEMbasis_STKDE)
  evaluation_STKDE <- eval.FEM(FEM_STKDE, locations = mesh.eval$nodes)
  evaluation_STKDE[idx_remove] <- NA
  evaluation_STKDE <- evaluation_STKDE/sum(evaluation_STKDE, na.rm = TRUE)
  
  if(sum(abs(t[time_index] - t_discrete) < 1e-06)){
    time_index_discrete <- which(abs(t[time_index] - t_discrete) < 1e-06)
    
    # STKDE-discrete
    coeff <- c(solution_STKDE_discrete$dens[,,time_index_discrete])
    FEM_STKDE_discrete <- FEM(coeff, FEMbasis_STKDE_discrete)
    evaluation_STKDE_discrete <- eval.FEM(FEM_STKDE_discrete, locations = mesh.eval$nodes)
    evaluation_STKDE_discrete[idx_remove] <- NA
    evaluation_STKDE_discrete <- evaluation_STKDE_discrete/sum(evaluation_STKDE_discrete, na.rm = TRUE)
    
    # LGCP
    evaluation_LGCP <- intensity_LGCP[["grid"]][[time_index_discrete]]
    fftgr <- discreteWindow(solution_LGCP)
    # NA Values Outside the Spatial Domain (LGCP enlarges the window)
    for(i in 1:dim(evaluation_LGCP)[1]){
      for(j in 1:dim(evaluation_LGCP)[1]){
        if(!fftgr[i,j])
          evaluation_LGCP[i,j] <- NA
      }
    }
    FEM_LGCP <- FEM(c(evaluation_LGCP), FEMbasis_LGCP)
    evaluation_LGCP <- eval.FEM(FEM_LGCP, locations = mesh.eval$nodes)
    evaluation_LGCP[idx_remove] <- NA
    evaluation_LGCP <- evaluation_LGCP/sum(evaluation_LGCP, na.rm = TRUE)
    
  }  
  
  # STKDE-separable
  marginal_time_STKDE <- f_time_STKDE_separable$estimate[findInterval(t[time_index], f_time_STKDE_separable$eval.points)]
  solution_STKDE_separable <- marginal_time_STKDE * f_space_STKDE_separable
  solution_STKDE_separable[idx_remove] <- NA
  solution_STKDE_separable <- solution_STKDE_separable/sum(solution_STKDE_separable, na.rm = TRUE)
  
  # KNNSTDE-separable
  marginal_time_KNNSTDE_separable <- f_time_KNNSTDE_separable[findInterval(t[time_index], seq(0, 1, length.out=25))]
  solution_KNNSTDE_separable <- marginal_time_KNNSTDE_separable * f_space_KNNSTDE_separable
  solution_KNNSTDE_separable[idx_remove] <- NA
  solution_KNNSTDE_separable <- solution_KNNSTDE_separable/sum(solution_KNNSTDE_separable, na.rm = TRUE)
  
  mean_sol_STDEPDE[,time_index] <- mapply(sum, mean_sol_STDEPDE[,time_index], evaluation_STDEPDE, na.rm = TRUE)
  
  if(sum(abs(t[time_index] - mesh_time) < 1e-06)){
    
    time_index_INLA <- which(abs(t[time_index] - mesh_time) < 1e-06)
    
    mean_sol_INLA[,time_index_INLA] <- mapply(sum, mean_sol_INLA[,time_index_INLA], evaluation_INLA, na.rm = TRUE)
    mean_sol_INLA_aggregated[,time_index_INLA] <- mapply(sum, mean_sol_INLA_aggregated[,time_index_INLA], evaluation_INLA_aggregated, na.rm = TRUE)
  
  }
  
  mean_sol_STKDE[,time_index] <- mapply(sum, mean_sol_STKDE[,time_index], evaluation_STKDE, na.rm = TRUE)
  
  if(sum(abs(t[time_index] - t_discrete) < 1e-06)){
    
    mean_sol_STKDE_discrete[,time_index_discrete] <- mapply(sum, mean_sol_STKDE_discrete[,time_index_discrete], evaluation_STKDE_discrete, na.rm = TRUE)
    mean_sol_LGCP[,time_index_discrete] <- mapply(sum, mean_sol_LGCP[,time_index_discrete], evaluation_LGCP, na.rm = TRUE)
    
  }
  
  mean_sol_STKDE_separable[,time_index] <- mapply(sum, mean_sol_STKDE_separable[,time_index], solution_STKDE_separable, na.rm = TRUE)
  mean_sol_KNNSTDE_separable[,time_index] <- mapply(sum, mean_sol_KNNSTDE_separable[,time_index], solution_KNNSTDE_separable, na.rm = TRUE)
    
}

#'############################################################################'#
#'############################################################################'#

## EXPORT RESULTS --------------------------------------------------------------
# Export estimates
dir.create(file.path(getwd(), "output"), showWarnings = FALSE)

estimates <- list(mean_sol_STDEPDE, mean_sol_LGCP, mean_sol_INLA,
                  mean_sol_INLA_aggregated, mean_sol_STKDE, mean_sol_STKDE_discrete, mean_sol_STKDE_separable,
                  mean_sol_KNNSTDE_separable)
save(estimates, file = paste0("output/estimates.rda"))

base::save.image("output/workspace.RData")

# Remove the content of the LGCP directory
unlink(paste0("LGCP/*"), recursive = T)

# Alternative: Remove LGCP files from directory
# for(proc in 1:processes){
#   unlink(paste0("LGCP/LGCP_",proc,"/*"))
# }

#'############################################################################'#
#'############################################################################'#

## VISUALIZATION ---------------------------------------------------------------
{
  dir.create(file.path(getwd(), "pictures"), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/sample")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/sample_history")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/STDEPDE")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/LGCP")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/INLA")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/INLA_aggregated")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/STKDE")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/STKDE_discrete")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/STKDE_separable")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/KNNSTDE_separable")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/time_bar")), showWarnings = FALSE)
  
  ### MAP PREPROCESSING --------------------------------------------------------
  # Convert boundary_nodes back to Original Coordinates (Northing - Easting)
  boundary_nodes <- boundary_nodes * 1000
  data <- cbind(xyt$x, xyt$y) * 1000
  mesh$nodes <- mesh$nodes * 1000

  # From British National Grid CRS (EPSG:27700) to latitude/longitude CRS (EPSG:4326 - WSG 84)
  boundary_df <- as.data.frame(boundary_nodes[,2:1])
  colnames(boundary_df) <- c("Easting", "Northing")
  boundary_sp <- boundary_df
  coordinates(boundary_sp)<-~Northing+Easting
  proj4string(boundary_sp)<-CRS("+init=epsg:27700")
  boundary_sp <- spTransform(boundary_sp, CRS("+init=epsg:4326"))
  
  data_df <- as.data.frame(data[,2:1])
  colnames(data_df) <- c("Easting", "Northing")
  data_sp <- data_df
  coordinates(data_sp)<-~Northing+Easting
  proj4string(data_sp)<-CRS("+init=epsg:27700")
  data_sp <- spTransform(data_sp, CRS("+init=epsg:4326"))
  
  boundary_latlong <- as.data.frame(boundary_sp@coords)
  colnames(boundary_latlong) <- c("long","lat")
  data_latlong <- as.data.frame(data_sp@coords)
  colnames(data_latlong) <- c("long","lat")
  
  # Center of the Southampton Region
  center_x <- round((max(boundary_nodes[,1]) + min(boundary_nodes[,1]))/2)
  center_y <- round((max(boundary_nodes[,2]) + min(boundary_nodes[,2]))/2)
  
  center_df <- as.data.frame(cbind(center_y, center_x))
  colnames(center_df) <- c("Easting", "Northing")
  center_sp <- center_df
  coordinates(center_sp)<-~Northing+Easting
  proj4string(center_sp)<-CRS("+init=epsg:27700")
  center_sp <- spTransform(center_sp, CRS("+init=epsg:4326"))
  
  center_latlong <- as.data.frame(center_sp@coords)
  colnames(center_latlong) <- c("long","lat")
  
  zoom <- min(MaxZoom(range(boundary_latlong$lat), range(boundary_latlong$long)))
  
  # Zoom of the Map
  # SouthamptonZoom <- GetMap(center = c(lat = center_latlong$lat, lon = center_latlong$long), destfile = "map/mapTiles/Southampton.z9.png", 
  #                           GRAYSCALE = TRUE, zoom = zoom, SCALE = 2, maptype = "terrain", extraURL = "&scale=2")
  
  # SouthamptonZoom <- GetMap.bbox(lonR = c(min(boundary_latlong$long), max(boundary_latlong$long)),
  #                                latR = c(min(boundary_latlong$lat), max(boundary_latlong$lat)),
  #                                center = c(lat = center_latlong$lat, lon = center_latlong$long),
  #                                size = c(640, 640), destfile = "map/mapTiles/MyTile.png", MINIMUMSIZE = FALSE,
  #                                RETURNIMAGE = TRUE, GRAYSCALE = FALSE, NEWMAP = TRUE, verbose = 0,
  #                                SCALE = 1, type = "terrain", extraURL = "&scale=2", zoom = 9)
  
  # Map
  SouthamptonZoom <- ReadMapTile(destfile = "map/mapTiles/Southampton.z9.png")
  #SouthamptonZoom$myTile <- RGB2GRAY(SouthamptonZoom$myTile)
  #PlotOnStaticMap(SouthamptonZoom)
  
  ## MESH [FIGURE 3] -----------------------------------------------------------
  string <- paste0("pictures/app1_mesh.png")
  png(string, width = 1280, height = 1280)
  Figure.MeshOnMap(mesh = mesh,
                   boundary = boundary_latlong, map = SouthamptonZoom,
                   col = "black", cex = 1.2)
  
  s <- sample(1:nrow(locations), 1000)
  
  PlotOnStaticMap(SouthamptonZoom,
                  lat = data_latlong$lat[s],
                  lon = data_latlong$long[s],
                  pch = 19, col = rgb(190, 0, 0, max = 255), add = TRUE, TrueProj = TRUE, cex = 1.5, size = c(1280, 1280))
  
  dev.off()
  img <- magick::image_read(string)
  img_cropped <- magick::image_crop(img, geometry=paste0(1000,"x",900,"+0+0"),
                                    gravity = "Center", repage = TRUE)
  magick::image_write(img_cropped, path = string)
  
  ## DATA ON MAP [FIGURE 1] ----------------------------------------------------
  # Time half width
  h <- 1/6
  
  # Time instants at which the solution is evaluated (midpoints of each month)
  t <- seq(from = h, to = 10-h, by = 2*h)
  
  # Discrete time instants at which the solution is evaluated (midpoints of each 3 month-period)
  t_discrete <- seq(0.5, 9.5, by = 1)
  
  for(time_index in 1:length(t)){
    # First Plot: Data with Past Locations
    string <- paste0("pictures/sample_history/app1_sample_history_", time_index,".png")
    png(string, width = 1280, height = 1280)
    if(time_index != 1){
      Figure.DatasetOnMap(data = data_latlong[which(times < t[time_index]),],
                          boundary = boundary_latlong, map = SouthamptonZoom,
                          col = "black", cex = 1.5)
      
      PlotOnStaticMap(SouthamptonZoom,
                      lat = data_latlong$lat[which(abs(times-(t[time_index]+h)) < h)],
                      lon = data_latlong$long[which(abs(times-(t[time_index]+h)) < h)],
                      pch = 19, col = rgb(190, 0, 0, max = 255), add = TRUE, TrueProj = TRUE, cex = 1.5, size = c(1280, 1280))
    } else {
      Figure.DatasetOnMap(data = data_latlong[which(abs(times-(t[time_index]+h)) < h),],
                          boundary = boundary_latlong, map = SouthamptonZoom,
                          col = rgb(190, 0, 0, max = 255), cex = 1.5)
    }
    dev.off()
    img <- magick::image_read(string)
    img_cropped <- magick::image_crop(img, geometry = paste0(1000,"x",900,"+0+0"),
                                      gravity = "Center", repage = TRUE)
    magick::image_write(img_cropped, path = string)
    
    # Second Plot: Only Current Data at Each Given Time Instant
    string <- paste0("pictures/sample/app1_sample_", time_index,".png")
    png(string, width = 1280, height = 1280)
    Figure.DatasetOnMap(data = data_latlong[which(abs(times-t[time_index]) < h),],
                        boundary = boundary_latlong, map = SouthamptonZoom,
                        col = rgb(190, 0, 0, max = 255), cex = 1.5)
    dev.off()
    img <- magick::image_read(string)
    img_cropped <- magick::image_crop(img, geometry = paste0(1000,"x",900,"+0+0"),
                                      gravity = "Center", repage = TRUE)
    magick::image_write(img_cropped, path = string)
  }

  ## ESTIMATES [FIGURE 9] ------------------------------------------------------
  # Time half width
  h <- 1/6
  
  # Time instants at which the solution is evaluated (midpoints of each month)
  t <- seq(from = h, to = 10-h, by = 2*h)
  
  # Discrete time instants at which the solution is evaluated (midpoints of each 3 month-period)
  t_discrete <- seq(0.5, 9.5, by = 1)
  
  # Plots
  M <- max(max(mean_sol_STDEPDE, na.rm = TRUE),
           max(mean_sol_LGCP, na.rm = TRUE),
           max(mean_sol_INLA, na.rm = TRUE),
           max(mean_sol_INLA_aggregated, na.rm = TRUE),
           max(mean_sol_STKDE, na.rm = TRUE),
           max(mean_sol_STKDE_discrete, na.rm = TRUE),
           max(mean_sol_STKDE_separable, na.rm = TRUE),
           max(mean_sol_KNNSTDE_separable, na.rm = TRUE), na.rm = TRUE)
  
  for(time_index in 1:length(t)) {
    t_img <- proc.time()
    
    # Third Plot: STDE-PDE Estimated Density at t[time_index]
    string <- paste0("pictures/STDEPDE/app1_STDEPDE_",time_index,".png")
    png(string, width = 975, height = 950)
    
    mean_sol_STDEPDE[idx_remove,time_index] <- NA
    Figure.DensityOnMap(evaluation = matrix(mean_sol_STDEPDE[,time_index], n, n), grid_nodes = mesh.eval$nodes, M = M)
    
    dev.off()
    
    # Fourth Plot: STKDE Estimated Density at t[time_index]
    string <- paste0("pictures/STKDE/app1_STKDE_",time_index,".png")
    png(string, width = 1280, height = 1280)
    
    mean_sol_STKDE[idx_remove,time_index] <- NA
    Figure.DensityOnMap(evaluation = matrix(mean_sol_STKDE[,time_index], n, n), grid_nodes = mesh.eval$nodes, M = M)
    
    dev.off()
    
    if(sum(abs(t[time_index] - t_discrete) < 1e-06)){
      time_index_discrete <- which(abs(t[time_index] - t_discrete) < 1e-06)
      
      # Fifth Plot: STKDE-discrete Estimated Density at t[time_index]
      string <- paste0("pictures/STKDE_discrete/app1_STKDE_discrete_",time_index,".png")
      png(string, width = 1280, height = 1280)
      
      mean_sol_STKDE_discrete[idx_remove,time_index_discrete] <- NA
      Figure.DensityOnMap(evaluation = matrix(mean_sol_STKDE_discrete[,time_index_discrete], n, n), grid_nodes = mesh.eval$nodes, M = M)
      
      dev.off()
      
      # Sixth Plot: LGCP Estimated Density at t[time_index]
      string <- paste0("pictures/LGCP/app1_LGCP_",time_index,".png")
      png(string, width = 1280, height = 1280)
      
      mean_sol_LGCP[idx_remove,time_index_discrete] <- NA
      Figure.DensityOnMap(evaluation = matrix(mean_sol_LGCP[,time_index_discrete], n, n), grid_nodes = mesh.eval$nodes, M = M)
      
      dev.off()
      
    }
    
    # Seventh Plot: STKDE-separable Estimated (Separable) Density at t[time_index]
    string <- paste0("pictures/STKDE_separable/app1_STKDE_separable_",time_index,".png")
    png(string, width = 1280, height = 1280)
    
    mean_sol_STKDE_separable[idx_remove,time_index] <- NA
    Figure.DensityOnMap(evaluation = matrix(mean_sol_STKDE_separable[,time_index], n, n), grid_nodes = mesh.eval$nodes, M = M)
  
    dev.off()
    
    # Eight Plot: KNNSTDE-separable Estimated (Separable) Density at t[time_index]
    string <- paste0("pictures/KNNSTDE_separable/app1_KNNSTDE_separable_",time_index,".png")
    png(string, width = 1280, height = 1280)
    
    mean_sol_KNNSTDE_separable[idx_remove,time_index] <- NA
    Figure.DensityOnMap(evaluation = matrix(mean_sol_KNNSTDE_separable[,time_index], n, n), grid_nodes = mesh.eval$nodes, M = M)
    
    dev.off()
    
    if(sum(abs(t[time_index] - mesh_time) < 1e-06)){
      time_index_INLA <- which(abs(t[time_index] - mesh_time) < 1e-06)
      
      # Ninth Plot: INLA Estimated Density at t[time_index]
      string <- paste0("pictures/INLA/app1_INLA_",time_index,".png")
      png(string, width = 1280, height = 1280)
      
      mean_sol_INLA[idx_remove,time_index_INLA] <- NA
      Figure.DensityOnMap(evaluation = matrix(mean_sol_INLA[,time_index_INLA], n, n), grid_nodes = mesh.eval$nodes, M = M)
      
      dev.off()
      
      # Tenth Plot: INLA-aggregated Estimated Density at t[time_index]
      string <- paste0("pictures/INLA_aggregated/app1_INLA_aggregated_",time_index,".png")
      png(string, width = 1280, height = 1280)
      
      mean_sol_INLA_aggregated[idx_remove,time_index_INLA] <- NA
      Figure.DensityOnMap(evaluation = matrix(mean_sol_INLA_aggregated[,time_index_INLA], n, n), grid_nodes = mesh.eval$nodes, M = M)
      
      dev.off()
      
    }
    
    print(paste("Images for t =", round(t[time_index],3), "done in", round(as.numeric(proc.time() - t_img)[3]), "seconds."))
    
  }
  
  # for(time_index in c(5,11)) {
  #   t_img <- proc.time()
  #   
  #   # Third Plot: STDE-PDE Estimated Density at t[time_index]
  #   string <- paste0("pictures/STDEPDE/app1_STDEPDE_",time_index,".png")
  #   png(string, width = 1280, height = 1280)
  # 
  #   Figure.DensityOnMap2(mesh = mesh.eval, boundary = boundary_latlong, map = SouthamptonZoom,
  #                        evaluation = mean_sol_STDEPDE[,time_index], col = "black", cex = 1.2)
  # 
  #   dev.off()
  #   img <- magick::image_read(string)
  #   img_cropped <- magick::image_crop(img, geometry = paste0(1000,"x",900,"+0+0"),
  #                                     gravity = "Center", repage = TRUE)
  #   magick::image_write(img_cropped, path = string)
  # 
  #   # Fourth Plot: STKDE Estimated Density at t[time_index]
  #   string <- paste0("pictures/STKDE/app1_STKDE_",time_index,".png")
  #   png(string, width = 1280, height = 1280)
  #   
  #   Figure.DensityOnMap2(mesh = mesh.eval, boundary = boundary_latlong, map = SouthamptonZoom,
  #                        evaluation = mean_sol_STKDE[,time_index], col = "black", cex = 1.2)
  #   dev.off()
  #   img <- magick::image_read(string)
  #   img_cropped <- magick::image_crop(img, geometry = paste0(1000,"x",900,"+0+0"),
  #                                     gravity = "Center", repage = TRUE)
  #   magick::image_write(img_cropped, path = string)
  #   
  #   if(sum(abs(t[time_index] - t_discrete) < 1e-06)){
  #     time_index_discrete <- which(abs(t[time_index] - t_discrete) < 1e-06)
  #     
  #     # Fifth Plot: STKDE-discrete Estimated Density at t[time_index]
  #     string <- paste0("pictures/STKDE_discrete/app1_STKDE_discrete_",time_index,".png")
  #     png(string, width = 1280, height = 1280)
  #     
  #     Figure.DensityOnMap2(mesh = mesh.eval, boundary = boundary_latlong, map = SouthamptonZoom,
  #                          evaluation = mean_sol_STKDE_discrete[,time_index_discrete], col = "black", cex = 1.2)
  #     dev.off()
  #     img <- magick::image_read(string)
  #     img_cropped <- magick::image_crop(img, geometry = paste0(1000,"x",900,"+0+0"),
  #                                       gravity = "Center", repage = TRUE)
  #     magick::image_write(img_cropped, path = string)
  #     
  #     # Sixth Plot: LGCP Estimated Density at t[time_index]
  #     string <- paste0("pictures/LGCP/app1_LGCP_",time_index,".png")
  #     png(string, width = 1280, height = 1280)
  #     
  #     Figure.DensityOnMap2(mesh = mesh.eval, boundary = boundary_latlong, map = SouthamptonZoom,
  #                          evaluation = mean_sol_LGCP[,time_index_discrete], col = "black", cex = 1.2)
  #     dev.off()
  #     img <- magick::image_read(string)
  #     img_cropped <- magick::image_crop(img, geometry = paste0(1000,"x",900,"+0+0"),
  #                                       gravity = "Center", repage = TRUE)
  #     magick::image_write(img_cropped, path = string)
  #     
  #   }
  #   
  #   # Seventh Plot: STKDE-separable Estimated (Separable) Density at t[time_index]
  #   string <- paste0("pictures/STKDE_separable/app1_STKDE_separable_",time_index,".png")
  #   png(string, width = 1280, height = 1280)
  #   
  #   Figure.DensityOnMap2(mesh = mesh.eval, boundary = boundary_latlong, map = SouthamptonZoom,
  #                        evaluation = mean_sol_STKDE_separable[,time_index], col = "black", cex = 1.2)
  #   dev.off()
  #   img <- magick::image_read(string)
  #   img_cropped <- magick::image_crop(img, geometry = paste0(1000,"x",900,"+0+0"),
  #                                     gravity = "Center", repage = TRUE)
  #   magick::image_write(img_cropped, path = string)
  #   
  #   # Eight Plot: KNNSTDE-separable Estimated (Separable) Density at t[time_index]
  #   string <- paste0("pictures/KNNSTDE_separable/app1_KNNSTDE_separable_",time_index,".png")
  #   png(string, width = 1280, height = 1280)
  #   
  #   Figure.DensityOnMap2(mesh = mesh.eval, boundary = boundary_latlong, map = SouthamptonZoom,
  #                        evaluation = mean_sol_KNNSTDE_separable[,time_index], col = "black", cex = 1.2)
  #   dev.off()
  #   img <- magick::image_read(string)
  #   img_cropped <- magick::image_crop(img, geometry = paste0(1000,"x",900,"+0+0"),
  #                                     gravity = "Center", repage = TRUE)
  #   magick::image_write(img_cropped, path = string)
  #   
  #   if(sum(abs(t[time_index] - mesh_time) < 1e-06)){
  #     time_index_INLA <- which(abs(t[time_index] - mesh_time) < 1e-06)
  #     
  #     # Ninth Plot: INLA Estimated Density at t[time_index]
  #     string <- paste0("pictures/INLA/app1_INLA_",time_index,".png")
  #     png(string, width = 1280, height = 1280)
  #     
  #     Figure.DensityOnMap2(mesh = mesh.eval, boundary = boundary_latlong, map = SouthamptonZoom,
  #                          evaluation = mean_sol_INLA[,time_index_INLA], col = "black", cex = 1.2)
  #     dev.off()
  #     img <- magick::image_read(string)
  #     img_cropped <- magick::image_crop(img, geometry = paste0(1000,"x",900,"+0+0"),
  #                                       gravity = "Center", repage = TRUE)
  #     magick::image_write(img_cropped, path = string)
  #     
  #     # Tenth Plot: INLA-aggregated Estimated Density at t[time_index]
  #     string <- paste0("pictures/INLA_aggregated/app1_INLA_aggregated_",time_index,".png")
  #     png(string, width = 1280, height = 1280)
  #     
  #     Figure.DensityOnMap2(mesh = mesh.eval, boundary = boundary_latlong, map = SouthamptonZoom,
  #                          evaluation = mean_sol_INLA_aggregated[,time_index_INLA], col = "black", cex = 1.2)
  #     dev.off()
  #     img <- magick::image_read(string)
  #     img_cropped <- magick::image_crop(img, geometry = paste0(1000,"x",900,"+0+0"),
  #                                       gravity = "Center", repage = TRUE)
  #     magick::image_write(img_cropped, path = string)
  # 
  #   }
  #   
  #   print(paste("Images for t =", round(t[time_index],3), "done in", round(as.numeric(proc.time() - t_img)[3]), "seconds."))
  #   
  # }

  ## LEGEND [FIGURE 9] ---------------------------------------------------------
  p <- colorRampPalette(c("#0E1E44","#3E6DD8","#68D061","#ECAF53", "#EB5F5F","#E11F1C"))(1000)
  palette(p)
  
  pdf("pictures/app1_legend.pdf", family = "serif", width = 11, height = 3)
  #x11(width = 11, height = 3)
  par(mai = c(1,0.75,0,0))
  plot(c(0, 112.5), c(0, 15), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame.plot = F)
  gradient.rect(0, 0, 100, 5, col = p, border = "black")
  axis(1,at = c(0,33.33,66.66,100), labels = c("0", "2", "4", "6"),
       lwd.ticks = 2, cex.axis = 2, lwd = 2)
  text(107,2, TeX("$\\times 10^{-4}$"), cex = 2)
  dev.off()
  
  print(paste0("M = ", round(M, 4)))
  
  ## TIME BAR --------------------------------------------------------------------
  for(time_index in 1:length(t)) {
    string <- paste0("pictures/time_bar/app1_time_bar_",time_index,".pdf")
    pdf(file = string, family = "serif", width = 7, height = 1.5)
    #x11(width = 7, height = 1.5)
    par(mar = c(0.01, 2, 0.01, 2))
    bar <- rbind(c(t[time_index]/10,0),c(t[time_index]/10,0.075))
    plot(bar, type = "l", bty = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "",
         pch = 3, lwd = 5, xlim = c(0,1), asp = 1)
    mtext("Time", line = -1.5, font = 2, cex = 1.5)
    axis(side = 1, c(0, 0.5, 1), pos = 0, labels = c("0", "5", "10"))
    dev.off()
    
  }
  
}

#'############################################################################'#
#'############################################################################'#

# CROSS VALIDATION ERROR -------------------------------------------------------
boundary_nodes <- boundary_nodes / 1000
mesh$nodes <- mesh$nodes / 1000

# Grid Size for Estimation
n <- 64

# Fine Grid for Evaluation
X <- seq(min(boundary_nodes[,1]), max(boundary_nodes[,1]), length.out = n)
Y <- seq(min(boundary_nodes[,2]), max(boundary_nodes[,2]), length.out = n)
grid <- expand.grid(X, Y)
mesh.eval <- create.mesh.2D(grid)
FEMbasis.eval <- create.FEM.basis(mesh.eval)

# Clean Grid on the Southampton Region
grid_clean <- NULL
idx_remove <- NULL
for(p in 1:nrow(mesh.eval$nodes)){
  if(point.in.polygon(mesh.eval$nodes[p,1], mesh.eval$nodes[p,2], boundary_nodes[,1], boundary_nodes[,2])>0){
    grid_clean <- rbind(grid_clean, mesh.eval$nodes[p,])
  }
  else{
    idx_remove <- c(idx_remove, p)
  }
}

# Number of Folds
K <- 10

# Subsample Size in Each Fold
n_k <- floor(N/K)

# Indices
folds <- sample(1:N, N, replace = FALSE)

# Cross-Validation Error Containers
CV_error_STDEPDE <- rep(0,K)
CV_error_LGCP <- rep(0,K)
CV_error_INLA <- rep(0,K)
CV_error_INLA_aggregated <- rep(0,K)
CV_error_STKDE <- rep(0,K)
CV_error_STKDE_discrete <- rep(0,K)
CV_error_STKDE_separable <- rep(0,K)
CV_error_KNNSTDE_separable <- rep(0,K)

for(k in 1:K){
  
  ### PREPROCESSING ------------------------------------------------------------
  if(k != K){
    
    test_idx <- ((k-1)*n_k+1):(k*n_k)
    train_idx <- (1:N)[!((1:N) %in% test_idx)]
    
    locations_train_k <- solution_STDEPDE$data[folds[train_idx],]
    times_train_k <- solution_STDEPDE$data_time[folds[train_idx]]
    
    locations_test_k <- solution_STDEPDE$data[folds[test_idx],]
    times_test_k <- solution_STDEPDE$data_time[folds[test_idx]]
    
  } else {
    
    test_idx <- ((k-1)*n_k+1):N
    train_idx <- (1:N)[!((1:N) %in% test_idx)]
    
    locations_train_k <- solution_STDEPDE$data[folds[train_idx],]
    times_train_k <- solution_STDEPDE$data_time[folds[train_idx]]
    
    locations_test_k <- solution_STDEPDE$data[folds[test_idx],]
    times_test_k <- solution_STDEPDE$data_time[folds[test_idx]]
    
  }
  
  discrete_times_train_k <- times_train_k %/% 1
  discrete_times_test_k <- times_test_k %/% 1
  
  # Plot
  # x11()
  # plot(boundary_df$Northing, boundary_df$Easting, type = "l", main = paste0("Training and test sets, iter = ",k))
  # legend("topleft", legend = c(paste0("Train: ", nrow(locations_train_k)), paste0("Test: ", nrow(locations_test_k))), fill = c("red", "green"))
  # points(locations_train_k, col = "red", pch = 19)
  # points(locations_test_k, col = "green", pch = 19)
  
  ## STDE-PDE: Spatio-Temporal Density Estimation with PDE Regulariz. ----------
  lambda <- 0.1
  # alternative: lambda <- 10^seq(from = -2, to = 0, by = 1)
  
  lambda_time <- 0.001
  # alternative: lambda_time <- 10^seq(from = -4, to = -2, by = 1)
  
  # Solution
  # [If lambda and/or lambda_time are vectors, to select the best proposals by
  # 10-folds CV, please specify: preprocess_method = "RightCV"]
  solution_STDEPDE_train <- DE.FEM.time(data = locations_train_k, data_time = times_train_k, FEMbasis = FEMbasis,
                                        mesh_time = mesh_time, lambda = lambda, tol1 = 1e-7, 
                                        lambda_time = lambda_time, fvec = NULL, heatStep = 0.1,
                                        heatIter = 10, print = F, nfolds = 10, nsimulations = 10000,
                                        step_method = "Wolfe_Method", direction_method = "L-BFGS5",
                                        preprocess_method = "NoCrossValidation")
  
  FEMfunction_STDEPDE_train <- FEM.time(solution_STDEPDE_train$g, mesh_time, FEMbasis, FLAG_PARABOLIC = FALSE)
  
  ## LGCP: Log-Gaussian Cox Process --------------------------------------------
  # Spatio-Temporal Planar Point Process
  win <- xyt$window
  tlim <- c(0, 10)
  xy_LGCP <- ppp(x = locations_train_k[,1], y = locations_train_k[,2], window = win)
  xyt_LGCP <- stppp(xy_LGCP, times_train_k, tlim = tlim)
  xyt_LGCP[["n"]] <- N
  xyt_LGCP <- integerise(xyt_LGCP)
  
  # Kernel Smoothed Intensity Function
  den <- density.ppp(xyt_LGCP, sigma = 1)
  sar <- spatialAtRisk(den)
  mut <- muEst(xyt_LGCP)
  
  exceed <- exceedProbs(c(1.5,2,3))
  
  # Temporary Folder
  dir.create(file.path(getwd(), "LGCP"), showWarnings = FALSE)
  tmpdr <- paste0(getwd(), "/LGCP")
  
  # Grid Size
  n_LGCP <- n/2
  
  # Solution
  solution_LGCP_train <- lgcpPredict(xyt = xyt_LGCP, T = as.integer(9), laglength = as.integer(9),
                                     model.parameters = lgcppars(sigma = 1.6, phi = 1.9, theta = 1.4),
                                     gridsize = n_LGCP, spatial.intensity = sar, temporal.intensity = mut,
                                     mcmc.control = mcmcpars(mala.length = 120000, burnin = 20000, retain = 100,
                                                           adaptivescheme = andrieuthomsh(inith = 1, alpha = 0.5, C = 1, targetacceptance = 0.574)),
                                     output.control = setoutput(gridfunction = dump2dir(dirname = tmpdr, forceSave = TRUE),
                                                              gridmeans = MonteCarloAverage("exceed"))) # to be run locally
  intensity_LGCP_train <- intens(solution_LGCP_train)
  grid_LGCP_train <- expand.grid(intensity_LGCP_train$xvals, intensity_LGCP_train$yvals)
  FEMbasis_LGCP_train <- create.FEM.basis(create.mesh.2D(grid_LGCP_train))
  
  ## INLA: SPDE Spatio-Temporal Model ------------------------------------
  # Temporal Discretization
  k_INLA <- length(mesh_time)
  mesh_time_INLA <- INLA::inla.mesh.1d(mesh_time)
  
  # Spatial Discretization
  min_bdy_x <- min(boundary_nodes[,1])
  max_bdy_x <- max(boundary_nodes[,1])
  min_bdy_y <- min(boundary_nodes[,2])
  max_bdy_y <- max(boundary_nodes[,2])
  
  domain_INLA <- SpatialPolygons(list(Polygons(list(Polygon(boundary_nodes)),"0")))
  mesh_INLA <- INLA::inla.mesh.2d(boundary = INLA::inla.sp2segment(domain_INLA),
                                  max.edge = c(6,6), cutoff = 0)
  # alternative: pass locations as additional argument to build the mesh through R-INLA
  
  # Plot of the Spatial Mesh
  # x11()
  # plot(mesh_INLA)
  # points(boundary_nodes, col = "red", pch = 19)
  
  # SPDE Model
  spde <- INLA::inla.spde2.pcmatern(mesh = mesh_INLA, prior.range = c(5, 0.01),
                                    prior.sigma = c(1, 0.01))
  m <- spde$n.spde
  
  # Spatio-Temporal Projection Matrix: Kronecker Product between the Spatial
  # Projector Matrix and the Group Projector one, i.e., the Temporal Dimension
  Ast <- INLA::inla.spde.make.A(mesh = mesh_INLA, loc = as.matrix(locations_train_k),
                                n.group = length(mesh_time_INLA$n),
                                group = times_train_k, group.mesh = mesh_time_INLA)
  
  # Space-Time Index Set
  idx_INLA <- INLA::inla.spde.make.index("s", spde$n.spde, n.group = mesh_time_INLA$n)
  
  # Dual Mesh
  dmesh <- book.mesh.dual(mesh_INLA)
  
  # Intersection with each Polygon from the Dual Mesh
  w <- sapply(1:length(dmesh), function(i) {
    if (gIntersects(dmesh[i,], domain_INLA))
      return(gArea(gIntersection(dmesh[i,], domain_INLA)))
    else return(0)
  })
  
  #gArea(domain_inla)
  
  # Spatio-Temporal Volume
  st.vol <- rep(w, k_INLA) * rep(diag(INLA::inla.mesh.fem(mesh_time_INLA)$c0), m)
  
  # Data Stack
  y <- rep(0:1, c(mesh_time_INLA$n * m, nrow(locations_train_k)))
  expected <- c(st.vol, rep(0, nrow(locations_train_k)))
  stk <- INLA::inla.stack(
    data = list(y = y, expect = expected), 
    A = list(rbind(Diagonal(n = k_INLA * m), Ast), 1), 
    effects = list(idx_INLA, list(a0 = rep(1, k_INLA * m + nrow(locations_train_k)))))
  
  # Model Fitting with the Gaussian Approximation
  pcrho <- list(prior = "pccor1", param = c(0.7, 0.7))
  form <- y ~ 0 + a0 + f(s, model = spde, group = s.group, 
                         control.group = list(model = "ar1",
                                              hyper = list(theta = pcrho)))
  
  res <- INLA::inla(form, family = "poisson", 
                    data = INLA::inla.stack.data(stk), E = expect,
                    control.predictor = list(A = INLA::inla.stack.A(stk)),
                    control.inla = list(strategy = "adaptive"))
  
  # NOTE: the exponential of the intercept plus the random effect at each
  #       space-time integration point is the relative risk at each of these
  #       points. This relative risk times the space-time volume will give the
  #       expected number of points (E(n)) at each one of these space-time
  #       locations. Summing over them will give a value that approaches the
  #       number of observations
  
  eta.at.integration.points <- res$summary.fix[1,1] + res$summary.ran$s$mean
  #c(n = nrow(locations_train_k), "E(n)" = sum(st.vol * exp(eta.at.integration.points)))
  
  # Grid Size
  n_INLA <- n
  
  # Projection over a Grid for each Time Knot
  r0 <- diff(range(boundary_nodes[, 1])) / diff(range(boundary_nodes[, 2]))
  prj <- INLA::inla.mesh.projector(mesh_INLA, xlim = range(boundary_nodes[, 1]),
                                   ylim = range(boundary_nodes[, 2]), dims = c(n_INLA, n_INLA)) 
  ov <- over(SpatialPoints(prj$lattice$loc), domain_INLA)
  m.prj <- lapply(1:k_INLA, function(j) {
    r <- INLA::inla.mesh.project(prj, res$summary.ran$s$mean[1:m + (j - 1) * m])
    # alternative: r <- INLA::inla.mesh.project(prj, exp(res_INLA$summary.fix$mean + res_INLA$summary.ran$s$mean[1:m + (j - 1) * m]) * st.vol[1:m + (j - 1) * m] / N)
    r[is.na(ov)] <- NA
    return(r) 
  })
  
  # Plot of the Fitted Latent Field at each Time Knot
  # igr <- apply(abs(outer(times_train_k, mesh_time_INLA$loc, "-")), 1, which.min)
  # zlm <- range(unlist(m.prj), na.rm = TRUE)
  # x11()
  # par(mfrow = c(2, 3), mar = c(0, 0, 0.9, 0))
  # for (j in 1:k_INLA) {
  # image2D(x = prj$x, y = prj$y, z = m.prj[[j]], asp = 1,
  #         xlab = "", zlim = zlm, axes = FALSE, col = heat.colors(100),
  #         main = paste0("Time: ", j))
  # #points(locations_train_k[igr == j, 1:2], pch = 19)
  # }
  
  # Solution
  x_INLA <- prj$x
  y_INLA <- prj$y
  grid_INLA <- expand.grid(x_INLA, y_INLA)
  solution_INLA_train <- m.prj
  FEMbasis_INLA_train <- create.FEM.basis(create.mesh.2D(grid_INLA))
  
  ## INLA-aggregated: SPDE Spatio-Temporal Aggregated Model (Large Dataset) ---------
  # Temporal Discretization
  k_INLA <- length(mesh_time)
  mesh_time_INLA <- INLA::inla.mesh.1d(mesh_time)
  
  # Spatial Discretization
  domain_INLA <- SpatialPolygons(list(Polygons(list(Polygon(boundary_nodes)),"0")))
  mesh_INLA <- INLA::inla.mesh.2d(boundary = INLA::inla.sp2segment(domain_INLA),
                                  max.edge = c(6,6), cutoff = 0)
  # alternative: pass locations as additional argument to build the mesh through R-INLA
  
  # Spatio-Temporal Planar Point Pattern
  xyt_INLA <- stppp(xy_LGCP, times_train_k, tlim = c(as.integer(min(times_train_k)), as.integer(max(times_train_k))))
  
  # Space-Time Aggregation
  dd <- deldir(mesh_INLA$loc[, 1], mesh_INLA$loc[, 2])
  tiles <- tile.list(dd)
  polys <- SpatialPolygons(lapply(1:length(tiles), function(i) {
    p <- cbind(tiles[[i]]$x, tiles[[i]]$y)
    n <- nrow(p)
    Polygons(list(Polygon(p[c(1:n, 1), ])), i)
  }))
  area <- factor(over(SpatialPoints(cbind(xyt_INLA$x, xyt_INLA$y)), polys),
                 levels = 1:length(polys))
  
  t.breaks <- sort(c(mesh_time_INLA$loc[c(1, k_INLA)],
                     mesh_time_INLA$loc[2:k_INLA - 1] / 2 + mesh_time_INLA$loc[2:k_INLA] / 2))
  time <- factor(findInterval(xyt_INLA$t, t.breaks), levels = 1:(length(t.breaks) - 1))
  
  # Distribution of Data Points on the Time Knots
  # table(time)
  
  agg.dat <- as.data.frame(table(area, time))
  for(j in 1:2){
    agg.dat[[j]] <- as.integer(as.character(agg.dat[[j]])) 
  }
  
  # Areas of the Polygons
  w.areas <- sapply(1:length(tiles), function(i) {
    p <- cbind(tiles[[i]]$x, tiles[[i]]$y)
    n <- nrow(p)
    pl <- SpatialPolygons(
      list(Polygons(list(Polygon(p[c(1:n, 1),])), i)))
    if (gIntersects(pl, domain_INLA))
      return(gArea(gIntersection(pl, domain_INLA)))
    else return(0)
  })
  
  # Time Width
  w.t <- diag(INLA::inla.mesh.fem(mesh_time_INLA)$c0)
  
  # Space-Time Volume (Area Unit per Time Unit) at each Polygon and Time Knot
  e0 <- w.areas[agg.dat$area] * (w.t[agg.dat$time])
  
  # Projector Matrix
  A.st <- INLA::inla.spde.make.A(mesh_INLA, mesh_INLA$loc[agg.dat$area, ],
                                 group = agg.dat$time, mesh.group = mesh_time_INLA)
  
  # SPDE Model Object
  spde <- INLA::inla.spde2.pcmatern(mesh_INLA, prior.sigma = c(1,0.01),
                                    prior.range = c(0.05,0.01))
  
  # Space-Time Index Set
  idx_INLA <- INLA::inla.spde.make.index("s", spde$n.spde, n.group = k_INLA)
  
  # Data Stack
  stk <- INLA::inla.stack(data = list(y = agg.dat$Freq, exposure = e0),
                          A = list(A.st, 1),
                          effects = list(idx_INLA, list(b0 = rep(1, nrow(agg.dat)))))
  
  # PC Prior on Correlation
  pcrho <- list(theta = list(prior = "pccor1", param = c(0.7, 0.7)))
  
  # Model Formula: Intercept + Spatial Effect + Temporal Effect
  formula <- y ~ 0 + b0 + f(s, model = spde, group = s.group,
                            control.group = list(model = "ar1", hyper = pcrho))
  
  # Model Fitting
  res <- INLA::inla(formula, family = "poisson",
                    data = INLA::inla.stack.data(stk), E = exposure, 
                    control.predictor = list(A = INLA::inla.stack.A(stk)),
                    control.inla = list(strategy ="adaptive"))
  
  # Value of mu and Summary for the Intercept
  # cbind(True = mu, res$summary.fixed[, 1:6])
  
  # Summary for the Hyperparameters
  # cbind(True = c(range, sigma, rho), res$summary.hyperpar[, c(1, 2, 3, 5)])
  
  # Grid Size
  n_INLA_aggregated <- n
  
  # Projection over a Grid for each Time Knot
  r0 <- diff(range(boundary_nodes[, 1])) / diff(range(boundary_nodes[, 2]))
  prj <- INLA::inla.mesh.projector(mesh_INLA, xlim = c(min_bdy_x, max_bdy_x), 
                                   ylim = c(min_bdy_y, max_bdy_y),
                                   dims = c(n_INLA_aggregated, n_INLA_aggregated))
  g.no.in <- is.na(over(SpatialPoints(prj$lattice$loc), domain_INLA))
  t.mean <- lapply(1:k_INLA, function(j) {
    z.j <- res$summary.ran$s$mean[idx_INLA$s.group == j]
    # alternative: z.j <- exp(res_INLA_aggregated$summary.fix$mean + res_INLA_aggregated$summary.ran$s$mean[idx_INLA$s.group == j]) * w.areas / N
    z <- INLA::inla.mesh.project(prj, z.j)
    z[g.no.in] <- NA
    return(z)
  })
  
  # Plot of the Spatial Surface at each Time Knot
  # zlims <- range(unlist(t.mean), na.rm = TRUE)
  # x11()
  # par(mfrow = c(2, 3), mar = c(0.1, 0.1, 1, 0.1))
  # for (j in 1:k_INLA) {
  #   image(prj$x, prj$y, t.mean[[j]], axes = FALSE, zlim = zlims,
  #         col = heat.colors(100), main = paste0("Time knot: ", j))
  #   #points(xyt$x[time == j], xyt$y[time == j], cex = 0.1, cex.main = 0.95)
  # }
  
  # Solution
  x_INLA_aggregated <- prj$x
  y_INLA_aggregated <- prj$y
  grid_INLA_aggregated <- expand.grid(x_INLA_aggregated, y_INLA_aggregated)
  solution_INLA_aggregated_train <- t.mean
  FEMbasis_INLA_aggregated_train <- create.FEM.basis(create.mesh.2D(grid_INLA_aggregated))
  
  ## STKDE: Spatio-Temporal Kernel Density Estimation --------------------------
  # Spatio-Temporal Planar Point Pattern
  xy_STKDE <- ppp(x = locations_train_k[,1], y = locations_train_k[,2], window = win)
  
  # Grid Size
  n_STKDE <- n
  
  # Time half width
  h <- 1/6
  t <- seq(from = h, to = 10-h, by = 2*h)
  
  # Solution
  solution_STKDE_train <- spattemp.density(pp = xy_STKDE, h = NULL, tt = times_train_k,
                                           lambda = NULL, tlim = c(min(t)-h, max(t)+h),
                                           sedge = "none", tedge = "none", 
                                           sres = n_STKDE, tres = length(t))
  
  x_STKDE <- solution_STKDE_train$spatial.z$xcol
  y_STKDE <- solution_STKDE_train$spatial.z$yrow
  grid_STKDE <- expand.grid(x_STKDE, y_STKDE)
  FEMbasis_STKDE_train <- create.FEM.basis(create.mesh.2D(grid_STKDE))
  
  ## STKDE-discrete: Spatio-Temporal Kernel Density Estimation --------------------------
  # Grid Size
  n_STKDE_discrete <- n
  
  # Solution
  x11(width = 25, height = 14)
  idx <- which(locations_train_k[,1] >= min_bdy_x & locations_train_k[,1] <= max_bdy_x & locations_train_k[,2] >= min_bdy_y & locations_train_k[,2] <= max_bdy_y)
  solution_STKDE_discrete_train <- stkde(xlong = locations_train_k[idx,1], ylat = locations_train_k[idx,2], ztime = discrete_times_train_k[idx],
                          xgrids = n_STKDE_discrete, ygrids = n_STKDE_discrete, breaks = 0.05, alpha = 0.05, nrowspar = 5, bwmethod = "normal-reference") # alternative: bwmethod = "cv.ml"
  x_STKDE_discrete <- seq(min(locations_train_k[idx,1]), max(locations_train_k[idx,1]), length.out = n_STKDE_discrete)
  y_STKDE_discrete <- seq(min(locations_train_k[idx,2]), max(locations_train_k[idx,2]), length.out = n_STKDE_discrete)
  grid_STKDE_discrete <- expand.grid(x_STKDE_discrete, y_STKDE_discrete)
  FEMbasis_STKDE_discrete_train <- create.FEM.basis(create.mesh.2D(grid_STKDE_discrete))
  dev.off()
  
  ## STKDE-separable: STKDE for 1st-Order Separable Spatio-Temporal Point Processes -----
  # Spatial Component
  f_space_STKDE_separable_train <- TDA::kde(locations_train_k, mesh.eval$nodes, 0.25)
  idx <- which(f_space_STKDE_separable_train == 0)
  f_space_STKDE_separable_train[idx] <- min(f_space_STKDE_separable_train[-idx])
  
  # Temporal Component
  f_time_STKDE_separable_train <- ks::kde(x = as.matrix(times_train_k), gridsize = 50)
  
  ## KNNSTDE-separable: k-Nearest Neighbors Spatio-Temporal Density Estimation ----------
  ##  for 1st-Order Separable Spatio-Temporal Point Processes
  # Spatial Component
  f_space_KNNSTDE_separable_train <- knnDE(as.matrix(locations_train_k), as.matrix(mesh.eval$nodes), 250)
  
  # Temporal Component
  f_time_KNNSTDE_separable_train <- knnDE(as.matrix(times_train_k), as.matrix(seq(-0.25, 10+0.25, length.out = 25)), 50)
  
  ## ERRORS --------------------------------------------------------------------
  # Time half width
  h <- 1/6
  
  # Time instants at which the solution is evaluated (midpoints of each month)
  t <- seq(from = h, to = 10-h, by = 2*h)
  
  # Discrete time instants at which the solution is evaluated (midpoints of each 3 month-period)
  t_discrete <- seq(0.5, 9.5, by = 1)
  
  for(time_index in 1:length(t)){
    
    # STDEPDE
    evaluation_STDEPDE_train <- eval.FEM.time(FEM.time = FEMfunction_STDEPDE_train, locations = mesh.eval$nodes, time.instants = t[time_index])
    evaluation_STDEPDE_train <- exp(evaluation_STDEPDE_train)
    evaluation_STDEPDE_train[idx_remove] <- NA
    evaluation_STDEPDE_train <- evaluation_STDEPDE_train/sum(evaluation_STDEPDE_train, na.rm = TRUE)
    
    evaluation_STDEPDE_test <- eval.FEM(FEM = FEM(evaluation_STDEPDE_train, FEMbasis.eval), locations = locations_test_k[which(abs(t[time_index]-times_test_k)<h),])
    
    cardinality_k <- 1 / (sum(!is.na(evaluation_STDEPDE_test)) / nrow(locations_test_k[which(abs(t[time_index]-times_test_k)<h),]))
    CV_error_STDEPDE[k] <- CV_error_STDEPDE[k] + sum(evaluation_STDEPDE_train^2, na.rm = TRUE)/sum(!is.na(evaluation_STDEPDE_train)) -2 * sum(evaluation_STDEPDE_test, na.rm = TRUE) * cardinality_k
    
    if(sum(abs(t[time_index] - mesh_time) < 1e-06)){
      
      time_index_INLA <- which(abs(t[time_index] - mesh_time) < 1e-06)
      
      # INLA
      coeff <- c(exp(solution_INLA_train[[time_index_INLA]]))
      FEM_INLA_train <- FEM(coeff, FEMbasis_INLA_train)
      evaluation_INLA_train <- eval.FEM(FEM_INLA_train, locations = mesh.eval$nodes)
      evaluation_INLA_train[idx_remove] <- NA
      evaluation_INLA_train <- evaluation_INLA_train/sum(evaluation_INLA_train, na.rm = TRUE)
      
      FEM_INLA_test <- FEM(evaluation_INLA_train, FEMbasis.eval)
      evaluation_INLA_test <- eval.FEM(FEM_INLA_test, locations = locations_test_k[which(abs(t[time_index]-times_test_k)<h),])
      
      cardinality_k <- 1 / (sum(!is.na(evaluation_INLA_test)) / nrow(locations_test_k[which(abs(t[time_index]-times_test_k)<h),]))
      CV_error_INLA[k] <- CV_error_INLA[k] + sum(evaluation_INLA_train^2, na.rm = TRUE)/sum(!is.na(evaluation_INLA_train)) -2 * sum(evaluation_INLA_test, na.rm = TRUE) * cardinality_k
      
      # INLA-aggregated
      coeff <- c(exp(solution_INLA_aggregated_train[[time_index_INLA]]))
      FEM_INLA_aggregated_train <- FEM(coeff, FEMbasis_INLA_aggregated_train)
      evaluation_INLA_aggregated_train <- eval.FEM(FEM_INLA_aggregated_train, locations = mesh.eval$nodes)
      evaluation_INLA_aggregated_train[idx_remove] <- NA
      evaluation_INLA_aggregated_train <- evaluation_INLA_aggregated_train/sum(evaluation_INLA_aggregated_train, na.rm = TRUE)

      FEM_INLA_aggregated_test <- FEM(evaluation_INLA_aggregated_train, FEMbasis.eval)
      evaluation_INLA_aggregated_test <- eval.FEM(FEM_INLA_aggregated_test, locations = locations_test_k[which(abs(t[time_index]-times_test_k)<h),])
      
      cardinality_k <- 1 / (sum(!is.na(evaluation_INLA_aggregated_test)) / nrow(locations_test_k[which(abs(t[time_index]-times_test_k)<h),]))
      CV_error_INLA_aggregated[k] <- CV_error_INLA_aggregated[k] + sum(evaluation_INLA_aggregated_train^2, na.rm = TRUE)/sum(!is.na(evaluation_INLA_aggregated_train)) -2 * sum(evaluation_INLA_aggregated_test, na.rm = TRUE) * cardinality_k
      
    }
    
    # STKDE
    coeff <- c(solution_STKDE_train$z[[time_index]]$v)
    FEM_STKDE_train <- FEM(coeff, FEMbasis_STKDE_train)
    evaluation_STKDE_train <- eval.FEM(FEM_STKDE_train, locations = mesh.eval$nodes)
    evaluation_STKDE_train[idx_remove] <- NA
    evaluation_STKDE_train <- evaluation_STKDE_train/sum(evaluation_STKDE_train, na.rm = TRUE)
    
    FEM_STKDE_test <- FEM(evaluation_STKDE_train, FEMbasis.eval)
    evaluation_STKDE_test <- eval.FEM(FEM_STKDE_test, locations = locations_test_k[which(abs(t[time_index]-times_test_k)<h),])
    
    cardinality_k <- 1 / (sum(!is.na(evaluation_STKDE_test)) / nrow(locations_test_k[which(abs(t[time_index]-times_test_k)<h),]))
    CV_error_STKDE[k] <- CV_error_STKDE[k] + sum(evaluation_STKDE_train^2, na.rm = TRUE)/sum(!is.na(evaluation_STKDE_train)) -2 * sum(evaluation_STKDE_test, na.rm = TRUE) * cardinality_k
    
    if(sum(abs(t[time_index] - t_discrete) < 1e-06)){
      
      time_index_discrete <- which(abs(t[time_index] - t_discrete) < 1e-06)[1]
      
      # STKDE-discrete
      coeff <- c(solution_STKDE_discrete_train$dens[,,time_index_discrete])
      FEM_STKDE_discrete_train <- FEM(coeff, FEMbasis_STKDE_discrete_train)
      evaluation_STKDE_discrete_train <- eval.FEM(FEM_STKDE_discrete_train, locations = mesh.eval$nodes)
      evaluation_STKDE_discrete_train[idx_remove] <- NA
      evaluation_STKDE_discrete_train <- evaluation_STKDE_discrete_train/sum(evaluation_STKDE_discrete_train, na.rm = TRUE)
      
      FEM_STKDE_discrete_test <- FEM(evaluation_STKDE_discrete_train, FEMbasis.eval)
      evaluation_STKDE_discrete_test <- eval.FEM(FEM_STKDE_discrete_test, locations = locations_test_k[which(abs(t[time_index]-times_test_k)<h),])
      
      cardinality_k <- 1 / (sum(!is.na(evaluation_STKDE_discrete_test)) / nrow(locations_test_k[which(abs(t[time_index]-times_test_k)<h),]))
      CV_error_STKDE_discrete[k] <- CV_error_STKDE_discrete[k] + sum(evaluation_STKDE_discrete_train^2, na.rm = TRUE)/sum(!is.na(evaluation_STKDE_discrete_train)) -2 * sum(evaluation_STKDE_discrete_test, na.rm = TRUE) * cardinality_k
      
      # LGCP
      evaluation_LGCP_train <- intensity_LGCP_train[["grid"]][[time_index_discrete]]
      fftgr <- discreteWindow(solution_LGCP_train)
      # NA Values Outside the Spatial Domain (LGCP enlarges the window)
      for(i in 1:dim(evaluation_LGCP_train)[1]){
        for(j in 1:dim(evaluation_LGCP_train)[1]){
          if(!fftgr[i,j])
            evaluation_LGCP_train[i,j] <- NA
        }
      }
      FEM_LGCP_train <- FEM(c(evaluation_LGCP_train), FEMbasis_LGCP_train)
      evaluation_LGCP_train <- eval.FEM(FEM_LGCP_train, locations = mesh.eval$nodes)
      evaluation_LGCP_train[idx_remove] <- NA
      evaluation_LGCP_train <- evaluation_LGCP_train/sum(evaluation_LGCP_train, na.rm = TRUE)
      
      FEM_LGCP_test <- FEM(evaluation_LGCP_train, FEMbasis.eval)
      evaluation_LGCP_test <- eval.FEM(FEM_LGCP_test, locations = locations_test_k[which(abs(t[time_index]-times_test_k)<h),])
      
      cardinality_k <- 1 / (sum(!is.na(evaluation_LGCP_test)) / nrow(locations_test_k[which(abs(t[time_index]-times_test_k)<h),]))
      CV_error_LGCP[k] <- CV_error_LGCP[k] + sum(evaluation_LGCP_train^2, na.rm = TRUE)/sum(!is.na(evaluation_LGCP_train)) -2 * sum(evaluation_LGCP_test, na.rm = TRUE) * cardinality_k
      
    }  
    
    # STKDE-separable
    marginal_time_STKDE_train <- f_time_STKDE_separable_train$estimate[findInterval(t[time_index], f_time_STKDE_separable_train$eval.points)]
    solution_STKDE_separable_train <- marginal_time_STKDE_train * f_space_STKDE_separable_train
    solution_STKDE_separable_train[idx_remove] <- NA
    solution_STKDE_separable_train <- solution_STKDE_separable_train/sum(solution_STKDE_separable_train, na.rm = TRUE)
    
    FEM_STKDE_separable_test <- FEM(solution_STKDE_separable_train, FEMbasis.eval)
    solution_STKDE_separable_test <- eval.FEM(FEM_STKDE_separable_test, locations = locations_test_k[which(abs(t[time_index]-times_test_k)<h),])
    
    cardinality_k <- 1 / (sum(!is.na(solution_STKDE_separable_test)) / nrow(locations_test_k[which(abs(t[time_index]-times_test_k)<h),]))
    CV_error_STKDE_separable[k] <- CV_error_STKDE_separable[k] + sum(solution_STKDE_separable_train^2, na.rm = TRUE)/sum(!is.na(solution_STKDE_separable_train)) -2 * sum(solution_STKDE_separable_test, na.rm = TRUE) * cardinality_k
    
    # KNNSTDE-separable
    marginal_time_KNNSTDE_separable_train <- f_time_KNNSTDE_separable_train[findInterval(t[time_index], seq(-0.25, 10.25, length.out=25))]
    solution_KNNSTDE_separable_train <- marginal_time_KNNSTDE_separable_train * f_space_KNNSTDE_separable_train
    solution_KNNSTDE_separable_train[idx_remove] <- NA
    solution_KNNSTDE_separable_train <- solution_KNNSTDE_separable_train/sum(solution_KNNSTDE_separable_train, na.rm = TRUE)
    
    FEM_KNNSTDE_separable_test <- FEM(solution_KNNSTDE_separable_train, FEMbasis.eval)
    solution_KNNSTDE_separable_test <- eval.FEM(FEM_KNNSTDE_separable_test, locations = locations_test_k[which(abs(t[time_index]-times_test_k)<h),])
    
    cardinality_k <- 1 / (sum(!is.na(solution_KNNSTDE_separable_test)) / nrow(locations_test_k[which(abs(t[time_index]-times_test_k)<h),]))
    CV_error_KNNSTDE_separable[k] <- CV_error_KNNSTDE_separable[k] + sum(solution_KNNSTDE_separable_train^2, na.rm = TRUE)/sum(!is.na(solution_KNNSTDE_separable_train)) -2 * sum(solution_KNNSTDE_separable_test, na.rm = TRUE) * cardinality_k
    
  }

  CV_error_STDEPDE[k] <- CV_error_STDEPDE[k]  / length(t)
  CV_error_LGCP[k] <- CV_error_LGCP[k] / length(t_discrete)
  CV_error_INLA[k] <- CV_error_INLA[k] / length(t)
  CV_error_INLA_aggregated[k] <- CV_error_INLA_aggregated[k] / length(t)
  CV_error_STKDE[k] <- CV_error_STKDE[k] / length(t)
  CV_error_STKDE_discrete[k] <- CV_error_STKDE_discrete[k] / length(t_discrete)
  CV_error_STKDE_separable[k] <- CV_error_STKDE_separable[k] / length(t)
  CV_error_KNNSTDE_separable[k] <- CV_error_KNNSTDE_separable[k] / length(t)
  
  print(paste("CV iteration", k, "done."))
  
}

## BOXPLOT [FIGURE 10] ---------------------------------------------------------
{
  blue_polimi <- rgb(114, 143, 165, max = 255, alpha = 65, names = "blue_polimi25") # #728FA5
  
  pdf(paste0("pictures/app1_cv_errors.pdf"), family = "serif", width = 12, height = 5.3)
  
  #x11(width = 12, height = 5.3)
  par(mfrow = c(1,3), mai = c(1.75,0.35,0.5,0.15))
  plot.new()
  boxplot(CV_error_STDEPDE,CV_error_LGCP, CV_error_INLA, CV_error_INLA_aggregated,
          CV_error_STKDE, CV_error_STKDE_discrete, CV_error_STKDE_separable, CV_error_KNNSTDE_separable, asp = 1, xaxt = "n", yaxt = "n")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = blue_polimi)
  grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 1)
  par(new = TRUE)
  boxplot(CV_error_STDEPDE, CV_error_LGCP, CV_error_INLA, CV_error_INLA_aggregated,
          CV_error_STKDE, CV_error_STKDE_discrete, CV_error_STKDE_separable, CV_error_KNNSTDE_separable,
          names = c("STDE-PDE", "LGCP", "INLA", "INLA-agg", "STKDE",
                    "STKDE-dis", "STKDE-sep", "KNNSTDE-sep"), main = "CV error",
          col = brewer.pal(8, "YlGnBu"), cex.lab = 2, cex.axis = 2, cex.main = 2, las = 3)
  
  dev.off()
}
  