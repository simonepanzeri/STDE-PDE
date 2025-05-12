#'############################################################################'#
#'##################### SPATIO-TEMPORAL DENSITY ESTIMATION ###################'# 
#'############# SIMULATION 5 - MIXTURE OF 4 GAUSSIAN DISTRIBUTIONS ###########'#
#'############################################################################'#

## LIBRARIES AND FUNCTIONS -----------------------------------------------------
source("libraries_5.R")

## SPATIO-TEMPORAL DISCRETIZATION ----------------------------------------------
# Domain Boundaries
x_lim <- c(-6,6)
y_lim <- c(-6,6)
Ti <- 0
Tf <- 1

# Domain Area
domain_area <- diff(x_lim) * diff(y_lim)

# Spatial 2D Mesh over [-6,6]x[-6,6] for Estimation
Xbound <- seq(x_lim[1], x_lim[2], length.out = 11)
Ybound <- seq(y_lim[1], y_lim[2], length.out = 11)
grid_XY <- expand.grid(Xbound, Ybound)
Bounds <- grid_XY[(grid_XY$Var1 %in% x_lim) | (grid_XY$Var2 %in% y_lim), ]
mesh <- create.mesh.2D(nodes = Bounds, order = 1)
mesh <- refine.mesh.2D(mesh, maximum_area = 0.125, minimum_angle = 30)
FEMbasis <- create.FEM.basis(mesh)

x11()
plot(mesh)
dev.off()

# Fine Spatial 2D Mesh over [-6,6]x[-6,6] for Evaluation
# Xbound.eval <- seq(x_lim[1], x_lim[2], length.out = 21)
# Ybound.eval <- seq(y_lim[1], y_lim[2], length.out = 21)
# grid_XY.eval <- expand.grid(Xbound.eval, Ybound.eval)
# Bounds.eval <- grid_XY.eval[(grid_XY.eval$Var1 %in% x_lim) | (grid_XY.eval$Var2 %in% y_lim), ]
# mesh.eval <- create.mesh.2D(nodes = Bounds.eval, order = 1)
# mesh.eval <- refine.mesh.2D(mesh.eval, maximum_area = 0.05, minimum_angle = 30) # for gridsize = 32
# 
# vertices.eval <- mesh.eval$nodes
# FEMbasis.eval <- create.FEM.basis(mesh.eval)

# x11()
# plot(mesh.eval)
# dev.off()

# Grid Size for Estimation
n <- 32

# Fine Grid for Evaluation
X <- seq(x_lim[1], x_lim[2], length.out = n)
Y <- seq(y_lim[1], y_lim[2], length.out = n)
grid <- expand.grid(X, Y)
mesh.eval <- create.mesh.2D(grid)
FEMbasis.eval <- create.FEM.basis(mesh.eval)

x11()
plot(mesh.eval)
dev.off()

# Temporal 1D Mesh over [0,1]
mesh_time <- c(0, seq(from = 0.1, to = 0.9, by = 0.2), 1)

## ESTIMATION PROCEDURE --------------------------------------------------------
# Number of Processes
processes <- 30

# Sample Size
NN <- 5000
# alternative: NN <- c(1000, 2500, 5000, 10000)

# Compute the Integral of the True Density over the Domain of Interest
t0 <- proc.time()

integral <- 0.9608377
#integral <- compute_integral()

print(paste0("Computation of the integral done in ", round(as.numeric((proc.time()-t0)[3])), " seconds."))

for (N in NN) {
  
  # Containers for Errors
  mise_STDEPDE <- rep(0, processes)
  mise_LGCP <- rep(0, processes)
  mise_INLA <- rep(0, processes)
  mise_INLA_aggregated <- rep(0, processes)
  mise_STKDE <- rep(0, processes)
  mise_STKDE_discrete <- rep(0, processes)
 
  KLdist_STDEPDE <- rep(0, processes)
  KLdist_LGCP <- rep(0, processes)
  KLdist_INLA <- rep(0, processes)
  KLdist_INLA_aggregated <- rep(0, processes)
  KLdist_STKDE <- rep(0, processes)
  KLdist_STKDE_discrete <- rep(0, processes)
  
  Hdist_STDEPDE <- rep(0, processes)
  Hdist_LGCP <- rep(0, processes)
  Hdist_INLA <- rep(0, processes)
  Hdist_INLA_aggregated <- rep(0, processes)
  Hdist_STKDE <- rep(0, processes)
  Hdist_STKDE_discrete <- rep(0, processes)
  
  # Containers for CPU Times
  CPUtimes_STDEPDE <- rep(0, processes)
  CPUtimes_LGCP <- rep(0, processes)
  CPUtimes_INLA <- rep(0, processes)
  CPUtimes_INLA_aggregated <- rep(0, processes)
  CPUtimes_STKDE <- rep(0, processes)
  CPUtimes_STKDE_discrete <- rep(0, processes)
  
  for(proc in 1:processes){
    # Time for One Process
    t_proc <- proc.time()
    
    ### DATA -------------------------------------------------------------------
    # Generate the Data
    generate.data(N, proc) # run only once to generate the data
    
    # Read the Data
    data <- read.table(paste0("data/[-6,6]x[-6,6]/",N,"data_",proc,".txt"))
    
    # Locations
    locations <- data[,1:2]
    
    # Times
    times <- data[,3]
    
    # Discrete Times
    discrete_times <- (times-0.05) %/% 0.1
    discrete_times <- discrete_times + 1
    #discrete_times[discrete_times == -1] <- 0
    #discrete_times[discrete_times == 9] <- 8

    ## STDE-PDE: Spatio-Temporal Density Estimation with PDE Regulariz. --------
    t0 <- proc.time()
    
    # Smoothing Parameters
    lambda <- 0.075
    # alternative: lambda <- 10^seq(from = -3, to = -1, by = 1)
    
    lambda_time <- 0.01
    # alternative: lambda_time <- 10^seq(from = -3, to = -1, by = 1)
    
    # Solution
    # [If lambda and/or lambda_time are vectors, to select the best proposals by
    # 10-folds CV, please specify: preprocess_method = "RightCV"]
    solution_STDEPDE <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                    mesh_time = mesh_time, lambda = lambda, tol1 = 1e-7,
                                    lambda_time = lambda_time, fvec = NULL, heatStep = 0.1,
                                    heatIter = 10, print = F, nfolds = 10, nsimulations = 10000,
                                    step_method = "Wolfe_Method", direction_method = "L-BFGS5",
                                    preprocess_method = "NoCrossValidation") 
    
    FEMfunction_STDEPDE <- FEM.time(solution_STDEPDE$g, mesh_time, FEMbasis, FLAG_PARABOLIC = F)
    
    # CPU Time
    CPUtime <- proc.time() - t0
    CPUtimes_STDEPDE[proc] <- as.numeric(CPUtime[3])
    
    ## LGCP: Log-Gaussian Cox Process ------------------------------------------
    t0 <- proc.time()
    
    # Spatial Domain
    win <- owin(xrange = x_lim, yrange = y_lim)
    win[["bdry"]][[1]][["x"]] <- c(x_lim[1], x_lim[2], x_lim[2],x_lim[1])
    win[["bdry"]][[1]][["y"]] <- c(y_lim[1],y_lim[1], y_lim[2], y_lim[2])
    win[["bdry"]][[1]][["area"]] <- diff(x_lim)*diff(y_lim)
    #save(win, file = "window")
    
    # Spatio-Temporal Planar Point Process
    discrete_times <- as.integer(discrete_times)
    idx <- which(locations[,1] >= x_lim[1] & locations[,1] <= x_lim[2] & locations[,2] >= y_lim[1] & locations[,2] <= y_lim[2])
    xy_LGCP <- ppp(x = locations[idx,1], y = locations[idx,2], window = win)
    xyt_LGCP <- stppp(xy_LGCP, discrete_times[idx], tlim = c(as.integer(min(discrete_times[idx])), as.integer(max(discrete_times[idx]))))
    xyt_LGCP[["n"]] <- length(idx)
    xyt_LGCP <- integerise(xyt_LGCP)
    
    # Kernel Smoothed Intensity Function
    den <- density.ppp(xyt_LGCP, sigma = 1)
    sar <- spatialAtRisk(den)
    
    #mut <- constantInTime(xyt_LGCP)
    mut <- muEst(xyt_LGCP)
    
    #gin <- ginhomAverage(xyt_LGCP,spatial.intensity = sar,temporal.intensity = mut)
    #kin <- KinhomAverage(xyt_LGCP,spatial.intensity = sar,temporal.intensity = mut) 
    
    exceed <- exceedProbs(c(1.5,2,3))
    
    # Temporary Folder
    if(proc == 1){
      dir.create(file.path(getwd(), "LGCP"), showWarnings = FALSE)
    }
    tmpdr <- paste0(getwd(), "/LGCP/LGCP_", proc)
    #tmpdr <- tempdir()
    
    # Grid Size
    n_LGCP <- 16
    
    # Solution
    solution_LGCP <- lgcpPredict(xyt = xyt_LGCP, T = as.integer(10), laglength = as.integer(10),
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
    
    # CPU Time
    CPUtime <- proc.time() - t0
    CPUtimes_LGCP[proc] <- as.numeric(CPUtime[3])

    ## INLA: sPDE Spatio-Temporal Model ----------------------------------------
    t0 <- proc.time()
    
    # Temporal Discretization
    #t <- seq(from = 0, to = 1, by = 0.1)
    t <- mesh_time
    k <- length(t)
    mesh_time_INLA <- INLA::inla.mesh.1d(t)
    
    # Spatial Discretization
    m_temp <- mesh$nodes[which(mesh$nodes[,2] == y_lim[1]),]
    bdy <- m_temp[order(m_temp[,1], decreasing = FALSE),]
    m_temp <- mesh$nodes[which(mesh$nodes[,1] == x_lim[2]),]
    bdy <- rbind(bdy, m_temp[order(m_temp[,2], decreasing = FALSE),])
    m_temp <- mesh$nodes[which(mesh$nodes[,2] == y_lim[2]),]
    bdy <- rbind(bdy, m_temp[order(m_temp[,1], decreasing = TRUE),])
    m_temp <- mesh$nodes[which(mesh$nodes[,1] == x_lim[1]),]
    bdy <- rbind(bdy, m_temp[order(m_temp[,2], decreasing = TRUE),])
    
    domain_INLA <- SpatialPolygons(list(Polygons(list(Polygon(bdy)),"0")))
    mesh_INLA <- INLA::inla.mesh.2d(boundary = INLA::inla.sp2segment(domain_INLA),
                                    max.edge = 0.7, cutoff = 1)
    # alternative: pass locations as additional argument to build the mesh through R-INLA
    
    # Plot of the Spatial Mesh
    # x11()
    # plot(mesh_INLA)
    # points(bdy, col = "red", pch = 19)
    
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
    
    #gArea(domain_inla)
    
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
    
    res_INLA <- INLA::inla(form, family = "poisson", 
                           data = INLA::inla.stack.data(stk), E = expect,
                           control.predictor = list(A = INLA::inla.stack.A(stk)),
                           control.inla = list(strategy = "adaptive"))
    
    # NOTE: the exponential of the intercept plus the random effect at each
    #       space-time integration point is the relative risk at each of these
    #       points. This relative risk times the space-time volume will give the
    #       expected number of points (E(n)) at each one of these space-time
    #       locations. Summing over them will give a value that approaches the
    #       number of observations
    
    eta.at.integration.points <- res_INLA$summary.fix[1,1] + res_INLA$summary.ran$s$mean
    #c(n = n, "E(n)" = sum(st.vol * exp(eta.at.integration.points)))
    
    # Grid Size
    n_INLA <- n
    
    # Projection over a Grid for each Time Knot
    r0 <- diff(range(bdy[, 1])) / diff(range(bdy[, 2]))
    prj <- INLA::inla.mesh.projector(mesh_INLA, xlim = range(bdy[, 1]),
                                     ylim = range(bdy[, 2]), dims = c(n_INLA, n_INLA)) 
    ov <- over(SpatialPoints(prj$lattice$loc), domain_INLA)
    m.prj <- lapply(1:k, function(j) {
      #E[N]: r <- INLA::inla.mesh.project(prj, exp(res_INLA$summary.fix$mean + res_INLA$summary.ran$s$mean[1:m + (j - 1) * m]) * st.vol[1:m + (j - 1) * m] / N)
      r <- INLA::inla.mesh.project(prj, exp(res_INLA$summary.fix$mean + res_INLA$summary.ran$s$mean[1:m + (j - 1) * m]) / N)
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
    
    # CPU Time
    CPUtime <- proc.time() - t0
    CPUtimes_INLA[proc] <- as.numeric(CPUtime[3])

    ## INLA-aggregated: sPDE Spatio-Temporal Aggregated Model (Large Dataset) ----
    t0 <- proc.time()
    
    # Temporal Discretization
    #t <- seq(from = 0, to = 1, by = 0.1)
    t <- mesh_time
    k <- length(t)
    mesh_time_INLA <- INLA::inla.mesh.1d(t)

    # Spatial Discretization
    m_temp <- mesh$nodes[which(mesh$nodes[,2] == y_lim[1]),]
    bdy <- m_temp[order(m_temp[,1], decreasing = FALSE),]
    m_temp <- mesh$nodes[which(mesh$nodes[,1] == x_lim[2]),]
    bdy <- rbind(bdy, m_temp[order(m_temp[,2], decreasing = FALSE),])
    m_temp <- mesh$nodes[which(mesh$nodes[,2] == y_lim[2]),]
    bdy <- rbind(bdy, m_temp[order(m_temp[,1], decreasing = TRUE),])
    m_temp <- mesh$nodes[which(mesh$nodes[,1] == x_lim[1]),]
    bdy <- rbind(bdy, m_temp[order(m_temp[,2], decreasing = TRUE),])
    
    domain_INLA <- SpatialPolygons(list(Polygons(list(Polygon(bdy)),"0")))
    mesh_INLA <- INLA::inla.mesh.2d(boundary = INLA::inla.sp2segment(domain_INLA),
                                    max.edge = 0.7, cutoff = 1)
    # alternative: pass locations as additional argument to build the mesh through R-INLA
    
    # Spatio-Temporal Planar Point Pattern
    xyt_INLA <- stppp(xy_LGCP, times[idx], tlim = c(0, 1))
    
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
    res_INLA_aggregated <- INLA::inla(formula, family = "poisson",
                                      data = INLA::inla.stack.data(stk), E = exposure, 
                                      control.predictor = list(A = INLA::inla.stack.A(stk)),
                                      control.inla = list(strategy = "adaptive"))
    
    # Value of mu and Summary for the Intercept
    # cbind(True = mu, res_INLA_aggregated$summary.fixed[, 1:6])
    
    # Summary for the Hyperparameters
    # cbind(True = c(range, sigma, rho), res_INLA_aggregated$summary.hyperpar[, c(1, 2, 3, 5)])
    
    # Grid Size
    n_INLA_aggregated <- n
    
    # Projection over a Grid for each Time Knot
    r0 <- diff(range(bdy[, 1])) / diff(range(bdy[, 2]))
    prj <- INLA::inla.mesh.projector(mesh_INLA, xlim = bbox(domain_INLA)[1, ], 
                                     ylim = bbox(domain_INLA)[2, ],
                                     dims = c(n_INLA_aggregated, n_INLA_aggregated))
    g.no.in <- is.na(over(SpatialPoints(prj$lattice$loc), domain_INLA))
    t.mean <- lapply(1:k, function(j) {
      # E[N]: z.j <- exp(res_INLA_aggregated$summary.fix$mean + res_INLA_aggregated$summary.ran$s$mean[idx_INLA$s.group == j]) * w.areas * w.t / N
      z.j <- exp(res_INLA_aggregated$summary.fix$mean + res_INLA_aggregated$summary.ran$s$mean[idx_INLA$s.group == j]) / N
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
    
    CPUtime <- proc.time() - t0
    CPUtimes_INLA_aggregated[proc] <- as.numeric(CPUtime[3])
    
    ## STKDE: Spatio-Temporal Kernel Density Estimation ------------------------
    t0 <- proc.time()
    
    # Spatial Domain
    win <- owin(xrange = x_lim, yrange = y_lim)
    win[["bdry"]][[1]][["x"]] <- c(x_lim[1], x_lim[2], x_lim[2],x_lim[1])
    win[["bdry"]][[1]][["y"]] <- c(y_lim[1],y_lim[1], y_lim[2], y_lim[2])
    win[["bdry"]][[1]][["area"]] <- diff(x_lim)*diff(y_lim)
    
    # Spatio-Temporal Planar Point Pattern
    idx <- which(locations[,1] >= x_lim[1] & locations[,1] <= x_lim[2] & locations[,2] >= y_lim[1] & locations[,2] <= y_lim[2])
    xy_STKDE <- ppp(x = locations[idx,1], y = locations[idx,2], window = win)
    
    # Grid Size
    n_STKDE <- n
    
    # Time Refinement
    t <- seq(from = 0, to = 1, by = 0.1)
    hh <- 1/(length(t)-1)
    
    # Solution
    solution_STKDE <- spattemp.density(pp = xy_STKDE, h = NULL, tt = times[idx],
                                       lambda = NULL, tlim = c(min(t)-hh/2, max(t)+hh/2),
                                       sedge = "none", tedge = "none", 
                                       sres = n_STKDE, tres = length(t))
    
    x_STKDE <- seq(x_lim[1], x_lim[2], length.out = n_STKDE)
    y_STKDE <- seq(y_lim[1], y_lim[2], length.out = n_STKDE)
    grid_STKDE <- expand.grid(x_STKDE, y_STKDE)
    FEMbasis_STKDE <- create.FEM.basis(create.mesh.2D(grid_STKDE))
    
    CPUtime <- proc.time() - t0
    CPUtimes_STKDE[proc] <- as.numeric(CPUtime[3])

    ## STKDE-discrete: Spatio-Temporal Kernel Density Estimation ------------------------
    t0 <- proc.time()
    
    # Grid Size
    n_STKDE_discrete <- n
    
    # Solution
    x11(width = 25, height = 14)
    idx <- which(abs(locations[,1]) <= x_lim[2] & abs(locations[,2]) <= y_lim[2])
    solution_STKDE_discrete <- stkde(xlong = locations[idx,1], ylat = locations[idx,2], ztime = discrete_times[idx],
                                     xgrids = n_STKDE_discrete, ygrids = n_STKDE_discrete, breaks = 0.05, alpha = 0.05, nrowspar = 5, bwmethod = "normal-reference") # "cv.ml"
    x_STKDE_discrete <- seq(min(locations[idx,1]), max(locations[idx,1]), length.out = n_STKDE_discrete)
    y_STKDE_discrete <- seq(min(locations[idx,2]), max(locations[idx,2]), length.out = n_STKDE_discrete)
    grid_STKDE_discrete <- expand.grid(x_STKDE_discrete, y_STKDE_discrete)
    FEMbasis_STKDE_discrete <- create.FEM.basis(create.mesh.2D(grid_STKDE_discrete))
    dev.off()
    
    CPUtime <- proc.time() - t0
    CPUtimes_STKDE_discrete[proc] <- as.numeric(CPUtime[3])
  
    ## ERRORS ------------------------------------------------------------------
    # Time instants at which the solution is evaluated
    #t <- seq(from = 0, to = 1, by = 0.1)
    t <- mesh_time
    
    # Mean Solutions over the Processes
    if (proc == 1){
      
      mean_sol_STDEPDE <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
      mean_sol_LGCP <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
      mean_sol_INLA <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
      mean_sol_INLA_aggregated <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
      mean_sol_STKDE <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
      mean_sol_STKDE_discrete <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
  
      true_density <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
        
    }

    for (time_index in 1:length(t)) {
      # True Density
      if(proc == 1){
        
        true_density[,time_index] <- dens.func(mesh.eval$nodes, t[time_index])/integral
        true_density[,time_index] <- true_density[,time_index]/sum(true_density[,time_index], na.rm = TRUE)
      
      }
       
      # STDE-PDE
      evaluation_STDEPDE <- eval.FEM.time(FEM.time = FEMfunction_STDEPDE, locations = mesh.eval$nodes, time.instants = t[time_index])
      evaluation_STDEPDE <- exp(evaluation_STDEPDE)
      evaluation_STDEPDE <- evaluation_STDEPDE/sum(evaluation_STDEPDE, na.rm = TRUE)
      
      # INLA
      #coeff <- c(exp(solution_INLA[[time_index]])/N)
      coeff <- c(solution_INLA[[time_index]])
      FEM_INLA <- FEM(coeff, FEMbasis_INLA)
      evaluation_INLA <- eval.FEM(FEM_INLA, mesh.eval$nodes)
      evaluation_INLA <- evaluation_INLA/sum(evaluation_INLA, na.rm = TRUE)
      
      # INLA-aggregated
      #coeff <- c(exp(solution_INLA_aggregated[[time_index]])/N)
      coeff <- c(solution_INLA_aggregated[[time_index]])
      FEM_INLA_aggregated <- FEM(coeff, FEMbasis_INLA_aggregated)
      evaluation_INLA_aggregated <- eval.FEM(FEM_INLA_aggregated, mesh.eval$nodes)
      evaluation_INLA_aggregated <- evaluation_INLA_aggregated/sum(evaluation_INLA_aggregated, na.rm = TRUE)
      
      # STKDE
      time_index_STKDE <- which(abs(t[time_index] - seq(from = 0, to = 1, by = 0.1)) < 1e-6)
      FEM_STKDE <- FEM(c(t(solution_STKDE$z[[time_index_STKDE]]$v)), FEMbasis_STKDE)
      evaluation_STKDE <- eval.FEM(FEM_STKDE, locations = mesh.eval$nodes)
      evaluation_STKDE <- evaluation_STKDE/sum(evaluation_STKDE, na.rm = TRUE)
      
      # STKDE-discrete
      time_index_STKDE_discrete <- which(abs(t[time_index] - seq(from = 0, to = 1, by = 0.1)) < 1e-6)
      FEM_STKDE_discrete <- FEM(c(solution_STKDE_discrete$dens[,,time_index_STKDE_discrete]), FEMbasis_STKDE_discrete)
      evaluation_STKDE_discrete <- eval.FEM(FEM_STKDE_discrete, locations = mesh.eval$nodes)
      evaluation_STKDE_discrete <- evaluation_STKDE_discrete/sum(evaluation_STKDE_discrete, na.rm = TRUE)
        
      # LGCP
      time_index_LGCP <- which(abs(t[time_index] - seq(from = 0, to = 1, by = 0.1)) < 1e-6)
      evaluation_LGCP <- intensity_LGCP[["grid"]][[time_index_LGCP]]
      fftgr <- discreteWindow(solution_LGCP)
      # NA Values Outside the Spatial Domain (LGCP enlarges the window)
      for(i in 1:dim(evaluation_LGCP)[1]){
        for(j in 1:dim(evaluation_LGCP)[1]){
          if(!fftgr[i,j])
            evaluation_LGCP[i,j] <- NA
        }
      }
      FEM_LGCP <- FEM(c(evaluation_LGCP)/N, FEMbasis_LGCP)
      evaluation_LGCP <- eval.FEM(FEM_LGCP, locations = mesh.eval$nodes)
      evaluation_LGCP[which(is.na(evaluation_LGCP))] <- min(evaluation_LGCP, na.rm = T)
      evaluation_LGCP <- evaluation_LGCP/sum(evaluation_LGCP, na.rm = TRUE)
      
      # Estimates and Errors Updates
      
      # STDE-PDE
      mise_STDEPDE[proc] <- mise_STDEPDE[proc] + 1*domain_area*mean((evaluation_STDEPDE - true_density[,time_index])^2, na.rm = TRUE)/length(t)
      for(p in 1:nrow(mesh.eval$nodes)){
        if(!is.na(evaluation_STDEPDE[p]) & evaluation_STDEPDE[p] != 0 & true_density[p,time_index] != 0){
          KLdist_STDEPDE[proc] <- KLdist_STDEPDE[proc] + (evaluation_STDEPDE[p]*log(evaluation_STDEPDE[p]/true_density[p,time_index]))/length(t)
        }
      }
      Hdist_STDEPDE[proc] <- Hdist_STDEPDE[proc] + domain_area*hellinger.distance(evaluation_STDEPDE, true_density[,time_index])/length(t)
        
      mean_sol_STDEPDE[,time_index] <- mapply(sum, mean_sol_STDEPDE[,time_index], evaluation_STDEPDE, na.rm = TRUE)
      
      # STKDE
      mise_STKDE[proc] <- mise_STKDE[proc] + 1*domain_area*mean((evaluation_STKDE - true_density[,time_index])^2, na.rm = TRUE)/length(t)
      for(p in 1:nrow(mesh.eval$nodes)){
        if(!is.na(evaluation_STKDE[p]) & evaluation_STKDE[p] != 0 & true_density[p,time_index] != 0){
          KLdist_STKDE[proc] <- KLdist_STKDE[proc] + (evaluation_STKDE[p]*log(evaluation_STKDE[p]/true_density[p,time_index]))/length(t)
        }
      }
      Hdist_STKDE[proc] <- Hdist_STKDE[proc] + domain_area*hellinger.distance(evaluation_STKDE, true_density[,time_index])/length(t)
      
      mean_sol_STKDE[,time_index] <- mapply(sum, mean_sol_STKDE[,time_index], evaluation_STKDE, na.rm = TRUE)
      
      # INLA
      mise_INLA[proc] <- mise_INLA[proc] + 1*domain_area*mean((evaluation_INLA - true_density[,time_index])^2, na.rm = TRUE)/length(t)
      for(p in 1:nrow(mesh.eval$nodes)){
        if(!is.na(evaluation_INLA[p]) & evaluation_INLA[p] != 0 & true_density[p,time_index] != 0){
          KLdist_INLA[proc] <- KLdist_INLA[proc] + (evaluation_INLA[p]*log(evaluation_INLA[p]/true_density[p,time_index]))/length(t)
        }
      }
      Hdist_INLA[proc] <- Hdist_INLA[proc] + domain_area*hellinger.distance(evaluation_INLA, true_density[,time_index])/length(t)
        
      mean_sol_INLA[,time_index] <- mapply(sum, mean_sol_INLA[,time_index], evaluation_INLA, na.rm = TRUE)
        
      # INLA-aggregated
      mise_INLA_aggregated[proc] <- mise_INLA_aggregated[proc] + 1*domain_area*mean((evaluation_INLA_aggregated - true_density[,time_index])^2, na.rm = TRUE)/length(t)
      for(p in 1:nrow(mesh.eval$nodes)){
        if(!is.na(evaluation_INLA_aggregated[p]) & evaluation_INLA_aggregated[p] != 0 & true_density[p,time_index] != 0){
          KLdist_INLA_aggregated[proc] <- KLdist_INLA_aggregated[proc] + (evaluation_INLA_aggregated[p]*log(evaluation_INLA_aggregated[p]/true_density[p,time_index]))/length(t)
        }
      }
      Hdist_INLA_aggregated[proc] <- Hdist_INLA_aggregated[proc] + domain_area*hellinger.distance(evaluation_INLA_aggregated, true_density[,time_index])/length(t)
        
      mean_sol_INLA_aggregated[,time_index] <- mapply(sum, mean_sol_INLA_aggregated[,time_index], evaluation_INLA_aggregated, na.rm = TRUE)
      
      # STKDE-discrete
      mise_STKDE_discrete[proc] <- mise_STKDE_discrete[proc] + 1*domain_area*mean((evaluation_STKDE_discrete - true_density[,time_index])^2, na.rm = TRUE)/length(t)
      for(p in 1:nrow(mesh.eval$nodes)){
        if(!is.na(evaluation_STKDE_discrete[p]) & evaluation_STKDE_discrete[p] != 0 & true_density[p,time_index] != 0){
          KLdist_STKDE_discrete[proc] <- KLdist_STKDE_discrete[proc] + (evaluation_STKDE_discrete[p]*log(evaluation_STKDE_discrete[p]/true_density[p,time_index]))/length(t)
        }
      }
      Hdist_STKDE_discrete[proc] <- Hdist_STKDE_discrete[proc] + domain_area*hellinger.distance(evaluation_STKDE_discrete, true_density[,time_index])/length(t)
        
      mean_sol_STKDE_discrete[,time_index] <- mapply(sum, mean_sol_STKDE_discrete[,time_index], evaluation_STKDE_discrete, na.rm = TRUE)
        
      # LGCP
      mise_LGCP[proc] <- mise_LGCP[proc] + 1*domain_area*mean((evaluation_LGCP - true_density[,time_index])^2, na.rm = TRUE)/length(unique(discrete_times))
      for(p in 1:nrow(mesh.eval$nodes)){
        if(!is.na(evaluation_LGCP[p]) & evaluation_LGCP[p] != 0 & true_density[p,time_index] != 0){
          KLdist_LGCP[proc] <- KLdist_LGCP[proc] + (evaluation_LGCP[p]*log(evaluation_LGCP[p]/true_density[p,time_index]))/length(unique(discrete_times))
        }
      }
      Hdist_LGCP[proc] <- Hdist_LGCP[proc] + domain_area*hellinger.distance(evaluation_LGCP, true_density[,time_index])/length(unique(discrete_times))
        
      mean_sol_LGCP[,time_index] <- mapply(sum, mean_sol_LGCP[,time_index], evaluation_LGCP, na.rm = TRUE)
      
      print(paste("Errors computed for time:", t[time_index]))
    }

    print(paste("Process", proc, "done in", round(as.numeric(proc.time() - t_proc)[3]), "seconds."))
    unlink(paste0(tmpdr,"/*"))
    
  }
  
  mean_sol_STDEPDE <- mean_sol_STDEPDE/processes
  mean_sol_LGCP <- mean_sol_LGCP/processes
  mean_sol_INLA <- mean_sol_INLA/processes
  mean_sol_INLA_aggregated <- mean_sol_INLA_aggregated/processes
  mean_sol_STKDE <- mean_sol_STKDE/processes
  mean_sol_STKDE_discrete <- mean_sol_STKDE_discrete/processes

  #'##########################################################################'#
  #'##########################################################################'#

  ## EXPORT RESULTS ------------------------------------------------------------
  # Export estimates and errors
  dir.create(file.path(getwd(), "output"), showWarnings = FALSE)
  
  estimates <- list(mean_sol_STDEPDE, mean_sol_LGCP, mean_sol_INLA,
                    mean_sol_INLA_aggregated, mean_sol_STKDE, mean_sol_STKDE_discrete)
  save(estimates, file = paste0("output/estimates_", N, "data.rda"))
  
  mise <- list(mise_STDEPDE, mise_LGCP, mise_INLA, mise_INLA_aggregated,
               mise_STKDE, mise_STKDE_discrete)
  save(mise, file = paste0("output/mise_", N, "data.rda"))
  
  KLdist <- list(KLdist_STDEPDE, KLdist_LGCP, KLdist_INLA,
                  KLdist_INLA_aggregated, KLdist_STKDE, KLdist_STKDE_discrete)
  save(KLdist, file = paste0("output/KLdist_", N, "data.rda"))
  
  Hdist <- list(Hdist_STDEPDE, Hdist_LGCP, Hdist_INLA, Hdist_INLA_aggregated,
               Hdist_STKDE, Hdist_STKDE_discrete)
  save(Hdist, file = paste0("output/Hdist_", N, "data.rda"))
  
  CPUtimes <- list(CPUtimes_STDEPDE, CPUtimes_LGCP, CPUtimes_INLA,
                   CPUtimes_INLA_aggregated, CPUtimes_STKDE, CPUtimes_STKDE_discrete)
  save(CPUtimes, file = paste0("output/CPUtimes_", N, "data.rda"))

  base::save.image(paste0("output/workspace_", N, "data.RData"))
  
  # Remove the content of the LGCP directory
  unlink(paste0("LGCP/*"), recursive = T)
  
  # Alternative: Remove LGCP files from directory
  # for(proc in 1:processes){
  #   unlink(paste0("LGCP/LGCP_",proc,"/*"))
  # }

  #'##########################################################################'#
  #'##########################################################################'#
  
  ## VISUALIZATION -------------------------------------------------------------
  
  ### ESTIMATES [FIGURE G.15] -----------------------------------------------------
  # Import estimates and errors
  # load(file = paste0("output/estimates_", N, "data.rda"))
  # mean_sol_STDEPDE <- estimates[[1]]
  # mean_sol_LGCP <- estimates[[2]]
  # mean_sol_INLA <- estimates[[3]]
  # mean_sol_INLA_aggregated <- estimates[[4]]
  # mean_sol_STKDE <- estimates[[5]]
  # mean_sol_STKDE_discrete <- estimates[[6]]
    
  # Time Instants for Evaluation
  #t <- seq(from = 0, to = 1, by = 0.1)
  t <- mesh_time
  h <- 0.1/2
    
  # Plots
  m <- min(min(true_density, na.rm = TRUE),
           min(mean_sol_STDEPDE, na.rm = TRUE),
           min(mean_sol_LGCP[,2:5], na.rm = TRUE),
           min(mean_sol_INLA, na.rm = TRUE),
           min(mean_sol_INLA_aggregated, na.rm = TRUE),
           min(mean_sol_STKDE, na.rm = TRUE),
           min(mean_sol_STKDE_discrete, na.rm = TRUE), na.rm = TRUE)
  
  M <- max(max(true_density, na.rm = TRUE),
           max(mean_sol_STDEPDE, na.rm = TRUE),
           max(mean_sol_LGCP[,2:5], na.rm = TRUE),
           max(mean_sol_INLA, na.rm = TRUE),
           max(mean_sol_INLA_aggregated, na.rm = TRUE),
           max(mean_sol_STKDE, na.rm = TRUE),
           max(mean_sol_STKDE_discrete, na.rm = TRUE), na.rm = TRUE)
            
  dir.create(file.path(getwd(), "pictures"), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/", N, "_data")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/", N, "_data/sample")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/", N, "_data/true_density")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/", N, "_data/STDEPDE")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/", N, "_data/LGCP")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/", N, "_data/INLA")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/", N, "_data/INLA_aggregated")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/", N, "_data/STKDE")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/", N, "_data/STKDE_discrete")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/", N, "_data/time_bar")), showWarnings = FALSE)
  
  for(time_index in 1:length(t)){
    t_img <- proc.time()
    
    # First Plot: Sample
    name <- paste0("pictures/", N, "_data/sample/sim5_sample_", N, "data_", t[time_index])
    plot.sample(locations[abs(times-t[time_index])<h,], name)
    
    # Second Plot: True Density at t[time_index]
    name <- paste0("pictures/", N, "_data/true_density/sim5_true_density_", N, "data_", t[time_index])
    plot.density(mesh.eval$nodes[,1], mesh.eval$nodes[,2], true_density[,time_index],
                 max_range = M,
                 filename = name, colorscale = "Jet")
    #plot.interactive.density(FEM(true_density[,time_index], FEMbasis.eval))
    
    # Third Plot: STDE-PDE Estimated Density at t[time_index]
    name <- paste0("pictures/", N, "_data/STDEPDE/sim5_STDEPDE_", N, "data_", t[time_index])
    plot.density(mesh.eval$nodes[,1], mesh.eval$nodes[,2], mean_sol_STDEPDE[,time_index],
                 max_range = M,
                 filename = name, colorscale = "Jet")
    #plot.interactive.density(FEM(mean_sol_STDEPDE[,time_index], FEMbasis.eval))
    
    # Fourth Plot: STKDE Estimated Density at t[time_index]
    name <- paste0("pictures/", N, "_data/STKDE/sim5_STKDE_", N, "data_", t[time_index])
    plot.density(mesh.eval$nodes[,1], mesh.eval$nodes[,2], mean_sol_STKDE[,time_index],
                 max_range = M,
                 filename = name, colorscale = "Jet")
    #plot.interactive.density(FEM(mean_sol_STKDE[,time_index], FEMbasis.eval))
  
    # Fifth Plot: STKDE-discrete Estimated Density at t[time_index]
    name <- paste0("pictures/", N, "_data/STKDE_discrete/sim5_STKDE_discrete_", N, "data_", t[time_index])
    plot.density(mesh.eval$nodes[,1], mesh.eval$nodes[,2], mean_sol_STKDE_discrete[,time_index],
                 max_range = M,
                 filename = name, colorscale = "Jet")
    #plot.interactive.density(FEM(mean_sol_STKDE_discrete[,time_index], FEMbasis.eval))
      
    # Sixth Plot: LGCP Estimated Density at t[time_index]
    name <- paste0("pictures/", N, "_data/LGCP/sim5_LGCP_", N, "data_", t[time_index])
    plot.density(mesh.eval$nodes[,1], mesh.eval$nodes[,2], mean_sol_LGCP[,time_index],
                 max_range = M,
                 filename = name, colorscale = "Jet")
    #plot.interactive.density(FEM(mean_sol_LGCP[,time_index], FEMbasis.eval))
      
    # Seventh Plot: INLA Estimated Density at t[time_index]
    name <- paste0("pictures/", N, "_data/INLA/sim5_INLA_", N, "data_", t[time_index])
    plot.density(mesh.eval$nodes[,1], mesh.eval$nodes[,2], mean_sol_INLA[,time_index],
                 max_range = M,
                 filename = name, colorscale = "Jet")
    #plot.interactive.density(FEM(mean_sol_INLA[,time_index], FEMbasis.eval))
      
    # Eight Plot: INLA-aggregated Estimated Density at t[time_index]
    name <- paste0("pictures/", N, "_data/INLA_aggregated/sim5_INLA_aggregated_", N, "data_", t[time_index])
    plot.density(mesh.eval$nodes[,1], mesh.eval$nodes[,2], mean_sol_INLA_aggregated[,time_index],
                 max_range = M,
                 filename = name, colorscale = "Jet")
    #plot.interactive.density(FEM(mean_sol_INLA_aggregated[,time_index], FEMbasis.eval))
      
    print(paste("Images for t =", t[time_index], "done in", round(as.numeric(proc.time() - t_img)[3]), "seconds."))
    
  }
  
  ### LEGEND [FIGURE G.15] ------------------------------------------------------
  pdf(paste0("pictures/", N, "_data/new_sim5_legend.pdf"), family = "serif", width = 11, height = 3)
  #x11(width = 11, height = 3)
  par(mai = c(1,0.75,0,0))
  plot(c(0, 112.5), c(0, 15), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame.plot = F)
  gradient.rect(0, 0, 100, 5, col = jet.col(1000), border = "black")
  # segments(seq(0,100,by=100/175), 5*rep(0,length=176),
  #          seq(0,100,by=100/175), 5*rep(1,length=176),
  #          col = "white", lwd = 0.05)
  # segments(c(0, 100, 100, 0), c(0, 0, 5, 5), c(100, 100, 0, 0), c(0, 5, 5, 0),
  #          col = "black", lwd = 2)
  axis(1, at = c(0,20,40,60,80,100), labels = c("0", "1.0", "2.0", "3.0", "4.0", "5.0"),
       lwd.ticks = 2, cex.axis = 2, lwd = 2)
  text(107,2, TeX("$\\times 10^{-3}$"), cex = 2)
  
  dev.off()
  
  print(paste0("M = ", round(M, 4)))
  
  ### TIME BAR ---------------------------------------------------------------
  for(time_index in 1:length(t)) {
    string <- paste0(paste0("pictures/", N, "_data/time_bar/sim5_time_bar_",time_index,".pdf"))
    pdf(file=string, family = "serif", width = 7, height = 1.5)
    #x11(width = 7, height = 1.5)
    par(mar = c(0.01, 2, 0.01, 2))
    bar <- rbind(c(t[time_index],0),c(t[time_index],0.075))
    plot(bar, type = "l", bty = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "",
         pch = 3, lwd = 5, xlim = c(0,1), asp = 1)
    mtext("Time", line = -1.5, font = 2, cex = 2)
    axis(side = 1, c(0, 0.25, 0.5, 0.75, 1), pos = 0,
         labels = c("0", "0.25", "0.50", "0.75", "1"), cex.axis = 2)
    dev.off()
  }
  
  ### ERRORS [FIGURE G.16] --------------------------------------------------------
  {
    blue_polimi <- rgb(114, 143, 165, max = 255, alpha = 65, names = "blue_polimi25") # #728FA5
    pdf(paste0("pictures/", N, "_data/sim5_errors.pdf"), family = "serif", width = 12, height = 5.3)
    #x11(width = 12, height = 5.3)
    par(mfrow = c(1,3), mai = c(1.75,0.35,0.5,0.15))
    boxplot(mise_STDEPDE, mise_INLA, mise_LGCP, mise_STKDE, asp = 1, xaxt = "n", yaxt = "n")
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = blue_polimi)
    grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 1)
    par(new = TRUE)
    boxplot(mise_STDEPDE, mise_INLA, mise_LGCP, mise_STKDE,
            names = c("STDE-PDE", "INLA", "LGCP", "STKDE"), main = TeX("$err_{L^2}$", bold = T),
            col = brewer.pal(5, "YlGnBu")[1:4], cex.lab = 2, cex.axis = 2, cex.main = 2, las = 3)
    
    boxplot(KLdist_STDEPDE, KLdist_INLA, KLdist_LGCP, KLdist_STKDE, asp = 1, xaxt = "n", yaxt = "n")
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = blue_polimi)
    grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 1)
    par(new = TRUE)
    boxplot(KLdist_STDEPDE, KLdist_INLA, KLdist_LGCP, KLdist_STKDE,
            names = c("STDE-PDE", "INLA", "LGCP", "STKDE"), main = "KL",
            col = brewer.pal(5, "YlGnBu")[1:4], cex.lab = 2, cex.axis = 2, cex.main = 2, las = 3)
    
    boxplot(Hdist_STDEPDE, Hdist_INLA, Hdist_LGCP, Hdist_STKDE, asp = 1, xaxt = "n", yaxt = "n")
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = blue_polimi)
    grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 1)
    par(new = TRUE)
    boxplot(Hdist_STDEPDE, Hdist_INLA, Hdist_LGCP, Hdist_STKDE,
            names = c("STDE-PDE", "INLA", "LGCP", "STKDE"), main = "H",
            col = brewer.pal(5, "YlGnBu")[1:4], cex.lab = 2, cex.axis = 2, cex.main = 2, las = 3)
    dev.off()
  } 
  
  ### CPU TIMES VISUALIZATION --------------------------------------------------
  pdf(paste0("pictures/", N, "_data/sim5_CPUtimes.pdf"), family = "serif", width = 12, height = 5.3)
  #x11(width = 12, height = 5.3)
  par(mfrow = c(1,3), mai = c(1.75,0.35,0.5,0.15))
  plot.new()
  boxplot(CPUtimes_STDEPDE, CPUtimes_INLA, CPUtimes_LGCP, CPUtimes_STKDE,
          asp = 1, xaxt = "n", yaxt = "n")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = blue_polimi)
  grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 1)
  par(new = TRUE)
  boxplot(CPUtimes_STDEPDE, CPUtimes_INLA, CPUtimes_LGCP, CPUtimes_STKDE,
          names = c("STDE-PDE", "INLA", "LGCP", "STKDE"), main = "CPU TIMES [seconds]",
          col = brewer.pal(5, "YlGnBu")[1:4], cex.lab = 2, cex.axis = 2, cex.main = 2, las = 3)
  dev.off()

  ### MESH VISUALIZATION -------------------------------------------------------
  name <- paste0("pictures/sim5_mesh")
  plot.mesh(mesh, name)
  
  name <- paste0("pictures/sim5_mesh_eval")
  plot.mesh(mesh.eval, name)
  
}


## SMOOTHING PARAMETERS SELECTION ----------------------------------------------
cat("Selecting smoothing parameters...\n")

lambda <- c(1, 1e-1, 1e-2, 1e-3, 1e-4)
lambda <- c(0.1, 0.075, 0.05, 0.025, 0.01)
lambda_time <- 0.01 # try other values of lambda_time

for(i in 1:length(lambda)){
  
  # Solution
  # [If lambda and/or lambda_time are vectors, to select the best proposals by
  # 10-folds CV, please specify: preprocess_method = "RightCV"]
  invisible(capture.output(solution_STDEPDE <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                           mesh_time = mesh_time, lambda = lambda[i], tol1 = 1e-7,
                                                           lambda_time = lambda_time, fvec = NULL, heatStep = 0.1,
                                                           heatIter = 10, print = TRUE, nfolds = 10, nsimulations = 10000,
                                                           step_method = "Wolfe_Method", direction_method = "L-BFGS5",
                                                           preprocess_method = "NoCrossValidation")))
  
  FEMfunction_STDEPDE <- FEM.time(solution_STDEPDE$g, mesh_time, FEMbasis,
                                  FLAG_PARABOLIC = FALSE)
  
  # Time instants at which the solution is evaluated
  t <- seq(from = Ti, to = Tf, length.out = 31)
  
  mise_STDEPDE <- 0
  mean_sol_STDEPDE <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
  true_density <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
  
  for(time_index in 1:length(t)) {
    
    true_density[,time_index] <- dens.func(mesh.eval$nodes, t[time_index])/integral
    true_density[,time_index] <- true_density[,time_index]/sum(true_density[,time_index], na.rm = TRUE)
    
    # STDE-PDE
    evaluation_STDEPDE <- eval.FEM.time(FEM.time = FEMfunction_STDEPDE, locations = mesh.eval$nodes, time.instants = t[time_index])
    evaluation_STDEPDE <- exp(evaluation_STDEPDE)
    evaluation_STDEPDE <- evaluation_STDEPDE/sum(evaluation_STDEPDE, na.rm = TRUE)
    
    # Estimates and Errors Updates
    
    # STDE-PDE
    mise_STDEPDE <- mise_STDEPDE + 1*domain_area*mean((evaluation_STDEPDE - true_density[,time_index])^2, na.rm = TRUE)/length(t)
    mean_sol_STDEPDE[,time_index] <- mapply(sum, mean_sol_STDEPDE[,time_index], evaluation_STDEPDE, na.rm = TRUE)
    
  }
  
  cat(paste0("lambda = ", lambda[i], "; L2 error = ", mise_STDEPDE, "\n"))

}



