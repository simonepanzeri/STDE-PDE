#'############################################################################'#
#'##################### SPATIO-TEMPORAL DENSITY ESTIMATION ###################'# 
#'############# SIMULATION 1 - MIXTURE OF 4 GAUSSIAN DISTRIBUTIONS ###########'#
#'############################################################################'#

## LIBRARIES AND FUNCTIONS -----------------------------------------------------
source("libraries_1.R")

## SPATIO-TEMPORAL DISCRETIZATION ----------------------------------------------
# Domain Boundaries
x_lim <- c(-6,6)
y_lim <- c(-6,6)

# Domain Area
domain_area <- diff(x_lim) * diff(y_lim)

# Spatial 2D Mesh over [-6,6]x[-6,6] for Estimation
Xbound <- seq(x_lim[1], x_lim[2], length.out = 11)
Ybound <- seq(y_lim[1], y_lim[2], length.out = 11)
grid_XY <- expand.grid(Xbound, Ybound)
Bounds <- grid_XY[(grid_XY$Var1 %in% x_lim) | (grid_XY$Var2 %in% y_lim), ]
mesh <- create.mesh.2D(nodes = Bounds, order = 1)
mesh <- refine.mesh.2D(mesh, maximum_area = 0.15, minimum_angle = 30)
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

for (N in NN) {
  
  # Containers for CPU Times
  CPUtimes_STDEPDE <- rep(0, processes)
  CPUtimes_INLA <- rep(0, processes)
  CPUtimes_INLA_aggregated <- rep(0, processes)
  
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
    
    ## STDE-PDE: Spatio-Temporal Density Estimation with PDE Regulariz. --------
    t0 <- proc.time()
    
    # Smoothing Parameters
    lambda <- 0.01
    # alternative: lambda <- 10^seq(from = -3, to = -1, by = 1)
    
    lambda_time <- 0.01
    # alternative: lambda_time <- 10^seq(from = -3, to = -1, by = 1)
    
    # Scaling Factor
    scaling_factor <- select.scaling(scaling = N^c(from = 0.4, to = 0.6, by = 0.1),
                                     data = data, mesh = mesh, mesh_time = mesh_time)
    
    # Solution
    # [If lambda and/or lambda_time are vectors, to select the best proposals by
    # 10-folds CV, please specify: preprocess_method = "RightCV"]
    solution_STDEPDE <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                    mesh_time = mesh_time, scaling = scaling_factor, lambda = lambda, tol1 = 1e-7,
                                    lambda_time = lambda_time, fvec = NULL, heatStep = 0.1,
                                    heatIter = 10, print = F, nfolds = 10, nsimulations = 10000,
                                    step_method = "Wolfe_Method", direction_method = "L-BFGS5",
                                    preprocess_method = "NoCrossValidation") 
    
    FEMfunction_STDEPDE <- FEM.time(solution_STDEPDE$g, mesh_time, FEMbasis, FLAG_PARABOLIC = FALSE)
    FEMfunction_STDEPDE_lwb <- FEM.time(solution_STDEPDE$g_CI_L, mesh_time, FEMbasis, FLAG_PARABOLIC = FALSE)
    FEMfunction_STDEPDE_upb <- FEM.time(solution_STDEPDE$g_CI_U, mesh_time, FEMbasis, FLAG_PARABOLIC = FALSE)
    
    # CPU Time
    CPUtime <- proc.time() - t0
    CPUtimes_STDEPDE[proc] <- as.numeric(CPUtime[3])

    ## INLA: sPDE Spatio-Temporal Model ----------------------------------------
    t0 <- proc.time()
    
    # Temporal Discretization
    k <- length(mesh_time)
    mesh_time_INLA <- INLA::inla.mesh.1d(mesh_time)
    
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
      #r <- INLA::inla.mesh.project(prj, res_INLA$summary.ran$s$mean[1:m + (j - 1) * m])
      r <- INLA::inla.mesh.project(prj, exp(res_INLA$summary.fix$mean + res_INLA$summary.ran$s$mean[1:m + (j - 1) * m]) * st.vol[1:m + (j - 1) * m] / N)
      r[is.na(ov)] <- NA
      return(r) 
    })
    
    IC.lwb.prj <- lapply(1:k, function(j) {
      #r <- INLA::inla.mesh.project(prj, res_INLA$summary.ran$s$`0.025quant`[1:m + (j - 1) * m])
      r <- INLA::inla.mesh.project(prj, exp(res_INLA$summary.fix$`0.025quant` + res_INLA$summary.ran$s$`0.025quant`[1:m + (j - 1) * m]) * st.vol[1:m + (j - 1) * m] / N)
      r[is.na(ov)] <- NA
      return(r) 
    })
    
    IC.upb.prj <- lapply(1:k, function(j) {
      #r <- INLA::inla.mesh.project(prj, res_INLA$summary.ran$s$`0.975quant`[1:m + (j - 1) * m])
      r <- INLA::inla.mesh.project(prj, exp(res_INLA$summary.fix$`0.975quant` + res_INLA$summary.ran$s$`0.975quant`[1:m + (j - 1) * m]) * st.vol[1:m + (j - 1) * m] / N)
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
    IC_lwb_INLA <- IC.lwb.prj
    IC_upb_INLA <- IC.upb.prj
    for(i in 1:k){
      factor <- floor(log(min(IC.lwb.prj[[i]]) / min(solution_INLA[[i]]), 10))+2
      IC_lwb_INLA[[i]] <- IC.lwb.prj[[i]] * 10^{-factor-1}
      IC_upb_INLA[[i]] <- IC.upb.prj[[i]] * 10^{factor+1}
    }
    FEMbasis_INLA <- create.FEM.basis(create.mesh.2D(grid_INLA))
    
    # CPU Time
    CPUtime <- proc.time() - t0
    CPUtimes_INLA[proc] <- as.numeric(CPUtime[3])
    
    ## INLA-aggregated: sPDE Spatio-Temporal Aggregated Model (Large Dataset) ----
    t0 <- proc.time()
    
    # Temporal Discretization
    k <- length(mesh_time)
    mesh_time_INLA <- INLA::inla.mesh.1d(mesh_time)
    
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
    
    # Spatial Domain
    win <- owin(xrange = x_lim, yrange = y_lim)
    win[["bdry"]][[1]][["x"]] <- c(x_lim[1], x_lim[2], x_lim[2],x_lim[1])
    win[["bdry"]][[1]][["y"]] <- c(y_lim[1],y_lim[1], y_lim[2], y_lim[2])
    win[["bdry"]][[1]][["area"]] <- diff(x_lim)*diff(y_lim)
    #save(win, file = "window")
    
    # Spatio-Temporal Planar Point Pattern
    idx <- which(locations[,1] >= x_lim[1] & locations[,1] <= x_lim[2] & locations[,2] >= y_lim[1] & locations[,2] <= y_lim[2])
    xy_INLA <- ppp(x = locations[idx,1], y = locations[idx,2], window = win)
    xyt_INLA <- stppp(xy_INLA, times[idx], tlim = c(0, 1))
    
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
      #z.j <- res_INLA_aggregated$summary.ran$s$mean[idx_INLA$s.group == j]
      z.j <- exp(res_INLA_aggregated$summary.fix$mean + res_INLA_aggregated$summary.ran$s$mean[idx_INLA$s.group == j]) * w.areas / N
      z <- INLA::inla.mesh.project(prj, z.j)
      z[g.no.in] <- NA
      return(z)
    })
    
    IC.lwb.prj <- lapply(1:k, function(j) {
      #r <- INLA::inla.mesh.project(prj, res_INLA_aggregated$summary.ran$s$`0.025quant`[idx_INLA$s.group == j])
      r <- INLA::inla.mesh.project(prj, exp(res_INLA_aggregated$summary.fix$`0.025quant` + res_INLA_aggregated$summary.ran$s$`0.025quant`[idx_INLA$s.group == j]) * w.areas / N)
      r[is.na(ov)] <- NA
      return(r) 
    })
    
    IC.upb.prj <- lapply(1:k, function(j) {
      #r <- INLA::inla.mesh.project(prj, res_INLA_aggregated$summary.ran$s$`0.975quant`[idx_INLA$s.group == j])
      r <- INLA::inla.mesh.project(prj, exp(res_INLA_aggregated$summary.fix$`0.975quant` + res_INLA_aggregated$summary.ran$s$`0.975quant`[idx_INLA$s.group == j]) * w.areas / N)
      r[is.na(ov)] <- NA
      return(r) 
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
    for(i in 1:k){
      solution_INLA_aggregated[[i]] <- solution_INLA_aggregated[[i]] * w.t[i]
    }
    IC_lwb_INLA_aggregated <- IC.lwb.prj
    for(i in 1:k){
      IC_lwb_INLA_aggregated[[i]] <- IC_lwb_INLA_aggregated[[i]] * w.t[i]
    }
    IC_upb_INLA_aggregated <- IC.upb.prj
    for(i in 1:k){
      IC_upb_INLA_aggregated[[i]] <- IC_upb_INLA_aggregated[[i]] * w.t[i]
    }
    FEMbasis_INLA_aggregated <- create.FEM.basis(create.mesh.2D(grid_INLA_aggregated))
    
    CPUtime <- proc.time() - t0
    CPUtimes_INLA_aggregated[proc] <- as.numeric(CPUtime[3])
    
    ## EVALUATION --------------------------------------------------------------
    # Time instants at which the solution is evaluated
    t <- mesh_time

    # Mean Solutions over the Processes
    if (proc == 1){
      
      mean_sol_STDEPDE <- matrix(rep(0, nrow(mesh.eval$nodes)*length(t)), nrow = nrow(mesh.eval$nodes), ncol = length(t))
      mean_sol_INLA <- matrix(rep(0, nrow(mesh.eval$nodes)*length(t)), nrow = nrow(mesh.eval$nodes), ncol = length(t))
      mean_sol_INLA_aggregated <- matrix(rep(0, nrow(mesh.eval$nodes)*length(t)), nrow = nrow(mesh.eval$nodes), ncol = length(t))

      mean_IC_coverage_STDEPDE <- matrix(rep(0, nrow(mesh.eval$nodes)*length(t)), nrow = nrow(mesh.eval$nodes), ncol = length(t))
      mean_IC_coverage_INLA <- matrix(rep(0, nrow(mesh.eval$nodes)*length(t)), nrow = nrow(mesh.eval$nodes), ncol = length(t))
      mean_IC_coverage_INLA_aggregated <- matrix(rep(0, nrow(mesh.eval$nodes)*length(t)), nrow = nrow(mesh.eval$nodes), ncol = length(t))
      
      true_instant <- matrix(rep(0, nrow(mesh.eval$nodes)*length(t)), nrow = nrow(mesh.eval$nodes), ncol = length(t))
      
    }
    
    for (time_index in 1:length(t)) {
      # True Density
      if(proc == 1){
        
        true_instant[,time_index] <- dens.func(mesh.eval$nodes, t[time_index])
        true_instant[,time_index] <- true_instant[,time_index]/sum(true_instant[,time_index], na.rm = TRUE)
        
      }
      
      # STDE-PDE
      evaluation_STDEPDE <- eval.FEM.time(FEM.time = FEMfunction_STDEPDE, locations = mesh.eval$nodes, time.instants = t[time_index])
      evaluation_STDEPDE <- exp(evaluation_STDEPDE)
      evaluation_STDEPDE <- evaluation_STDEPDE/sum(evaluation_STDEPDE, na.rm = TRUE)
      
      mean_sol_STDEPDE[,time_index] <- mapply(sum, mean_sol_STDEPDE[,time_index], evaluation_STDEPDE, na.rm = TRUE)
      
      evaluation_STDEPDE_lwb <- eval.FEM.time(FEM.time = FEMfunction_STDEPDE_lwb, locations = mesh.eval$nodes, time.instants = t[time_index])
      evaluation_STDEPDE_lwb <- exp(evaluation_STDEPDE_lwb)
      evaluation_STDEPDE_lwb <- evaluation_STDEPDE_lwb/sum(evaluation_STDEPDE_lwb, na.rm = TRUE)
      
      evaluation_STDEPDE_upb <- eval.FEM.time(FEM.time = evaluation_STDEPDE_upb, locations = mesh.eval$nodes, time.instants = t[time_index])
      evaluation_STDEPDE_upb <- exp(evaluation_STDEPDE_upb)
      evaluation_STDEPDE_upb <- evaluation_STDEPDE_upb/sum(evaluation_STDEPDE_upb, na.rm = TRUE)
      
      true <- true_instant[,time_index]
      condition <- (!is.na(true) & !is.na(evaluation_STDEPDE_lwb) & !is.na(evaluation_STDEPDE_upb) &
                      true >= evaluation_STDEPDE_lwb & true <= evaluation_STDEPDE_upb)
      
      mean_IC_coverage_STDEPDE[,time_index] <- mapply(sum, mean_IC_coverage_STDEPDE[,time_index], as.numeric(condition), na.rm = TRUE)
      
      # INLA
      coeff <- c(solution_INLA[[time_index]])
      FEM_INLA <- FEM(coeff, FEMbasis_INLA)
      evaluation_INLA <- eval.FEM(FEM_INLA, mesh.eval$nodes)
      evaluation_INLA <- evaluation_INLA/sum(evaluation_INLA, na.rm = TRUE)
      
      mean_sol_INLA[,time_index] <- mapply(sum, mean_sol_INLA[,time_index], evaluation_INLA, na.rm = TRUE)
      
      coeff <- c(IC_lwb_INLA[[time_index]])
      FEM_IC_lwb_INLA <- FEM(coeff, FEMbasis_INLA)
      evaluation_IC_lwb_INLA <- eval.FEM(FEM_IC_lwb_INLA, mesh.eval$nodes)
      
      coeff <- c(IC_upb_INLA[[time_index]])
      FEM_IC_upb_INLA <- FEM(coeff, FEMbasis_INLA)
      evaluation_IC_upb_INLA <- eval.FEM(FEM_IC_upb_INLA, mesh.eval$nodes)
      
      true <- true_instant[,time_index]
      condition <- (!is.na(true) & !is.na(evaluation_IC_lwb_INLA) & !is.na(evaluation_IC_upb_INLA) &
                      true >= evaluation_IC_lwb_INLA & true <= evaluation_IC_upb_INLA)
      mean_IC_coverage_INLA[,time_index] <- mapply(sum, mean_IC_coverage_INLA[,time_index], as.numeric(condition), na.rm = TRUE)
      
      # INLA-aggregated
      coeff <- c(solution_INLA_aggregated[[time_index]])
      FEM_INLA_aggregated <- FEM(coeff, FEMbasis_INLA_aggregated)
      evaluation_INLA_aggregated <- eval.FEM(FEM_INLA_aggregated, mesh.eval$nodes)
      evaluation_INLA_aggregated <- evaluation_INLA_aggregated/sum(evaluation_INLA_aggregated, na.rm = TRUE)
      
      mean_sol_INLA_aggregated[,time_index] <- mapply(sum, mean_sol_INLA_aggregated[,time_index], evaluation_INLA_aggregated, na.rm = TRUE)
      
      coeff <- c(IC_lwb_INLA_aggregated[[time_index]])
      FEM_IC_lwb_INLA_aggregated <- FEM(coeff, FEMbasis_INLA_aggregated)
      evaluation_IC_lwb_INLA_aggregated <- eval.FEM(FEM_IC_lwb_INLA_aggregated, mesh.eval$nodes)
      
      coeff <- c(IC_upb_INLA_aggregated[[time_index]])
      FEM_IC_upb_INLA_aggregated <- FEM(coeff, FEMbasis_INLA_aggregated)
      evaluation_IC_upb_INLA_aggregated <- eval.FEM(FEM_IC_upb_INLA_aggregated, mesh.eval$nodes)
      
      true <- true_instant[,time_index]
      condition <- (!is.na(true) & !is.na(evaluation_IC_lwb_INLA_aggregated) & !is.na(evaluation_IC_upb_INLA_aggregated) &
                      true >= evaluation_IC_lwb_INLA_aggregated & true <= evaluation_IC_upb_INLA_aggregated)
      mean_IC_coverage_INLA_aggregated[,time_index] <- mapply(sum, mean_IC_coverage_INLA_aggregated[,time_index], as.numeric(condition), na.rm = TRUE)
    
    }
    
    print(paste("Process", proc, "done in", round(as.numeric(proc.time() - t_proc)[3]), "seconds."))
    
  }
    
    mean_sol_INLA <- mean_sol_INLA/processes
    mean_sol_INLA_aggregated <- mean_sol_INLA_aggregated/processes
    
    mean_IC_coverage_INLA <- mean_IC_coverage_INLA/processes
    mean_IC_coverage_INLA_aggregated <- mean_IC_coverage_INLA_aggregated/processes
    
}

colMeans(mean_IC_coverage_STDEPDE[,2:6])
colMeans(mean_IC_coverage_INLA[,2:6])
colMeans(mean_IC_coverage_INLA_aggregated[,2:6])

mean_IC_coverage <- list(mean_IC_coverage_STDEPDE, mean_IC_coverage_INLA, mean_IC_coverage_INLA_aggregated)
save(mean_IC_coverage, file = paste0("output/mean_IC_coverage", N, "data.rda"))

## COVERAGE PROBABILITY [FIGURE 6] ---------------------------------------------
mean_IC_coverage_STDEPDE <- as.data.frame(mean_IC_coverage_STDEPDE[,2:6])
mean_IC_coverage_INLA <- as.data.frame(mean_IC_coverage_INLA[,2:6])
mean_IC_coverage_INLA_aggregated <- as.data.frame(mean_IC_coverage_INLA_aggregated[,2:6])

# Grid
n <- 32
X <- seq(-6, 6, length.out = n)
Y <- seq(-6, 6, length.out = n)
grid <- expand.grid(X, Y)

t <- seq(from = 0.1, to = 0.9, by = 0.2)

# Inference at Locations where True Density is Far from 0
true_log_density <- NULL

idx <- list()

for(i in 1:length(t)){
  true_log_density[[i]] <- log(dens.func(grid, t[i]))
  idx[[i]] <- which(true_g[[i]]>-6)
}

selected_values_STDEPDE <- lapply(1:5, function(i) {
  mean_IC_coverage_STDEPDE[idx[[i]],i]
})

selected_values_INLA <- lapply(1:5, function(i) {
  mean_IC_coverage_INLA[idx[[i]],i]
})

selected_values_INLA_aggregated <- lapply(1:5, function(i) {
  mean_IC_coverage_INLA_aggregated[idx[[i]],i]
})

blue_polimi <- rgb(114, 143, 165, max = 255, alpha = 65, names = "blue_polimi25") # #728FA5

pdf(paste0("pictures/", N, "_data/sim1_coverage_probability.pdf"), family = "serif", width = 12, height = 5.3)
#x11(width = 12, height = 5.3)
par(mfrow = c(1,3), mai = c(1.5,0.35,0.75,0.15))
boxplot(selected_values_STDEPDE, outline = FALSE, ylim = c(0,1),
        asp = 1, xaxt = "n", yaxt = "n")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = blue_polimi)
grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 1)
par(new = TRUE)
boxplot(selected_values_STDEPDE, outline = FALSE, ylim = c(0,1),
        names = c("0.10", "0.30", "0.50", "0.70", "0.90"), main = "\nSTDE-PDE",
        col = brewer.pal(8, "YlGnBu")[1], cex.lab = 2, cex.axis = 2, cex.main = 2, xlab = "Time")
abline(h = 0.95, lwd = 2, col = "red")

boxplot(selected_values_INLA, outline = FALSE, ylim = c(0,1),
        asp = 1, xaxt = "n", yaxt = "n")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = blue_polimi)
grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 1)
par(new = TRUE)
boxplot(selected_values_INLA, outline = FALSE, ylim = c(0,1),
        names = c("0.10", "0.30", "0.50", "0.70", "0.90"), main = "Coverage probability\nINLA",
        col = brewer.pal(8, "YlGnBu")[3], cex.lab = 2, cex.axis = 2, cex.main = 2, xlab = "Time")
abline(h = 0.95, lwd = 2, col = "red")

boxplot(selected_values_INLA_aggregated, outline = FALSE, ylim = c(0,1),
        asp = 1, xaxt = "n", yaxt = "n")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = blue_polimi)
grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 1)
par(new = TRUE)
boxplot(selected_values_INLA_aggregated, outline = FALSE, ylim = c(0,1),
        names = c("0.10", "0.30", "0.50", "0.70", "0.90"), main = "\nINLA-agg",
        col = brewer.pal(8, "YlGnBu")[4], cex.lab = 2, cex.axis = 2, cex.main = 2, xlab = "Time")
abline(h = 0.95, lwd = 2, col = "red")
dev.off()
