#'############################################################################'#
#'##################### SPATIO-TEMPORAL DENSITY ESTIMATION ###################'# 
#'############### SIMULATION 2 - MIXTURE OF 5 KENT DISTRIBUTIONS #############'#
#'############################################################################'#

## LIBRARIES AND FUNCTIONS -----------------------------------------------------
source("libraries_2.R")

## SPATIO-TEMPORAL DISCRETIZATION ----------------------------------------------
# Domain Area
domain_area <- 4*pi

# Spatial 2D Mesh over the Unitary Sphere for Estimation
vertices <- read.table("mesh/sphere.606.vertices.txt", quote = "\"", comment.char = "")
triangles <- read.table("mesh/sphere.606.triangles.txt", quote = "\"", comment.char = "")
mesh <- create.mesh.2.5D(nodes = vertices[,1:3], triangles = triangles[,1:3])
FEMbasis <- create.FEM.basis(mesh)

plot.mesh(mesh)
rgl.close()

# Fine Spatial Mesh over the Unitary Sphere for Evaluation
vertices.eval <- read.table("mesh/sphere.5016.vertices.txt", row.names = 1)
vertices.eval.proj <- projection.points.2.5D(mesh, vertices.eval)
triangles.eval <- read.table("mesh/sphere.5016.triangles.txt", row.names = 1)
mesh.eval <- create.mesh.2.5D(nodes = vertices.eval, triangles = triangles.eval)
FEMbasis.eval <- create.FEM.basis(mesh.eval)

plot.mesh(mesh.eval)
rgl.close()

# Fine Grid Size for Evaluation
n <- 32

# Temporal 1D Mesh Over [0,1]
mesh_time <- seq(from = 0, to = 1, by = 0.1)

## ESTIMATION PROCEDURE --------------------------------------------------------
# Number of Processes
processes <- 30

# Sample Size
NN <- 10000
# alternative: NN <- c(1000, 2500, 5000, 10000)

for(N in NN){
  
  # Containers for Errors
  mise_STDEPDE <- rep(0, processes)
  mise_STKDE_discrete <- rep(0, processes)
  mise_STKDE_separable <- rep(0, processes)
  mise_KNNSTDE_separable <- rep(0, processes)
  
  KLdist_STDEPDE <- rep(0, processes)
  KLdist_STKDE_discrete <- rep(0, processes)
  KLdist_STKDE_separable <- rep(0, processes)
  KLdist_KNNSTDE_separable <- rep(0, processes)
  
  Hdist_STDEPDE <- rep(0, processes)
  Hdist_STKDE_discrete <- rep(0, processes)
  Hdist_STKDE_separable <- rep(0, processes)
  Hdist_KNNSTDE_separable <- rep(0, processes)
  
  # Containers for CPU Times
  CPUtimes_STDEPDE <- rep(0, processes)
  CPUtimes_STKDE_discrete <- rep(0, processes)
  CPUtimes_STKDE_separable <- rep(0, processes)
  CPUtimes_KNNSTDE_separable <- rep(0, processes)
  
  for(proc in 1:processes){
    # Time for One Process
    t_proc <- proc.time()
    
    ### DATA -------------------------------------------------------------------
    # Generate the Data
    generate.data(N, proc) # run only once to generate the data
    
    # Read the Data
    data <- read.table(paste0("data/",N,"data_",proc,".txt"))
    
    # Locations
    locations <- data[,1:3]
    
    # Times
    times <- data[,4]
    
    # Discrete Times
    discrete_times <- (times-0.05) %/% 0.1
    discrete_times[discrete_times == -1] <- 0
    discrete_times[discrete_times == 9] <- 8
  
    ## STDE-PDE: Spatio-Temporal Density Estimation with PDE Regulariz. --------
    t0 <- proc.time()
    
    # Smoothing Parameters
    lambda <- 0.01
    # alternative: lambda <- 10^seq(from = -3, to = -1, by = 1)
    
    lambda_time <- 0.001
    # alternative: lambda_time <- 10^seq(from = -4, to = -2, by = 1)
    
    # Solution
    # [If lambda and/or lambda_time are vectors, to select the best proposals by
    # 10-folds CV, please specify: preprocess_method = "RightCV"]
    solution_STDEPDE <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                    mesh_time = mesh_time, lambda = lambda, tol1 = 1e-7,
                                    lambda_time = lambda_time, fvec = NULL, heatStep = 0.1,
                                    heatIter = 10, print = F, nfolds = 10, nsimulations = 10000,
                                    step_method = "Wolfe_Method", direction_method = "L-BFGS5",
                                    preprocess_method = "NoCrossValidation")
    
    FEMfunction_STDEPDE <- FEM.time(solution_STDEPDE$g, mesh_time, FEMbasis, FLAG_PARABOLIC=F)
    
    # CPU Time
    CPUtime <- proc.time() - t0
    CPUtimes_STDEPDE[proc] <- CPUtime[3]
    
    ## STKDE-discrete: Spatio-Temporal Kernel Density Estimation ---------------
    t0 <- proc.time()
    
    # Radius
    R <- 1
    
    # Grid Size
    n_STKDE_discrete <- n
    
    # Latitudes and Longitudes of the Observed Locations
    lat <- asin(locations[,3] / R) * 180 / pi
    long <- atan2(locations[,2], locations[,1]) * 180 / pi
    
    # Latitudes and Longitudes of Mesh Nodes for Evaluation
    mesh.eval.nodes.lat <- asin(FEMbasis.eval$mesh$nodes[,3] / R) * 180 / pi
    mesh.eval.nodes.long <- atan2(FEMbasis.eval$mesh$nodes[,2], FEMbasis.eval$mesh$nodes[,1]) * 180 / pi
    
    # Grid
    long_STKDE_discrete <- seq(from = min(long), to = max(long), length.out = n_STKDE_discrete)
    lat_STKDE_discrete <- seq(from = min(lat), to = max(lat), length.out = n_STKDE_discrete)
    grid_STKDE_discrete_latlong <- expand.grid(long_STKDE_discrete, lat_STKDE_discrete)
    names(grid_STKDE_discrete_latlong) <- c("long", "lat")
    FEMbasis_STKDE_discrete <- create.FEM.basis(create.mesh.2D(grid_STKDE_discrete_latlong))
    
    # Solution
    x11(width = 25, height = 14)
    solution_STKDE_discrete <- stkde(xlong = long, ylat = lat, ztime = discrete_times,
                                     xgrids = n_STKDE_discrete, ygrids = n_STKDE_discrete, breaks = 0.05,
                                     alpha = 0.05, nrowspar = 5, bwmethod = "normal-reference") # alternative: bwmethod = "cv.ml"
    dev.off()
  
    # CPU Time
    CPUtime <- proc.time() - t0
    CPUtimes_STKDE_discrete[proc] <- CPUtime[3]
    
    ## STKDE-separable: STKDE for 1st-Order Separable Spatio-Temporal Point Processes ----
    t0 <- proc.time()
    
    # Spatial Component
    locations_spherical <- sphere.cart.to.sphere.spherical(as.matrix(locations))
    kappa <- kde.compute.concentration(locations_spherical)
    f_space_STKDE_separable <- kde.fhat.cart(vertices.eval, as.matrix(locations), kappa = kappa)
    
    # Temporal Component
    f_time_STKDE_separable <- ks::kde(x = as.matrix(times), gridsize = 50)
    
    # CPU Time
    CPUtime <- proc.time() - t0
    CPUtimes_STKDE_separable[proc] <- CPUtime[3]
    
    ## KNNSTDE-separable: k-Nearest Neighbors Spatio-Temporal Density Estimation ---------
    ## for 1st-Order Separable Spatio-Temporal Point Processes
    t0 <- proc.time()
    
    # Spatial Component
    f_space_KNNSTDE_separable <- knnDE(as.matrix(locations), as.matrix(vertices.eval[,1:3]), 500)
    
    # Temporal Component
    f_time_KNNSTDE_separable <- knnDE(as.matrix(times), as.matrix(seq(from = 0, to = 10, length.out = 25)), 50)
    
    # CPU Time
    CPUtime <- proc.time() - t0
    CPUtimes_KNNSTDE_separable[proc] <- CPUtime[3]
    
    ## ERRORS ------------------------------------------------------------------
    # Time instants at which the solution is evaluated
    t <- mesh_time
    t_discrete <- seq(from = 0.1, to = 0.9, by = 0.1)
    
    # Time Half Width
    h <- (t_discrete[2]-t_discrete[1])/2
    
    # Mean Solutions over the Processes
    if (proc == 1){
      
      mean_sol_STDEPDE <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
      mean_sol_STKDE_discrete <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t_discrete))
      mean_sol_STKDE_separable <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
      mean_sol_KNNSTDE_separable <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
      
      true_instant <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
      true_discrete <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t_discrete))
        
    }
    
    for(time_index in 1:length(t)) {
      
      # True Density
      if(proc == 1){
        
        true_instant[,time_index] <- dens.func(vertices.eval.proj, rep(t[time_index], nrow(vertices.eval.proj)))
        true_instant[,time_index] <- true_instant[,time_index] / sum(true_instant[,time_index], na.rm = TRUE)
        
        if(sum(abs(t[time_index] - t_discrete) < 1e-06)){
          time_index_discrete <- which(abs(t[time_index] - t_discrete) < 1e-06)
        
          instants <- 11
          true_subinterval <- rep(0,nrow(vertices.eval.proj))
          subinterval <- seq(t[time_index]-h, t[time_index]+h, length.out = instants)
          for(time_subinterval in subinterval){
            coeff <- dens.func(vertices.eval.proj, rep(time_subinterval, nrow(vertices.eval.proj)))
            coeff <- coeff / sum(coeff)
            true_subinterval <- true_subinterval + (coeff)*(subinterval[2]-subinterval[1])
          }
          
          true_discrete[,time_index_discrete] <- true_subinterval / sum(true_subinterval)
          
        }
      }
      
      # STDE-PDE
      evaluation_STDEPDE <- eval.FEM.time(FEM.time = FEMfunction_STDEPDE, locations = vertices.eval.proj, time.instants = t[time_index])
      evaluation_STDEPDE <- exp(evaluation_STDEPDE)
      evaluation_STDEPDE <- evaluation_STDEPDE / sum(evaluation_STDEPDE, na.rm = TRUE)
      
      if(sum(abs(t[time_index] - t_discrete) < 1e-06)){
        time_index_discrete <- which(abs(t[time_index] - t_discrete) < 1e-06)
        
        # STKDE-discrete
        FEM_STKDE_discrete <- FEM(c(solution_STKDE_discrete$dens[,,time_index_discrete]), FEMbasis_STKDE_discrete)
        evaluation_STKDE_discrete <- eval.FEM(FEM_STKDE_discrete, locations = cbind(mesh.eval.nodes.long, mesh.eval.nodes.lat))
        evaluation_STKDE_discrete[is.na(evaluation_STKDE_discrete)] <- 0.0
        evaluation_STKDE_discrete <- evaluation_STKDE_discrete / sum(evaluation_STKDE_discrete, na.rm = TRUE)
        
      }  
      
      # STKDE-separable
      marginal_time_STKDE_separable <- f_time_STKDE_separable$estimate[findInterval(t[time_index], f_time_STKDE_separable$eval.points)]
      solution_STKDE_separable <- marginal_time_STKDE_separable * f_space_STKDE_separable
      solution_STKDE_separable <- solution_STKDE_separable / sum(solution_STKDE_separable, na.rm = TRUE)
      
      # KNNSTDE-separable
      marginal_time_KNNSTDE_separable <- f_time_KNNSTDE_separable[findInterval(t[time_index], seq(from = 0, to = 1, length.out = 25))]
      solution_KNNSTDE_separable <- marginal_time_KNNSTDE_separable * f_space_KNNSTDE_separable
      solution_KNNSTDE_separable <- solution_KNNSTDE_separable / sum(solution_KNNSTDE_separable, na.rm = TRUE)
      
      # Estimates and Errors Updates
      
      # STDE-PDE
      mise_STDEPDE[proc] <- mise_STDEPDE[proc] + domain_area*mean((evaluation_STDEPDE - true_instant[,time_index])^2, na.rm = TRUE)/length(t)
      KLdist_STDEPDE[proc] <- KLdist_STDEPDE[proc] + abs(kullback_leibler_distance(evaluation_STDEPDE, true_instant[,time_index], testNA = F, unit = "log", epsilon = 1e-4))/length(t)
      Hdist_STDEPDE[proc] <- Hdist_STDEPDE[proc] + domain_area*hellinger.distance(evaluation_STDEPDE, true_instant[,time_index])/length(t)
      
      mean_sol_STDEPDE[,time_index] <- mapply(sum, mean_sol_STDEPDE[,time_index], evaluation_STDEPDE, na.rm = TRUE)
      
      if(sum(abs(t[time_index] - t_discrete) < 1e-06)){
        # STKDE-discrete
        time_index_discrete <- which(abs(t[time_index] - t_discrete) < 1e-06)
        
        mise_STKDE_discrete[proc] <- mise_STKDE_discrete[proc] + domain_area*mean((evaluation_STKDE_discrete - true_discrete[,time_index_discrete])^2, na.rm = TRUE)/length(t_discrete)
        KLdist_STKDE_discrete[proc] <- KLdist_STKDE_discrete[proc] + abs(kullback_leibler_distance(evaluation_STKDE_discrete, true_discrete[,time_index_discrete], testNA = F, unit = "log", epsilon = 1e-4))/length(t_discrete)
        Hdist_STKDE_discrete[proc] <- Hdist_STKDE_discrete[proc] + domain_area*hellinger.distance(evaluation_STKDE_discrete, true_discrete[,time_index_discrete])/length(t_discrete)
        
        mean_sol_STKDE_discrete[,time_index_discrete] <- mapply(sum, mean_sol_STKDE_discrete[,time_index_discrete], evaluation_STKDE_discrete, na.rm = TRUE)
      }
      
      # STKDE-separable
      mise_STKDE_separable[proc] <- mise_STKDE_separable[proc] + domain_area*mean((solution_STKDE_separable - true_instant[,time_index])^2, na.rm = TRUE)/length(t)
      KLdist_STKDE_separable[proc] <- KLdist_STKDE_separable[proc] + abs(kullback_leibler_distance(solution_STKDE_separable, true_instant[,time_index], testNA = F, unit = "log", epsilon = 1e-4))/length(t)
      Hdist_STKDE_separable[proc] <- Hdist_STKDE_separable[proc] + domain_area*hellinger.distance(solution_STKDE_separable, true_instant[,time_index])/length(t)
      
      mean_sol_STKDE_separable[,time_index] <- mapply(sum, mean_sol_STKDE_separable[,time_index], solution_STKDE_separable, na.rm = TRUE)
      
      # KNNSTDE-separable
      mise_KNNSTDE_separable[proc] <- mise_KNNSTDE_separable[proc] + domain_area*mean((solution_KNNSTDE_separable - true_instant[,time_index])^2, na.rm = TRUE)/length(t)
      KLdist_KNNSTDE_separable[proc] <- KLdist_KNNSTDE_separable[proc] + abs(kullback_leibler_distance(solution_KNNSTDE_separable, true_instant[,time_index], testNA = F, unit = "log", epsilon = 1e-4))/length(t)
      Hdist_KNNSTDE_separable[proc] <- Hdist_KNNSTDE_separable[proc] + domain_area*hellinger.distance(solution_KNNSTDE_separable, true_instant[,time_index])/length(t)
      
      mean_sol_KNNSTDE_separable[,time_index] <- mapply(sum, mean_sol_KNNSTDE_separable[,time_index], solution_KNNSTDE_separable, na.rm = TRUE)
      
    } 
    
    print(paste("Process", proc, "done in", round(as.numeric(proc.time() - t_proc)[3]), "seconds."))
    
  }
  
  mean_sol_STDEPDE <- mean_sol_STDEPDE/processes
  mean_sol_STKDE_discrete <- mean_sol_STKDE_discrete/processes
  mean_sol_STKDE_separable <- mean_sol_STKDE_separable/processes
  mean_sol_KNNSTDE_separable <- mean_sol_KNNSTDE_separable/processes

  #'##########################################################################'#
  #'##########################################################################'#
  
  ## EXPORT RESULTS ------------------------------------------------------------
  # Export estimates and errors
  dir.create(file.path(getwd(), "output"), showWarnings = FALSE)
  
  estimates <- list(mean_sol_STDEPDE, mean_sol_STKDE_discrete, mean_sol_STKDE_separable, mean_sol_KNNSTDE_separable)
  save(estimates, file = paste0("output/estimates_", N, "data.rda"))
  
  mise <- list(mise_STDEPDE, mise_STKDE_discrete, mise_STKDE_separable, mise_KNNSTDE_separable)
  save(mise, file = paste0("output/mise_", N, "data.rda"))
  
  KLdist <- list(KLdist_STDEPDE, KLdist_STKDE_discrete, KLdist_STKDE_separable, KLdist_KNNSTDE_separable)
  save(KLdist, file = paste0("output/KLdist_", N, "data.rda"))
  
  Hdist <- list(Hdist_STDEPDE, Hdist_STKDE_discrete, Hdist_STKDE_separable, Hdist_KNNSTDE_separable)
  save(Hdist, file = paste0("output/Hdist_", N, "data.rda"))
  
  CPUtimes <- list(CPUtimes_STDEPDE, CPUtimes_STKDE_discrete, CPUtimes_STKDE_separable, CPUtimes_KNNSTDE_separable)
  save(CPUtimes, file = paste0("output/CPUtimes_", N, "data.rda"))
  
  base::save.image(paste0("output/workspace_", N, "data.RData"))
  
  #'##########################################################################'#
  #'##########################################################################'#

  ## VISUALIZATION -------------------------------------------------------------
  {
    ### ESTIMATES [FIGURE 8] ---------------------------------------------------
    # Import estimates and errors
    # load(file = paste0("output/estimates_", N, "data.rda"))
    # mean_sol_STDEPDE <- estimates[[1]]
    # mean_sol_STKDE_discrete <- estimates[[2]]
    # mean_sol_STKDE_separable <- estimates[[3]]
    # mean_sol_KNNSTDE_separable <- estimates[[4]]
    
    # Time instants at which the solution is evaluated
    t <- mesh_time
    t_discrete <- seq(from = 0.1, to = 0.9, by = 0.1)
    
    # Time Half Width
    h <- (t_discrete[2]-t_discrete[1])/2
    
    # Plots
    M <- max(max(true_instant, na.rm = TRUE),
             # max(mean_sol_STDEPDE, na.rm = TRUE),
             # max(mean_sol_STKDE_discrete*instants, na.rm = TRUE),
             # max(mean_sol_STKDE_separable, na.rm = TRUE),
             # max(mean_sol_KNNSTDE_separable, na.rm = TRUE),
             na.rm = TRUE)
    
    dir.create(file.path(getwd(), "pictures"), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data")), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data/sample")), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data/true_instant")), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data/STDEPDE")), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data/STKDE_discrete")), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data/STKDE_separable")), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data/KNNSTDE_separable")), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data/time_bar")), showWarnings = FALSE)
    
    for(time_index in 1:length(t)){
      
      t_img <- proc.time()
      
      # First Plot: Sample
      plot.sample(FEM(rep(0, dim(vertices.eval.proj)[1]), FEMbasis.eval), locations[abs(times-t[time_index])<h,])
      string <- paste0("pictures/", N, "_data/sample/sim2_sample_", N, "data_", t[time_index], ".png")
      snapshot3d(string, fmt = "png", width = 750, height = 750, webshot = rgl.useNULL())
      rgl.close()
      
      # Second Plot: True Density at t[time_index]
      plot.FEM(FEM(true_instant[,time_index], FEMbasis.eval), m = 0, M = M)
      string <- paste0("pictures/", N, "_data/true_instant/sim2_true_instant_", N, "data_", t[time_index], ".png")
      snapshot3d(string, fmt = "png", width = 750, height = 750, webshot = rgl.useNULL())
      rgl.close()
      
      # Third Plot: STDE-PDE Estimated Density at t[time_index]
      plot.FEM(FEM(mean_sol_STDEPDE[,time_index], FEMbasis.eval), m = 0, M = M)
      string <- paste0("pictures/", N, "_data/STDEPDE/sim2_STDEPDE_", N, "data_", t[time_index], ".png")
      snapshot3d(string, fmt = "png", width = 750, height = 750, webshot = rgl.useNULL())
      rgl.close()
      
      if(sum(abs(t[time_index] - t_discrete) < 1e-06)){
        time_index_discrete <- which(abs(t[time_index] - t_discrete) < 1e-06)
        
        # Fourth Plot: STKDE-discrete Estimated Density at t[time_index]
        plot.FEM(FEM(mean_sol_STKDE_discrete[,time_index_discrete], FEMbasis.eval), m = 0, M = M)
        string <- paste0("pictures/", N, "_data/STKDE_discrete/sim2_STKDE_discrete_", N, "data_", t[time_index], ".png")
        snapshot3d(string, fmt = "png", width = 750, height = 750, webshot = rgl.useNULL())
        rgl.close()
      }
      
      # Fifth Plot: STKDE-separable Estimated (Separable) Density at t[time_index]
      plot.FEM(FEM(mean_sol_STKDE_separable[,time_index], FEMbasis.eval), m = 0, M = M)
      string <- paste0("pictures/", N, "_data/STKDE_separable/sim2_STKDE_separable_", N, "data_", t[time_index], ".png")
      snapshot3d(string, fmt = "png", width = 750, height = 750, webshot = rgl.useNULL())
      rgl.close()
      
      # Sixth Plot: KNNSTDE-separable Estimated (Separable) Density at t[time_index]
      plot.FEM(FEM(mean_sol_KNNSTDE_separable[,time_index], FEMbasis.eval), m = 0, M = M)
      string <- paste0("pictures/", N, "_data/KNNSTDE_separable/sim2_KNNSTDE_separable_", N, "data_", t[time_index], ".png")
      snapshot3d(string, fmt = "png", width = 750, height = 750, webshot = rgl.useNULL())
      rgl.close()
      
      print(paste("Images for t =", t[time_index], "done in", round(as.numeric(proc.time() - t_img)[3]), "seconds."))
      
    }
  
    ## LEGEND [FIGURE 8] -------------------------------------------------------
    pdf(paste0("pictures/", N, "_data/sim2_legend.pdf"), family = "serif", width = 11, height = 3)
    #x11(width = 11, height = 3)
    par(mai = c(1,0.75,0,0))
    plot(c(0, 112.5), c(0, 15), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame.plot = F)
    gradient.rect(0, 0, 100, 5, col = jet.col(1000), border = "black")
    axis(1, at = c(0,33.33,66.66,100), labels = c("0", "1", "2", "3"),
         lwd.ticks = 2, cex.axis = 2, lwd = 2)
    text(107, 2, TeX("$\\times 10^{-3}$"), cex = 2)
    
    dev.off()
    
    print(paste0("M = ", round(M, 4)))
    
    ## ERRORS [FIGURE 7] -------------------------------------------------------
    {
      blue_polimi <- rgb(114, 143, 165, max = 255, alpha = 65, names = "blue_polimi25") # #728FA5
      
      pdf(paste0("pictures/", N, "_data/sim2_errors.pdf"), family = "serif", width = 12, height = 5.3)
      #x11(width = 12, height = 5.3)
      par(mfrow = c(1,3), mai = c(1.75,0.35,0.5,0.15))
      boxplot(mise_STDEPDE, mise_STKDE_discrete, mise_STKDE_separable, mise_KNNSTDE_separable, asp = 1, xaxt = "n", yaxt = "n")
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = blue_polimi)
      grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 1)
      par(new = TRUE)
      boxplot(mise_STDEPDE, mise_STKDE_discrete, mise_STKDE_separable, mise_KNNSTDE_separable,
              names = c("STDE-PDE", "STKDE-dis", "STKDE-sep", "KNNSTDE-sep"), main = TeX("$err_{L^2}$", bold = T),
              col = brewer.pal(8, "YlGnBu")[c(1,6,7,8)], cex.lab = 2, cex.axis = 2, cex.main = 2, las = 3)
      
      boxplot(KLdist_STDEPDE, KLdist_STKDE_discrete, KLdist_STKDE_separable, KLdist_KNNSTDE_separable, asp = 1, xaxt = "n", yaxt = "n")
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = blue_polimi)
      grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 1)
      par(new = TRUE)
      boxplot(KLdist_STDEPDE, KLdist_STKDE_discrete, KLdist_STKDE_separable, KLdist_KNNSTDE_separable,
              names = c("STDE-PDE", "STKDE-dis", "STKDE-sep", "KNNSTDE-sep"), main = "KL",
              col = brewer.pal(8, "YlGnBu")[c(1,6,7,8)], cex.lab = 2, cex.axis = 2, cex.main = 2, las = 3)
      
      boxplot(Hdist_STDEPDE, Hdist_STKDE_discrete, Hdist_STKDE_separable, Hdist_KNNSTDE_separable, asp = 1, xaxt = "n", yaxt = "n")
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = blue_polimi)
      grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 1)
      par(new = TRUE)
      boxplot(Hdist_STDEPDE, Hdist_STKDE_discrete, Hdist_STKDE_separable, Hdist_KNNSTDE_separable,
              names = c("STDE-PDE", "STKDE-dis", "STKDE-sep", "KNNSTDE-sep"), main = "H",
              col = brewer.pal(8, "YlGnBu")[c(1,6,7,8)], cex.lab = 2, cex.axis = 2, cex.main = 2, las = 3)
      dev.off()
    }
    
    ## CPU TIMES ---------------------------------------------------------------
    pdf(paste0("pictures/", N, "_data/sim2_CPUtimes.pdf"), family = "serif", width = 12, height = 5.3)
    #x11(width = 12, height = 5.3)
    par(mfrow = c(1,3), mai = c(1.75,0.35,0.5,0.15))
    plot.new()
    boxplot(CPUtimes_STDEPDE, CPUtimes_STKDE_discrete, CPUtimes_STKDE_separable, CPUtimes_KNNSTDE_separable,
            asp = 1, xaxt = "n", yaxt = "n")
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = blue_polimi)
    grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 1)
    par(new = TRUE)
    boxplot(CPUtimes_STDEPDE, CPUtimes_STKDE_discrete, CPUtimes_STKDE_separable, CPUtimes_KNNSTDE_separable,
            names = c("STDE-PDE", "STKDE-dis", "STKDE-sep", "KNNSTDE-sep"), main = "CPU Times [seconds]",
            col = brewer.pal(8, "YlGnBu")[c(1,6,7,8)], cex.lab = 2, cex.axis = 2, cex.main = 2, las = 3)
    dev.off()
    
    ## MESH --------------------------------------------------------------------
    string <- paste0("pictures/", N, "_data/sim2_mesh.png")
    plot.mesh(mesh)
    rgl.points(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3], col = "black")
    snapshot3d(string, fmt = "png", width = 750, height = 750, webshot = rgl.useNULL())
    rgl.close()
    
    string <- paste0("pictures/", N, "_data/sim2_mesh_eval.png")
    plot.mesh(mesh.eval)
    rgl.points(mesh.eval$nodes[,1], mesh.eval$nodes[,2], mesh.eval$nodes[,3], col = "black")
    snapshot3d(string, fmt = "png", width = 750, height = 750, webshot = rgl.useNULL())
    rgl.close()
    
    ## TIME BAR ----------------------------------------------------------------
    for(time_index in 1:length(t)){
    
      string <- paste0("pictures/", N, "_data/time_bar/sim2_time_bar", time_index, ".pdf")
      pdf(file = string, family = "serif", width = 7, height = 1.5)
      #x11(width = 7, height = 1.5)
      par(mar = c(0.01, 2, 0.01, 2))
      bar <- rbind(c(t[time_index],0), c(t[time_index],0.075))
      plot(bar, type = "l", bty = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "",
           pch = 3, lwd = 5, xlim = c(0,1), asp = 1)
      mtext("Time", line = -1.5, font = 2, cex = 1.5)
      axis(side = 1, c(0, 0.5, 1), pos = 0, labels = c("0","0.5","1"))
      dev.off()
      
    }
    
  }
}
