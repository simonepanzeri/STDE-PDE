#'############################################################################'#
#'########## SIMULATION 5 - INHOMOGENEOUS POINT PROCESS ON easynet ###########'#
#'############################################################################'#

graphics.off()
rm(list=ls())
options(warn = -1)

set.seed(23)

## LIBRARIES AND UTILITIES -----------------------------------------------------
if(!require(pacman)) install.packages("pacman")
pacman::p_load("rstudioapi")

# Setting the Working Directory 
setwd(dirname(getActiveDocumentContext()$path))

# Loading Packages and Auxiliary Functions
file.sources <- list.files("utils/", pattern = "*.R", full.names = TRUE)
sapply(file.sources, source, .GlobalEnv)

# Create the Current Directory
dir.create(paste0(getwd(),"/simulation_5"))
setwd(paste0(getwd(),"/simulation_5"))
foldername <- getwd()
dir.create(paste0(getwd(),"/data"))

## AVAILABLE DOMAINS -----------------------------------------------------------
domains <- c("estevan", "ontario", "london", "munster", "chicago", "valencia",
             "eastbourne", "medellin", "easynet", "simplenet") 

## SPATIO-TEMPORAL DISCRETIZATION ----------------------------------------------
d <- which(domains == "easynet")

# Spatial 1.5D Mesh over the Domain for Estimation
sett <- setting(domains[d]) 
mesh <- sett$mesh
mesh <- refine.mesh.1.5D(mesh, 0.025)
FEMbasis <- create.FEM.basis(mesh)
nnodes <- nrow(FEMbasis$mesh$nodes)
spatstat.linnet <- as.linnet(mesh)

# Fine Spatial 1.5D Mesh over the Domain for Evaluation
mesh.eval <- refine.mesh.1.5D(mesh, 0.023)
FEMbasis.eval <- create.FEM.basis(mesh.eval)
nnodes.eval <- nrow(FEMbasis.eval$mesh$nodes)
spatstat.linnet.eval <- as.linnet(mesh.eval)

# Grid Size for Estimation
n <- 512

# Domain Area
domain_area <- sum(fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis) %*% rep(1, nrow(FEMbasis$mesh$nodes)))

# Temporal 1D Mesh Over [0,1]
Ti <- 0
Tf <- 1
mesh_time <- seq(from = Ti, to = Tf, by = 0.2)

## ESTIMATION PROCEDURE --------------------------------------------------------
# Number of Processes
processes <- 50

# Sample Size
NN <- 500

for(N in NN){
  setwd(foldername)
  
  cat(paste0("[N = ", N, "]\n"))
  
  # Generate the Data (Inhomogeneous Poisson Process)
  cat("Computing the weight for lambda...\n")
  w <- compute.weight(fun = inhomogeneous_intensity4, f = FEMbasis.eval, N = N)
  
  # Generate the Data (Inhomogeneous Poisson Process)
  cat("Generating the data...\n")
  if(!file.exists(paste0(getwd(),"/data/dataset_",processes,"x",N,".RData"))){
    
    dataset <- rpoistlpp(lambda = inhomogeneous_intensity4, lmax = NULL,
                         a = Ti, b = Tf, L = spatstat.linnet, check = TRUE,
                         nsim = processes)
    
    save(dataset, file = paste0(getwd(),"/data/dataset_",processes,"x",N,".RData"))
    
  } else {
    
    load(paste0(getwd(),"/data/dataset_",processes,"x",N,".RData"))
    
  }
  
  # Containers for Errors
  mise_STDEPDE <- rep(0, processes)
  mise_STLNPP <- rep(0, processes)
  mise_HEAT <- rep(0, processes)
  mise_EQUAL_SPLIT <- rep(0, processes)
  mise_VORONOI <- rep(0, processes)
  mise_VORONOI_SEP <- rep(0, processes)
  
  # Containers for CPU Times
  CPUtimes_STDEPDE <- rep(0, processes)
  CPUtimes_STLNPP <- rep(0, processes)
  CPUtimes_HEAT <- rep(0, processes)
  CPUtimes_EQUAL_SPLIT <- rep(0, processes)
  CPUtimes_VORONOI <- rep(0, processes)
  CPUtimes_VORONOI_SEP <- rep(0, processes)
  
  cat("###############################################\n")
  
  for(proc in 1:processes){
    cat(paste0("########## PROCESS ", proc, "/", processes, " FOR N = ", N, " ##########\n"))
    
    # Time for One Process
    t_proc <- proc.time()
    
    ### DATA -------------------------------------------------------------------
    # Read the Data
    if(processes == 1){
      observations <- dataset[[proc]]
    } else{
      observations <- dataset[[proc]]$data
    }
    
    # Locations
    locations <- cbind(observations$x, observations$y)
    locations <- projection.points.1.5D(mesh, locations)
    
    # Times
    times <- observations$t
    
    ## STDE-PDE: Spatio-Temporal Intensity Estimation with PDE Regulariz. ------
    t0 <- proc.time()
    
    # Smoothing Parameters
    lambda <- 0.05
    lambda_time <- 0.0125
    
    # Solution
    # [If lambda and/or lambda_time are vectors, to select the best proposals by
    # 10-folds CV, please specify: preprocess_method = "RightCV"]
    invisible(capture.output(solution_STDEPDE <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                             mesh_time = mesh_time, lambda = lambda, tol1 = 1e-7,
                                                             lambda_time = lambda_time, fvec = NULL, heatStep = 0.1,
                                                             heatIter = 10, print = FALSE, nfolds = 10, nsimulations = 10000,
                                                             step_method = "Wolfe_Method", direction_method = "L-BFGS5",
                                                             preprocess_method = "NoCrossValidation")))
    
    FEMfunction_STDEPDE <- FEM.time(solution_STDEPDE$g, mesh_time, FEMbasis,
                                    FLAG_PARABOLIC = FALSE)
    
    # CPU Time
    CPUtime <- as.numeric(proc.time() - t0)[3]
    CPUtimes_STDEPDE[proc] <- CPUtime
    
    cat(paste0("STDE-PDE done in ", round(CPUtime, 2), " seconds.\n"))
    
    ## STLNPP: Spatio-Temporal Kernel Intensity Estimation ---------------------
    t0 <- proc.time()
    
    observations_stlpp <- as.stlpp(locations[,1], locations[,2], times, spatstat.linnet)
    Xs <- as.lpp.stlpp(x = observations_stlpp)
    Xt <- lpp(X=cbind(times,rep(0,n)), 
              L=linnet_interval(startp=0, endp=1))
    
    solution_STLNPP <- density(observations_stlpp,
                               lbw = as.numeric(bw.lppl(X=Xs)),
                               tbw = as.numeric(bw.lppl(X=Xt)),
                               dimyx = n, dimt = n,
                               diggle = TRUE, at = "pixels", verbose = FALSE)
    
    tgrid_STLNPP <- attr(solution_STLNPP, "tempden")$x
    
    # CPU Time
    CPUtime <- as.numeric(proc.time() - t0)[3]
    CPUtimes_STLNPP[proc] <- CPUtime
    
    cat(paste0("STLNPP done in ", round(CPUtime, 2), " seconds.\n"))
    
    ## HEAT: Spatio-Temporal Separable Kernel Intensity Estimation based on the Heat Equation -----
    t0 <- proc.time()

    observations_stlpp <- as.stlpp(locations[,1], locations[,2], times, spatstat.linnet)
    solution_HEAT <- density1Dkernel(observations_stlpp, at = "pixels",
                                     leaveoneout = FALSE, iterMax = 1e+06,
                                     verbose = FALSE)

    tgrid_HEAT <- attr(solution_HEAT, "tgrid")

    # CPU Time
    CPUtime <- as.numeric(proc.time() - t0)[3]
    CPUtimes_HEAT[proc] <- CPUtime

    cat(paste0("HEAT done in ", round(CPUtime, 2), " seconds.\n"))

    ## EQUAL-SPLIT: Spatio-Temporal Separable Kernel Intensity Estimation based on the the Okabe-Sugihara Equal-Split Algorithm -----
    t0 <- proc.time()

    observations_stlpp <- as.stlpp(locations[,1], locations[,2], times, spatstat.linnet)
    solution_EQUAL_SPLIT <- density1Dkernel(observations_stlpp, at = "pixels",
                                            leaveoneout = FALSE, iterMax = 1e+06,
                                            kernel = "epanechnikov", verbose = FALSE)

    tgrid_EQUAL_SPLIT <- attr(solution_EQUAL_SPLIT, "tgrid")

    # CPU Time
    CPUtime <- as.numeric(proc.time() - t0)[3]
    CPUtimes_EQUAL_SPLIT[proc] <- CPUtime

    cat(paste0("EQUAL-SPLIT done in ", round(CPUtime, 2), " seconds.\n"))
    
    ## VORONOI: Spatio-Temporal Pseudo-Separable Intensity Estimation based on the Voronoi-Dirichlet Tessellation -----
    t0 <- proc.time()
    
    observations_stlpp <- as.stlpp(locations[,1], locations[,2], times, spatstat.linnet)
    invisible(capture.output(solution_VORONOI <- densityVoronoi(observations_stlpp, f = 0.9, nrep = 10, separable = FALSE,
                                                                dimt = n, at = "pixels")))
    
    tgrid_VORONOI <- attr(solution_VORONOI, "tgrid")
    
    # CPU Time
    CPUtime <- as.numeric(proc.time() - t0)[3]
    CPUtimes_VORONOI[proc] <- CPUtime
    
    cat(paste0("VORONOI done in ", round(CPUtime, 2), " seconds.\n"))
    
    ## VORONOI-SEP: Spatio-Temporal Separable Intensity Estimation based on the Voronoi-Dirichlet Tessellation -----
    t0 <- proc.time()
    
    observations_stlpp <- as.stlpp(locations[,1], locations[,2], times, spatstat.linnet)
    invisible(capture.output(solution_VORONOI_SEP <- densityVoronoi(observations_stlpp, f = 1, nrep = 10, separable = TRUE,
                                                                    dimt = n, at = "pixels")))
    
    tgrid_VORONOI_SEP <- attr(solution_VORONOI_SEP, "tgrid")
    
    # CPU Time
    CPUtime <- as.numeric(proc.time() - t0)[3]
    CPUtimes_VORONOI_SEP[proc] <- CPUtime
    
    cat(paste0("VORONOI-SEP done in ", round(CPUtime, 2), " seconds.\n"))
    
    ## ERRORS ------------------------------------------------------------------
    cat(paste0("Computing the errors...\n"))
    
    # Time instants at which the solution is evaluated
    t <- seq(from = Ti, to = Tf, length.out = 31)
    
    # Mean Solutions over the Processes
    if (proc == 1){
      
      mean_sol_STDEPDE <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
      mean_sol_STLNPP <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
      mean_sol_HEAT <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
      mean_sol_EQUAL_SPLIT <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
      mean_sol_VORONOI <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
      mean_sol_VORONOI_SEP <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
      
      true_intensity <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
      
    }
    
    for(time_index in 1:length(t)) {
      
      # True Density
      if(proc == 1){
        
        true_intensity[,time_index] <- inhomogeneous_intensity4(mesh.eval$nodes[,1], mesh.eval$nodes[,2], t[time_index])
        
      }
      
      # STDE-PDE
      evaluation_STDEPDE <- eval.FEM.time(FEM.time = FEMfunction_STDEPDE, locations = mesh.eval$nodes, time.instants = t[time_index])
      evaluation_STDEPDE <- exp(evaluation_STDEPDE)*nrow(locations)
      
      # STLNPP
      idx <- which.min(abs(tgrid_STLNPP - t[time_index]))
      evaluation_STLNPP <- as.linfun(solution_STLNPP[[idx]])(mesh.eval$nodes[,1], mesh.eval$nodes[,2])
      
      # HEAT
      idx <- which.min(abs(tgrid_HEAT - t[time_index]))
      evaluation_HEAT <- as.linfun(solution_HEAT[[idx]])(mesh.eval$nodes[,1], mesh.eval$nodes[,2])

      # EQUAL-SPLIT
      idx <- which.min(abs(tgrid_EQUAL_SPLIT - t[time_index]))
      evaluation_EQUAL_SPLIT <- as.linfun(solution_EQUAL_SPLIT[[idx]])(mesh.eval$nodes[,1], mesh.eval$nodes[,2])
      
      # VORONOI
      idx <- which.min(abs(tgrid_VORONOI - t[time_index]))
      evaluation_VORONOI <- as.linfun(solution_VORONOI[[idx]])(mesh.eval$nodes[,1], mesh.eval$nodes[,2])
      
      # VORONOI-SEP
      idx <- which.min(abs(tgrid_VORONOI_SEP - t[time_index]))
      evaluation_VORONOI_SEP <- as.linfun(solution_VORONOI_SEP[[idx]])(mesh.eval$nodes[,1], mesh.eval$nodes[,2])
      
      # Estimates and Errors Updates
      
      # STDE-PDE
      #mise_STDEPDE[proc] <- mise_STDEPDE[proc] + domain_area*mean((evaluation_STDEPDE - true_intensity[,time_index])^2, na.rm = TRUE)/length(t)
      mise_STDEPDE[proc] <- mise_STDEPDE[proc] + compute_err(FEMbasis.eval, evaluation_STDEPDE, true_intensity[,time_index]) / length(t)
      mean_sol_STDEPDE[,time_index] <- mapply(sum, mean_sol_STDEPDE[,time_index], evaluation_STDEPDE, na.rm = TRUE)
      
      # STLNPP
      #mise_STLNPP[proc] <- mise_STLNPP[proc] + domain_area*mean((evaluation_STLNPP - true_intensity[,time_index])^2, na.rm = TRUE)/length(t)
      mise_STLNPP[proc] <- mise_STLNPP[proc] + compute_err(FEMbasis.eval, evaluation_STLNPP, true_intensity[,time_index]) / length(t)
      mean_sol_STLNPP[,time_index] <- mapply(sum, mean_sol_STLNPP[,time_index], evaluation_STLNPP, na.rm = TRUE)
      
      # HEAT
      #mise_HEAT[proc] <- mise_HEAT[proc] + domain_area*mean((evaluation_HEAT - true_intensity[,time_index])^2, na.rm = TRUE)/length(t)
      mise_HEAT[proc] <- mise_HEAT[proc] + compute_err(FEMbasis.eval, evaluation_HEAT, true_intensity[,time_index]) / length(t)
      mean_sol_HEAT[,time_index] <- mapply(sum, mean_sol_HEAT[,time_index], evaluation_HEAT, na.rm = TRUE)

      # EQUAL-SPLIT
      #mise_EQUAL_SPLIT[proc] <- mise_EQUAL_SPLIT[proc] + domain_area*mean((evaluation_EQUAL_SPLIT - true_intensity[,time_index])^2, na.rm = TRUE)/length(t)
      mise_EQUAL_SPLIT[proc] <- mise_EQUAL_SPLIT[proc] + compute_err(FEMbasis.eval, evaluation_EQUAL_SPLIT, true_intensity[,time_index]) / length(t)
      mean_sol_EQUAL_SPLIT[,time_index] <- mapply(sum, mean_sol_EQUAL_SPLIT[,time_index], evaluation_EQUAL_SPLIT, na.rm = TRUE)
      
      # VORONOI
      #mise_VORONOI[proc] <- mise_VORONOI[proc] + domain_area*mean((evaluation_VORONOI - true_intensity[,time_index])^2, na.rm = TRUE)/length(t)
      mise_VORONOI[proc] <- mise_VORONOI[proc] + compute_err(FEMbasis.eval, evaluation_VORONOI, true_intensity[,time_index]) / length(t)
      mean_sol_VORONOI[,time_index] <- mapply(sum, mean_sol_VORONOI[,time_index], evaluation_VORONOI, na.rm = TRUE)
      
      # VORONOI-SEP
      #mise_VORONOI_SEP[proc] <- mise_VORONOI_SEP[proc] + domain_area*mean((evaluation_VORONOI_SEP - true_intensity[,time_index])^2, na.rm = TRUE)/length(t)
      mise_VORONOI_SEP[proc] <- mise_VORONOI_SEP[proc] + compute_err(FEMbasis.eval, evaluation_VORONOI_SEP, true_intensity[,time_index]) / length(t)
      mean_sol_VORONOI_SEP[,time_index] <- mapply(sum, mean_sol_VORONOI_SEP[,time_index], evaluation_VORONOI_SEP, na.rm = TRUE)
      
    } 
    
    cat(paste("Process", proc, "done in", round(as.numeric(proc.time() - t_proc)[3]), "seconds.\n"))
    cat("###############################################\n")
    
  }
  
  mean_sol_STDEPDE <- mean_sol_STDEPDE/processes
  mean_sol_STLNPP <- mean_sol_STLNPP/processes
  mean_sol_HEAT <- mean_sol_HEAT/processes
  mean_sol_EQUAL_SPLIT <- mean_sol_EQUAL_SPLIT/processes
  mean_sol_VORONOI <- mean_sol_VORONOI/processes
  mean_sol_VORONOI_SEP <- mean_sol_VORONOI_SEP/processes
  
  #'##########################################################################'#
  #'##########################################################################'#
  
  cat("[RESULTS]\n")
  
  # Create the folder where to save the results
  current_datetime <- Sys.time()
  formatted_datetime <- format(current_datetime, "[%Y-%m-%d]-[%H-%M-%S]")
  
  path <- paste0(getwd(),"/results-",formatted_datetime)
  dir.create(path, showWarnings = FALSE)
  setwd(path)
  
  ## EXPORT RESULTS ------------------------------------------------------------
  cat("Saving estimates, errors and CPU times...\n")
  
  # Export estimates and errors
  dir.create(file.path(getwd(), "output"), showWarnings = FALSE)
  
  estimates <- list(mean_sol_STDEPDE, mean_sol_STLNPP, mean_sol_HEAT, mean_sol_EQUAL_SPLIT,
                    mean_sol_VORONOI, mean_sol_VORONOI_SEP)
  save(estimates, file = paste0("output/estimates_", N, "data.RData"))
  
  mise <- list(mise_STDEPDE, mise_STLNPP, mise_HEAT, mise_EQUAL_SPLIT,
               mise_VORONOI, mise_VORONOI_SEP)
  save(mise, file = paste0("output/mise_", N, "data.RData"))
  
  CPUtimes <- list(CPUtimes_STDEPDE, CPUtimes_STLNPP, CPUtimes_HEAT, CPUtimes_EQUAL_SPLIT,
                   CPUtimes_VORONOI, CPUtimes_VORONOI_SEP)
  save(CPUtimes, file = paste0("output/CPUtimes_", N, "data.RData"))
  
  #base::save.image(paste0("output/workspace_", N, "data.RData"))
  
  #'##########################################################################'#
  #'##########################################################################'#
  
  cat("Saving pictures...\n")
  
  ## VISUALIZATION -------------------------------------------------------------
  {
    ### ESTIMATES --------------------------------------------------------------
    # Import estimates and errors
    # load(file = paste0("output/estimates_", N, "data.rda"))
    # mean_sol_STDEPDE <- estimates[[1]]
    # mean_sol_STLNPP <- estimates[[2]]
    # mean_sol_HEAT <- estimates[[3]]
    # mean_sol_EQUAL_SPLIT <- estimates[[4]]
    # mean_sol_VORONOI <- estimates[[5]]
    # mean_sol_VORONOI_SEP <- estimates[[6]]
    
    # Time instants at which the solution is evaluated
    t <- seq(from = Ti, to = Tf, length.out = 31)
    
    # Time Half Width
    h <- (mesh_time[2]-mesh_time[1])/2
    
    # Plots
    m <- min(min(true_intensity, na.rm = TRUE),
             min(mean_sol_STDEPDE, na.rm = TRUE),
             min(mean_sol_STLNPP, na.rm = TRUE),
             min(mean_sol_HEAT, na.rm = TRUE),
             min(mean_sol_EQUAL_SPLIT, na.rm = TRUE),
             # min(mean_sol_VORONOI, na.rm = TRUE),
             # min(mean_sol_VORONOI_SEP, na.rm = TRUE),
             na.rm = TRUE)
    
    M <- max(max(true_intensity, na.rm = TRUE),
             max(mean_sol_STDEPDE, na.rm = TRUE),
             max(mean_sol_STLNPP, na.rm = TRUE),
             max(mean_sol_HEAT, na.rm = TRUE),
             max(mean_sol_EQUAL_SPLIT, na.rm = TRUE),
             # max(mean_sol_VORONOI, na.rm = TRUE),
             # max(mean_sol_VORONOI_SEP, na.rm = TRUE),
             na.rm = TRUE)
    
    dir.create(file.path(getwd(), "pictures"), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data")), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data/TIME_BAR")), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data/SAMPLE")), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data/TRUE_INTENSITY")), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data/STDEPDE")), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data/STLNPP")), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data/HEAT")), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data/EQUAL_SPLIT")), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data/VORONOI")), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data/VORONOI_SEP")), showWarnings = FALSE)
    
    color_palette <- c("jet.col", "viridis", "magma", "plasma", "inferno")
    
    pb <- txtProgressBar(min = 1, max = length(t), style = 3, char = "=")
    
    for(time_index in 1:length(t)){
      
      # First Plot: Sample
      idx <- which(abs(times - t[time_index]) < h)
      plot(mesh, linewidth = 1) + geom_point(data = data.frame(x = locations[idx,1], y = locations[idx,2]),
                                             aes(x = x, y = y), color = "red3", size = 5)
      plotname <- paste0("pictures/", N, "_data/SAMPLE/sim5_SAMPLE_", N, "data_", time_index, ".pdf")
      ggsave(filename = plotname, plot = last_plot(), width = 8, height = 8)
    
      for(col in color_palette) {
        
        # Second Plot: True Density at t[time_index]
        plotname <- paste0("pictures/", N, "_data/TRUE_INTENSITY/sim5_TRUE_INTENSITY_", N, "data_", col, "_", time_index)
        plot(FEM(true_intensity[,time_index], FEMbasis.eval), m = m, M = M,
             colormap = col, filename = plotname, showLegend = FALSE, linewidth = 2.5)
      
        # Third Plot: STDE-PDE Estimated Density at t[time_index]
        plotname <- paste0("pictures/", N, "_data/STDEPDE/sim5_STDEPDE_", N, "data_", col, "_", time_index)
        plot(FEM(mean_sol_STDEPDE[,time_index], FEMbasis.eval), m = m, M = M,
             colormap = col, filename = plotname, showLegend = FALSE, linewidth = 2.5)
        
        # Fourth Plot: STLNPP Estimated Density at t[time_index]
        plotname <- paste0("pictures/", N, "_data/STLNPP/sim5_STLNPP_", N, "data_", col, "_", time_index)
        plot(FEM(mean_sol_STLNPP[,time_index], FEMbasis.eval), m = m, M = M,
             colormap = col, filename = plotname, showLegend = FALSE, linewidth = 2.5)
        
        # Fifth Plot: HEAT Estimated Density at t[time_index]
        plotname <- paste0("pictures/", N, "_data/HEAT/sim5_HEAT_", N, "data_", col, "_", time_index)
        plot(FEM(mean_sol_HEAT[,time_index], FEMbasis.eval), m = m, M = M,
             colormap = col, filename = plotname, showLegend = FALSE, linewidth = 2.5)
        
        # Sixth Plot: EQUAL-SPLIT Estimated Density at t[time_index]
        plotname <- paste0("pictures/", N, "_data/EQUAL_SPLIT/sim5_EQUAL_SPLIT_", N, "data_", col, "_", time_index)
        plot(FEM(mean_sol_EQUAL_SPLIT[,time_index], FEMbasis.eval), m = m, M = M,
             colormap = col, filename = plotname, showLegend = FALSE, linewidth = 2.5)
        
        # Seventh Plot: VORONOI Estimated Density at t[time_index]
        plotname <- paste0("pictures/", N, "_data/VORONOI/sim5_VORONOI_", N, "data_", col, "_", time_index)
        plot(FEM(mean_sol_VORONOI[,time_index], FEMbasis.eval), m = m, M = M,
             colormap = col, filename = plotname, showLegend = FALSE, linewidth = 2.5)
        
        # Eight Plot: VORONOI-SEP Estimated Density at t[time_index]
        plotname <- paste0("pictures/", N, "_data/VORONOI_SEP/sim5_VORONOI_SEP_", N, "data_", col, "_", time_index)
        plot(FEM(mean_sol_VORONOI_SEP[,time_index], FEMbasis.eval), m = m, M = M,
             colormap = col, filename = plotname, showLegend = FALSE, linewidth = 2.5)
      
      }
      
      setTxtProgressBar(pb, time_index)
      
    }
    
    cat("\n")
    
    ## LEGEND ------------------------------------------------------------------
    for(col in color_palette) {
      
      pdf(paste0("pictures/", N, "_data/sim5_legend_", col, ".pdf"), family = "serif", width = 11, height = 3)
      #x11(width = 11, height = 3)
      par(mai = c(1,0.75,0,0))
      plot(c(0, 112.5), c(0, 15), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame.plot = F)
      gradient.rect(0, 0, 100, 5, col = match.fun(col)(100), border = "black")
      axis(1, at = c(0,25,50,75,100), labels = c("0", "0.25", "0.50", "0.75", "1"),
           lwd.ticks = 2, cex.axis = 2, lwd = 2)
      text(107, 2, TeX("$\\times 10^{2}$"), cex = 2)
      
      dev.off()
      
    }
    
    cat(paste0("m = ", round(m, 4), "; M = ", round(M, 4), "\n"))
    
    ## ERRORS ------------------------------------------------------------------
    {
      blue_polimi <- rgb(114, 143, 165, max = 255, alpha = 65, names = "blue_polimi25") # #728FA5
      
      pdf(paste0("pictures/", N, "_data/sim5_errors.pdf"), family = "serif", width = 11.75, height = 5.25)
      #x11(width = 11.75, height = 5.25)
      par(mfrow = c(1,3), mai = c(1.75,0.75,0.5,0.15))
      plot.new()
      boxplot(mise_STDEPDE, mise_EQUAL_SPLIT, mise_HEAT, mise_STLNPP,
              mise_VORONOI, mise_VORONOI_SEP,
              asp = 1, xaxt = "n", yaxt = "n")
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = blue_polimi)
      grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 1)
      par(new = TRUE)
      boxplot(mise_STDEPDE, mise_EQUAL_SPLIT, mise_HEAT, mise_STLNPP,
              mise_VORONOI, mise_VORONOI_SEP,
              names = c("STDE-PDE", "STKDE-EPAN", "STKDE-HEAT", "STKDE-QUICK", "STVDE", "STVDE-SEP"),
              ylab = TeX("$err_{L^2}$", bold = T),
              main = "Simulation 4",
              col = brewer.pal(6, "YlGnBu"), cex.lab = 2, cex.axis = 2, cex.main = 2, las = 3)
      plot.new()
      dev.off()
    }
    
    {
      blue_polimi <- rgb(114, 143, 165, max = 255, alpha = 65, names = "blue_polimi25") # #728FA5
      
      pdf(paste0("pictures/", N, "_data/sim5_errors_zoom.pdf"), family = "serif", width = 11.75, height = 5.25)
      #x11(width = 11.75, height = 5.25)
      par(mfrow = c(1,3), mai = c(1.75,0.75,0.5,0.15))
      plot.new()
      boxplot(mise_STDEPDE, mise_EQUAL_SPLIT, mise_HEAT, mise_STLNPP,
              asp = 1, xaxt = "n", yaxt = "n")
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = blue_polimi)
      grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 1)
      par(new = TRUE)
      boxplot(mise_STDEPDE, mise_EQUAL_SPLIT, mise_HEAT, mise_STLNPP,
              names = c("STDE-PDE", "STKDE-EPAN", "STKDE-HEAT", "STKDE-QUICK"),
              ylab = TeX("$err_{L^2}$", bold = T),
              main = "Simulation 4",
              col = brewer.pal(6, "YlGnBu")[1:4], cex.lab = 2, cex.axis = 2, cex.main = 2, las = 3)
      plot.new()
      dev.off()
    }
    
    ## CPU TIMES ---------------------------------------------------------------
    pdf(paste0("pictures/", N, "_data/sim5_CPUtimes.pdf"), family = "serif", width = 11.75, height = 5.25)
    #x11(width = 11.75, height = 5.25)
    par(mfrow = c(1,3), mai = c(1.75,0.75,0.5,0.15))
    plot.new()
    boxplot(CPUtimes_STDEPDE, CPUtimes_EQUAL_SPLIT, CPUtimes_HEAT, CPUtimes_STLNPP,
            CPUtimes_VORONOI, CPUtimes_VORONOI_SEP,
            asp = 1, xaxt = "n", yaxt = "n")
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = blue_polimi)
    grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 1)
    par(new = TRUE)
    boxplot(CPUtimes_STDEPDE, CPUtimes_EQUAL_SPLIT, CPUtimes_HEAT, CPUtimes_STLNPP,
            CPUtimes_VORONOI, CPUtimes_VORONOI_SEP,
            names = c("STDE-PDE", "STKDE-EPAN", "STKDE-HEAT", "STKDE-QUICK", "STVDE", "STVDE-SEP"),
            yalb = "CPU Times [seconds]",
            main = "Simulation 4",
            col = brewer.pal(6, "YlGnBu"), cex.lab = 2, cex.axis = 2, cex.main = 2, las = 3)
    dev.off()
    
    ## MESH --------------------------------------------------------------------
    plot(mesh, linewidth = 1) + geom_point(data = data.frame(x = mesh$nodes[,1], y = mesh$nodes[,2]),
                                           aes(x = x, y = y), color = "red3", size = 2)
    ggsave(filename = paste0("pictures/", N, "_data/sim5_mesh.pdf"),
           plot = last_plot(), width = 8, height = 8)
    
    plot(mesh.eval, linewidth = 1) + geom_point(data = data.frame(x = mesh.eval$nodes[,1], y = mesh.eval$nodes[,2]),
                                                aes(x = x, y = y), color = "red3", size = 2)
    ggsave(filename = paste0("pictures/", N, "_data/sim5_mesh_eval.pdf"),
           plot = last_plot(), width = 8, height = 8)
    
    ## TIME BAR ----------------------------------------------------------------
    for(time_index in 1:length(t)) {
      string <- paste0(paste0("pictures/", N, "_data/TIME_BAR/sim5_time_bar_",time_index,".pdf"))
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
    
  }
  
  cat("###############################################\n")
  
}


## SMOOTHING PARAMETERS SELECTION ----------------------------------------------
cat("Selecting smoothing parameters...\n")
options(OutDec = ",")

lambda_space <- c(10, 1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5)
lambda_time <- c(10, 1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5)

# lambda_space <- c(0.15, 0.125, 0.075, 0.05, 0.025, 0.01, 0.005)
# lambda_time <- c(0.02, 0.015, 0.0125, 0.0075, 0.005, 0.0025, 0.0001)

for(i in 1:length(lambda_space)){
  for(j in 1:length(lambda_time)){
    
    # Solution
    # [If lambda and/or lambda_time are vectors, to select the best proposals by
    # 10-folds CV, please specify: preprocess_method = "RightCV"]
    invisible(capture.output(solution_STDEPDE <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                             mesh_time = mesh_time, lambda = lambda_space[i], tol1 = 1e-4,
                                                             lambda_time = lambda_time[j], fvec = NULL, heatStep = 0.1,
                                                             heatIter = 10, print = TRUE, nfolds = 10, nsimulations = 5000,
                                                             step_method = "Wolfe_Method", direction_method = "L-BFGS5",
                                                             preprocess_method = "NoCrossValidation")))
    
    FEMfunction_STDEPDE <- FEM.time(solution_STDEPDE$g, mesh_time, FEMbasis,
                                    FLAG_PARABOLIC = FALSE)
    
    # Time instants at which the solution is evaluated
    t <- seq(from = Ti, to = Tf, length.out = 31)
    
    mise_STDEPDE <- 0
    mean_sol_STDEPDE <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
    true_intensity <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
    
    for(time_index in 1:length(t)) {
      
      true_intensity[,time_index] <- inhomogeneous_intensity4(mesh.eval$nodes[,1], mesh.eval$nodes[,2], t[time_index])
      
      # STDE-PDE
      evaluation_STDEPDE <- eval.FEM.time(FEM.time = FEMfunction_STDEPDE, locations = mesh.eval$nodes, time.instants = t[time_index])
      evaluation_STDEPDE <- exp(evaluation_STDEPDE)*nrow(locations)
      
      # Estimates and Errors Updates
      
      # STDE-PDE
      #mise_STDEPDE <- mise_STDEPDE + domain_area*mean((evaluation_STDEPDE - true_intensity[,time_index])^2, na.rm = TRUE)/length(t)
      mise_STDEPDE <- mise_STDEPDE + compute_err(FEMbasis.eval, evaluation_STDEPDE, true_intensity[,time_index]) / length(t)
      mean_sol_STDEPDE[,time_index] <- mapply(sum, mean_sol_STDEPDE[,time_index], evaluation_STDEPDE, na.rm = TRUE)
      
    }
    
    cat(paste0("lambda_space = ", lambda_space[i], "; lambda_time = ", lambda_time[j],
               "; L2 error = ", round(mise_STDEPDE, digits = 5), "\n"))
  }
}
