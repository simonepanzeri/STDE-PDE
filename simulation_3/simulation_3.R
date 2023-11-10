#'############################################################################'#
#'##################### SPATIO-TEMPORAL DENSITY ESTIMATION ###################'# 
#'########### SIMULATION 3 - MIXTURE OF DENSITIES OVER 2.5D DOMAIN ###########'#
#'############################################################################'#

## LIBRARIES AND FUNCTIONS -----------------------------------------------------
source("libraries_3.R")

## SPATIO-TEMPORAL DISCRETIZATION ----------------------------------------------
# Domain Area
domain_area <- MeshSurface(mesh)

# Spatial 2.5D Mesh over the Domain for Estimation
vertices <- read.table("mesh/simulation3.surface.vertices.txt")
triangles <- read.table("mesh/simulation3.surface.triangles.txt")
mesh <- create.mesh.2.5D(nodes = vertices, triangles = triangles[,1:3])
FEMbasis <- create.FEM.basis(mesh)

plot.mesh(mesh)
rgl.close()

# Fine Spatial Mesh over the Domain for Evaluation
mesh.eval <- refine.by.splitting.mesh.2.5D(mesh)
FEMbasis.eval <- create.FEM.basis(mesh.eval)
vertices.eval <- mesh.eval$nodes
vertices.eval.proj <- vertices.eval

plot.mesh(mesh.eval)
rgl.close()

# Super Fine Spatial Mesh over the Domain for Visualization
vertices.visual <- read.table("mesh/fine.surface.vertices.txt")
triangles.visual <- read.table("mesh/fine.surface.triangles.txt")
mesh.visual <- create.mesh.2.5D(nodes = vertices.visual[,1:3], triangles = triangles.visual[,1:3])
FEMbasis.visual <- create.FEM.basis(mesh.visual)
vertices.visual.proj <- mesh.visual$nodes

plot.mesh(mesh.visual)
rgl.close()

# Temporal 1D Mesh over [0,1]
mesh_time <- seq(from = 0, to = 1, length.out = 6)

## ESTIMATION PROCEDURE --------------------------------------------------------
# Number of processes
processes <- 30

# Sample Size
NN <- 10000
# alternative: NN <- c(1000, 2500, 5000, 10000)

for(N in NN){
  
  # Containers for Errors
  mise_STDEPDE <- rep(0, processes)
  mise_KNNSTDE_separable <- rep(0, processes)
  
  KLdist_STDEPDE <- rep(0, processes)
  KLdist_KNNSTDE_separable <- rep(0, processes)
  
  Hdist_STDEPDE <- rep(0, processes)
  Hdist_KNNSTDE_separable <- rep(0, processes)
  
  # Containers for CPU Times
  CPUtimes_STDEPDE <- rep(0, processes)
  CPUtimes_KNNSTDE_separable <- rep(0, processes)
  
  for(proc in 1:processes){
    # Time for One Process
    t_proc <- proc.time()
    
    ### DATA -------------------------------------------------------------------
    # Generate the data
    generate.data(N, proc) # run only once to generate the data
    
    # Read the Data
    data <- read.table(paste0("data/",N,"data_",proc,".txt"))
    
    # Locations
    locations <- data[,1:3]
    
    # Times
    times <- data[,4]

    ## STDE-PDE: Spatio-Temporal Density Estimation with PDE Regulariz. --------
    t0 <- proc.time()
    
    # Smoothing Parameters
    lambda <- 0.4
    # alternative: lambda <- 10^seq(from = -1, to = 0, by = 0.5)
    
    lambda_time <- 0.4
    # alternative: lambda_time <- 10^seq(from = -1, to = 0, by = 0.5)
    
    # Solution
    # [If lambda and/or lambda_time are vectors, to select the best proposals by
    # 10-folds CV, please specify: preprocess_method = "RightCV"]
    solution_STDEPDE <- DE.FEM.time(data = locations, data_time = times,
                                    FEMbasis = FEMbasis, mesh_time = mesh_time,
                                    lambda = lambda, lambda_time = lambda_time,
                                    fvec = NULL, heatIter = 10, print = FALSE,
                                    direction_method = "L-BFGS5", tol1 = 1e-7,
                                    preprocess_method = "NoCrossValidation", nfolds = 10,
                                    nsimulations = 10000)
    
    FEMfunction_STDEPDE <- FEM.time(solution_STDEPDE$g, mesh_time, FEMbasis, FLAG_PARABOLIC = FALSE)
    
    # CPU Time
    CPUtime <- proc.time() - t0
    CPUtimes_STDEPDE[proc] <- as.numeric(CPUtime[3])
    
    ## KNNSTDE-separable: kNN Spatio-Temporal Density Estimation ---------------
    t0 <- proc.time()
    
    # Spatial Component
    f_space_KNNSTDE_separable <- knnDE(as.matrix(locations), as.matrix(vertices[,1:3]), 150)
    
    # Temporal Component
    f_time_KNNSTDE_separable <- knnDE(as.matrix(times), as.matrix(seq(from = 0, to = 1, length.out = 25)), 50)
  
    # CPU Time
    CPUtime <- proc.time() - t0
    CPUtimes_KNNSTDE_separable[proc] <- as.numeric(CPUtime[3])
    
    ## ERRORS ------------------------------------------------------------------
    # Time instants at which the solution is evaluated
    t <- seq(from = 0, to = 1, by = 0.05)
    
    # Temporal Mesh Step Size
    h <- (mesh_time[2] - mesh_time[1])/2
    
    # Mean Solutions over the Processes
    if (proc == 1){
      
      mean_sol_STDEPDE <- matrix(nrow = nrow(vertices.eval), ncol = length(t))
      mean_sol_KNNSTDE_separable <- matrix(nrow = nrow(vertices.eval), ncol = length(t))
      
      true_instant <-  matrix(nrow = nrow(vertices.eval), ncol = length(t))
      
    }
    
    for (time_index in 1:length(t)) {
      
      # True Density
      if(proc == 1){

        true_instant[,time_index] <- dens.func(mesh.eval$nodes) / sum(dens.func(mesh.eval$nodes))
        
      }
      
      # STDE-PDE
      evaluation_STDEPDE <- eval.FEM.time(FEM.time = FEMfunction_STDEPDE, locations = vertices.eval.proj, time.instants = t[time_index])
      evaluation_STDEPDE <- exp(evaluation_STDEPDE)
      evaluation_STDEPDE <- evaluation_STDEPDE / sum(evaluation_STDEPDE)
      
      # KNNSTDE-separable
      marginal_time_KNNSTDE_separable <- f_time_KNNSTDE_separable[findInterval(t[time_index], seq(from = 0, to = 1, length.out = 25))]
      solution_KNNSTDE_separable <- marginal_time_KNNSTDE_separable * f_space_KNNSTDE_separable
      solution_KNNSTDE_separable <- eval.FEM(FEM(solution_KNNSTDE_separable, FEMbasis), locations = vertices.eval.proj)
      solution_KNNSTDE_separable <- solution_KNNSTDE_separable / sum(solution_KNNSTDE_separable)
      
      # Estimates and Errors Updates
      
      # STDE-PDE
      mise_STDEPDE[proc] <- mise_STDEPDE[proc] + domain_area*mean((evaluation_STDEPDE - true_instant[,time_index])^2, na.rm = TRUE)/length(t)
      KLdist_STDEPDE[proc] <- KLdist_STDEPDE[proc] + abs(kullback_leibler_distance(evaluation_STDEPDE, true_instant[,time_index], testNA = FALSE, unit = "log", epsilon = 1e-4))/nrow(vertices.eval)
      Hdist_STDEPDE[proc] <- Hdist_STDEPDE[proc] + domain_area*hellinger.distance(evaluation_STDEPDE, true_instant[,time_index])/length(t)
      
      mean_sol_STDEPDE[,time_index] <- mapply(sum, mean_sol_STDEPDE[,time_index], evaluation_STDEPDE, na.rm = TRUE)
        
      # KNNSTDE-separable
      mise_KNNSTDE_separable[proc] <- mise_KNNSTDE_separable[proc] + domain_area*mean((solution_KNNSTDE_separable - true_instant[,time_index])^2, na.rm = TRUE)/length(t)
      KLdist_KNNSTDE_separable[proc] <- KLdist_KNNSTDE_separable[proc] + abs(kullback_leibler_distance(solution_KNNSTDE_separable, true_instant[,time_index], testNA = FALSE, unit = "log", epsilon = 1e-4))/nrow(vertices.eval)
      Hdist_KNNSTDE_separable[proc] <- Hdist_KNNSTDE_separable[proc] + domain_area*hellinger.distance(solution_KNNSTDE_separable, true_instant[,time_index])/length(t)
      
      mean_sol_KNNSTDE_separable[,time_index] <- mapply(sum, mean_sol_KNNSTDE_separable[,time_index], solution_KNNSTDE_separable, na.rm = TRUE)
      
    } 
    
    print(paste("Process", proc, "done in", round(as.numeric(proc.time() - t_proc)[3]), "seconds."))
    
  }
  
  mean_sol_STDEPDE <- mean_sol_STDEPDE/processes
  mean_sol_KNNSTDE_separable <- mean_sol_KNNSTDE_separable/processes
  
  #'##########################################################################'#
  #'##########################################################################'#
  
  ## EXPORT RESULTS ------------------------------------------------------------
  # Export estimates and errors
  dir.create(file.path(getwd(), "output"), showWarnings = FALSE)
  
  estimates <- list(mean_sol_STDEPDE, mean_sol_KNNSTDE_separable)
  save(estimates, file = paste0("output/estimates_", N, "data.rda"))
  
  mise <- list(mise_STDEPDE, mise_KNNSTDE_separable)
  save(mise, file = paste0("output/estimates_", N, "data.rda"))
  
  KLdist <- list(KLdist_STDEPDE, KLdist_KNNSTDE_separable)
  save(KLdist, file = paste0("output/KLdist_", N, "data.rda"))
  
  Hdist <- list(Hdist_STDEPDE, Hdist_KNNSTDE_separable)
  save(Hdist, file = paste0("output/Hdist_", N, "data.rda"))
  
  CPUtimes <- list(CPUtimes_STDEPDE, CPUtimes_KNNSTDE_separable)
  save(CPUtimes, file = paste0("output/CPUtimes_", N, "data.rda"))
  
  base::save.image(paste0("output/workspace_", N, "data.RData"))
  
  #'##########################################################################'#
  #'##########################################################################'#
  
  ## VISUALIZATION -------------------------------------------------------------
  {
    ### ESTIMATES [FIGURE F.12 and F.13] ---------------------------------------
    # Time Instants for Evaluation
    t <- seq(from = 0, to = 1, by = 0.05)
    
    # Temporal Mesh Step Size
    h <- (mesh_time[2] - mesh_time[1])/2
    
    # Plots
    m <- 0
    M <- max(max(true_instant, na.rm = TRUE),
             max(mean_sol_STDEPDE, na.rm = TRUE),
             max(mean_sol_KNNSTDE_separable, na.rm = TRUE),
             na.rm = TRUE)
    
    dir.create(file.path(getwd(), "pictures"), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data")), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data/sample")), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data/true_instant")), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data/STDEPDE")), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data/KNNSTDE-separable")), showWarnings = FALSE)
    dir.create(file.path(getwd(), paste0("pictures/", N, "_data/time_bar")), showWarnings = FALSE)
    
    for(time_index in 1:length(t)){
      
      t_img <- proc.time()
        
      # First Plot: Sample
      string <- paste0("pictures/", N, "_data/sample/sim3_sample_", t[time_index], ".png")
      plot.sample(FEM(rep(0,nrow(mesh.visual$nodes)),FEMbasis.visual), locations[abs(times-t[time_index])<h,])
      snapshot3d(string, fmt = "png", width = 600, height = 400, webshot = rgl.useNULL())
      rgl.close()
        
      # Second Plot: True Density at t[time_index]
      string <- paste0("pictures/", N, "_data/true_instant/sim3_true_instant_", t[time_index], ".png")
      plot.FEM(FEM(true_instant[,time_index], FEMbasis.eval), M)
      snapshot3d(string, fmt = "png", width = 600, height = 400, webshot = rgl.useNULL())
      rgl.close()
        
      # Third Plot: STDE-PDE Estimated Density at t[time_index]
      string <- paste0("pictures/", N, "_data/STDEPDE/sim3_STDEPDE_", t[time_index], ".png")
      plot.FEM(FEM(mean_sol_STDEPDE[,time_index], FEMbasis.eval), M)
      snapshot3d(string, fmt = "png", width = 600, height = 400, webshot = rgl.useNULL())
      rgl.close()
        
      # Fourth Plot: KNNSTDE Estimated Density at t[time_index]
      string <- paste0("pictures/", N, "_data/KNNSTDE-separable/sim3_KNNSTDE_separable_", t[time_index], ".png")
      plot.FEM(FEM(mean_sol_KNNSTDE_separable[,time_index], FEMbasis.eval), M)
      snapshot3d(string, fmt = "png", width = 600, height = 400, webshot = rgl.useNULL())
      rgl.close()
        
      print(paste("Images for t =", t[time_index], "done in", round(as.numeric(proc.time() - t_img)[3]), "seconds."))
        
    }
    
    ## LEGEND [FIGURE F.13] ----------------------------------------------------
    pdf(paste0("pictures/", N, "_data/sim3_legend.pdf"), family = "serif", width = 11, height = 3)
    par(mai = c(1,0.75,0,0))
    plot(c(0, 112.5), c(0, 15), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame.plot = F)
    gradient.rect(0, 0, 100, 5, col = jet.col(1000), border = "black")
    axis(1, at = c(0,33.33,66.66,100), labels = c("0", "0.5", "1.0", "1.5"),
         lwd.ticks = 2, cex.axis = 2, lwd = 2)
    text(107, 2, TeX("$\\times 10^{-3}$"), cex = 2)
    dev.off()
    
    print(paste0("M = ", round(M, 4)))
  
    ## ERRORS [FIGURE F.14] ----------------------------------------------------
    {
      
      blue_polimi <- rgb(114, 143, 165, max = 255, alpha = 65, names = "blue_polimi25") # #728FA5
      
      pdf(paste0("pictures/", N, "_data/sim3_errors.pdf"), family = "serif", width = 12, height = 5.3)
      #x11(width = 12, height = 5.3)
      par(mfrow = c(1,3), mai = c(1.75,0.35,0.5,0.15))
      boxplot(mise_STDEPDE, mise_KNNSTDE_separable, asp = 1, xaxt = "n", yaxt = "n")
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = blue_polimi)
      grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 1)
      par(new = TRUE)
      boxplot(mise_STDEPDE, mise_KNNSTDE_separable,
              names = c("STDE-PDE", "KNNSTDE-sep"), main = TeX("$err_{L^2}$", bold = T),
              col = brewer.pal(8, "YlGnBu")[c(1,8)], cex.lab = 2, cex.axis = 2, cex.main = 2, las = 3)
      axis(side = 2, at = c(0.0005,0.0015,0.0025), cex.lab = 2, cex.axis = 2)
      boxplot(KLdist_STDEPDE, KLdist_KNNSTDE_separable, asp = 1, xaxt="n", yaxt = "n")
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = blue_polimi)
      grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 1)
      par(new = TRUE)
      boxplot(KLdist_STDEPDE, KLdist_KNNSTDE_separable,
              names=c("STDE-PDE", "KNNSTDE-sep"), main = "KL",
              col = brewer.pal(8, "YlGnBu")[c(1,8)], cex.lab = 2, cex.axis = 2, cex.main = 2, las = 3)
      axis(side = 2, at = c(0,0.0010,0.0020), cex.lab = 2, cex.axis = 2)
      boxplot(Hdist_STDEPDE, Hdist_KNNSTDE_separable, asp = 1, xaxt = "n", yaxt = "n")
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = blue_polimi)
      grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 1)
      par(new = TRUE)
      boxplot(Hdist_STDEPDE, Hdist_KNNSTDE_separable,
              names = c("STDE-PDE", "KNNSTDE-sep"), main = "H",
              col = brewer.pal(8, "YlGnBu")[c(1,8)], cex.lab = 2, cex.axis = 2, cex.main = 2, las = 3)
      dev.off()
      
    }
    
    ## CPU TIMES ---------------------------------------------------------------
    pdf(paste0("pictures/", N, "_data/sim3_CPUtimes.pdf"), family = "serif", width = 12, height = 5.3)
    #x11(width = 12, height = 5.3)
    par(mfrow = c(1,3), mai = c(1.75,0.35,0.5,0.15))
    plot.new()
    boxplot(CPUtimes_STDEPDE, CPUtimes_KNNSTDE_separable, asp = 1, xaxt = "n", yaxt = "n")
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = blue_polimi)
    grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 1)
    par(new = TRUE)
    boxplot(CPUtimes_STDEPDE, CPUtimes_KNNSTDE_separable,
            names = c("STDE-PDE", "KNNSTDE-sep"), main = "CPU TIMES [seconds]",
            col = brewer.pal(8, "YlGnBu")[c(1,8)], cex.lab = 2, cex.axis = 2, cex.main = 2, las = 3)
    dev.off()
      
    ## MESH --------------------------------------------------------------------
    string <- paste0("pictures/", N, "_data/sim3_mesh.png")
    plot.mesh(mesh)
    snapshot3d(string, fmt = "png", width = 600, height = 400, webshot = rgl.useNULL())
    rgl.close()
    
    string <- paste0("pictures/", N, "_data/sim3_mesh.eval.png")
    plot.mesh(mesh.eval)
    snapshot3d(string, fmt = "png", width = 600, height = 400, webshot = rgl.useNULL())
    rgl.close()
      
    ## DOMAIN ------------------------------------------------------------------
    string <- paste0("pictures/", N, "_data/sim3_domain.png")
    plot.FEM(FEM(rep(0, nrow(mesh.visual$nodes)), FEMbasis.visual))
    snapshot3d(string, fmt = "png", width = 600, height = 400, webshot = rgl.useNULL())
    rgl.close()
      
    ## TIME BAR ----------------------------------------------------------------
    for(time_index in 1:length(t)){
      
      string <- paste0("pictures/", N, "_data/time_bar/sim3_time_bar", time_index, ".pdf")
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
