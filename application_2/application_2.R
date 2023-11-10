#'############################################################################'#
#'##################### SPATIO-TEMPORAL DENSITY ESTIMATION ###################'# 
#'######################## APPLICATION 2 - EARTHQUAKES #######################'#
#'############################################################################'#

## LIBRARIES AND FUNCTIONS -----------------------------------------------------
source("libraries_2.R")

## IMPORT DATA -----------------------------------------------------------------
# Load the Data (Earthquakes with Magnitude > 4.5)
earthquakes_df <- read.csv(file = "data/query_4.75.csv")
times <- earthquakes_df[,1]

# Number of Observations
N <- nrow(earthquakes_df)

# Convert Dataframe for Visualization
earthquakes <- tibble::as_tibble(earthquakes_df)
earthquakes_sf <- sf::st_as_sf(earthquakes, coords = c("longitude", "latitude"), crs = 4326)

# Convert Time
starting_time <- min(earthquakes_sf$time)
times <- as.numeric(ymd_hms(earthquakes_df[,1])) - as.numeric(ymd_hms(starting_time)) # number of seconds elapsed from starting_time
times <- times / (60*60*24) # number of days elapsed from starting_time
times <- times / 365 # number of years elapsed from starting_time

# Discrete Times
times <- times/max(times)
discrete_times <- times %/% 0.1

# Convert Coordinates from (lat,lon) to (x,y,z)
R <- 1 # radius of the unit sphere
x <- R * cos(earthquakes_df$latitude * pi / 180) * cos(earthquakes_df$longitude * pi / 180)
y <- R * cos(earthquakes_df$latitude * pi / 180) * sin(earthquakes_df$longitude * pi / 180)
z <- R * sin(earthquakes_df$latitude * pi / 180)
locations <- cbind(x,y,z)
locations <- as.data.frame(locations)

## SPATIO-TEMPORAL DISCRETIZATION ----------------------------------------------
# Spatial 2.5D Mesh over the Unit Sphere for Estimation
vertices <- read.table("mesh/sphere.3097.vertices.txt", quote = "\"", comment.char = "")
triangles <- read.table("mesh/sphere.3097.triangles.txt", quote = "\"", comment.char = "")
mesh <- create.mesh.2.5D(nodes = vertices[,2:4], triangles = triangles[,2:4])
FEMbasis <- create.FEM.basis(mesh)
 
# (Simplified) Spatial 2.5D Mesh over the Unit Sphere for Estimation
vertices <- read.table("mesh/sphere.simplified.2000.vertices.txt", quote = "\"", comment.char = "")
triangles <- read.table("mesh/sphere.simplified.2000.triangles.txt", quote = "\"", comment.char = "")
mesh <- create.mesh.2.5D(nodes = vertices[,2:4], triangles = triangles[,2:4])
FEMbasis <- create.FEM.basis(mesh)

# Fine Spatial 2.5D Mesh over the Unit Sphere for Evaluation
vertices.eval <- read.table("mesh/sphere.5016.vertices.txt", row.names = 1)
vertices.eval.proj <- projection.points.2.5D(mesh, vertices.eval)
triangles.eval <- read.table("mesh/sphere.5016.triangles.txt", row.names = 1)
mesh.eval <- create.mesh.2.5D(nodes = vertices.eval, triangles = triangles.eval)
FEMbasis.eval <- create.FEM.basis(mesh.eval)

vertices <- vertices.eval
triangles <- triangles.eval
mesh <- mesh.eval
FEMbasis <- FEMbasis.eval

# Fine Grid Size for Evaluation
n <- 64

# Temporal 1D Mesh Over [0,1]
mesh_time <- c(0, seq(from = 0.15, to = 0.95, by = 0.2), 1)

## ESTIMATION PROCEDURE --------------------------------------------------------

### STDE-PDE: Spatio-Temporal Density Estimation with PDE Regulariz. -----------
# Smoothing Parameters
lambda <- 0.01
# alternative: lambda <- 10^seq(from = -3, to = -1, by = 1)

lambda_time <- 0.001
# alternative: lambda_time <- 10^seq(from = -4, to = -2, by = 1)

# Solution
# [If lambda and/or lambda_time are vectors, to select the best proposals by
# 10-folds CV, please specify: preprocess_method = "RightCV"]
solution_STDEPDE <- DE.FEM.time(data = locations, data_time = times,
                                FEMbasis = FEMbasis, mesh_time = mesh_time,
                                lambda = lambda, lambda_time = lambda_time,
                                fvec = NULL, heatIter = 10, print = TRUE,
                                direction_method = "L-BFGS5", tol1 = 1e-7,
                                preprocess_method = "NoCrossValidation", nfolds = 10,
                                nsimulations = 10000)

FEMfunction_STDEPDE <- FEM.time(solution_STDEPDE$g, mesh_time, FEMbasis, FLAG_PARABOLIC = FALSE)

### STKDE-discrete: Spatio-Temporal Kernel Density Estimation ------------------
# Radius
R <- 1

# Grid Size
n_STKDE_discrete <- n

# Latitudes and Longitudes of the Observed Locations
lat <- asin(locations[,3] / R)*180/pi
long <- atan2(locations[,2], locations[,1])*180/pi

# Grid
x_STKDE_discrete <- seq(from = min(long), to = max(long), length.out = n_STKDE_discrete)
y_STKDE_discrete <- seq(from = min(lat), to = max(lat), length.out = n_STKDE_discrete)
grid_STKDE_discrete <- expand.grid(x_STKDE_discrete, y_STKDE_discrete)
FEMbasis_STKDE_discrete <- create.FEM.basis(create.mesh.2D(grid_STKDE_discrete))

# Latitudes and Longitudes of Mesh Nodes for Evaluation
mesh.eval.nodes.lat <- asin(FEMbasis.eval$mesh$nodes[,3] / R) * 180 / pi
mesh.eval.nodes.long <- atan2(FEMbasis.eval$mesh$nodes[,2], FEMbasis.eval$mesh$nodes[,1]) * 180 / pi

# Solution
x11(width = 25, height = 14)
solution_STKDE_discrete <- stkde(xlong = long, ylat = lat, ztime = discrete_times,
                                 xgrids = n_STKDE_discrete, ygrids = n_STKDE_discrete, breaks = 0.05,
                                 alpha = 0.05, nrowspar = 5, bwmethod = "cv.ml") # alternative: bwmethod = "normal-reference"
dev.off()

### STKDE-separable: STKDE for 1st-Order Separable Spatio-Temporal Point Processes ----
# Spatial Component
locations_spherical <- sphere.cart.to.sphere.spherical(as.matrix(locations))
kappa <- kde.compute.concentration(locations_spherical)
f_space_STKDE_separable <- kde.fhat.cart(vertices.eval, as.matrix(locations), kappa = kappa)

# Temporal Component
f_time_STKDE_separable <- ks::kde(x = as.matrix(times), gridsize = 50)

### KNNSTDE-separable: k-Nearest Neighbors Spatio-Temporal Density Estimation -----------
###  for 1st-Order Separable Spatio-Temporal Point Processes
# Spatial Component
f_space_KNNSTDE_separable <- knnDE(as.matrix(locations), as.matrix(vertices.eval[,1:3]), 25)

# Temporal Component
f_time_KNNSTDE_separable <- knnDE(as.matrix(times), as.matrix(seq(from = 0, to = 1, length.out = 25)), 25)

### ESTIMATES ------------------------------------------------------------------
# Time half width
h <- 1/60

# Time instants at which the solution is evaluated (midpoints of each month)
t <- seq(from = h, to = 1-h, by = 2*h)

# Discrete time instants at which the solution is evaluated (midpoints of each 3 month-period)
t_discrete <- seq(from = 0.05, to = 0.95, by = 0.1)

mean_sol_STDEPDE <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
mean_sol_STKDE_discrete <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t_discrete))
mean_sol_STKDE_separable <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
mean_sol_KNNSTDE_separable <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))

for(time_index in 1:length(t)) {
  
  # STDE-PDE
  evaluation_STDEPDE <- eval.FEM.time(FEM.time = FEMfunction_STDEPDE, locations = vertices.eval.proj, time.instants = t[time_index])
  evaluation_STDEPDE <- exp(evaluation_STDEPDE)
  evaluation_STDEPDE <- evaluation_STDEPDE / sum(evaluation_STDEPDE, na.rm = TRUE)
  
  mean_sol_STDEPDE[,time_index] <- mapply(sum, mean_sol_STDEPDE[,time_index], evaluation_STDEPDE, na.rm = TRUE)
  
  if(sum(abs(t[time_index] - t_discrete) < 1e-06)){
    time_index_discrete <- which(abs(t[time_index] - t_discrete) < 1e-06)
    
    # STKDE-discrete
    FEM_STKDE_discrete <- FEM(c(solution_STKDE_discrete$dens[,,time_index_discrete]), FEMbasis_STKDE_discrete)
    evaluation_STKDE_discrete <- eval.FEM(FEM_STKDE_discrete, locations = cbind(mesh.eval.nodes.long, mesh.eval.nodes.lat))
    evaluation_STKDE_discrete[is.na(evaluation_STKDE_discrete)] <- 0.0
    evaluation_STKDE_discrete <- evaluation_STKDE_discrete / sum(evaluation_STKDE_discrete, na.rm = TRUE)
    
    mean_sol_STKDE_discrete[,time_index_discrete] <- mapply(sum, mean_sol_STKDE_discrete[,time_index_discrete], evaluation_STKDE_discrete, na.rm = TRUE)
  }  
  
  # STKDE-separable
  marginal_time_STKDE_separable <- f_time_STKDE_separable$estimate[findInterval(t[time_index], f_time_STKDE_separable$eval.points)]
  solution_STKDE_separable <- marginal_time_STKDE_separable * f_space_STKDE_separable
  solution_STKDE_separable <- solution_STKDE_separable / sum(solution_STKDE_separable, na.rm = TRUE)
  
  mean_sol_STKDE_separable[,time_index] <- mapply(sum, mean_sol_STKDE_separable[,time_index], solution_STKDE_separable, na.rm = TRUE)
  
  # KNNSTDE-separable
  marginal_time_KNNSTDE_separable <- f_time_KNNSTDE_separable[findInterval(t[time_index], seq(from = 0, to = 1, length.out = 25))]
  solution_KNNSTDE_separable <- marginal_time_KNNSTDE_separable * f_space_KNNSTDE_separable
  solution_KNNSTDE_separable <- solution_KNNSTDE_separable / sum(solution_KNNSTDE_separable, na.rm = TRUE)
  
  mean_sol_KNNSTDE_separable[,time_index] <- mapply(sum, mean_sol_KNNSTDE_separable[,time_index], solution_KNNSTDE_separable, na.rm = TRUE)
  
} 

save.image("aaaa.RData")

#'############################################################################'#
#'############################################################################'#

## VISUALIZATION ---------------------------------------------------------------
{
  ### ESTIMATES [FIGURE 2 and 11] ----------------------------------------------
  # Time half width
  h <- 1/60
  
  # Time instants at which the solution is evaluated (midpoints of each month)
  t <- seq(from = h, to = 1-h, by = 2*h)
  
  # Discrete time instants at which the solution is evaluated (midpoints of each 3 month-period)
  t_discrete <- seq(from = 0.05, to = 0.95, by = 0.1)
  
  # Color Palette
  color_palette <- colorRampPalette(c("#0E1E44", "#3E6DD8", "#68D061", "#ECAF53", "#EB5F5F", "#E11F1C"))(1000)
  # alternative: color_palette = jet.col(n = 1000, alpha = 0.8)
  
  palette(color_palette)
  
  # Plots
  M <- max(max(mean_sol_STDEPDE, na.rm = TRUE),
           max(mean_sol_STKDE_discrete, na.rm = TRUE),
           max(mean_sol_STKDE_separable, na.rm = TRUE),
           max(mean_sol_KNNSTDE_separable, na.rm = TRUE), na.rm = TRUE)
  
  dir.create(file.path(getwd(), "pictures"), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/sample")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/STDEPDE")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/STKDE_discrete")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/STKDE_separable")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/KNNSTDE_separable")), showWarnings = FALSE)
  dir.create(file.path(getwd(), paste0("pictures/time_bar")), showWarnings = FALSE)
  
  view <- c("userMatrix1", "userMatrix2", "userMatrix3", "userMatrix4")
  l <- c("a","b","c","d")

  for(time_index in 1:length(t)){
    
    M <- max(mean_sol_STDEPDE[,time_index], na.rm = TRUE)
    
    for(v in view){
      
      t_img <- proc.time()
      
      idx <- which(v == view)
      ll <- l[idx]
      
      # First Plot: Sample
      plot.sample(FEM(rep(0, dim(vertices.eval.proj)[1]), FEMbasis.eval), locations[abs(times-t[time_index])<1.5*h,], v)
      string <- paste0("pictures/sample/app2_sample_", time_index, "_", ll, ".png")
      snapshot3d(string, fmt = "png", width = 750, height = 750, webshot = rgl.useNULL())
      rgl.close()
      
      # Second Plot: STDE-PDE Estimated Density at t[time_index]
      plot.FEM(FEM(mean_sol_STDEPDE[,time_index], FEMbasis.eval), v, m = 0, M = M)
      string <- paste0("pictures/STDEPDE/app2_STDEPDE_", time_index, "_", ll, ".png")
      snapshot3d(string, fmt = "png", width = 750, height = 750, webshot = rgl.useNULL())
      rgl.close()
      
      if(sum(abs(t[time_index] - t_discrete) < 1e-06)){
        time_index_discrete <- which(abs(t[time_index] - t_discrete) < 1e-06)
        
        # Third Plot: STKDE-discrete Estimated Density at t[time_index]
        plot.FEM(FEM(mean_sol_STKDE_discrete[,time_index_discrete], FEMbasis.eval), v, m = 0, M = M*2)
        string <- paste0("pictures/STKDE_discrete/app2_STKDE_discrete_", time_index, "_", ll, ".png")
        snapshot3d(string, fmt = "png", width = 750, height = 750, webshot = rgl.useNULL())
        rgl.close()
      }
      
      # Fourth Plot: STKDE-separable Estimated (Separable) Density at t[time_index]
      plot.FEM(FEM(mean_sol_STKDE_separable[,time_index], FEMbasis.eval), v, m = 0, M = M)
      string <- paste0("pictures/STKDE_separable/app2_STKDE_separable_", time_index, "_", ll, ".png")
      snapshot3d(string, fmt = "png", width = 750, height = 750, webshot = rgl.useNULL())
      rgl.close()
      
      # Fifth Plot: KNNSTDE-separable Estimated (Separable) Density at t[time_index]
      plot.FEM(FEM(mean_sol_KNNSTDE_separable[,time_index], FEMbasis.eval), v, m = 0, M = M)
      string <- paste0("pictures/KNNSTDE_separable/app2_KNNSTDE_separable_", time_index, "_", ll, ".png")
      snapshot3d(string, fmt = "png", width = 750, height = 750, webshot = rgl.useNULL())
      rgl.close()
      
    }
    
    print(paste("Images for t =", round(t[time_index],2), "done in", round(as.numeric(proc.time() - t_img)[3]), "seconds."))
    
  }
  
  ## TIME BAR ------------------------------------------------------------------
  for(time_index in 1:length(t)) {
    
    string <- paste0("pictures/time_bar/app2_time_bar_",time_index,".pdf")
    pdf(file = string, family = "serif", width = 7, height = 1.5)
    
    #x11(width = 7, height = 1.5)
    par(mar = c(0.01, 2, 0.01, 2))
    bar <- rbind(c(t[time_index],0),c(t[time_index],0.075))
    plot(bar, type = "l", bty = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "",
         pch = 3, lwd = 5, xlim = c(0,1), asp = 1)
    mtext("Time", line = -1.5, font = 2, cex = 1.5)
    axis(side = 1, c(0, 0.5, 1), pos = 0, labels = c("0", "0.5", "1"))
    
    dev.off()
    
  }
  
  ## LEGEND [FIGURE 11] --------------------------------------------------------
  pdf(paste0("pictures/app2_legend.pdf"), family = "serif", width = 11, height = 3)
  #x11(width = 11,height = 3)
  par(mai = c(1,0.75,0,0))
  plot(c(0, 112.5), c(0, 15), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame.plot = F)
  gradient.rect(0, 0, 100, 5, col = color_palette, border = "black")
  axis(1,at = c(0,25,50,75,100), labels = c("0", "0.5", "1", "1.5", "2"),
       lwd.ticks = 2, cex.axis = 2, lwd = 2)
  text(107,2, TeX("$\\times 10^{-2}$"), cex = 2)
  dev.off()
  
  print(paste0("M = ", round(M, 4)))
  
  ## MESH [FIGURE 3] -----------------------------------------------------------
  string <- paste0("pictures/app2_mesh.png")
  plot.mesh(mesh, "userMatrix2")
  s <- sample(1:N, 1000)
  rgl.spheres(x = locations[s,1], y = locations[s,2], z = locations[s,3], radius = 0.005, color = rgb(190, 0, 0, max = 255), aspect = T, add = T)
  snapshot3d(string, fmt = "png", width = 750, height = 750, webshot = rgl.useNULL())
  rgl.close()
}

#'############################################################################'#
#'############################################################################'#

# CROSS VALIDATION ERROR -------------------------------------------------------
# Number of Folds
K <- 10

# Subsample Size in Each Fold
n_k <- floor(N/K)

# Indices
folds <- sample(1:N, N, replace = FALSE)

# Cross-Validation Error Containers
CV_error_STDEPDE <- rep(0,K)
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
  
  discrete_times_train_k <- times_train_k %/% 0.1
  discrete_times_test_k <- times_test_k %/% 0.1
  
  # Plot
  # x11()
  # plot(boundary_df$Northing, boundary_df$Easting, type = "l", main = paste0("Training and test sets, iter = ",k))
  # legend("topleft", legend = c(paste0("Train: ", nrow(locations_train_k)), paste0("Test: ", nrow(locations_test_k))), fill = c("red", "green"))
  # points(locations_train_k, col = "red", pch = 19)
  # points(locations_test_k, col = "green", pch = 19)
  
  ## STDE-PDE: Spatio-Temporal Density Estimation with PDE Regulariz. ----------
  # Smoothing Parameters
  lambda <- 0.01
  # alternative: lambda <- 10^seq(from = -3, to = -1, by = 1)
  
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
  
  ## STKDE-discrete: Spatio-Temporal Kernel Density Estimation -----------------
  # Grid Size
  n_STKDE_discrete <- n
  
  # Latitudes and Longitudes of the Observed Locations
  lat <- asin(locations_train_k[,3] / R)*180/pi
  long <- atan2(locations_train_k[,2], locations_train_k[,1])*180/pi
  
  # Grid
  x_STKDE_discrete <- seq(min(long), max(long), length.out = n_STKDE_discrete)
  y_STKDE_discrete <- seq(min(lat), max(lat), length.out = n_STKDE_discrete)
  grid_STKDE_discrete <- expand.grid(x_STKDE_discrete, y_STKDE_discrete)
  FEMbasis_STKDE_discrete_train <- create.FEM.basis(create.mesh.2D(grid_STKDE_discrete))
  
  # Latitudes and Longitudes of Mesh Nodes for Evaluation
  mesh.eval.nodes.lat <- asin(FEMbasis.eval$mesh$nodes[,3] / R) * 180 / pi
  mesh.eval.nodes.long <- atan2(FEMbasis.eval$mesh$nodes[,2], FEMbasis.eval$mesh$nodes[,1]) * 180 / pi
  
  # Solution
  x11(width = 25, height = 14)
  solution_STKDE_discrete_train <- stkde(xlong = long, ylat = lat, ztime = discrete_times_train_k,
                                         xgrids = n_STKDE_discrete, ygrids = n_STKDE_discrete, breaks = 0.05,
                                         alpha = 0.05, nrowspar = 5, bwmethod = "normal-reference") # alternative: bwmethod = "cv.ml"
  dev.off()
  
  ## STKDE-separable: STKDE for 1st-Order Separable Spatio-Temporal Point Processes ----
  # Spatial Component
  locations_spherical <- sphere.cart.to.sphere.spherical(as.matrix(locations_train_k))
  kappa <- kde.compute.concentration(locations_spherical)
  f_space_STKDE_separable_train <- kde.fhat.cart(vertices.eval, as.matrix(locations_train_k), kappa = kappa)
  
  # Temporal Component
  f_time_STKDE_separable_train <- ks::kde(x = as.matrix(times_train_k), gridsize = 50)
  
  ## KNNSTDE-separable: k-Nearest Neighbors Spatio-Temporal Density Estimation -----------
  ###  for 1st-Order Separable Spatio-Temporal Point Processes
  # Spatial Component
  f_space_KNNSTDE_separable_train <- knnDE(as.matrix(locations_train_k), as.matrix(vertices.eval[,1:3]), 1000)
  
  # Temporal Component
  f_time_KNNSTDE_separable_train <- knnDE(as.matrix(times_train_k), as.matrix(seq(0, 1, length.out = 25)), 25)
  
  ## ERRORS --------------------------------------------------------------------
  # Time half width
  h <- 1/60
  
  # Time instants at which the solution is evaluated (midpoints of each month)
  t <- seq(from = h, to = 1-h, by = 2*h)
  
  # Discrete time instants at which the solution is evaluated (midpoints of each 3 month-period)
  t_discrete <- seq(0.05, 0.95, by = 0.1)
  
  for(time_index in 1:length(t)){
    
    # STDEPDE
    evaluation_STDEPDE_train <- eval.FEM.time(FEM.time = FEMfunction_STDEPDE_train, locations = vertices.eval.proj, time.instants = t[time_index])
    evaluation_STDEPDE_train <- exp(evaluation_STDEPDE_train)
    evaluation_STDEPDE_train <- evaluation_STDEPDE_train / sum(evaluation_STDEPDE_train, na.rm = TRUE)
    
    evaluation_STDEPDE_test <- eval.FEM(FEM = FEM(evaluation_STDEPDE_train, FEMbasis.eval), locations = projection.points.2.5D(mesh.eval, locations_test_k[which(abs(t[time_index]-times_test_k)<h),]))
    
    cardinality_k <- 1 / (sum(!is.na(evaluation_STDEPDE_test)) / nrow(locations_test_k[which(abs(t[time_index]-times_test_k)<h),]))
    CV_error_STDEPDE[k] <- CV_error_STDEPDE[k] + sum(evaluation_STDEPDE_train^2, na.rm = TRUE) / sum(!is.na(evaluation_STDEPDE_train)) -2 * sum(evaluation_STDEPDE_test, na.rm = TRUE) * cardinality_k
    
    if(sum(abs(t[time_index] - t_discrete) < 1e-06)){
      
      time_index_discrete <- which(abs(t[time_index] - t_discrete) < 1e-06)[1]
      
      # STKDE-discrete
      coeff <- c(solution_STKDE_discrete_train$dens[,,time_index_discrete])
      FEM_STKDE_discrete_train <- FEM(coeff, FEMbasis_STKDE_discrete_train)
      evaluation_STKDE_discrete_train <- eval.FEM(FEM_STKDE_discrete_train, locations = cbind(mesh.eval.nodes.long, mesh.eval.nodes.lat))
      evaluation_STKDE_discrete_train <- evaluation_STKDE_discrete_train / sum(evaluation_STKDE_discrete_train, na.rm = TRUE)
      
      FEM_STKDE_discrete_test <- FEM(evaluation_STKDE_discrete_train, FEMbasis.eval)
      evaluation_STKDE_discrete_test <- eval.FEM(FEM_STKDE_discrete_test, locations = projection.points.2.5D(mesh.eval, locations_test_k[which(abs(t[time_index]-times_test_k)<h),]))
      
      cardinality_k <- 1 / (sum(!is.na(evaluation_STKDE_discrete_test)) / nrow(locations_test_k[which(abs(t[time_index]-times_test_k)<h),]))
      CV_error_STKDE_discrete[k] <- CV_error_STKDE_discrete[k] + sum(evaluation_STKDE_discrete_train^2, na.rm = TRUE) / sum(!is.na(evaluation_STKDE_discrete_train)) -2 * sum(evaluation_STKDE_discrete_test, na.rm = TRUE) * cardinality_k
      
    }  
    
    # STKDE-separable
    marginal_time_STKDE_train <- f_time_STKDE_separable_train$estimate[findInterval(t[time_index], f_time_STKDE_separable_train$eval.points)]
    solution_STKDE_separable_train <- marginal_time_STKDE_train * f_space_STKDE_separable_train
    solution_STKDE_separable_train <- solution_STKDE_separable_train / sum(solution_STKDE_separable_train, na.rm = TRUE)
    
    FEM_STKDE_separable_test <- FEM(solution_STKDE_separable_train, FEMbasis.eval)
    solution_STKDE_separable_test <- eval.FEM(FEM_STKDE_separable_test, locations = projection.points.2.5D(mesh.eval, locations_test_k[which(abs(t[time_index]-times_test_k)<h),]))
    
    cardinality_k <- 1 / (sum(!is.na(solution_STKDE_separable_test)) / nrow(locations_test_k[which(abs(t[time_index]-times_test_k)<h),]))
    CV_error_STKDE_separable[k] <- CV_error_STKDE_separable[k] + sum(solution_STKDE_separable_train^2, na.rm = TRUE) / sum(!is.na(solution_STKDE_separable_train)) -2 * sum(solution_STKDE_separable_test, na.rm = TRUE) * cardinality_k
    
    # KNNSTDE-separable
    marginal_time_KNNSTDE_separable_train <- f_time_KNNSTDE_separable_train[findInterval(t[time_index], seq(-0.25, 10.25, length.out = 25))]
    solution_KNNSTDE_separable_train <- marginal_time_KNNSTDE_separable_train * f_space_KNNSTDE_separable_train
    solution_KNNSTDE_separable_train <- solution_KNNSTDE_separable_train / sum(solution_KNNSTDE_separable_train, na.rm = TRUE)
    
    FEM_KNNSTDE_separable_test <- FEM(solution_KNNSTDE_separable_train, FEMbasis.eval)
    solution_KNNSTDE_separable_test <- eval.FEM(FEM_KNNSTDE_separable_test, locations = projection.points.2.5D(mesh.eval, locations_test_k[which(abs(t[time_index]-times_test_k)<h),]))
    
    cardinality_k <- 1 / (sum(!is.na(solution_KNNSTDE_separable_test)) / nrow(locations_test_k[which(abs(t[time_index]-times_test_k)<h),]))
    CV_error_KNNSTDE_separable[k] <- CV_error_KNNSTDE_separable[k] + sum(solution_KNNSTDE_separable_train^2, na.rm = TRUE) / sum(!is.na(solution_KNNSTDE_separable_train)) -2 * sum(solution_KNNSTDE_separable_test, na.rm = TRUE) * cardinality_k
    
  }
  
  CV_error_STDEPDE[k] <- CV_error_STDEPDE[k] / length(t)
  CV_error_STKDE_discrete[k] <- CV_error_STKDE_discrete[k] / length(t_discrete)
  CV_error_STKDE_separable[k] <- CV_error_STKDE_separable[k] / length(t)
  CV_error_KNNSTDE_separable[k] <- CV_error_KNNSTDE_separable[k] / length(t)
  
  print(paste("CV iteration", k, "done."))
  
}

## BOXPLOT [FIGURE 10] ---------------------------------------------------------
{
  blue_polimi <- rgb(114, 143, 165, max = 255, alpha = 65, names = "blue_polimi25") # #728FA5
  
  pdf(paste0("pictures/app2_cv_errors.pdf"), family = "serif", width = 12, height = 5.3)
  #x11(width = 12,height = 5.3)
  par(mfrow = c(1,3), mai = c(1.75,0.35,0.5,0.15))
  plot.new()
  boxplot(CV_error_STDEPDE, CV_error_STKDE_discrete, CV_error_STKDE_separable, CV_error_KNNSTDE_separable, asp = 1, xaxt = "n", yaxt = "n")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = blue_polimi)
  grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 1)
  par(new = TRUE)
  boxplot(CV_error_STDEPDE, CV_error_STKDE_discrete, CV_error_STKDE_separable, CV_error_KNNSTDE_separable,
          names = c("STDE-PDE", "STKDE-dis", "STKDE-sep", "KNNSTDE-sep"), main = "CV error",
          col = brewer.pal(8, "YlGnBu")[c(1,6,7,8)], cex.lab = 2, cex.axis = 2, cex.main = 2, las = 3)
  dev.off()
}

save.image("output/cv_errors.RData")



