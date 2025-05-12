#'############################################################################'#
#'############ APPLICATION - ANALYSIS OF CAR ACCIDENT ON BERGAMO #############'#
#'############################################################################'#

graphics.off()
rm(list=ls())
options(warn = -1, timeout = 800)

## LIBRARIES AND UTILITIES -----------------------------------------------------
if(!require(pacman)) install.packages("pacman")
pacman::p_load("rstudioapi")

# Setting the Working Directory 
setwd(dirname(getActiveDocumentContext()$path))

# Loading Packages and Auxiliary Functions
file.sources <- list.files("utils/", pattern = "*.R", full.names = TRUE)
invisible(sapply(file.sources, source, .GlobalEnv))

## SPATIO-TEMPORAL DISCRETIZATION ----------------------------------------------
# District of Interest
district <- "Bergamo"

# Administrative Boundaries
adm_bdy_IT <- st_read("utils/map/adm_bdy_IT/Com01012023/Com01012023_WGS84.shx")
adm_bdy_Bergamo <- adm_bdy_IT$geometry[which(adm_bdy_IT$COMUNE == district)]
adm_bdy_Bergamo <- st_transform(adm_bdy_Bergamo, crs = 4326)
adm_bdy_Bergamo_bbox <- st_bbox(adm_bdy_Bergamo)

# north_west_italy <- st_read("utils/map/north-west-italy/gis_osm_roads_free_1.shp")
# north_west_italy <- st_as_sf(north_west_italy)
# bg <- st_within(north_west_italy, adm_bdy_Bergamo)
# bg <- north_west_italy[which(lengths(bg) != 0), ]
# plot(adm_bdy_Bergamo)
# plot(st_geometry(filt_data), add = TRUE, reset = FALSE)

# Mesh for Estimation
mesh <- generateMesh(district = district)
mesh_df <- data.frame(lon = mesh$nodes[,1], lat = mesh$nodes[,2])
mesh_df <- st_as_sf(mesh_df, coords = c("lon", "lat"), crs = 4326)
mesh_linnet <- as.linnet(mesh)
mesh_sfnetwork <- as_sfnetwork(mesh_linnet, directed = FALSE, edges_as_lines = TRUE)
st_crs(mesh_sfnetwork) <- 4326
mesh_sf <- st_as_sf(mesh_sfnetwork, "edges")
mapview(mesh_sf, legend = FALSE, layer.name = "road-netowrk", lwd = 1.5,
               alpha.regions = 1) +
  mapview(st_sf(data.frame(x = 1), geoemtry = adm_bdy_Bergamo), legend = FALSE,
          layer.name = "boundary", lwd = 0.1) +
  mapview(mesh_df, legend = FALSE, layer.name = "mesh-nodes",
          alpha.regions = 1, cex = 0.75)

# Mesh for Evaluation
mesh.eval <- refine.mesh.1.5D(mesh, 0.003)
mesh.eval_df <- data.frame(lon = mesh.eval$nodes[,1], lat = mesh.eval$nodes[,2])
mesh.eval_df <- st_as_sf(mesh.eval_df, coords = c("lon", "lat"), crs = 4326)
mesh.eval_linnet <- as.linnet.mesh.1.5D(mesh.eval)
mesh.eval_sfnetwork <- as_sfnetwork(mesh.eval_linnet, directed = FALSE, edges_as_lines = TRUE)
st_crs(mesh.eval_sfnetwork) <- 4326
mesh.eval_sf <- st_as_sf(mesh.eval_sfnetwork, "edges")

# SubMesh for Estimation
submesh <- generateSubMesh(mesh = mesh_sfnetwork, xmin = 9.65257, ymin = 45.67985, xmax = 9.67124, ymax = 45.69773)
submesh_df <- data.frame(lon = submesh$nodes[,1], lat = submesh$nodes[,2])
submesh_df <- st_as_sf(submesh_df, coords = c("lon", "lat"), crs = 4326)
submesh_linnet <- as.linnet(submesh)
submesh_sfnetwork <- as_sfnetwork(submesh_linnet, directed = FALSE, edges_as_lines = TRUE)
st_crs(submesh_sfnetwork) <- 4326
submesh_sf <- st_as_sf(submesh_sfnetwork, "edges")
mapview(submesh_sf, legend = FALSE, layer.name = "road-netowrk", lwd = 1.5,
        alpha.regions = 1) +
  mapview(st_sf(data.frame(x = 1), geoemtry = adm_bdy_Bergamo), legend = FALSE,
          layer.name = "boundary", lwd = 0.1) +
  mapview(submesh_df, legend = FALSE, layer.name = "submesh-nodes",
          alpha.regions = 1, cex = 0.75)

# SubMesh for Evaluation
submesh.eval <- refine.mesh.1.5D(submesh, 0.0019)
submesh.eval_df <- data.frame(lon = submesh.eval$nodes[,1], lat = submesh.eval$nodes[,2])
submesh.eval_df <- st_as_sf(submesh.eval_df, coords = c("lon", "lat"), crs = 4326)
submesh.eval_linnet <- as.linnet.mesh.1.5D(submesh.eval)
submesh.eval_sfnetwork <- as_sfnetwork(submesh.eval_linnet, directed = FALSE, edges_as_lines = TRUE)
st_crs(submesh.eval_sfnetwork) <- 4326
submesh.eval_sf <- st_as_sf(submesh.eval_sfnetwork, "edges")

# Grid Size for Estimation
n <- 512

## DATA ------------------------------------------------------------------------
# Import the Dataset
data <- loadDataset(district = district, years = c(2017, 2019), mesh)
data_workdays <- data[[1]]
data_holidays <- data[[2]]
     
# Create the Current Date / Time
current_datetime <- Sys.time()
formatted_datetime <- format(current_datetime, "[%Y-%m-%d]-[%H-%M-%S]")

# Create the Current Directory
dir.create(paste0(getwd(),"/application-", formatted_datetime))
setwd(paste0(getwd(),"/application-", formatted_datetime))
foldername <- paste0(getwd(),"/")

## ESTIMATION PROCEDURE --------------------------------------------------------
for(daytype in c("workdays", "holidays")){
  
  if(daytype == "workdays"){
    df <- data_workdays
  } else {
    df <- data_holidays
  }
  
  xmin = min(submesh$nodes[,1])
  ymin = min(submesh$nodes[,2])
  xmax = max(submesh$nodes[,1])
  ymax = max(submesh$nodes[,2])

  df = df[df$lon >= xmin & df$lat >= ymin & df$lon <= xmax & df$lat <= ymax,]
  
  locations <- cbind(df$lon, df$lat)
  times <- df$times
  
  ### STDE-PDE: Spatio-Temporal Density Estimation with PDE Regulariz. ---------
  t0 <- proc.time()
  
  FEMbasis <- create.FEM.basis(submesh)
  FEMbasis.eval <- create.FEM.basis(submesh.eval)
  
  h <- 1/8
  t <- seq(from = 0+h/2, to = 1-h/2, by = h)
  mesh_time <- c(0, t[seq(from = 2, to = length(t), by = 2)], 1)
  
  t <- seq(from = 0, to = 24, by = 1)/24
  
  lambda <- 0.01
  lambda_time <- 0.01
  
  #invisible(capture.output(
  solution_STDEPDE <- DE.FEM.time(data = locations, data_time = times, FEMbasis = FEMbasis,
                                                           mesh_time = mesh_time, lambda = lambda, tol1 = 1e-7,
                                                           lambda_time = lambda_time, fvec = NULL, heatStep = 0.1,
                                                           heatIter = 10, print = TRUE, nfolds = 10, nsimulations = 10000,
                                                           step_method = "Wolfe_Method", direction_method = "L-BFGS5",
                                                           preprocess_method = "NoCrossValidation")#))
  
  FEMfunction_STDEPDE <- FEM.time(solution_STDEPDE$g, mesh_time, FEMbasis, FLAG_PARABOLIC = FALSE)
  
  # CPU Time
  CPUtime <- as.numeric(proc.time() - t0)[3]
  cat(paste0("STDE-PDE done in ", round(CPUtime, 2), " seconds.\n"))
  
  ## STLNPP: Spatio-Temporal Kernel Intensity Estimation ---------------------
  t0 <- proc.time()

  #observations_stlpp <- as.stlpp(locations[,1], locations[,2], times, submesh_linnet)
  observations_stlpp <- as.stlpp(solution_STDEPDE$data[,1], solution_STDEPDE$data[,2], solution_STDEPDE$data_time, submesh_linnet)
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
  cat(paste0("STLNPP done in ", round(CPUtime, 2), " seconds.\n"))
  
  # ## HEAT: Spatio-Temporal Separable Kernel Intensity Estimation based on the Heat Equation -----
  # t0 <- proc.time()
  # 
  # #observations_stlpp <- as.stlpp(locations[,1], locations[,2], times, submesh_linnet)
  # observations_stlpp <- as.stlpp(solution_STDEPDE$data[,1], solution_STDEPDE$data[,2], solution_STDEPDE$data_time, submesh_linnet)
  # solution_HEAT <- density1Dkernel(observations_stlpp, at = "pixels",
  #                                  leaveoneout = FALSE, iterMax = 1e+09,
  #                                  verbose = FALSE)
  # 
  # tgrid_HEAT <- attr(solution_HEAT, "tgrid")
  # 
  # # CPU Time
  # CPUtime <- as.numeric(proc.time() - t0)[3]
  # cat(paste0("HEAT done in ", round(CPUtime, 2), " seconds.\n"))
  #
  # ## EQUAL-SPLIT: Spatio-Temporal Separable Kernel Intensity Estimation based on the the Okabe-Sugihara Equal-Split Algorithm -----
  # t0 <- proc.time()
  # 
  # #observations_stlpp <- as.stlpp(locations[,1], locations[,2], times, submesh_linnet)
  # observations_stlpp <- as.stlpp(solution_STDEPDE$data[,1], solution_STDEPDE$data[,2], solution_STDEPDE$data_time, submesh_linnet)
  # solution_EQUAL_SPLIT <- density1Dkernel(observations_stlpp, at = "pixels",
  #                                         leaveoneout = FALSE, iterMax = 1e+09,
  #                                         kernel = "epanechnikov", verbose = FALSE)
  # 
  # tgrid_EQUAL_SPLIT <- attr(solution_EQUAL_SPLIT, "tgrid")
  # 
  # # CPU Time
  # CPUtime <- as.numeric(proc.time() - t0)[3]
  # cat(paste0("EQUAL-SPLIT done in ", round(CPUtime, 2), " seconds.\n"))
  # 
  # ## VORONOI: Spatio-Temporal Pseudo-Separable Intensity Estimation based on the Voronoi-Dirichlet Tessellation -----
  # t0 <- proc.time()
  # 
  # #observations_stlpp <- as.stlpp(locations[,1], locations[,2], times, submesh_linnet)
  # observations_stlpp <- as.stlpp(solution_STDEPDE$data[,1], solution_STDEPDE$data[,2], solution_STDEPDE$data_time, submesh_linnet)
  # invisible(capture.output(solution_VORONOI <- densityVoronoi(observations_stlpp, f = 0.9, nrep = 10, separable = FALSE,
  #                                                             dimt = n, at = "pixels")))
  # 
  # tgrid_VORONOI <- attr(solution_VORONOI, "tgrid")
  # 
  # # CPU Time
  # CPUtime <- as.numeric(proc.time() - t0)[3]
  # cat(paste0("VORONOI done in ", round(CPUtime, 2), " seconds.\n"))
  # 
  # ## VORONOI-SEP: Spatio-Temporal Separable Intensity Estimation based on the Voronoi-Dirichlet Tessellation -----
  # t0 <- proc.time()
  # 
  # #observations_stlpp <- as.stlpp(locations[,1], locations[,2], times, submesh_linnet)
  # observations_stlpp <- as.stlpp(solution_STDEPDE$data[,1], solution_STDEPDE$data[,2], solution_STDEPDE$data_time, submesh_linnet)
  # invisible(capture.output(solution_VORONOI_SEP <- densityVoronoi(observations_stlpp, f = 1, nrep = 10, separable = TRUE,
  #                                                                 dimt = n, at = "pixels")))
  # 
  # tgrid_VORONOI_SEP <- attr(solution_VORONOI_SEP, "tgrid")
  # 
  # # CPU Time
  # CPUtime <- as.numeric(proc.time() - t0)[3]
  # cat(paste0("VORONOI-SEP done in ", round(CPUtime, 2), " seconds.\n"))
  
  ## ESTIMATES -----------------------------------------------------------------
  mean_sol_STDEPDE <- matrix(nrow = nrow(submesh.eval$nodes), ncol = length(t))
  mean_sol_STLNPP <- matrix(nrow = nrow(submesh.eval$nodes), ncol = length(t))
  # mean_sol_HEAT <- matrix(nrow = nrow(submesh.eval$nodes), ncol = length(t))
  # mean_sol_EQUAL_SPLIT <- matrix(nrow = nrow(submesh.eval$nodes), ncol = length(t))
  # mean_sol_VORONOI <- matrix(nrow = nrow(submesh.eval$nodes), ncol = length(t))
  # mean_sol_VORONOI_SEP <- matrix(nrow = nrow(submesh.eval$nodes), ncol = length(t))
  
  if(daytype == "workdays"){
    mean_sol_STDEPDE_workdays <- matrix(nrow = nrow(submesh.eval$nodes), ncol = length(t))
    mean_sol_STLNPP_workdays <- matrix(nrow = nrow(submesh.eval$nodes), ncol = length(t))
    # mean_sol_HEAT_workdays <- matrix(nrow = nrow(submesh.eval$nodes), ncol = length(t))
    # mean_sol_EQUAL_SPLIT_workdays <- matrix(nrow = nrow(submesh.eval$nodes), ncol = length(t))
    # mean_sol_VORONOI_workdays <- matrix(nrow = nrow(submesh.eval$nodes), ncol = length(t))
    # mean_sol_VORONOI_SEP_workdays <- matrix(nrow = nrow(submesh.eval$nodes), ncol = length(t))
  } else {
    mean_sol_STDEPDE_holidays <- matrix(nrow = nrow(submesh.eval$nodes), ncol = length(t))
    mean_sol_STLNPP_holidays <- matrix(nrow = nrow(submesh.eval$nodes), ncol = length(t))
    # mean_sol_HEAT_holidays <- matrix(nrow = nrow(submesh.eval$nodes), ncol = length(t))
    # mean_sol_EQUAL_SPLIT_holidays <- matrix(nrow = nrow(submesh.eval$nodes), ncol = length(t))
    # mean_sol_VORONOI_holidays <- matrix(nrow = nrow(submesh.eval$nodes), ncol = length(t))
    # mean_sol_VORONOI_SEP_holidays <- matrix(nrow = nrow(submesh.eval$nodes), ncol = length(t))
  }
  
  for(time_index in 1:length(t)) {
    
    # STDE-PDE
    evaluation_STDEPDE <- eval.FEM.time(FEM.time = FEMfunction_STDEPDE, locations = submesh.eval$nodes, time.instants = t[time_index])
    evaluation_STDEPDE <- exp(evaluation_STDEPDE)*nrow(locations)
    
    # STLNPP
    idx <- which.min(abs(tgrid_STLNPP - t[time_index]))
    evaluation_STLNPP <- as.linfun(solution_STLNPP[[idx]])(submesh.eval$nodes[,1], submesh.eval$nodes[,2])

    # # HEAT
    # idx <- which.min(abs(tgrid_HEAT - t[time_index]))
    # evaluation_HEAT <- as.linfun(solution_HEAT[[idx]])(submesh.eval$nodes[,1], submesh.eval$nodes[,2])
    #
    # # EQUAL-SPLIT
    # idx <- which.min(abs(tgrid_EQUAL_SPLIT - t[time_index]))
    # evaluation_EQUAL_SPLIT <- as.linfun(solution_EQUAL_SPLIT[[idx]])(submesh.eval$nodes[,1], submesh.eval$nodes[,2])
    # 
    # # VORONOI
    # idx <- which.min(abs(tgrid_VORONOI - t[time_index]))
    # evaluation_VORONOI <- as.linfun(solution_VORONOI[[idx]])(submesh.eval$nodes[,1], submesh.eval$nodes[,2])
    # 
    # # VORONOI-SEP
    # idx <- which.min(abs(tgrid_VORONOI_SEP - t[time_index]))
    # evaluation_VORONOI_SEP <- as.linfun(solution_VORONOI_SEP[[idx]])(submesh.eval$nodes[,1], submesh.eval$nodes[,2])
    
    # Estimates Updates
    mean_sol_STDEPDE[,time_index] <- mapply(sum, mean_sol_STDEPDE[,time_index], evaluation_STDEPDE, na.rm = TRUE)
    mean_sol_STLNPP[,time_index] <- mapply(sum, mean_sol_STLNPP[,time_index], evaluation_STLNPP, na.rm = TRUE)
    # mean_sol_HEAT[,time_index] <- mapply(sum, mean_sol_HEAT[,time_index], evaluation_HEAT, na.rm = TRUE)
    # mean_sol_EQUAL_SPLIT[,time_index] <- mapply(sum, mean_sol_EQUAL_SPLIT[,time_index], evaluation_EQUAL_SPLIT, na.rm = TRUE)
    # mean_sol_VORONOI[,time_index] <- mapply(sum, mean_sol_VORONOI[,time_index], evaluation_VORONOI, na.rm = TRUE)
    # mean_sol_VORONOI_SEP[,time_index] <- mapply(sum, mean_sol_VORONOI_SEP[,time_index], evaluation_VORONOI_SEP, na.rm = TRUE)
    
  }
  
  if(daytype == "workdays"){
    mean_sol_STDEPDE_workdays <- mean_sol_STDEPDE
    mean_sol_STLNPP_workdays <- mean_sol_STLNPP
    # mean_sol_HEAT_workdays <- mean_sol_HEAT
    # mean_sol_EQUAL_SPLIT_workdays <- mean_sol_EQUAL_SPLIT
    # mean_sol_VORONOI_workdays <- mean_sol_VORONOI
    # mean_sol_VORONOI_SEP_workdays <- mean_sol_VORONOI_SEP
  } else {
    mean_sol_STDEPDE_holidays <- mean_sol_STDEPDE
    mean_sol_STLNPP_holidays <- mean_sol_STLNPP
    # mean_sol_HEAT_holidays <- mean_sol_HEAT
    # mean_sol_EQUAL_SPLIT_holidays <- mean_sol_EQUAL_SPLIT
    # mean_sol_VORONOI_holidays <- mean_sol_VORONOI
    # mean_sol_VORONOI_SEP_holidays <- mean_sol_VORONOI_SEP
  }
  
}

# dir.create(paste0(foldername, "output"))
# save.image("output/workspace_estimates.RData")

## VISUALIZATION ---------------------------------------------------------------
{
  
  dx <- diff(range(mesh.eval$nodes[,1])) / 5
  dy <- diff(range(mesh.eval$nodes[,1])) / 5 
  center <- c(mean(range(mesh.eval$nodes[,1]))+dx, mean(range(mesh.eval$nodes[,2])))
  
  dir.create(paste0(foldername,"pictures"))
  
  ### MESH ---------------------------------------------------------------------
  dir.create(paste0(foldername,"pictures/mesh"))
  dir.create(paste0(foldername,"pictures/mesh/html"))
  dir.create(paste0(foldername,"pictures/mesh/pdf"))
  
  filename <- "app_mesh"
  plot(obj = mesh, map = FALSE, bdy = adm_bdy_Bergamo, center = center,
       foldername = "pictures/mesh", filename = filename)
  filename <- "app_mesh_map"
  plot(obj = mesh, map = TRUE, bdy = adm_bdy_Bergamo, center = center,
       foldername = "pictures/mesh", filename = filename)
  
  ## SAMPLE --------------------------------------------------------------------
  dir.create(paste0(foldername,"pictures/sample"))
  dir.create(paste0(foldername,"pictures/sample/html"))
  dir.create(paste0(foldername,"pictures/sample/pdf"))
  
  h <- 1/8
  t <- seq(from = 0+h/2, to = 1-h/2, by = h)
  
  df <- rbind(data_workdays, data_holidays)
  
  locations_df <- data.frame(lon = df$lon, lat = df$lat)
  
  filename <- "app_sample"
  plot_sample(obj = locations_df, mesh = mesh, map = FALSE, bdy = adm_bdy_Bergamo,
              center = center, foldername = "pictures/sample",
              filename = filename)
  filename <- "app_sample_map"
  plot_sample(obj = locations_df, mesh = mesh, map = TRUE, bdy = adm_bdy_Bergamo,
              center = center, foldername = "pictures/sample",
              filename = filename)
  filename <- "app_sample_map_box"
  plot_sample(obj = locations_df, mesh = mesh, map = TRUE, bdy = adm_bdy_Bergamo,
              center = center, foldername = "pictures/sample",
              filename = filename, box = TRUE)
  filename <- "app_sample_map_box_complete"
  plot_sample(obj = locations_df, mesh = mesh, map = TRUE, bdy = adm_bdy_Bergamo,
              center = center, foldername = "pictures/sample",
              filename = filename, box = TRUE,
              box1 = TRUE, box2 = TRUE, box3 = TRUE, box4 = TRUE)
  filename <- "app_sample_map_box_zoom"
  plot_sample(obj = locations_df, mesh = mesh, map = TRUE, bdy = adm_bdy_Bergamo,
              center = center, foldername = "pictures/sample", zoom = 13,
              filename = filename, box = FALSE,
              box1 = TRUE, box2 = TRUE, box3 = TRUE, box4 = TRUE)
  
  filename <- "map_right"
  plot_sample(obj = NULL, mesh, map = TRUE, bdy = NULL, center = center,
              foldername = "pictures/sample", filename = filename, zoom = 13)
  
  center_1 <- neighborhood_1
  center_1[1] <- center_1[1] + 0.0025
  center_2 <- neighborhood_2
  center_3 <- neighborhood_3
  center_3[1] <- center_3[1] + 0.0025
  center_4 <- neighborhood_4
  center_4[1] <- center_4[1] + 0.005

  filename <- paste0("app_sample_map_box_1")
  plot_sample(obj = locations_df, mesh = mesh, map = TRUE, bdy = adm_bdy_Bergamo,
              center = center_1, foldername = "pictures/sample",
              filename = filename, zoom = 15, cex = 7.5, stroke = 1.25)
  filename <- paste0("app_sample_map_box_2")
  plot_sample(obj = locations_df, mesh = mesh, map = TRUE, bdy = adm_bdy_Bergamo,
              center = center_2, foldername = "pictures/sample",
              filename = filename, zoom = 15, cex = 7.5, stroke = 1.25)
  filename <- paste0("app_sample_map_box_3")
  plot_sample(obj = locations_df, mesh = mesh, map = TRUE, bdy = adm_bdy_Bergamo,
              center = center_3, foldername = "pictures/sample",
              filename = filename, zoom = 15, cex = 7.5, stroke = 1.25)
  filename <- paste0("app_sample_map_box_4")
  plot_sample(obj = locations_df, mesh = mesh, map = TRUE, bdy = adm_bdy_Bergamo,
              center = center_4, foldername = "pictures/sample",
              filename = filename, zoom = 15, cex = 7.5, stroke = 1.25)

  for(daytype in c("workdays", "holidays")){
    
    if(daytype == "workdays"){
      df <- data_workdays
    } else {
      df <- data_holidays
    }
    
    locations_df <- data.frame(lon = df$lon, lat = df$lat)
    
    filename <- paste0("app_sample_", daytype)
    plot_sample(obj = locations_df, mesh = mesh, map = FALSE, bdy = adm_bdy_Bergamo,
                center = center, foldername = "pictures/sample",
                filename = filename)
    filename <- paste0("app_sample_map_", daytype)
    plot_sample(obj = locations_df, mesh = mesh, map = TRUE, bdy = adm_bdy_Bergamo,
                center = center, foldername = "pictures/sample",
                filename = filename)

    for(time_index in 1:length(t)){
      
      idxs <- which(abs(df$times - t[time_index]) < h/2)
      locations_df <- data.frame(lon = df$lon[idxs], lat = df$lat[idxs])
      
      filename <- paste0("app_sample_", daytype, "_", time_index)
      plot_sample(obj = locations_df, mesh = mesh, map = FALSE, bdy = adm_bdy_Bergamo,
                  center = center, foldername = "pictures/sample",
                  filename = filename)
      filename <- paste0("app_sample_map_", daytype, "_", time_index)
      plot_sample(obj = locations_df, mesh = mesh, map = TRUE, bdy = adm_bdy_Bergamo,
                  center = center, foldername = "pictures/sample",
                  filename = filename, zoom = 14, cex = 7.5, stroke = 1.5)
    }
    
  }
  
  ## ESTIMATES -----------------------------------------------------------------
  dir.create(paste0(foldername,"pictures/estimates/"))
  
  color_palette <- c("jet.col", "viridis", "magma", "plasma", "inferno")
  color_palette <- "jet.col"
  
  for(col in color_palette){
    
    for(daytype in c("workdays", "holidays")){
      
      if(daytype == "workdays"){
        
        df <- mean_sol_STDEPDE_workdays
        m <- min(min(mean_sol_STDEPDE_workdays, na.rm = TRUE),
                 # min(mean_sol_STLNPP_workdays, na.rm = TRUE),
                 # min(mean_sol_HEAT_workdays, na.rm = TRUE),
                 # min(mean_sol_EQUAL_SPLIT_workdays, na.rm = TRUE),
                 # min(mean_sol_VORONOI_workdays, na.rm = TRUE),
                 # min(mean_sol_VORONOI_SEP_workdays, na.rm = TRUE),
                 na.rm = TRUE)
        
        M <- max(max(mean_sol_STDEPDE_workdays, na.rm = TRUE),
                 # max(mean_sol_STLNPP_workdays, na.rm = TRUE),
                 # max(mean_sol_HEAT_workdays, na.rm = TRUE),
                 # max(mean_sol_EQUAL_SPLIT_workdays, na.rm = TRUE),
                 # max(mean_sol_VORONOI_workdays, na.rm = TRUE),
                 # max(mean_sol_VORONOI_SEP_workdays, na.rm = TRUE),
                 na.rm = TRUE)
        M <- as.numeric(quantile(mean_sol_STDEPDE_workdays, 0.99975))
        
      } else {
        
        df <- mean_sol_STDEPDE_holidays
        m <- min(min(mean_sol_STDEPDE_holidays, na.rm = TRUE),
                 # min(mean_sol_STLNPP_holidays, na.rm = TRUE),
                 # min(mean_sol_HEAT_holidays, na.rm = TRUE),
                 # min(mean_sol_EQUAL_SPLIT_holidays, na.rm = TRUE),
                 # min(mean_sol_VORONOI_holidays, na.rm = TRUE),
                 # min(mean_sol_VORONOI_SEP_holidays, na.rm = TRUE),
                 na.rm = TRUE)
        
        M <- max(max(mean_sol_STDEPDE_holidays, na.rm = TRUE),
                 # max(mean_sol_STLNPP_holidays, na.rm = TRUE),
                 # max(mean_sol_HEAT_holidays, na.rm = TRUE),
                 # max(mean_sol_EQUAL_SPLIT_holidays, na.rm = TRUE),
                 # max(mean_sol_VORONOI_holidays, na.rm = TRUE),
                 # max(mean_sol_VORONOI_SEP_holidays, na.rm = TRUE),
                 na.rm = TRUE)
        M <- as.numeric(quantile(mean_sol_STDEPDE_holidays, 0.99975))
        
      }
      
      dir.create(paste0(foldername,"pictures/estimates/STDEPDE"))
      dir.create(paste0(foldername,"pictures/estimates/STDEPDE/html"))
      dir.create(paste0(foldername,"pictures/estimates/STDEPDE/pdf"))
      
      for(time_index in 1:length(t)){

        #filename <- paste0("app_STDEPDE_", daytype, "_", col, "_", time_index)
        filename <- paste0("app_STDEPDE_", daytype, "_", col, "_zoom_4_", time_index)
        plot(obj = FEM(df[,time_index], FEMbasis.eval),
             map = FALSE, bdy = adm_bdy_Bergamo, center = center,
             colormap = "jet.col", m = m, M = M, showLegend = FALSE,
             foldername = "pictures/estimates/STDEPDE",
             filename = filename, zoom = 14)
        # filename <- paste0("app_STDEPDE_map_", daytype, "_", col, "_", time_index)
        # plot(obj = FEM(df[,time_index], FEMbasis.eval),
        #      map = TRUE, bdy = adm_bdy_Bergamo, center = center,
        #      colormap = "jet.col", m = m, M = M, showLegend = FALSE,
        #      foldername = "pictures/estimates/STDEPDE",
        #      filename = filename)

      }

    }

    dir.create(paste0(foldername,"pictures/estimates/STLNPP"))
    dir.create(paste0(foldername,"pictures/estimates/STLNPP/html"))
    dir.create(paste0(foldername,"pictures/estimates/STLNPP/pdf"))

    for(daytype in c("workdays", "holidays")){

      if(daytype == "workdays"){

        df <- mean_sol_STLNPP_workdays

        m <- min(min(mean_sol_STDEPDE_workdays, na.rm = TRUE),
                 min(mean_sol_STLNPP_workdays, na.rm = TRUE),
                 # min(mean_sol_HEAT_workdays, na.rm = TRUE),
                 # min(mean_sol_EQUAL_SPLIT_workdays, na.rm = TRUE),
                 # min(mean_sol_VORONOI_workdays, na.rm = TRUE),
                 # min(mean_sol_VORONOI_SEP_workdays, na.rm = TRUE),
                 na.rm = TRUE)

        M <- max(max(mean_sol_STDEPDE_workdays, na.rm = TRUE),
                 # max(mean_sol_STLNPP_workdays, na.rm = TRUE),
                 # max(mean_sol_HEAT_workdays, na.rm = TRUE),
                 # max(mean_sol_EQUAL_SPLIT_workdays, na.rm = TRUE),
                 # max(mean_sol_VORONOI_workdays, na.rm = TRUE),
                 # max(mean_sol_VORONOI_SEP_workdays, na.rm = TRUE),
                 na.rm = TRUE)
        M <- as.numeric(quantile(mean_sol_STLNPP_workdays, 0.99))
        
      } else {

        df <- mean_sol_STLNPP_holidays

        m <- min(min(mean_sol_STDEPDE_holidays, na.rm = TRUE),
                 min(mean_sol_STLNPP_holidays, na.rm = TRUE),
                 # min(mean_sol_HEAT_holidays, na.rm = TRUE),
                 # min(mean_sol_EQUAL_SPLIT_holidays, na.rm = TRUE),
                 # min(mean_sol_VORONOI_holidays, na.rm = TRUE),
                 # min(mean_sol_VORONOI_SEP_holidays, na.rm = TRUE),
                 na.rm = TRUE)

        M <- max(max(mean_sol_STDEPDE_holidays, na.rm = TRUE),
                 # max(mean_sol_STLNPP_holidays, na.rm = TRUE),
                 # max(mean_sol_HEAT_holidays, na.rm = TRUE),
                 # max(mean_sol_EQUAL_SPLIT_holidays, na.rm = TRUE),
                 # max(mean_sol_VORONOI_holidays, na.rm = TRUE),
                 # max(mean_sol_VORONOI_SEP_holidays, na.rm = TRUE),
                 na.rm = TRUE)
        M <- as.numeric(quantile(mean_sol_STDEPDE_holidays, 0.999))

      }

      for(time_index in 1:length(t)){

        #filename <- paste0("app_STLNPP_", daytype, "_", col, "_zoom_", time_index)
        filename <- paste0("app_STLNPP_", daytype, "_", col, "_zoom_4_", time_index)
        
        plot(obj = FEM(df[,time_index], FEMbasis.eval),
             map = FALSE, bdy = adm_bdy_Bergamo, center = center,
             colormap = "jet.col", m = m, M = M, showLegend = FALSE,
             foldername = "pictures/estimates/STLNPP",
             filename = filename)
        # filename <- paste0("app_STLNPP_map_", daytype, "_", col, "_", time_index)
        # plot(obj = FEM(df[,time_index], FEMbasis.eval),
        #      map = TRUE, bdy = adm_bdy_Bergamo, center = center,
        #      colormap = "jet.col", m = m, M = M, showLegend = FALSE,
        #      foldername = "pictures/estimates/STLNPP",
        #      filename = filename)

      }

    }
  }
  
  for(neigh in 1:4){
    
    if(neigh == 1){
      center_ <- neighborhood_1
    }
    if(neigh == 2){
      center_ <- neighborhood_2
    }
    if(neigh == 3){
      center_ <- neighborhood_3
    }
    if(neigh == 4){
      center_ <- neighborhood_4
    }
    center_ <- projection.points.1.5D(submesh.eval, center_)
    center_ <- data.frame(lon = center_[1], lat = center_[2])
    
    for(daytype in c("workdays", "holidays")){

      if(daytype == "workdays"){

        df <- data_workdays
        mean_sol_STDEPDE <- mean_sol_STDEPDE_workdays
        mean_sol_STLNPP <- mean_sol_STLNPP_workdays

      } else {

        df <- data_holidays
        mean_sol_STDEPDE <- mean_sol_STDEPDE_holidays
        mean_sol_STLNPP <- mean_sol_STLNPP_holidays

      }

      foldername = "pictures/estimates/"
      filename <- paste0(foldername, "neighborhood_", daytype, "_", neigh, ".pdf")
      plot_neighborhood(df, 4*mean_sol_STDEPDE, 6*mean_sol_STLNPP, FEMbasis.eval,
                        center_, filename)

    }
  
  }

}
  
## CROSS-VALIDATION ERROR ------------------------------------------------------
for(daytype in c("workdays", "holidays")){
  
  if(daytype == "workdays"){
    df <- data_workdays
  } else {
    df <- data_holidays
  }
  
  locations <- cbind(df$lon, df$lat)
  times <- df$times
  
  # Sample Size
  N <- nrow(df)
  
  # Number of Folds
  K <- 10
  
  # Subsample Size in Each Fold
  n_k <- floor(N/K)
  
  # Indices
  folds <- sample(1:N, N, replace = FALSE)
  
  # Cross-Validation Error Containers
  CV_error_STDEPDE <- rep(0,K)
  CV_error_STLNPP <- rep(0,K)
  # CV_error_HEAT <- rep(0,K)
  # CV_error_EQUAL_SPLIT <- rep(0,K)
  # CV_error_VORONOI <- rep(0,K)
  # CV_error_VORONOI_SEP <- rep(0,K)
  
  if(daytype == "workdays"){
    CV_error_STDEPDE_workdays <- rep(0,K)
    CV_error_STLNPP_workdays <- rep(0,K)
    # CV_error_HEAT_workdays <- rep(0,K)
    # CV_error_EQUAL_SPLIT_workdays <- rep(0,K)
    # CV_error_VORONOI_workdays <- rep(0,K)
    # CV_error_VORONOI_SEP_workdays <- rep(0,K)
  } else {
    CV_error_STDEPDE_holidays <- rep(0,K)
    CV_error_STLNPP_holidays <- rep(0,K)
    # CV_error_HEAT_holidays <- rep(0,K)
    # CV_error_EQUAL_SPLIT_holidays <- rep(0,K)
    # CV_error_VORONOI_holidays <- rep(0,K)
    # CV_error_VORONOI_SEP_holidays <- rep(0,K)
  }
  
  for(k in 1:K){
    
    t_iter <- proc.time()
    
    ### PREPROCESSING ------------------------------------------------------------
    if(k != K){
      
      test_idx <- ((k-1)*n_k+1):(k*n_k)
      train_idx <- (1:N)[!((1:N) %in% test_idx)]
      
      locations_train_k <- locations[folds[train_idx],]
      times_train_k <- times[folds[train_idx]]
      
      locations_test_k <- locations[folds[test_idx],]
      times_test_k <- times[folds[test_idx]]
      
    } else {
      
      test_idx <- ((k-1)*n_k+1):N
      train_idx <- (1:N)[!((1:N) %in% test_idx)]
      
      locations_train_k <- locations[folds[train_idx],]
      times_train_k <- times[folds[train_idx]]
      
      locations_test_k <- locations[folds[test_idx],]
      times_test_k <- times[folds[test_idx]]
      
    }
    
    ## STDE-PDE: Spatio-Temporal Density Estimation with PDE Regulariz. ----------
    lambda <- 0.01
    # alternative: lambda <- 10^seq(from = -2, to = 0, by = 1)
    
    lambda_time <- 0.01
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
    
    ## STLNPP: Spatio-Temporal Kernel Intensity Estimation ---------------------
    observations_stlpp <- as.stlpp(locations_train_k[,1], locations_train_k[,2], times_train_k, mesh_linnet)
    solution_STLNPP_train <- density(observations_stlpp, dimyx = n, dimt = n,
                                     diggle = TRUE, at = "pixels", verbose = FALSE)
    
    tgrid_STLNPP_train <- attr(solution_STLNPP_train, "tempden")$x
    
    # ## HEAT: Spatio-Temporal Separable Kernel Intensity Estimation based on the Heat Equation -----
    # observations_stlpp <- as.stlpp(locations_train_k[,1], locations_train_k[,2], times_train_k, mesh_linnet)
    # solution_HEAT_train <- density1Dkernel(observations_stlpp, at = "pixels",
    #                                        leaveoneout = FALSE, iterMax = 1e+09,
    #                                        verbose = FALSE)
    # 
    # tgrid_HEAT_train <- attr(solution_HEAT_train, "tgrid")
    # 
    # ## EQUAL-SPLIT: Spatio-Temporal Separable Kernel Intensity Estimation based on the the Okabe-Sugihara Equal-Split Algorithm -----
    # observations_stlpp <- as.stlpp(locations_train_k[,1], locations_train_k[,2], times_train_k, mesh_linnet)
    # solution_EQUAL_SPLIT_train <- density1Dkernel(observations_stlpp, at = "pixels",
    #                                               leaveoneout = FALSE, iterMax = 1e+09,
    #                                               kernel = "epanechnikov", verbose = FALSE)
    # 
    # tgrid_EQUAL_SPLIT_train <- attr(solution_EQUAL_SPLIT_train, "tgrid")
    # 
    # ## VORONOI: Spatio-Temporal Pseudo-Separable Intensity Estimation based on the Voronoi-Dirichlet Tessellation -----
    # observations_stlpp <- as.stlpp(locations_train_k[,1], locations_train_k[,2], times_train_k, mesh_linnet)
    # invisible(capture.output(solution_VORONOI_train <- densityVoronoi(observations_stlpp, f = 0.9, nrep = 10, separable = FALSE,
    #                                                                   dimt = n, at = "pixels")))
    # 
    # tgrid_VORONOI_train <- attr(solution_VORONOI_train, "tgrid")
    # 
    # ## VORONOI-SEP: Spatio-Temporal Separable Intensity Estimation based on the Voronoi-Dirichlet Tessellation -----
    # observations_stlpp <- as.stlpp(locations_train_k[,1], locations_train_k[,2], times_train_k, mesh_linnet)
    # invisible(capture.output(solution_VORONOI_SEP_train <- densityVoronoi(observations_stlpp, f = 1, nrep = 10, separable = TRUE,
    #                                                                       dimt = n, at = "pixels")))
    # 
    # tgrid_VORONOI_SEP_train <- attr(solution_VORONOI_SEP_train, "tgrid")
    
    ## ERRORS --------------------------------------------------------------------
    # Time half width
    h <- 1/8
    
    # Time instants at which the solution is evaluated (midpoints of each month)
    t <- seq(from = 0+h/2, to = 1-h/2, by = h)
    
    for(time_index in 1:length(t)){
      
      # STDEPDE
      evaluation_STDEPDE_train <- eval.FEM.time(FEM.time = FEMfunction_STDEPDE_train, locations = mesh.eval$nodes, time.instants = t[time_index])
      evaluation_STDEPDE_train <- exp(evaluation_STDEPDE_train)*nrow(locations_train_k)
      evaluation_STDEPDE_train[idx_remove] <- NA
      evaluation_STDEPDE_train <- evaluation_STDEPDE_train/sum(evaluation_STDEPDE_train, na.rm = TRUE)
      
      evaluation_STDEPDE_test <- eval.FEM(FEM = FEM(evaluation_STDEPDE_train, FEMbasis.eval), locations = locations_test_k[which(abs(t[time_index]-times_test_k)<h),])
      
      cardinality_k <- 1 / (sum(!is.na(evaluation_STDEPDE_test)) / nrow(locations_test_k[which(abs(t[time_index]-times_test_k)<h),]))
      CV_error_STDEPDE[k] <- CV_error_STDEPDE[k] + sum(evaluation_STDEPDE_train^2, na.rm = TRUE)/sum(!is.na(evaluation_STDEPDE_train)) -2 * sum(evaluation_STDEPDE_test, na.rm = TRUE) * cardinality_k
      
      # STLNPP
      idx <- which.min(abs(tgrid_STLNPP_train - t[time_index]))
      evaluation_STLNPP_train <- as.linfun(solution_STLNPP_train[[idx]])(mesh.eval$nodes[,1], mesh.eval$nodes[,2])
      
      evaluation_STLNPP_test <- eval.FEM(FEM = FEM(evaluation_STLNPP_train, FEMbasis.eval), locations = locations_test_k[which(abs(t[time_index]-times_test_k)<h),])
      
      cardinality_k <- 1 / (sum(!is.na(evaluation_STLNPP_test)) / nrow(locations_test_k[which(abs(t[time_index]-times_test_k)<h),]))
      CV_error_STLNPP[k] <- CV_error_STLNPP[k] + sum(evaluation_STLNPP_train^2, na.rm = TRUE)/sum(!is.na(evaluation_STLNPP_train)) -2 * sum(evaluation_STLNPP_test, na.rm = TRUE) * cardinality_k
      
      # # HEAT
      # idx <- which.min(abs(tgrid_HEAT_train - t[time_index]))
      # evaluation_HEAT_train <- as.linfun(solution_HEAT_train[[idx]])(mesh.eval$nodes[,1], mesh.eval$nodes[,2])
      # 
      # evaluation_HEAT_test <- eval.FEM(FEM = FEM(evaluation_HEAT_train, FEMbasis.eval), locations = locations_test_k[which(abs(t[time_index]-times_test_k)<h),])
      # 
      # cardinality_k <- 1 / (sum(!is.na(evaluation_HEAT_test)) / nrow(locations_test_k[which(abs(t[time_index]-times_test_k)<h),]))
      # CV_error_HEAT[k] <- CV_error_HEAT[k] + sum(evaluation_HEAT_train^2, na.rm = TRUE)/sum(!is.na(evaluation_HEAT_train)) -2 * sum(evaluation_HEAT_test, na.rm = TRUE) * cardinality_k
      # 
      # # EQUAL-SPLIT
      # idx <- which.min(abs(tgrid_EQUAL_SPLIT_train - t[time_index]))
      # evaluation_EQUAL_SPLIT_train <- as.linfun(solution_EQUAL_SPLIT_train[[idx]])(mesh.eval$nodes[,1], mesh.eval$nodes[,2])
      # 
      # evaluation_EQUAL_SPLIT_test <- eval.FEM(FEM = FEM(evaluation_EQUAL_SPLIT_train, FEMbasis.eval), locations = locations_test_k[which(abs(t[time_index]-times_test_k)<h),])
      # 
      # cardinality_k <- 1 / (sum(!is.na(evaluation_EQUAL_SPLIT_test)) / nrow(locations_test_k[which(abs(t[time_index]-times_test_k)<h),]))
      # CV_error_EQUAL_SPLIT[k] <- CV_error_EQUAL_SPLIT[k] + sum(evaluation_EQUAL_SPLIT_train^2, na.rm = TRUE)/sum(!is.na(evaluation_EQUAL_SPLIT_train)) -2 * sum(evaluation_EQUAL_SPLIT_test, na.rm = TRUE) * cardinality_k
      # 
      # # VORONOI
      # idx <- which.min(abs(tgrid_VORONOI_train - t[time_index]))
      # evaluation_VORONOI_train <- as.linfun(solution_VORONOI_train[[idx]])(mesh.eval$nodes[,1], mesh.eval$nodes[,2])
      # 
      # evaluation_VORONOI_test <- eval.FEM(FEM = FEM(evaluation_VORONOI_train, FEMbasis.eval), locations = locations_test_k[which(abs(t[time_index]-times_test_k)<h),])
      # 
      # cardinality_k <- 1 / (sum(!is.na(evaluation_VORONOI_test)) / nrow(locations_test_k[which(abs(t[time_index]-times_test_k)<h),]))
      # CV_error_VORONOI[k] <- CV_error_VORONOI[k] + sum(evaluation_VORONOI_train^2, na.rm = TRUE)/sum(!is.na(evaluation_VORONOI_train)) -2 * sum(evaluation_VORNOI_test, na.rm = TRUE) * cardinality_k
      # 
      # # VORONOI-SEP
      # idx <- which.min(abs(tgrid_VORONOI_SEP_train - t[time_index]))
      # evaluation_VORONOI_SEP_train <- as.linfun(solution_VORONOI_SEP[[idx]])(mesh.eval$nodes[,1], mesh.eval$nodes[,2])
      # 
      # evaluation_VORONOI_SEP_test <- eval.FEM(FEM = FEM(evaluation_VORONOI_SEP_train, FEMbasis.eval), locations = locations_test_k[which(abs(t[time_index]-times_test_k)<h),])
      # 
      # cardinality_k <- 1 / (sum(!is.na(evaluation_VORONOI_SEP_test)) / nrow(locations_test_k[which(abs(t[time_index]-times_test_k)<h),]))
      # CV_error_VORONOI_SEP[k] <- CV_error_VORONOI_SEP[k] + sum(evaluation_VORONOI_SEP_train^2, na.rm = TRUE)/sum(!is.na(evaluation_VORONOI_SEP_train)) -2 * sum(evaluation_VORNOI_SEP_test, na.rm = TRUE) * cardinality_k
      
    }
    
    CV_error_STDEPDE[k] <- CV_error_STDEPDE[k]  / length(t)
    CV_error_STLNPP[k] <- CV_error_STLNPP[k] / length(t)
    # CV_error_HEAT[k] <- CV_error_HEAT[k] / length(t)
    # CV_error_EQUAL_SPLIT[k] <- CV_error_EQUAL_SPLIT[k] / length(t)
    # CV_error_VORONOI[k] <- CV_error_VORONOI[k] / length(t)
    # CV_error_VORONOI_SEP[k] <- CV_error_VORONOI_SEP[k] / length(t)
    
    print(paste("CV iteration", k, "done in", round(as.numeric(proc.time() - t_iter)[3]), "seconds."))
    
  }
  
  if(daytype == "workdays"){
    CV_error_STDEPDE_workdays <- CV_error_STDEPDE
    CV_error_STLNPP_workdays <- CV_error_STLNPP
    # CV_error_HEAT_workdays <- CV_error_HEAT
    # CV_error_EQUAL_SPLIT_workdays <- CV_error_EQUAL_SPLIT
    # CV_error_VORONOI_workdays <- CV_error_VORONOI
    # CV_error_VORONOI_SEP_workdays <- CV_error_VORONOI_SEP
  } else {
    CV_error_STDEPDE_holidays <- CV_error_STDEPDE
    CV_error_STLNPP_holidays <- CV_error_STLNPP
    # CV_error_HEAT_holidays <- CV_error_HEAT
    # CV_error_EQUAL_SPLIT_holidays <- CV_error_EQUAL_SPLIT
    # CV_error_VORONOI_holidays <- CV_error_VORONOI
    # CV_error_VORONOI_SEP_holidays <- CV_error_VORONOI_SEP
  }
  
}

# dir.create(paste0(foldername, "output"))
# save.image("output/workspace_cv_errors.RData")

## BOXPLOT ---------------------------------------------------------------------
# CV_error_STDEPDE_workdays <- 50+c(5, 7, 8, 9, 10, 12, 14, 15, 18, 20)
# CV_error_STLNPP_workdays <- 90+1.5*c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
# CV_error_HEAT_workdays <- 10+0.9*c(120, 105, 110, 115, 120, 125, 130, 125, 130, 135)
# CV_error_EQUAL_SPLIT_workdays <- 120+c(-10, -8, -6, -4, -2, 0, 2, 4, 6, 8)
# CV_error_VORONOI_workdays <- 100+c(50, 55, 60, 55, 70, 75, 80, 90, 100, 105)
# CV_error_VORONOI_SEP_workdays <- 100+c(25, 30, 40, 45, 50, 70, 75, 90, 100, 120)
# 
# CV_error_STDEPDE_holidays <- 60+c(-5, 7, 21, 9, 10, 12, 14, 15, 18, 20)
# CV_error_STLNPP_holidays <- 110+1.5*c(-1, 2, 13, 4, 5, 6, 7, 8, 9, 10)
# CV_error_HEAT_holidays <- 20+0.9*c(110, 105, 110, 115, 120, 125, 130, 125, 120, 125)
# CV_error_EQUAL_SPLIT_holidays <- 115+c(-8, -6, -4, -2, 0, 3, 2, 4, 6, 6)
# CV_error_VORONOI_holidays <- 115+c(80, 75, 80, 55, 70, 75, 80, 90, 90, 95)
# CV_error_VORONOI_SEP_holidays <- 120+c(15, 30, 40, 45, 50, 70, 75, 90, 100, 120)

{
  for(daytype in c("workdays", "holidays")){
    
    if(daytype == "workdays"){
      CV_error_STDEPDE <- CV_error_STDEPDE_workdays
      CV_error_STLNPP <- CV_error_STLNPP_workdays
      # CV_error_HEAT <- CV_error_STLNPP_workdays
      # CV_error_EQUAL_SPLIT <- CV_error_EQUAL_SPLIT_workdays
      # CV_error_VORONOI <- CV_error_VORONOI_workdays
      # CV_error_VORONOI_SEP <- CV_error_VORONOI_SEP_workdays
    } else {
      CV_error_STDEPDE <- CV_error_STDEPDE_holidays
      CV_error_STLNPP <- CV_error_STLNPP_holidays
      # CV_error_HEAT <- CV_error_HEAT_holidays
      # CV_error_EQUAL_SPLIT <- CV_error_EQUAL_SPLIT_holidays
      # CV_error_VORONOI <- CV_error_VORONOI_holidays
      # CV_error_VORONOI_SEP <- CV_error_VORONOI_SEP_holidays
    }
    
    blue_polimi <- rgb(114, 143, 165, max = 255, alpha = 65, names = "blue_polimi25") # #728FA5
    
    pdf(paste0("pictures/app_CV_errors_", daytype, ".pdf"), family = "serif", width = 12, height = 5.3)
    
    #x11(width = 12, height = 5.3)
    par(mfrow = c(1,3), mai = c(1.75,0.35,0.5,0.15))
    plot.new()
    boxplot(CV_error_STDEPDE, CV_error_STLNPP, #CV_error_HEAT,
            #CV_error_EQUAL_SPLIT, CV_error_VORONOI, CV_error_VORONOI_SEP,
            asp = 1, xaxt = "n", yaxt = "n")
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = blue_polimi)
    grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 1)
    par(new = TRUE)
    boxplot(CV_error_STDEPDE, CV_error_STLNPP, #CV_error_HEAT,
            #CV_error_EQUAL_SPLIT, CV_error_VORONOI, CV_error_VORONOI_SEP,
            names = c("STDE-PDE", "STKDE-QUICK", "STKDE-HEAT", "STKDE-EPAN", "STVDE", "STVDE-SEP"),
            main = ifelse(daytype == "workdays", "Working days", "Weekends and holidays"),
            col = brewer.pal(6, "YlGnBu")[c(1,4)], cex.lab = 2, cex.axis = 2, cex.main = 2, las = 3)
    
    dev.off()
    
  }
}

# Total length of the road network
sum(st_length(mesh.eval_sf))