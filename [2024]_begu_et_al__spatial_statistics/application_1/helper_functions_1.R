#'############################################################################'#
#'##################### SPATIO-TEMPORAL DENSITY ESTIMATION ###################'# 
#'##################### HELPER FUNCTIONS FOR APPLICATION 1 ###################'#
#'############################################################################'#

## R-INLA UTILITY --------------------------------------------------------------
book.mesh.dual <- function(mesh) {
  if (mesh$manifold == "R2") {
    ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
      colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
    library(parallel)
    pls <- mclapply(1:mesh$n, function(i) {
      p <- unique(Reduce("rbind", lapply(1:3, function(k) {
        j <- which(mesh$graph$tv[,k] == i)
        if (length(j) > 0) 
          return(rbind(ce[j, , drop = FALSE],
                       cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 1], 
                             mesh$loc[mesh$graph$tv[j, k], 2] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 2])/2))
        else return(ce[j, , drop = FALSE])
      })))
      j1 <- which(mesh$segm$bnd$idx[,1] == i)
      j2 <- which(mesh$segm$bnd$idx[,2] == i)
      if ((length(j1)>0) | (length(j2)>0)) {
        p <- unique(rbind(mesh$loc[i, 1:2], p,
                          mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2]/2, 
                          mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2]/2))
        yy <- p[,2]-mean(p[,2])/2-mesh$loc[i, 2]/2
        xx <- p[,1]-mean(p[,1])/2-mesh$loc[i, 1]/2
      }
      else {
        yy <- p[,2]-mesh$loc[i, 2]
        xx <- p[,1]-mesh$loc[i, 1]
      }
      Polygon(p[order(atan2(yy,xx)), ])
    })
    return(SpatialPolygons(lapply(1:mesh$n, function(i)
      Polygons(list(pls[[i]]), i))))
  }
  else stop("It only works for R2!")
}

## PLOT A MESH ON A MAP --------------------------------------------------------
Figure.MeshOnMap <- function(mesh, boundary, map, col, cex = 1){
  
  # Convert the Mesh
  mesh_df <- as.data.frame(mesh$nodes[,2:1])
  colnames(mesh_df) <- c("Easting", "Northing")
  mesh_sp <- mesh_df
  coordinates(mesh_sp)<-~Northing+Easting
  proj4string(mesh_sp)<-CRS("+init=epsg:27700")
  mesh_sp <- spTransform(mesh_sp, CRS("+init=epsg:4326"))
  mesh_latlong <- as.data.frame(mesh_sp@coords)
  colnames(mesh_latlong) <- c("long","lat")
  
  # Plot the Map
  par(mar = c(0, 0, 0, 0))
  PlotOnStaticMap(map, lat = 0, lon = 0, pch = 19, col = col,
                  add = FALSE, TrueProj = TRUE, cex = cex, size = c(1280, 1280))
  
  # Color the Region of Interest
  in_p <- sp::Polygon(rbind(boundary, boundary[1,]), hole = TRUE)
  in_sp <- sp::SpatialPolygons(list(sp::Polygons(list(in_p), ID = "IN")))
  
  out_p <- sp::Polygon(coords = rbind(map$BBOX$ll,
                                      cbind(map$BBOX$ll[1], map$BBOX$ur[2]),
                                      map$BBOX$ur,
                                      cbind(map$BBOX$ur[1], map$BBOX$ll[2]),
                                      map$BBOX$ll)[,2:1],
                       hole = FALSE)
  out_sp <- sp::SpatialPolygons(list(sp::Polygons(list(out_p), ID = "OUT")))
  
  sp <- rgeos::gDifference(out_sp, in_sp)
  
  #mycol <- rgb(114, 143, 165, max = 255, alpha=125)
  mycol <- rgb(110, 110, 110, max = 255, alpha=150)
  adjustcolor(mycol, alpha.f = 0.3)
  
  PlotPolysOnStaticMap(map, in_sp, border = NULL, col = mycol, 
                       verbose = 0, add = TRUE)
  
  # Overlay the Mesh
  ranPts <- 1:nrow(mesh_latlong)
  PlotOnStaticMap(map, lat = mesh_latlong$lat, lon = mesh_latlong$lon, pch = 19, col = col,
                  add = TRUE, TrueProj = TRUE, cex = cex, size = c(1280, 1280))
  
  PlotOnStaticMap(map, lat = c(boundary$lat, boundary$lat[1]), lon = c(boundary$long, boundary$long[1]), type = "l", lwd = 3, col = "black",
                  add = TRUE, TrueProj = TRUE, size = c(1280, 1280))
  
  for(i in 1:nrow(mesh$triangles)){
    PlotOnStaticMap(map,
                    lat = c(mesh_latlong$lat[mesh$triangles[i,1]], mesh_latlong$lat[mesh$triangles[i,2]]),
                    lon = c(mesh_latlong$long[mesh$triangles[i,1]], mesh_latlong$long[mesh$triangles[i,2]]),
                    type = "l", lwd = 1.75, col = "black",
                    add = TRUE, TrueProj = TRUE, size = c(1280, 1280))
    
    PlotOnStaticMap(map,
                    lat = c(mesh_latlong$lat[mesh$triangles[i,1]], mesh_latlong$lat[mesh$triangles[i,3]]),
                    lon = c(mesh_latlong$long[mesh$triangles[i,1]], mesh_latlong$long[mesh$triangles[i,3]]),
                    type = "l", lwd = 1.75, col = "black",
                    add = TRUE, TrueProj = TRUE, size = c(1280, 1280))
    
    PlotOnStaticMap(map,
                    lat = c(mesh_latlong$lat[mesh$triangles[i,2]], mesh_latlong$lat[mesh$triangles[i,3]]),
                    lon = c(mesh_latlong$long[mesh$triangles[i,2]], mesh_latlong$long[mesh$triangles[i,3]]),
                    type = "l", lwd = 1.75, col = "black",
                    add = TRUE, TrueProj = TRUE, size = c(1280, 1280))
  }
  
}

## PLOT THE DATASET ON A MAP ---------------------------------------------------
Figure.DatasetOnMap <- function(data, boundary, map, col, cex = 1){
  
  # Plot the Map
  par(mar = c(0, 0, 0, 0))
  PlotOnStaticMap(map, lat = 0, lon = 0, pch = 19, col = col,
                  add = FALSE, TrueProj = TRUE, cex = cex, size = c(1280, 1280))
  
  # Color the Region of Interest
  in_p <- sp::Polygon(rbind(boundary, boundary[1,]), hole = TRUE)
  in_sp <- sp::SpatialPolygons(list(sp::Polygons(list(in_p), ID = "IN")))
  
  out_p <- sp::Polygon(coords = rbind(map$BBOX$ll,
                                      cbind(map$BBOX$ll[1], map$BBOX$ur[2]),
                                      map$BBOX$ur,
                                      cbind(map$BBOX$ur[1], map$BBOX$ll[2]),
                                      map$BBOX$ll)[,2:1],
                       hole = FALSE)
  out_sp <- sp::SpatialPolygons(list(sp::Polygons(list(out_p), ID = "OUT")))
  
  sp <- rgeos::gDifference(out_sp, in_sp)
  
  #mycol <- rgb(114, 143, 165, max = 255, alpha=125)
  mycol <- rgb(110, 110, 110, max = 255, alpha=150)
  adjustcolor(mycol, alpha.f = 0.3)
  PlotPolysOnStaticMap(map, in_sp, border = NULL, col = mycol,
                       verbose = 0, add = T)
  
  # Add the Points in the Dataset
  ranPts = 1:nrow(data)
  PlotOnStaticMap(map, lat = data$lat, lon = data$long, pch = 19, col = col,
                  add = TRUE, TrueProj = TRUE, cex = cex, size = c(1280, 1280))
  
  # Add the Boundaries
  PlotOnStaticMap(map, lat = c(boundary$lat, boundary$lat[1]), lon = c(boundary$long, boundary$long[1]), type = "l", lwd = 3, col = "black",
                  add = TRUE, TrueProj = TRUE, size = c(1280, 1280))
  
}

## PLOT A DENSITY ON A MAP -----------------------------------------------------
Figure.DensityOnMap2 <- function(mesh, boundary, map, evaluation, max_range = NULL, col = "black", cex = 1){
  if (is.null(max_range)) { max_range = max(evaluation)}
  
  # Convert the Mesh
  mesh_df <- as.data.frame(mesh$nodes[,2:1])
  colnames(mesh_df) <- c("Easting", "Northing")
  mesh_sp <- mesh_df
  coordinates(mesh_sp)<-~Northing+Easting
  proj4string(mesh_sp)<-CRS("+init=epsg:27700")
  mesh_sp <- spTransform(mesh_sp, CRS("+init=epsg:4326"))
  mesh_latlong <- as.data.frame(mesh_sp@coords)
  colnames(mesh_latlong) <- c("long","lat")

  # Plot the Map
  par(mar = c(0, 0, 0, 0))
  PlotOnStaticMap(map, lat = 0, lon = 0, pch = 19, col = col,
                  add = FALSE, TrueProj = TRUE, cex = cex, size = c(1280, 1280))
  
  # Color the Region of Interest
  # in_p <- sp::Polygon(rbind(boundary, boundary[1,]), hole = TRUE)
  # in_sp <- sp::SpatialPolygons(list(sp::Polygons(list(in_p), ID = "IN")))
  # 
  # out_p <- sp::Polygon(coords = rbind(map$BBOX$ll,
  #                                     cbind(map$BBOX$ll[1], map$BBOX$ur[2]),
  #                                     map$BBOX$ur,
  #                                     cbind(map$BBOX$ur[1], map$BBOX$ll[2]),
  #                                     map$BBOX$ll)[,2:1],
  #                      hole = FALSE)
  # out_sp <- sp::SpatialPolygons(list(sp::Polygons(list(out_p), ID = "OUT")))
  # 
  # sp <- rgeos::gDifference(out_sp, in_sp)
  # 
  # #mycol <- rgb(114, 143, 165, max = 255, alpha = 125)
  # mycol <- rgb(110, 110, 110, max = 255, alpha = 150)
  # adjustcolor(mycol, alpha.f = 0.3)
  # 
  # PlotPolysOnStaticMap(map, in_sp, border = NULL, col = mycol, 
  #                      verbose = 0, add = TRUE)
  
  # Overlay the Mesh Nodes
  # ranPts = 1:nrow(mesh_latlong)
  # PlotOnStaticMap(map, lat = mesh_latlong$lat, lon = mesh_latlong$lon, pch = 19, col = col,
  #                 add = TRUE, TrueProj = TRUE, cex = cex, size = c(1280, 1280))
  # 
  
  # Add the Boundaries
  PlotOnStaticMap(map, lat = c(boundary$lat, boundary$lat[1]), lon = c(boundary$long, boundary$long[1]), type = "l", lwd = 3, col = "black",
                  add = TRUE, TrueProj = TRUE, size = c(1280, 1280))
  
  # Color Palette
  p <- colorRampPalette(c("#0E1E44", "#3E6DD8", "#68D061", "#ECAF53", "#EB5F5F", "#E11F1C"))(1000)
  palette(p)
  
  # Assign Colors to the Mesh Triangles based on the Values of the Mesh Nodes
  dens_col <- rep(0, nrow(mesh$triangles))
  mean_eval <- rep(0, nrow(mesh$triangles))
  
  for(i in 1:nrow(mesh$triangles)){
    mean_eval[i] <- (evaluation[mesh$triangles[i,1]] + evaluation[mesh$triangles[i,2]] + evaluation[mesh$triangles[i,3]])/3
  }
  
  for(i in 1:nrow(mesh$triangles)){
    n_col <- round(1000*(mean_eval[i]-0)/(max_range-0),0)+1
    if(n_col > 999){
      n_col <- 1000
    }
    dens_col[i] <- p[n_col]
  }
  
  # Plot the Density
  for(i in 1:nrow(mesh$triangles)){

    in_p <- sp::Polygon(rbind(cbind(mesh_latlong$long[mesh$triangles[i,1]], mesh_latlong$lat[mesh$triangles[i,1]]),
                              cbind(mesh_latlong$long[mesh$triangles[i,2]], mesh_latlong$lat[mesh$triangles[i,2]]),
                              cbind(mesh_latlong$long[mesh$triangles[i,3]], mesh_latlong$lat[mesh$triangles[i,3]]),
                              cbind(mesh_latlong$long[mesh$triangles[i,1]], mesh_latlong$lat[mesh$triangles[i,1]])))
    
    in_sp <- sp::SpatialPolygons(list(sp::Polygons(list(in_p), ID = "IN")))
    
    PlotPolysOnStaticMap(map, in_sp, col = dens_col[i], border = NA, add = TRUE,
                         TrueProj = TRUE, size = c(1280, 1280))
    
  }

}

## PLOT A DENSITY ON A MAP (ALTERNATIVE) ---------------------------------------
Figure.DensityOnMap <- function(evaluation, grid_nodes, M = NULL){
  if(is.null(M)) {M = max(evaluation, na.rm = TRUE)}
  
  # Convert the Mesh
  grid_nodes <- as.matrix(grid_nodes*1000)
  grid_df <- as.data.frame(grid_nodes[,2:1])
  colnames(grid_df) <- c("Easting", "Northing")
  grid_sp <- grid_df
  coordinates(grid_sp)<-~Northing+Easting
  proj4string(grid_sp)<-CRS("+init=epsg:27700")
  grid_sp <- spTransform(grid_sp, CRS("+init=epsg:4326"))
  
  grid_latlong <- as.data.frame(grid_sp@coords)
  colnames(grid_latlong) <- c("long","lat")
  
  # Get Hypsometric and Bathymetric Data from GEBCO Server
  bath <- readGEBCO.bathy("map/GEBCO/gebco_2021_n51.5_s50.625_w-1.98_e-0.69.nc")
  
  # Color Palettes for Sea and Land
  blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
  greys <- c(grey(0.6), grey(0.93), grey(0.99))
  # alternative: greys <- colorspace::terrain_hcl(12, c = c(65,0), power = c(1/3, 1.5))
  
  # Plot the Hypsometric/Bathymetric Map
  par(mar = c(0, 0, 0, 0))
  plot(bath, n = 1, land = T, im = T, lwd = 0.05, xlim = c(-1.98,-0.69), ylim = c(50.65, 51.45),
       axes = FALSE, xlab = "", ylab = "", asp = 1.475,
       bpal = list(c(0, max(bath), greys), c(min(bath), 0, blues)))
  #points(rbind(boundary_latlong, boundary_latlong[1,]), type = "l", lwd = 3)
  
  # Transform Data into a Bathy Object
  n <- sqrt(dim(grid_nodes))[1]
  rownames(evaluation) <- grid_latlong[1:n,1]
  colnames(evaluation) <- grid_latlong[seq(1,dim(grid_latlong)[1], by = n),2]
  class(evaluation) <- "bathy"
  
  # Color Palette
  p <- colorRampPalette(c("#0E1E44", "#3E6DD8", "#68D061", "#ECAF53", "#EB5F5F", "#E11F1C"))(1000)
  palette(p)
  
  # Plot the Density
  plot(evaluation, land = TRUE, image = TRUE, lwd = 1, lty = 1,
       bpal = col2alpha(palette(p), .9),
       add = TRUE, drawlabels = FALSE, deep = quantile(evaluation[!is.na(evaluation)], 0.9),
       shallow = M, step = quantile(evaluation[!is.na(evaluation)], 0.99),
       col = "black") # use deep, shallow, and step to adjust contour lines
  plot(outline.buffer(evaluation), add = TRUE, n = 1, lwd = 4) # Outline of the data
  
  # Add Cities Locations and Names
  # map.cities(country = "UK", label = TRUE, minpop = 50000, cex = 3)
  
}
