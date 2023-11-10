#'############################################################################'#
#'##################### SPATIO-TEMPORAL DENSITY ESTIMATION ###################'# 
#'##################### HELPER FUNCTIONS FOR APPLICATION 2 ###################'#
#'############################################################################'#

## GRAPHICAL SETTINGS ----------------------------------------------------------
zoom = 0.6
userMatrix = rbind(c( 0.96563137,  0.1774523, -0.1899119,    0),
                   c(-0.03294301,  0.8083354,  0.5877997,    0),
                   c( 0.25781897, -0.5613416,  0.7863999,    0),
                   c( 0.00000000,  0.0000000,  0.0000000 ,   1))
windowRect = c(150,  150, 420,  420)

#pp <- par3d(no.readonly = TRUE)

userMatrix1 = rbind(c ( 3.261129e-01,  0.9285929,    0.17710334,    0),
                    c ( 9.456277e-05, -0.1873772,    0.98228776,    0),
                    c ( 9.453306e-01, -0.3203200,   -0.06119408,    0),
                    c ( 0.000000e+00,  0.0000000,    0.00000000,    1))

userMatrix2 = rbind(c (-0.5708339,    -0.81879389,   0.06103384,    0),
                    c ( 0.1302170,    -0.01688673,   0.99134153,    0),
                    c (-0.8106737,     0.57383895,   0.11626063,    0),
                    c ( 0.0000000,     0.00000000,   0.00000000,    1))

userMatrix3 = rbind(c ( 0.9270199,     0.3670499,   -0.07686605,    0),
                    c ( 0.1287258,    -0.1189351,    0.98452210,    0),
                    c ( 0.3522265,    -0.9225665,   -0.15750417,    0),
                    c ( 0.0000000,     0.0000000,    0.00000000,    1))

userMatrix4 = rbind(c (-0.93542612,    0.3507093,    0.04450328,    0),
                    c (-0.07357776,   -0.3162677,    0.94581199,    0),
                    c ( 0.34578037,    0.8814633,    0.32164946,    0),
                    c ( 0.00000000,    0.0000000,    0.00000000,    1))

## PLOT A mesh OBJECT ----------------------------------------------------------
plot.mesh <- function(mesh, view = userMatrix, M = NULL, m = NULL, ...){
  if(!is.matrix(view)){
    userMatrix = eval(parse(text = view))
  }
  
  FEMbasis <- create.FEM.basis(mesh = mesh)
  FEM <- FEM(rep(0, nrow(mesh$nodes)), FEMbasis)
  
  if (is.null(m)) { m = min(FEM$coeff)}
  if (is.null(M)) { M = max(FEM$coeff)}
  triangles = c(t(FEM$FEMbasis$mesh$triangles))
  ntriangles = nrow(FEM$FEMbasis$mesh$triangles)
  order = FEM$FEMbasis$mesh$order
  nodes = FEM$FEMbasis$mesh$nodes
  edges = matrix(rep(0, 6*ntriangles), ncol = 2)
  for(i in 0:(ntriangles-1)){
    edges[3*i+1,] = c(triangles[3*order*i+1], triangles[3*order*i+2])
    edges[3*i+2,] = c(triangles[3*order*i+1], triangles[3*order*i+3])
    edges[3*i+3,] = c(triangles[3*order*i+2], triangles[3*order*i+3])
  }
  edges = edges[!duplicated(edges),]
  edges <- as.vector(t(edges))
  
  coeff = FEM$coeff
  
  FEMbasis = FEM$FEMbasis
  
  mesh = FEMbasis$mesh
  
  #p = jet.col(n = 1000, alpha = 0.8)
  p <- colorRampPalette(c("#0E1E44","#3E6DD8","#68D061","#ECAF53", "#EB5F5F","#E11F1C"))(1000)
  palette(p)
  
  ncolor = length(p)
  
  nsurf = dim(coeff)[[2]]
  for (isurf in 1:nsurf)
  {
    open3d(zoom = zoom, userMatrix = userMatrix, windowRect = windowRect)
    rgl.pop("lights") 
    light3d(specular = "black") 
    
    xyz <- shapefile_lines()
    rgl.linestrips(x = xyz[,1], y = xyz[,2], z = xyz[,3], color = "royalblue")
    
    diffrange = M - m
    
    col = coeff[triangles,isurf]
    col = (col - min(coeff[,isurf]))/diffrange*(ncolor-1)+1
    
    rgl.triangles(x = nodes[triangles ,1], y = nodes[triangles ,2],
                  z = nodes[triangles,3],
                  color = col,...)
    rgl.lines(x = nodes[edges ,1], y = nodes[edges ,2],
              z = nodes[edges,3],
              color = "black",...)
    aspect3d("iso")
    
    if (nsurf > 1 && isurf < nsurf)
    {readline("Press a button for the next plot...")}
  }
}


## PLOT A FEM OBJECT WITH JET COLORMAP WITH RANGE [m,M] ------------------------
plot.FEM <- function(FEM, view = userMatrix, M = NULL, m = NULL, ...){
  if(!is.matrix(view)){
    userMatrix = eval(parse(text = view))
  }
  
  if (is.null(m)) { m = min(FEM$coeff)}
  if (is.null(M)) { M = max(FEM$coeff)}
  triangles = c(t(FEM$FEMbasis$mesh$triangles))
  ntriangles = nrow(FEM$FEMbasis$mesh$triangles)
  order = FEM$FEMbasis$mesh$order
  nodes = FEM$FEMbasis$mesh$nodes
  edges = matrix(rep(0, 6*ntriangles), ncol = 2)
  for(i in 0:(ntriangles-1)){
    edges[3*i+1,] = c(triangles[3*order*i+1], triangles[3*order*i+2])
    edges[3*i+2,] = c(triangles[3*order*i+1], triangles[3*order*i+3])
    edges[3*i+3,] = c(triangles[3*order*i+2], triangles[3*order*i+3])
  }
  edges = edges[!duplicated(edges),]
  edges <- as.vector(t(edges))
  
  coeff = FEM$coeff
  
  FEMbasis = FEM$FEMbasis
  
  mesh = FEMbasis$mesh
  
  #p = jet.col(n = 1000, alpha = 0.8)
  p <- colorRampPalette(c("#0E1E44","#3E6DD8","#68D061","#ECAF53", "#EB5F5F","#E11F1C"))(1000)
  palette(p)
  
  ncolor = length(p)
  
  nsurf = dim(coeff)[[2]]
  for (isurf in 1:nsurf)
  {
    open3d(zoom = zoom, userMatrix = userMatrix, windowRect = windowRect)
    rgl.pop("lights") 
    light3d(specular = "black") 
    
    xyz <- shapefile_lines()
    rgl.linestrips(x = xyz[,1], y = xyz[,2], z = xyz[,3], color = "azure4")
    
    diffrange = M - m
    
    col = coeff[triangles,isurf]
    col = (col - min(coeff[,isurf]))/diffrange*(ncolor-1)+1
    
    rgl.triangles(x = nodes[triangles ,1], y = nodes[triangles ,2],
                  z = nodes[triangles,3],
                  color = col,...)
    # rgl.lines(x = nodes[edges ,1], y = nodes[edges ,2],
    #           z = nodes[edges,3],
    #           color = "black",...)
    aspect3d("iso")
    
    if (nsurf > 1 && isurf < nsurf)
    {readline("Press a button for the next plot...")}
  }
}

## PLOT THE SAMPLE ON THE EARTH AT A GIVEN TIME --------------------------------
plot.sample <- function(FEM, sample, view = userMatrix, M = NULL, m = NULL, ...){
  if(!is.matrix(view)){
    userMatrix = eval(parse(text = view))
  }
  
  if (is.null(m)) { m = min(FEM$coeff)}
  if (is.null(M)) { M = max(FEM$coeff)}
  triangles = c(t(FEM$FEMbasis$mesh$triangles))
  ntriangles = nrow(FEM$FEMbasis$mesh$triangles)
  order = FEM$FEMbasis$mesh$order
  nodes = FEM$FEMbasis$mesh$nodes
  edges = matrix(rep(0, 6*ntriangles), ncol = 2)
  for(i in 0:(ntriangles-1)){
    edges[3*i+1,] = c(triangles[3*order*i+1], triangles[3*order*i+2])
    edges[3*i+2,] = c(triangles[3*order*i+1], triangles[3*order*i+3])
    edges[3*i+3,] = c(triangles[3*order*i+2], triangles[3*order*i+3])
  }
  edges = edges[!duplicated(edges),]
  edges <- as.vector(t(edges))
  
  coeff = FEM$coeff
  
  FEMbasis = FEM$FEMbasis
  
  mesh = FEMbasis$mesh
  
  #p = jet.col(n = 1000, alpha = 0.8)
  p <- colorRampPalette(c("#0E1E44","#3E6DD8","#68D061","#ECAF53", "#EB5F5F","#E11F1C"))(1000)
  palette(p)
  
  ncolor = length(p)
  
  nsurf = dim(coeff)[[2]]
  for (isurf in 1:nsurf)
  {
    open3d(zoom = zoom, userMatrix = userMatrix, windowRect = windowRect)
    rgl.pop("lights") 
    light3d(specular = "black")
    
    xyz <- shapefile_lines()
    rgl.linestrips(x = xyz[,1], y = xyz[,2], z = xyz[,3], color = "grey10")
    
    diffrange = M - m
    
    col = coeff[triangles,isurf]
    col = (col - min(coeff[,isurf]))/diffrange*(ncolor-1)+1
    
    rgl.spheres(x = sample[,1], y = sample[,2], z = sample[,3], radius = 0.0125, color = rgb(190, 0, 0, max = 255), aspect = T, add = T)
    
    rgl.triangles(x = nodes[triangles ,1], y = nodes[triangles ,2],
                  z = nodes[triangles,3],
                  color = col,...)
    # rgl.lines(x = nodes[edges ,1], y = nodes[edges ,2],
    #           z = nodes[edges,3],
    #           color = "black",...)
    aspect3d("iso")
    
    if (nsurf > 1 && isurf < nsurf)
    {readline("Press a button for the next plot...")}
  }
}

## MAPVIEW INTERACTIVE PLOT ----------------------------------------------------
plot_mapview <- function(data_sf) {
  mapview(data_sf)
}

## 2D ggplot BETWEEN TWO DATES -------------------------------------------------
plot_ggplot <- function(initial_time, final_time, data_sf, t) {
  world <- st_read("data/geometry/ne_110m_admin_0_map_units/ne_110m_admin_0_map_units.shp") %>% 
           st_transform(crs = "+proj=longlat +datum=WGS84")
  
  starting_time <- paste0(sub("T.*", "", data_sf[1,1]$time),"T00:00:00Z")
  
  time_window <- c(initial_time, final_time)
  time_window <- (as.numeric(ymd_hms(time_window)) - as.numeric(ymd_hms(starting_time)))/(60*60*24)
  time_window <- time_window/365
  ggplot() +
    geom_sf(data = world, colour = "#ffffff20", fill = "#2d2d2d60", size = .5) +
    geom_sf(data = data_sf[t>=time_window[1] & t<=time_window[2],], size = .5, colour = "black") +
    coord_sf(crs = st_crs(world), datum = NA)
    #+ theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5)) +
    #labs(title = "EARTHQUAKES WITH MAGNITUDES > 5", x = NULL, y = NULL,
    #     caption = paste0("Data from: United States Geological Survey\nBetween ", sub("T.*", "", initial_time), " and ", sub("T.*", "", final_time)))
}


## 2D ggplot (DISCRETE TIMES) --------------------------------------------------

plot_ggplot_dt <- function(data_sf, t, idx) {
  world <- st_read("data/geometry/ne_110m_admin_0_map_units/ne_110m_admin_0_map_units.shp") %>% 
    st_transform(crs = "+proj=longlat +datum=WGS84")
  
  ggplot() +
    geom_sf(data = world, colour = "#ffffff20", fill = "#2d2d2d60", size = .5) +
    geom_sf(data = data_sf[t == idx,], size = 1, colour = "black") +
    coord_sf(crs = st_crs(world), datum = NA) +
    theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5)) +
    labs(title = "EARTHQUAKES WITH\nMAGNITUDES > 4.5", x = NULL, y = NULL,
         caption = paste0("Data from: United States Geological Survey\ntime_index = ", idx))
}

## CONVERT FROM REAL NUMBER TO DATE --------------------------------------------
convert_date <- function(starting_time, shift) {
  shift <- shift*60*60*24
  d <- as_datetime(shift, origin = ymd_hms(starting_time))
  
  return(d)
}

## GET COASTLINE BOUNDARIES  --------------------------------------------------
shapefile_lines <- function() {
  # save <- par3d(skipRedraw=TRUE)
  # on.exit(par3d(save))
  r <- 1 # earth radius in meters
  shp <- readShapeSpatial(fn = "data/geometry/ne_110m_coastline/ne_110m_coastline.shp")
  # Combine lines in a single matrix
  mat <- do.call(rbind, sapply(1:length(shp@lines), function(i) rbind(shp@lines[[i]]@Lines[[1]]@coords,c(NA,NA))))
  # Convert spherical to cartesian
  xyz <- rgl.sph2car(mat[,2], mat[,1], radius = r)
  # bg3d("black")
  # Draw sphere
  # rgl.sphere(ng = 100, radius = r, col = "gray50", specular = "black", add = T)
  # plot3d(xyz, col = "white", add = T, type = "l")
  
  return (xyz)
}

rgl.sphere <- function (x = 0, y = NULL, z = NULL, ng = 50, radius = 1, color = "white", add = F, set_normals = T, ...) {
  
  if(length(ng) == 1) ng <- c(ng, ng)
  nlon <- ng[1]; nlat <- ng[2]
  
  lat <- matrix(seq(90, -90, len = nlat)*pi/180, nlon, nlat, byrow = TRUE)
  long <- matrix(seq(-180, 180, len = nlon)*pi/180, nlon, nlat)
  
  vertex  <- rgl:::rgl.vertex(x, y, z)
  nvertex <- rgl:::rgl.nvertex(vertex)
  radius  <- rbind(vertex, rgl:::rgl.attr(radius, nvertex))[4,]
  color  <- rbind(vertex, rgl:::rgl.attr(color, nvertex))[4,]
  
  for(i in 1:nvertex) {
    add2 <- if(!add) i>1 else T
    x <- vertex[1,i] + radius[i]*cos(lat)*cos(long)
    y <- vertex[2,i] + radius[i]*cos(lat)*sin(long)
    z <- vertex[3,i] + radius[i]*sin(lat)
    if(set_normals)
      persp3d(x, y, z, add = add2, color = color[i], normal_x = x, normal_y = y, normal_z = z, ...)
    else
      persp3d(x, y, z, add = add2, color = color[i], ...)
  }
}

rgl.sph2car <- function(lat = 0, lon = 0, radius = 1, deg = T, precise = T) {
  if(deg) {
    if(precise){
      lat <- lat/180
      lon <- lon/180
      x <- radius*cospi(lat)*cospi(lon)
      y <- radius*cospi(lat)*sinpi(lon)
      z <- radius*sinpi(lat)
    }else{
      lat <- lat*pi/180
      lon <- lon*pi/180
      x <- radius*cos(lat)*cos(lon)
      y <- radius*cos(lat)*sin(lon)
      z <- radius*sin(lat)
    }
  }
  return(matrix(c(x,y,z), nrow = length(x), ncol = 3, dimnames = list(NULL, c("x","y","z"))))
}

## FIND THE ZONE FOR UTM COORDINATES -------------------------------------------
find_one_utm_zone <- function(longitude, latitude) {
  
  # Special zones for Svalbard
  if (latitude >= 72.0 && latitude <= 84.0 ) {
    if (longitude >= 0.0  && longitude <  9.0)
      return("31X");
    if (longitude >= 9.0  && longitude < 21.0)
      return("33X")
    if (longitude >= 21.0 && longitude < 33.0)
      return("35X")
    if (longitude >= 33.0 && longitude < 42.0)
      return("37X")
  }
  # Special zones for Norway
  if (latitude >= 56.0 && latitude < 64.0 ) {
    if (longitude >= 0.0  && longitude <  3.0)
      return("31V");
    if (longitude >= 3.0  && longitude < 12.0)
      return("32V")
  }
  
  # North + South Poles
  
  if (latitude > 84.0){
    if ((longitude+180)%%360-180 < 0) {return("Y")}
    if ((longitude+180)%%360-180 > 0) {return("Z")}
  } else if (latitude < -80.0){
    if ((longitude+180)%%360-180 < 0) {return("A")}
    if ((longitude+180)%%360-180 > 0) {return("B")}
  }
  
  # Everything in the middle
  
  if ( (latitude > -80.0) && (latitude <= 84.0) ){
    
    mid_zones <- LETTERS[c(3:8,10:14,16:24)] # C to X, skip I and O
    utm_letter <- mid_zones[ min(floor( (latitude + 80) / 8 )+1 , 20) ]
    utm_number <- (floor( (longitude + 180) / 6 ) %% 60) + 1 # modulo in case longitude is 0 to 360 instead of -180 to 180
    utm_zone <- paste0(utm_number, utm_letter)
    return(utm_zone)
    
  } else {
    stop("lat long not valid (or something else broke)")
  }
}

find_utm_zone <- function(lon, lat){
  purrr::map2_chr(.x = lon, .y = lat, .f = find_one_utm_zone)
}

## (long, lat) TO UTM COORDINATES ----------------------------------------------
LongLat2UTM <- function(x,y,zone) {
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  sp::proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")
  res <- NULL
  for(i in 1:length(x)) {
    res <- rbind(res, as.data.frame(sp::spTransform(xy[1,], CRS(paste("+proj=utm +zone=", zone[i], " ellps=WGS84", sep = "")))))
  }
  return(as.data.frame(res))
}
