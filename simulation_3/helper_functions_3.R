#'############################################################################'#
#'##################### SPATIO-TEMPORAL DENSITY ESTIMATION ###################'# 
#'###################### HELPER FUNCTIONS FOR SIMULATION 3 ###################'#
#'############################################################################'#

## SAMPLE n POINTS UNIFORMLY ON THE MANIFOLD -----------------------------------
uniform_surface <- function(n){
  load("mesh/simulation3.fullPoints.proj.RData")
  numFullPoints = nrow(fullPoints.proj)
  unifPointsIndex = sample(1:numFullPoints, n)
  unifPoints = fullPoints.proj[unifPointsIndex,]
  
  return (unifPoints)
}

## TRUE DENSITY FUNCTION -------------------------------------------------------
dens.func <- function(x){
  return (((sin(x[,1])+cos(x[,2])+x[,3]+5.159043)^3)/22188.75)
}

## SAMPLE n POINTS FROM THE DENSITY densityFEM ---------------------------------
vertices <- read.table("mesh/simulation3.surface.vertices.txt")
triangles <- read.table("mesh/simulation3.surface.triangles.txt")
mesh <- create.mesh.2.5D(nodes = vertices, triangles = triangles[,1:3])
FEMbasis <- create.FEM.basis(mesh)

true = dens.func(mesh$nodes)
densityFEM <- FEM(coeff = true, FEMbasis)

data_true_surface <- function(n, density.function = densityFEM){
  data_true <- c()
  # repeat until the desired sample dimension n
  while(length(data_true[,1]) < n){
    # sample uniformly
    data_unif <- uniform_surface(n)
    # evaluate the density function on the generated points
    f <- eval.FEM(density.function, locations = data_unif)
    f[is.na(f)] = max(f, na.rm = TRUE)
    # for each generated point, decide to keep it or discard it on the base of the 
    # criterium: runif(1,0,1) < 5*density.value
    for(i in 1:length(data_unif[,1])){
      if(runif(1,0,1) < 5*f[i]){
        data_true <- rbind(data_true, data_unif[i,])
      }
    }
  }
  data = data_true[1:n, ]
}

## DATA GENERATION -------------------------------------------------------------
generate.data <- function(N, proc){
  dir.create(file.path(getwd(), "data"), showWarnings = FALSE)
  
  set.seed(100+proc)
  
  locations <- data_true_surface(N)
  times <- rbeta(N,2,2)
  
  data <- cbind(locations, times)
  
  write.table(data, paste0("data/",N,"data_",proc,".txt"), row.names = FALSE, col.names = FALSE)
}

## HELLINGER DISTANCE ----------------------------------------------------------
hellinger.distance <- function(p, q){
  p <- as.numeric(p)
  q <- as.numeric(q)
  
  if(length(p) != length(q)){
    print("Different lengths")
  }
  
  return ( sqrt(0.5*mean((sqrt(p)-sqrt(q))^2)) )
}

## PLOT A SAMPLE ---------------------------------------------------------------
plot.sample <- function(FEM, sample, M = NULL, m = NULL, ...){
  
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
  
  p = jet.col(n = 1000, alpha = 0.8)
  # alternative color palette: p <- colorRampPalette(c("#0E1E44", "#3E6DD8", "#68D061", "#ECAF53", "#EB5F5F", "#E11F1C"))(1000)
  
  palette(p)
  
  ncolor = length(p)
  
  nsurf = dim(coeff)[[2]]
  for (isurf in 1:nsurf)
  {
    open3d(zoom = zoom, userMatrix = userMatrix, windowRect = windowRect)
    rgl.pop("lights") 
    light3d(specular="black") 
    
    diffrange = M - m
    
    col = coeff[triangles,isurf]
    col = (col - min(coeff[,isurf]))/diffrange*(ncolor-1)+1
    
    rgl.spheres(x = sample[,1], y = sample[,2], z = sample[,3], radius = 0.05, colo = rgb(190, 0, 0, max = 255))
    
    rgl.triangles(x = nodes[triangles ,1], y = nodes[triangles ,2],
                  z = nodes[triangles,3],
                  color = col,...)
    # rgl.lines(x = nodes[edges ,1], y = nodes[edges ,2],
    #           z = nodes[edges,3],
    #           color = "black",...)
    aspect3d("iso")
    
    if (nsurf > 1 && isurf<nsurf)
    {readline("Press a button for the next plot...")}
  }
}

## PLOT A mesh OBJECT ----------------------------------------------------------
plot.mesh <- function(mesh, M = NULL, m = NULL,...){
  
  FEM <- FEM(rep(0,nrow(mesh$nodes)), create.FEM.basis(mesh))
  
  if (is.null(m)) { m = min(FEM$coeff)}
  if (is.null(M)) { M = max(FEM$coeff)}
  triangles = c(t(FEM$FEMbasis$mesh$triangles))
  ntriangles = nrow(FEM$FEMbasis$mesh$triangles)
  order = FEM$FEMbasis$mesh$order
  nodes = FEM$FEMbasis$mesh$nodes
  edges = matrix(rep(0, 6*ntriangles), ncol = 2)
  for(i in 0:(ntriangles-1)){
    edges[3*i+1,]= c(triangles[3*order*i+1], triangles[3*order*i+2])
    edges[3*i+2,]= c(triangles[3*order*i+1], triangles[3*order*i+3])
    edges[3*i+3,]= c(triangles[3*order*i+2], triangles[3*order*i+3])
  }
  edges = edges[!duplicated(edges),]
  edges <- as.vector(t(edges))
  
  coeff = FEM$coeff
  
  FEMbasis = FEM$FEMbasis
  
  mesh = FEMbasis$mesh
  
  p = jet.col(n = 1000, alpha = 0.8)
  # alternative color palette: p <- colorRampPalette(c("#0E1E44", "#3E6DD8", "#68D061", "#ECAF53", "#EB5F5F", "#E11F1C"))(1000)
  
  palette(p)
  
  ncolor = length(p)
  
  nsurf = dim(coeff)[[2]]
  for (isurf in 1:nsurf)
  {
    open3d(zoom = zoom, userMatrix = userMatrix, windowRect = windowRect)
    rgl.pop("lights") 
    light3d(specular = "black") 
    
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
    
    if (nsurf > 1 && isurf<nsurf)
    {readline("Press a button for the next plot...")}
  }
}

## PLOT A FEM OBJECT WITH JET COLORMAP WITH RANGE [m,M] ------------------------
plot.FEM <- function(FEM, M = NULL, m = NULL, ...){
  
  if (is.null(m)) { m = min(FEM$coeff)}
  if (is.null(M)) { M = max(FEM$coeff)}
  triangles = c(t(FEM$FEMbasis$mesh$triangles))
  ntriangles = nrow(FEM$FEMbasis$mesh$triangles)
  order = FEM$FEMbasis$mesh$order
  nodes = FEM$FEMbasis$mesh$nodes
  edges = matrix(rep(0, 6*ntriangles), ncol = 2)
  for(i in 0:(ntriangles-1)){
    edges[3*i+1,]= c(triangles[3*order*i+1], triangles[3*order*i+2])
    edges[3*i+2,]= c(triangles[3*order*i+1], triangles[3*order*i+3])
    edges[3*i+3,]= c(triangles[3*order*i+2], triangles[3*order*i+3])
  }
  edges = edges[!duplicated(edges),]
  edges <- as.vector(t(edges))
  
  coeff = FEM$coeff
  
  FEMbasis = FEM$FEMbasis
  
  mesh = FEMbasis$mesh
  
  p = jet.col(n = 1000, alpha = 0.8)
  # alternative color palette: p <- colorRampPalette(c("#0E1E44", "#3E6DD8", "#68D061", "#ECAF53", "#EB5F5F", "#E11F1C"))(1000)
  
  palette(p)
  
  ncolor = length(p)
  
  nsurf = dim(coeff)[[2]]
  for (isurf in 1:nsurf)
  {
    open3d(zoom = zoom, userMatrix = userMatrix, windowRect = windowRect)
    rgl.pop("lights") 
    light3d(specular = "black") 
    
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
    
    if (nsurf > 1 && isurf<nsurf)
    {readline("Press a button for the next plot...")}
  }
}

#pp <- par3d(no.readonly = TRUE)

zoom = 0.4
userMatrix = rbind(c( 0.96563137,  0.1774523, -0.1899119,   0.8),
                   c( -0.03294301,  0.8083354,  0.5877997,  0.25),
                   c( 0.25781897, -0.5613416,  0.7863999,   0),
                   c(0.00000000,  0.0000000,  0.0000000 ,   1))
windowRect = c(20,30,600,400)

## COMPUTE THE AREA OF A SINGLE MESH TRIANGLE ----------------------------------
AreaOfTriangle <- function(p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z)
{
  ax = p2x - p1x
  ay = p2y - p1y
  az = p2z - p1z
  bx = p3x - p1x
  by = p3y - p1y
  bz = p3z - p1z
  cx = ay*bz - az*by
  cy = az*bx - ax*bz
  cz = ax*by - ay*bx
  
  return (0.5 * sqrt(cx*cx + cy*cy + cz*cz))
}    

## COMPUTE THE AREA OF A MESH ON A SURFACE -------------------------------------
MeshSurface <- function(mesh)
{
  Area <- 0
  
  for (i in 1:nrow(mesh$triangles)) {
    Area <- Area + AreaOfTriangle(mesh$nodes[mesh$triangles[i,1],1], mesh$nodes[mesh$triangles[i,1],2], mesh$nodes[mesh$triangles[i,1],3],
                                  mesh$nodes[mesh$triangles[i,2],1], mesh$nodes[mesh$triangles[i,2],2], mesh$nodes[mesh$triangles[i,2],3],
                                  mesh$nodes[mesh$triangles[i,3],1], mesh$nodes[mesh$triangles[i,3],2], mesh$nodes[mesh$triangles[i,3],3])
  }
  
  return (Area)
}
