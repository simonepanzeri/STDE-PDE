#'##############################################################################
#'##################### SPATIO-TEMPORAL DENSITY ESTIMATION ##################### 
#'###################### HELPER FUNCTIONS FOR SIMULATION 2 #####################
#'##############################################################################

## GRAPHICAL SETTINGS ----------------------------------------------------------
zoom = 0.6
userMatrix = rbind(c(  0.96563137,  0.1774523, -0.1899119,  0),
                   c( -0.03294301,  0.8083354,  0.5877997,  0),
                   c(  0.25781897, -0.5613416,  0.7863999,  0),
                   c(  0.00000000,  0.0000000,  0.0000000,  1))
windowRect = c(150,  150, 420,  420)

#pp <- par3d(no.readonly = TRUE)

## SAMPLE N POINTS UNIFORMLY ON THE UNITARY SPHERE -----------------------------
uniform.sphere <- function(N){
  theta = 2*pi*runif(N, 0, 1)
  phi = acos(1-2*runif(N, 0, 1))
  
  x = sin(phi) * cos(theta)
  y = sin(phi) * sin(theta)
  z = cos(phi)
  
  return(cbind(x,y,z))
}

## SAMPLE N POINTS FROM A DENSITY FUNCTION -------------------------------------
data.true.sphere <- function(N, density.function = dens.func){
  data_true <- c()
  # repeat until the desired sample dimension n
  while(length(data_true[,1]) < N){
    # Sample Uniformly
    data_unif <- uniform.sphere(N)
    t <- runif(N,0,1)
    
    # Evaluate the Density Function on the Generated Points
    f <- dens.func(data_unif, t)
    f[is.na(f)] = max(f, na.rm = TRUE)
    
    # For Each Generated Point, Decide to Keep it or Discard it on the Base of the 
    # Criterium: runif(1,0,1) < density.value
    for(i in 1:length(data_unif[,1])){
      if(runif(1,0,1) < f[i]){
        data_true <- rbind(data_true, c(data_unif[i,], t[i]))
      }
    }
  }
  data <- data_true[1:N, ]
  if(N == 1) {
    data <- t(data)
  }
  data <- as.data.frame(data)
  colnames(data) <- c('x','y','z','t')
  
  return (data)
}

## TRUE DENSITY FUNCTION -------------------------------------------------------
mixture.true <- function(x, mu1, k1, beta1, gamma11, gamma12, 
                         mu2, k2, beta2, gamma21, gamma22, 
                         mu3, k3, beta3, gamma31, gamma32, 
                         mu4, k4, beta4, gamma41, gamma42, 
                         mu5, k5, beta5, gamma51, gamma52){  
  G1 <- cbind(mu1, gamma11, gamma12)
  G2 <- cbind(mu2, gamma21, gamma22)
  G3 <- cbind(mu3, gamma31, gamma32)
  G4 <- cbind(mu4, gamma41, gamma42)
  G5 <- cbind(mu5, gamma51, gamma52)
  
  return (dkent(x, G1, param = c(k1, beta1)) / 5 +
          dkent(x, G2, param = c(k2, beta2)) / 5 + 
          dkent(x, G3, param = c(k3, beta3)) / 5 + 
          dkent(x, G4, param = c(k4, beta4)) / 5 +
          dkent(x, G5, param = c(k5, beta5)) / 5 )
}

## MIXTURE PARAMETERS (5 KENT DISTRIBUTIONS) -----------------------------------
get.mu <- function(time){
  mu1 <- c(-0.5+time, -0.5+time, 0.8+time) 
  mu1 <- mu1 / sqrt( sum(mu1^2) )
  
  mu2 <- c(-0.3, -0.3, 0.2)
  mu2 <- mu2 / sqrt( sum(mu2^2) )
  
  mu3 <- c(0.5, -0.5, 0.8)
  mu3 <- mu3 / sqrt( sum(mu3^2) )
  
  mu4 <- c(0.2, -1, 0)
  mu4 <- mu4 / sqrt( sum(mu4^2) )
  
  mu5 <- c(0.6+time, -0.5+time, 0.3+time)
  mu5 <- mu5 / sqrt( sum(mu5^2) )
  
  return (list(mu1,mu2,mu3,mu4,mu5))
}

{
  
  gamma11 <- c(-0.7789378, 0.6157424, 0.1188163)
  gamma12 <- c(-0.5695773, -0.6154000, -0.5448528)
  k1 = 18
  beta1 = 0
  
  gamma21 <- c(-0.8651146, 0.3803316, -0.3269933)
  gamma22 <- c(0.1482597, -0.4288975, -0.8911038)
  k2 = 15
  beta2 = 7
  
  gamma31 <- c(-0.66647307, -0.74323532,-0.05843723)
  gamma32 <- c(0.5753645, -0.4629244, -0.6742824)
  k3 = 20
  beta3 = 10
  
  gamma41 <- c( 0.6364658, -0.0303920, -0.7707059)
  gamma42 <-  c(-0.7545879, -0.2314437, -0.6140285)
  k4 = 20
  beta4 = 7
  
  gamma51 <- c( 0.6364658, -0.0303920, -0.7707059)
  gamma52 <- c( 0.7545879, -0.2314437, -0.6140285)
  k5 = 20
  beta5 = 4
  
}

## DENSITY FUNCTION ------------------------------------------------------------
dens.func <- function(x, time){
  if(length(time) == 1) {
    x <- t(x)
  }
  result <- NULL
  for (i in 1:length(time)) {
    param <- get.mu(time[i])
    mu1 <- param[[1]]
    mu2 <- param[[2]]
    mu3 <- param[[3]]
    mu4 <- param[[4]]
    mu5 <- param[[5]]
    
    current <- mixture.true(x[i,], mu1, k1, beta1, gamma11, gamma12, 
                            mu2, k2, beta2, gamma21, gamma22, 
                            mu3, k3, beta3, gamma31, gamma32, 
                            mu4, k4, beta4, gamma41, gamma42, 
                            mu5, k5, beta5, gamma51, gamma52)
    
    result <- c(result, current)
  }
  
  return (result)
}

## DATA GENERATION -------------------------------------------------------------
generate.data <- function(N, proc){
  dir.create(file.path(getwd(), "data"), showWarnings = FALSE)
  
  set.seed(100+proc)
  
  # Data
  data <- data.true.sphere(N)
  write.table(data, paste0("data/",N,"data_",proc,".txt"), row.names = F, col.names = F)
}

## HELLINGER DISTANCE ----------------------------------------------------------
hellinger.distance <- function(p, q){
  p <- as.numeric(p)
  q <- as.numeric(q)
  
  if(length(p) != length(q)){
    print("Different lengths")
  }
  
  return ( sqrt(0.5*mean((sqrt(p)-sqrt(q))^2, na.rm = TRUE)) )
}

## PLOT A mesh OBJECT ----------------------------------------------------------
plot.mesh <- function(mesh, M = NULL, m = NULL, ...){
  
  FEM = FEM(rep(0, nrow(mesh$nodes)), create.FEM.basis(mesh))
  
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

## PLOT A SAMPLE OF POINTS -----------------------------------------------------
plot.sample <- function(FEM, sample, M = NULL, m = NULL, ...){
  
  if (is.null(m)) { m = min(FEM$coeff)}
  if (is.null(M)) { M = max(FEM$coeff)}
  triangles = c(t(FEM$FEMbasis$mesh$triangles))
  ntriangles = nrow(FEM$FEMbasis$mesh$triangles)
  order = FEM$FEMbasis$mesh$order
  nodes = FEM$FEMbasis$mesh$nodes
  edges = matrix(rep(0,6*ntriangles), ncol = 2)
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
    light3d(specular = "black") 
    
    diffrange = M - m
    
    col = coeff[triangles,isurf]
    col = (col - min(coeff[,isurf]))/diffrange*(ncolor-1)+1
    
    rgl.spheres(x = sample[,1], y = sample[,2], z = sample[,3], radius = 0.015,
                color = rgb(190, 0, 0, max = 255))
    
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

## FEM PLOT (Colormap with Range [m,M]) ----------------------------------------
plot.FEM <- function(FEM, M = NULL, m = NULL, ...){
  
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
