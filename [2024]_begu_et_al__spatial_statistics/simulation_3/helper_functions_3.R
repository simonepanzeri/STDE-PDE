#'##############################################################################
#'##################### SPATIO-TEMPORAL DENSITY ESTIMATION ##################### 
#'###################### HELPER FUNCTIONS FOR SIMULATION 3 #####################
#'##############################################################################

## GRAPHICAL SETTINGS ----------------------------------------------------------
zoom = 0.6
userMatrix = rbind( c( -0.85268313,  0.4860866,  0.1914454,    0),
                    c( -0.04926674, -0.4396369,  0.8968235,    0),
                    c( 0.52010053,   0.7552745,  0.3988189,    0),
                    c( 0.00000000,   0.0000000,  0.0000000,    1))
windowRect = c(150, 150, 420, 420)

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
mixture.true <- function(x, mu1, k1, beta1, gamma11, gamma12, mu2, k2, beta2, gamma21,gamma22, 
                         mu3, k3, beta3, gamma31, gamma32, mu4, k4, beta4, gamma41, gamma42,  
                         mu5, k5, beta5, gamma51, gamma52, mu6, k6, beta6, gamma61, gamma62,
                         mu7, k7, beta7, gamma71, gamma72, mu8, k8, beta8, gamma81, gamma82,
                         mu9, k9, beta9, gamma91, gamma92, mu10, k10, beta10, gamma101, gamma102,
                         mu11, k11, beta11, gamma111, gamma112, mu12, k12, beta12, gamma121, gamma122,
                         mu13, k13, beta13, gamma131, gamma132, mu14, k14, beta14, gamma141, gamma142,
                         mu15, k15, beta15, gamma151, gamma152, mu16, k16, beta16, gamma161, gamma162,
                         mu17, k17, beta17, gamma171, gamma172, mu18, k18, beta18, gamma181, gamma182,
                         mu19, k19, beta19, gamma191, gamma192, mu20, k20, beta20, gamma201, gamma202,
                         mu21, k21, beta21, gamma211, gamma212, mu22, k22, beta22, gamma221, gamma222,
                         mu23, k23, beta23, gamma231, gamma232,
                         mu24, k24, beta24, gamma241, gamma242,
                         mu25, k25, beta25, gamma251, gamma252){  
  G1 <- cbind(mu1, gamma11, gamma12)
  G2 <- cbind(mu2, gamma21, gamma22)
  G3 <- cbind(mu3, gamma32, gamma31)
  G4 <- cbind(mu4, gamma41, gamma42)
  G5 <- cbind(mu5, gamma52, gamma51)
  G6 <- cbind(mu6, gamma62, gamma61)
  G7 <- cbind(mu7, gamma71, gamma72)
  G8 <- cbind(mu8, gamma82, gamma81)
  G9 <- cbind(mu9, gamma92, gamma91)
  G10 <- cbind(mu10, gamma102, gamma101)
  G11 <- cbind(mu11, gamma112, gamma111)
  G12 <- cbind(mu12, gamma122, gamma121)
  G13 <- cbind(mu13, gamma131, gamma132)
  G14 <- cbind(mu14, gamma141, gamma142)
  G15 <- cbind(mu15, gamma151, gamma152)
  G16 <- cbind(mu16, gamma161, gamma162)
  G17 <- cbind(mu17, gamma171, gamma172)
  G18 <- cbind(mu18, gamma181, gamma182)
  G19 <- cbind(mu19, gamma192, gamma191)
  G20 <- cbind(mu20, gamma202, gamma201)
  G21 <- cbind(mu21, gamma212, gamma211)
  G22 <- cbind(mu22, gamma222, gamma221)
  G23 <- cbind(mu23, gamma232, gamma231)
  G24 <- cbind(mu24, gamma242, gamma241)
  G25 <- cbind(mu25, gamma251, gamma252)
  
  return (dkent(x, G1, param=c(k1, beta1)) / 25 +
            dkent(x, G2, param=c(k2, beta2)) / 25 +
            dkent(x, G3, param=c(k3, beta3)) / 25 +
            dkent(x, G4, param=c(k4, beta4)) / 25 +
            dkent(x, G5, param=c(k5, beta5)) / 25 +
            dkent(x, G6, param=c(k6, beta6)) / 25 +
            dkent(x, G7, param=c(k7, beta7)) / 25 +
            dkent(x, G8, param=c(k8, beta8)) / 25 +
            dkent(x, G9, param=c(k9, beta9)) / 25 +
            dkent(x, G10, param=c(k10, beta10)) / 25 +
            dkent(x, G11, param=c(k11, beta11)) / 25 + 
            dkent(x, G12, param=c(k12, beta12)) / 25 +
            dkent(x, G13, param=c(k13, beta13)) / 25 +
            dkent(x, G14, param=c(k14, beta14)) / 25 +
            dkent(x, G15, param=c(k15, beta15)) / 25 +
            dkent(x, G16, param=c(k16, beta16)) / 25 +
            dkent(x, G17, param=c(k17, beta17)) / 25 +
            dkent(x, G18, param=c(k18, beta18)) / 25 +
            dkent(x, G19, param=c(k19, beta19)) / 25 +
            dkent(x, G20, param=c(k20, beta20)) / 25 +
            dkent(x, G21, param=c(k21, beta21)) / 25 +
            dkent(x, G22, param=c(k22, beta22)) / 25 +
            dkent(x, G23, param=c(k23, beta23)) / 25 +
            dkent(x, G24, param=c(k24, beta24)) / 25 +
            dkent(x, G25, param=c(k25, beta25)) / 25
  )
}

## MIXTURE PARAMETERS (25 KENT DISTRIBUTIONS) -----------------------------------
get.kb <- function(time){
  k <- 400
  b <- 125+100*abs(time-0.5)
  
  return (list(k,b))
}

{
  
  # Line 1
  mu1 <- c(1,0,0)
  mu2 <- c(0.8,-0.1,-0.1)
  mu2 <- mu2 / sqrt( sum(mu2^2) )
  mu3 <- c(0.6,-0.2,-0.2)
  mu3 <- mu3 / sqrt( sum(mu3^2) )
  mu4 <- c(0.4,-0.3,-0.3)
  mu4 <- mu4 / sqrt( sum(mu4^2) )
  
  gamma11 <- c(-0.001656101, -0.968065656,   0.250691329)
  gamma12 <- c(-0.0007877991, 0.2506928582, -0.9680663563)
  gamma21 <- c(-0.1628647,   -0.4337402,    -0.8861967)
  gamma22 <- c( 0.05557397,   0.89273019,   -0.44715136)
  gamma31 <- c(-0.34029681,  -0.93697252,   -0.07925005)
  gamma32 <- c(-0.2588086,    0.1743533,    -0.9500626)
  gamma41 <- c(-0.7321271,   -0.5268741,    -0.4317334)
  gamma42 <- c(-0.04951012,   0.67329713,   -0.73771249)
  
  # Line 2
  mu5 <- c(0,0,1)
  mu6 <- c(0.1,0.1,0.8)
  mu6 <- mu6 / sqrt( sum(mu6^2) )
  mu7 <- c(0.2,0.2,0.5)
  mu7 <- mu7 / sqrt( sum(mu7^2) )
  mu8 <- c(0.3,0.3,0.3)
  mu8 <- mu8 / sqrt( sum(mu8^2) )
  mu9 <- c(0.4,0.4,0.1)
  mu9 <- mu9 / sqrt( sum(mu9^2) )
  mu10 <- c(0.5,0.5,-0.1)
  mu10 <- mu10 / sqrt( sum(mu10^2) )
  mu11 <- c(0.15,0.15,0.1)
  mu11 <- mu11 / sqrt( sum(mu11^2) )
  
  gamma51 <- c( -7.510675e-01, 6.602255e-01, -3.442408e-05)
  gamma52 <- c( -0.660214721, -0.751054928,   0.005746117)
  gamma61 <- c( -0.9277247,    0.3666829,     0.0697894)
  gamma62 <- c( -0.3535596,   -0.9231982,     0.1506671)
  gamma71 <- c( -0.7111897,   -0.5093999,     0.4844801)
  gamma72 <- c(  0.6120213,   -0.7877217,     0.0701740)
  gamma81 <- c( -0.7466229,    0.6595548,     0.0868430)
  gamma82 <- c( -0.3262194,   -0.4767600,     0.8162603)
  gamma91 <- c( -0.6765086,    0.7167571,    -0.1691015)
  gamma92 <- c( -0.2423570628, 0.0001443617,  0.9701871125)
  gamma101 <- c(-0.6293300,    0.6982024,     0.3412582)
  gamma102 <- c( 0.3319063,   -0.1555797,     0.9303941)
  gamma111 <- c(-0.7680521,    0.5322785,     0.3560556)
  gamma112 <- c( 0.002049244, -0.553954592,   0.832544360)
  
  # Line 3
  mu12 <- c(0.6,-0.2,-0.1)
  mu12 <- mu12 / sqrt( sum(mu12^2) )
  
  gamma121 <-c( -0.3412487, -0.9139369,  -0.2197012)
  gamma122 <- c(-0.08342826, 0.26225792, -0.96138468)
  
  # Line 4
  mu13 <- c(-1,0,0)
  mu13 <- mu13 / sqrt( sum(mu13^2) )
  mu14 <- c(-0.8,0.2,0)
  mu14 <- mu14 / sqrt( sum(mu14^2) )
  mu15 <- c(-0.6,0.5,0)
  mu15 <- mu15 / sqrt( sum(mu15^2) )
  
  gamma131 <- c(-0.01162783,  -0.99930072,  -0.03553692)
  gamma132 <- c( 0.001223968, -0.035553520,  0.999367024)
  gamma141 <- c(-0.2338769,   -0.9648624,   -0.1197584)
  gamma142 <- c(-0.03279129,  -0.11527664,   0.99279204)
  gamma151 <- c(-0.6388464,   -0.7546472,   -0.1496089)
  gamma152 <- c(-0.0936074,   -0.1167740,    0.9887373)
  
  # Line 5
  mu16 <- c(-0.6,0,0.1)
  mu16 <- mu16 / sqrt( sum(mu16^2) )
  mu17 <- c(-0.6,-0.1,0.3)
  mu17 <- mu17 / sqrt( sum(mu17^2) )
  mu18 <- c(-0.6,0.2,-0.2)
  mu18 <- mu18 / sqrt( sum(mu18^2) )
  mu19 <- c(-0.9,0.1,-0.6)
  mu19 <- mu19 / sqrt( sum(mu19^2) )
  
  gamma161 <- c(-0.1539147,    0.5103210,   -0.8460985)
  gamma162 <- c(-0.08387069,  -0.85995926,  -0.50342405)
  gamma171 <- c(-0.4683750,    0.2761541,   -0.8392638)
  gamma172 <- c( 0.002319274, -0.949511193, -0.313724586)
  gamma181 <- c(-0.4246560,   -0.5672990,    0.7055771)
  gamma182 <- c( 0.03828766,   0.76738745,   0.64003949)
  gamma191 <- c(-0.4100383,   -0.7714573 ,   0.4865410)
  gamma192 <- c(-0.3877780,    0.6302894,    0.6725797)
  
  # Line 6
  mu20 <- c(0,-1,0.8)
  mu20 <- mu20 / sqrt( sum(mu20^2) )
  mu21 <- c(0,-1,0.5)
  mu21 <- mu21 / sqrt( sum(mu21^2) )
  mu22 <- c(0,-0.8,0.2)
  mu22 <- mu22 / sqrt( sum(mu22^2))
  mu23 <- c(0,-1,-0.1)
  mu23 <- mu23 / sqrt( sum(mu23^2) )
  mu24 = c(0.2,0.4,0.1)
  mu24 <- mu24 / sqrt( sum(mu24^2) )
  mu25 = c(0.2,0.7,0.1)
  mu25 <- mu25 / sqrt( sum(mu25^2) )
  
  gamma201 <- c( -0.8839177,   0.2926913,    0.3647207)
  gamma202 <- c( -0.4676421,  -0.5541359,   -0.6886540)
  gamma211 <- c( -0.9267158,  -0.1686490,   -0.3357905)
  gamma212 <- c(  0.3757628,  -0.4164846,   -0.8278544)
  gamma221 <- c( -0.8907620,   0.1103146,    0.4408784)
  gamma222 <- c( -0.4544509,  -0.2251213,   -0.8618554)
  gamma231 <- c( -0.84513209,  0.05302266,  -0.53192137)
  gamma232 <- c(  0.53454499,  0.09064177,  -0.84026527)
  gamma241 <- c( -0.7745285,   0.2406619,    0.5849679)
  gamma242 <- c(  0.4627985,  -0.4148049,    0.7834248)
  gamma251 <- c( -0.95816540,  0.28523697,  -0.02364172)
  gamma252 <- c( -0.06147375, -0.12441959,   0.99032356)
}

## DENSITY FUNCTION ------------------------------------------------------------
dens.func <- function(x, time)
{
  param <- get.kb(time)
  k <- param[[1]]
  b <- param[[2]]
  
  if(length(time) == 1) {
    x <- t(x)
  }
  result <- NULL
  for (i in 1:length(time)) {
    current <- mixture.true(x[i,], mu1, k, b, gamma11, gamma12, mu2, k, b, gamma21, gamma22, 
                            mu3, k, b, gamma31, gamma32, mu4, k, b, gamma41, gamma42, 
                            mu5, k, b, gamma51, gamma52, mu6, k, b, gamma61, gamma62,
                            mu7, k, b, gamma71, gamma72, mu8, k, b, gamma81, gamma82,
                            mu9, k, b, gamma91, gamma92, mu10, k, b, gamma101, gamma102,
                            mu11, k, b, gamma111, gamma112, mu12, k, b, gamma121, gamma122,  
                            mu13, k, b, gamma131, gamma132, mu14, k, b, gamma141, gamma142,
                            mu15, k, b, gamma151, gamma152,mu16, k, b, gamma161, gamma162,
                            mu17, k, b, gamma171, gamma172, mu18, k, b, gamma181, gamma182,
                            mu19, k, b, gamma191, gamma192, mu20, k, b, gamma201, gamma202,
                            mu21, k, b, gamma211, gamma212, mu22, k, b, gamma221, gamma222,
                            mu23, k, b, gamma231, gamma232,
                            mu24, k, b, gamma241, gamma242,
                            mu25, k, b, gamma251, gamma252)
    
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
