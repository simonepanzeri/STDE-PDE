#'############################################################################'#
#'##################### SPATIO-TEMPORAL DENSITY ESTIMATION ###################'# 
#'###################### HELPER FUNCTIONS FOR SIMULATION 5 ###################'#
#'############################################################################'#

# t <- seq(from = 0, to = 1, by = 0.1)
# integral <- compute_integral()
# true_density <- matrix(nrow = nrow(mesh.eval$nodes), ncol = length(t))
# 
# for (time_index in 1:length(t)) {
#   true_density[,time_index] <- dens.func(mesh.eval$nodes, t[time_index])/integral
#   true_density[,time_index] <- true_density[,time_index]/sum(true_density[,time_index], na.rm = TRUE)
# }
# 
# M <- max(true_density)
# 
# for (time_index in 1:length(t)) {
#   name <- paste0("pictures/", N, "_data/true_density/sim4_true_instant_", N, "data_", t[time_index])
#   plot.density(mesh.eval$nodes[,1], mesh.eval$nodes[,2], true_density[,time_index],
#                max_range = M,
#                filename = name, colorscale = "Jet")
# }


## MIXTURE WEIGHTS (NORMALIZED PRIORS) -----------------------------------------
get.priors <- function(time){
  var1 <- dgamma(10*time+2, 7.5, rate = 1)
  var2 <- dgamma(10*time+2, 7.5, rate = 1)
  var3 <- dgamma(10*time+2, 7.5, rate = 1)
  var4 <- dgamma(10*(0.5+sign(0.5-time)*(time-0.5))+2, 7.5, rate = 1)
  priors <- c(0.7*var1, 0.8*var2, 0.4*var3, 0.25*var4)
  priors <- priors/sum(priors)
  return (priors)
}

# p <- matrix(nrow = 4, ncol = length(t))
# for (time_index in 1:length(t)) {
#   p[,time_index] <- get.priors(t[time_index])
# }
# plot(cbind(t,p[1,]), type = 'l', ylim=range(p), col='red', lwd = 2)
# points(cbind(t,p[2,]), type = 'l', col='blue', lwd = 2)
# points(cbind(t,p[3,]), type = 'l', col='green', lwd = 2)
# points(cbind(t,p[4,]), type = 'l', col='orange', lwd = 2)

## DISTRIBUTION PARAMETERS -----------------------------------------------------
range <- c(-6, 6)

get.parameters <- function(d, time) {
  if (d == 1) {
    # First Gaussian distribution
    mu <- c(-3.5+1.5*time, -3.5+1.5*time) # Mean
    sigma <- 3*matrix(c(0.8, 0.4, 0.4, 0.8), 2, 2) # Covariance matrix
  }
  else if (d == 2) {
    # Second Gaussian distribution
    mu <- c(1.7+time, -2-1.4*time) # Mean
    sigma <- 3*matrix(c(1.5-0.5*time, -0.5*time, -0.5*time, 1.5-0.5*time), 2, 2) # Covariance matrix
  }
  else if (d == 3) {
    # Third Gaussian distribution
    mu <- c(-2-time, 1.5+1.7*time) # Mean
    sigma <- matrix(c(1.6, -1.2*time, -1.2*time, 1.6), 2, 2) # Covariance matrix
  }
  else if (d == 4) {
    ## Fourth Gaussian distribution
    mu <- c(2, 2-1.5*time) # Mean
    sigma <- 1.5*matrix(c(1, 0.9-0.3*time, 0.9-0.3*time, 1), 2, 2) # Covariance matrix
  }
  return (list(mu,sigma))
}

## DATA GENERATION -------------------------------------------------------------
generate.data <- function(N, proc){
  dir.create(file.path(getwd(), "data"), showWarnings = FALSE)
  dir.create(file.path(getwd(), "data/[-6,6]x[-6,6]"), showWarnings = FALSE)
  
  set.seed(100+proc)

  # Number of components of the mixture
  D <- 4 
  
  # Data
  times <- c()
  locations <- c()
  distribution <- c()
  n <- 1
  while (n <= N){
    t <- runif(1,0,1)
    p <- get.priors(t)
    d <- sample(c(1,2,3,4), 1, replace = T, prob = p)
    
    parameters <- get.parameters(d, t)
    l <- mvtnorm::rmvnorm(n = 1, mean = parameters[[1]], sigma = parameters[[2]], method = "svd")
    
    if(between(l[1], range[1], range[2]) & between(l[2], range[1], range[2])){
      
      distribution <- c(distribution, d)
      
      locations <- rbind(locations, l)
      times <- c(times, t)
      
      n <- n+1
    }
  }
  data <- cbind(locations, times)

  write.table(data, paste0("data/[-6,6]x[-6,6]/",N,"data_",proc,".txt"), row.names = F, col.names = F)
}

## PLOT SAMPLE -----------------------------------------------------------------
plot.sample <- function(sample){
  plot(sample, col = rgb(190, 0, 0, max = 255), pch = 20, xlim = range, ylim = range, xlab = "", ylab = "", asp = 1, xaxt = "n", yaxt = "n")
  points(mean(sample[,1]), mean(sample[,2]), pch = 20, cex = 1.5, lwd = 2, col = "black", xlim = range, ylim = range)
  #grid()
}

## DENSITY FUNCTION ------------------------------------------------------------
dens.func <- function(locations, time){
  p <- get.priors(time)
  p <- p/sum(p)
  dens <- p[1]*mvtnorm::dmvnorm(locations, mean = get.parameters(1,time)[[1]], sigma = get.parameters(1,time)[[2]]) +
          p[2]*mvtnorm::dmvnorm(locations, mean = get.parameters(2,time)[[1]], sigma = get.parameters(2,time)[[2]]) +
          p[3]*mvtnorm::dmvnorm(locations, mean = get.parameters(3,time)[[1]], sigma = get.parameters(3,time)[[2]]) +
          p[4]*mvtnorm::dmvnorm(locations, mean = get.parameters(4,time)[[1]], sigma = get.parameters(4,time)[[2]])
  
  return (dens)
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

## R-INLA UTILITY --------------------------------------------------------------
book.mesh.dual <- function(mesh) {
  if (mesh$manifold == "R2") {
    ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
      colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
    library(parallel)
    pls <- mclapply(1:mesh$n, function(i) {
      p <- unique(Reduce("rbind", lapply(1:3, function(k) {
        j <- which(mesh$graph$tv[,k] == i)
        if (length(j)>0) 
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

## PLOT MESH USING PLOTLY ------------------------------------------------------
plot.mesh <- function(x, filename, ...){ 
  edges <- x$edges
  p <- plot_ly(width = 1000, height = 1000) %>% layout(scene = list(
    aspectratio = list(
      x = 1,
      y = 1
    )),
    xaxis = list(
      title = "",
      showgrid = F,
      zeroline = F,
      showticklabels = F,
      ticks = ""),
    yaxis = list(
      title = "",
      showgrid = F,
      zeroline = F,
      showticklabels = F,
      ticks = ""),
    margin = list(
      b = 0,
      l = 0,
      r = 14,
      t = 13
    )) %>%
    add_markers(x = x$nodes[,1],
                y = x$nodes[,2], 
                color = I("black"),
                hoverinfo = "text",
                text = paste("</br><b> Coordinates:", round(x$nodes[,1],2),
                             round(x$nodes[,2],2)), 
                showlegend = T, 
                visible = T) %>%
    add_segments(x = x$nodes[edges[,1],1],
                 y = x$nodes[edges[,1],2],
                 xend = x$nodes[edges[,2],1],
                 yend = x$nodes[edges[,2],2], 
                 color = I("black"),
                 showlegend = F) 
  
  plotly::export(p, file = paste0(filename,".png")) 
  
  # Alternative to Export Images
  # saveWidget(p, paste0(filename, ".html"))
  # webshot(paste0(filename,".html"), paste0(filename,".png"), delay = 2)
}

## PLOT SAMPLE USING PLOTLY ----------------------------------------------------
plot.sample <- function(coordinates, filename, ...){
  
  DATA <- data.frame(x = coordinates[,1], y = coordinates[,2])
  
  ay <- list(
    showline = TRUE,
    mirror = "ticks",
    linecolor = toRGB("black"),
    linewidth = 2,
    range = range
  )
  
  ax <- list(
    showline = TRUE,
    mirror = "ticks",
    linecolor = toRGB("black"),
    linewidth = 2,
    range = range
  )
  
  p <- plot_ly(data = DATA, x = ~x, y = ~y, type = "scatter", mode = "markers",
               width = 1000, height = 1000,
               marker = list(size = 10,
                             #color = "black",
                             color = rgb(190, 0, 0, max = 255),
                             line = list(#color = "black",
                                         color = rgb(190, 0, 0, max = 255),
                                         width = 2))) %>%
    layout(scene = list(
      aspectmode = "data",
      aspectratio = list(
        x = 1,
        y = 1
      )),
      xaxis = list(
        title = "",
        showgrid = F,
        zeroline = F,
        showticklabels = F,
        ticks = ""),
      yaxis = list(
        title = "",
        showgrid = F,
        zeroline = F,
        showticklabels = F,
        ticks = ""),
      margin = list(
        b = 0,
        l = 0,
        r = 14,
        t = 13
      ))
  
  p <- p %>% layout(xaxis = ax, yaxis = ay)
  
  plotly::export(p, file = paste0(filename,".png"))
  
  # Alternative to Export Images
  # saveWidget(p, paste0(filename, ".html"))
  # webshot(paste0(filename,".html"), paste0(filename,".png"), delay = 2)
}

## PLOT DENSITY USING PLOTLY ---------------------------------------------------
plot.density <- function(X, Y, Z, max_range = NULL, filename,
                         iso = NULL, showlines = FALSE,...){
  
  if (is.null(iso)) {iso = 1000}
  if (is.null(max_range)) {max_range = max(Z)}
  
  DATA <- data.frame(x = X, y = Y, z = Z)
  
  ay <- list(
    showline = TRUE,
    mirror = "ticks",
    linecolor = toRGB("black"),
    linewidth = 2,
    range = range
  )
  
  ax <- list(
    showline = TRUE,
    mirror = "ticks",
    linecolor = toRGB("black"),
    linewidth = 2,
    range = range
  )
  
  if (showlines){
    # col = colorRamp(c("#00007F", "blue", "#007FFF", "#009EFF", "#00C5FF", "cyan", "#7FFFB9",
    #                   "#7FFF7F", "#B1FF7F", "yellow", "#FF7F00", "#FF4200", "red", "#7F0000"))
    
    p <- plot_ly(DATA, x = ~x, y = ~y, z = ~z, intensity = ~z, color = ~z, type = "contour",
                 width = 1000, height = 1000, showscale = F,
                 contours = list(
                   start = 0,
                   end = max_range,
                   size = (max_range-0)/iso,
                   showlines = showlines
                   #showlabels = T
                 ),
                 line = list(color = "white"),
                 #colors = col,
                 ...) 
  } else {
    p <- plot_ly(DATA, x = ~x, y = ~y, z = ~z, intensity = ~z, color = ~z, type = "contour",
                 width = 1000, height = 1000, showscale = F,
                 contours = list(
                   start = 0,
                   end = max_range,
                   size = (max_range-0)/iso,
                   showlines = FALSE
                   #showlabels = T
                 ), ...) 
  }
  
  p <- p %>%
    layout(scene = list(
      aspectmode = "data",
      aspectratio = list(
        x = 1,
        y = 1
      )),
      xaxis = list(
        title = "",
        showgrid = F,
        zeroline = F,
        showticklabels = F,
        ticks = ""),
      yaxis = list(
        title = "",
        showgrid = F,
        zeroline = F,
        showticklabels = F,
        ticks = ""),
      margin = list(
        b = 0,
        l = 0,
        r = 14,
        t = 13
      )
    )
  
  
  p <- p %>% layout(xaxis = ax, yaxis = ay)
  
  plotly::export(p, file = paste0(filename,".png"))
  
  # Alternative to Export Images
  # saveWidget(p, paste0(filename, ".html"))
  # webshot(paste0(filename,".html"), paste0(filename,".png"), delay = 2)
  
}

## PLOT INTERACTIVE DENSITY IN 3D ----------------------------------------------
plot.interactive.density <-  function(f, ...){
  plot_data <- data.frame(X = f$FEMbasis$mesh$nodes[,1], 
                          Y = f$FEMbasis$mesh$nodes[,2],
                          Z = f$coeff,
                          coeff = f$coeff)
  I = (f$FEMbasis$mesh$triangles[,1]-1); J = (f$FEMbasis$mesh$triangles[,2]-1); K = (f$FEMbasis$mesh$triangles[,3]-1)
  fig <- plot_ly(plot_data, type = "mesh3d", x = ~X, y = ~Y,  z = ~Z, 
                 i = I, j = J, k = K,
                 intensity = ~coeff, color = ~coeff,
                 contours = list(showlabels = TRUE),
                 colorbar = list(title = ""), ...) %>%
    layout(xaxis = list(title = ""),
           yaxis = list(title = ""),
           scene = list(
             camera = list(
               eye = list(x = 1.25, 
                          y = -1.25, 
                          z = 1.25))))
  fig
}

## COMPUTE INTEGRALS OVER THE SPATIO-TEMPORAL DOMAIN OF INTEREST ---------------
integral_time <- function(time){
    
  mixture <- function(x) {
    weights <- rep(0, 4)
    densities <- numeric(4)
    for (i in 1:4) {
      if (i == 1) {
        mu <- c(-3.5+1.5*time, -3.5+1.5*time) # Mean
        sigma <- 3*matrix(c(0.8, 0.4, 0.4, 0.8), 2, 2) # Covariance matrix
        weights[i] <- dgamma(10*time+2, 7.5, rate = 1)
      } else if (i == 2) {
        mu <- c(1.7+time, -2-1.4*time) # Mean
        sigma <- 3*matrix(c(1.5-0.5*time, -0.5*time, -0.5*time, 1.5-0.5*time), 2, 2) # Covariance matrix
        weights[i] <- dgamma(10*time+2, 7.5, rate = 1)
      } else if (i == 3) {
        mu <- c(-2-time, 1.5+1.7*time) # Mean
        sigma <- matrix(c(1.6, -1.2*time, -1.2*time, 1.6), 2, 2) # Covariance matrix
        weights[i] <- dgamma(10*time+2, 7.5, rate = 1)
      } else if (i == 4) {
        mu <- c(2, 2-1.5*time) # Mean
        sigma <- 1.5*matrix(c(1, 0.9-0.3*time, 0.9-0.3*time, 1), 2, 2) # Covariance matrix
        weights[i] <- dgamma(10*(0.5+sign(0.5-time)*(time-0.5))+2, 7.5, rate = 1)
      }
        
      # Calculate the Density of the selected Gaussian Distribution at (x, y)
      x_mu <- rep(0,2)
      x_mu[1] <- x[1] - mu[1]
      x_mu[2] <- x[2] - mu[2]
      inv_sigma <- solve(sigma)
      exponent <- -0.5 * t(x_mu) %*% inv_sigma %*% x_mu
      densities[i] <- exp(exponent) / sqrt((2 * pi) ^ 2 * det(sigma))
    }
      
    # Calculate the Mixture Density by Summing the Weighted Densities
    weights <- c(0.7*weights[1], 0.8*weights[2], 0.4*weights[3], 0.25*weights[4])
    weights <- weights/sum(weights)
    mixture_density <- sum(weights * densities)
      
    # Return the Mixture Density
    return(mixture_density)
  }
    
  int <- adaptIntegrate(mixture, lowerLimit = c(-6, -6), upperLimit = c(6, 6))
  return(int$integral)
}

compute_integral <- function(time_refinement = 0.01){
  integral_values <- sapply(seq(0, 1, by = time_refinement), integral_time)
  integral <- mean(integral_values)
  return (integral)
}