#'############################################################################'#
#'################## HELPER FUNCTIONS FOR SIMULATION STUDIES #################'#
#'############################################################################'#

set.seed(23)

## FIND THE PATH OF size CONSECUTIVE POINTS IN THE spatstat.linnet L STARTING FROM THE idx-th POINT -----
find.closest.nodes <- function(idx, L, size) {
  
  #pb <- txtProgressBar(min = 0, max = size-1, style = 3, char = "=")
  
  idxs <- idx
  
  nodes.lpp <- ppp(x = L$vertices$x, y = L$vertices$y, window = L$window)
  
  while(length(idxs) < size){
    
    PP <- ppp(x = L$vertices$x[idxs[length(idxs)]], y = L$vertices$y[idxs[length(idxs)]], window = nodes.lpp$window)
    ND <- crossdist.lpp(lpp(nodes.lpp, L), lpp(PP, L))
    ND <- cbind(ND, seq(1,nrow(ND),by=1))
    ND <- as.data.frame(ND)
    colnames(ND) <- c('dist','idxs')
    ND <- ND[order(ND$dist, decreasing = FALSE),]
    
    i <- 1
    while(ND$idxs[i] %in% idxs){
        i <- i+1
    }
    idxs <- c(idxs, ND$idxs[i])
    
    #setTxtProgressBar(pb, length(idxs))
    
  }
  
  return (idxs[1:size])
  
}

## COMPUTE THE SHORTEST PATH STARTING FROM NODE i TO NODE j ON THE spatstat.linnet L -----
shortestpath <- function(L, i, j) {
  L <- as.linnet(L, sparse=FALSE)
  d <- L$dpath
  m <- L$m
  to <- L$to
  from <- L$from
  path <- i
  leader <- i
  repeat {
    k <- setdiff(which(m[leader,]), path)
    leader <- k[which.min(d[i,k] + d[k, j])]
    path <- c(path, leader)
    if(leader == j) break
  }
  return(path)
}

joinpoints <- function(X, i, ...) {
  co <- coords(X[i])
  lines(co[,1], co[,2], ...)
}

## COMPUTE THE PATHS FOLLOWED BY THE SOURCES IN THE spatstat.linnet L OVER TIME -----
compute.sources2 <- function(L){
  
  # Eastbourne
  idxs1 <- shortestpath(L, 332, 38)
  idxs2 <- shortestpath(L, 62, 126)
  idxs3 <- shortestpath(L, 89, 116)
  idxs4 <- shortestpath(L, 319, 342)
    
  # x11()
  # plot(mesh, linewidth=1) +
  #   geom_point(data = data.frame(x = mesh$nodes[,1], y = mesh$nodes[,2]),
  #                                aes(x = x, y = y), color = "red3", size = 3) +
  #   geom_text(data = data.frame(x = mesh$nodes[,1], y = mesh$nodes[,2]),
  #             aes(x = x, y = y, label = 1:nrow(mesh$nodes)), vjust = -0.5, color = "blue")
  # 
  # plot(mesh$nodes, pch=19, col='black')
  # points(mesh$nodes[idxs1,1], mesh$nodes[idxs1,2], pch = 19, col = 'blue')
  # points(mesh$nodes[idxs2,1], mesh$nodes[idxs2,2], pch = 19, col = 'red')
  # points(mesh$nodes[idxs3,1], mesh$nodes[idxs3,2], pch = 19, col = 'green')
  # points(mesh$nodes[idxs4,1], mesh$nodes[idxs4,2], pch = 19, col = 'red')
  
  sources <- list(idxs1, idxs2, idxs3, idxs4)

  return (sources)
  
}

## COMPUTE THE PATHS FOLLOWED BY THE ADDITIONAL SOURCES IN THE spatstat.linnet L OVER TIME -----
compute.sources1 <- function(L){
  
  # Eastbourne
  idxs <- shortestpath(L, 20, 15)
  idxs <- c(idxs, shortestpath(L, 176, 57))
  idxs <- c(idxs, rev(shortestpath(L, 52, 233)))

  # x11()
  # plot(mesh, linewidth=1) +
  #   geom_point(data = data.frame(x = mesh$nodes[,1], y = mesh$nodes[,2]),
  #                                aes(x = x, y = y), color = "red3", size = 3) +
  #   geom_text(data = data.frame(x = mesh$nodes[,1], y = mesh$nodes[,2]),
  #             aes(x = x, y = y, label = 1:nrow(mesh$nodes)), vjust = -0.5, color = "blue")
  # 
  # plot(mesh$nodes, pch=19, col='black')
  # points(mesh$nodes[idxs,1], mesh$nodes[idxs,2], pch = 19, col = 'blue')
  
  sources <- list(idxs)
  
  return (sources)
  
}

## SELECT THE INDICES OF THE SOURCES, COMPUTING THE PATH COMPLETION PERCENTAGE ACCORDING TO t AND tlim -----
select.sources <- function(t, tlim, sources){
  
  idxs <- matrix(nrow = length(t), ncol = length(sources))
  for(time_index in 1:length(t)){
    
    path_percentage <- (t[time_index]-tlim[1])/(diff(tlim))
    
    for(i in 1:length(sources)){
      
      idx_source <- floor(path_percentage*length(sources[[i]]))
      idx_source <- ifelse(idx_source == 0, 1, idx_source)
      idxs[time_index,i] <- sources[[i]][idx_source]
      
    }
  }
  
  return (idxs)

}

## COMPUTE THE INTEGRAL OF A SPATIO-TEMPORAL INTENSITY FUNCTION -----
compute.integral.at.t <- function(fun, f, t){
  
  integral.at.t <- sum(fdaPDE:::CPP_get.FEM.Mass.Matrix(f) %*%
                         fun(x = f$mesh$nodes[,1], y = f$mesh$nodes[,2], t = t, weight = 1000))
    
  return (integral.at.t / 1000)
  
}

compute.integral <- function(fun, f, time_refinement = 0.01, tlim = c(0, 1)){
  
  times <- seq(from = tlim[1], to = tlim[2], by = time_refinement)
  
  pb <- txtProgressBar(min = tlim[1], max = tlim[2], style = 3, char = "=")
  integral_values <- lapply(times, function(i){
    a <- compute.integral.at.t(fun = fun, f = f, t = i)
    setTxtProgressBar(pb, i)
    return (a)
    }
  )
  cat("\n")
  
  integral <- mean(unlist(integral_values))
  
  return (integral)
}

compute.weight <- function(fun, f, N){
  
  return (N / compute.integral(fun = fun, f = f))
  
}

## INHOMOGENEOUS SPATIO-TEMPORAL INTENSITY FUNCTION [SIMULATION 1] -----
inhomogeneous_intensity1 <- function(x, y, t, seg, tp,
                                     weight = w, sigma = 0.3, tlim = c(0,1),
                                     L = spatstat.linnet){
  
  sources <- compute.sources1(L)
  current_sources <- select.sources(t = t, tlim = tlim, sources = sources)
  
  nodes.lpp <- ppp(x = L$vertices$x, y = L$vertices$y, window = L$window)
  
  PP <- ppp(x = x, y = y, window = nodes.lpp$window)
  ND <- crossdist.lpp(lpp(nodes.lpp, L), lpp(PP, L))
  
  coeff <- rep(0, length(t))
  
  if(length(t) == 1){
    
    return ( weight * 1/sqrt(2*pi*sigma^2) * exp(-ND[current_sources[,1],]^2/(2*sigma^2)) )
    
  } else {
    
    coeff <- rep(0, length(x))
    
    for(i in 1:length(x)){
      
      coeff[i] <- weight * 1/sqrt(2*pi*sigma^2) * exp(-ND[current_sources[i,1],i]^2/(2*sigma^2))
      
    }
    
    return ( coeff )
    
  }
  
}

## INHOMOGENEOUS SPATIO-TEMPORAL INTENSITY FUNCTION [SIMULATION 2] -----
inhomogeneous_intensity2 <- function(x, y, t, seg, tp,
                                     weight = w, sigma = 0.075, tlim = c(0,1),
                                     L = spatstat.linnet){
  
  sources <- compute.sources2(L)
  current_sources <- select.sources(t = t, tlim = tlim, sources = sources)
  
  nodes.lpp <- ppp(x = L$vertices$x, y = L$vertices$y, window = L$window)
  
  PP <- ppp(x = x, y = y, window = nodes.lpp$window)
  ND <- crossdist.lpp(lpp(nodes.lpp, L), lpp(PP, L))
  
  coeff <- rep(0, length(t))
  
  if(length(t) == 1){
    
    return ( weight * ( 0.25 * 1/sqrt(2*pi*sigma^2) * exp(-ND[current_sources[,1],]^2/(2*sigma^2)) + 
                          0.25 * 1/sqrt(2*pi*sigma^2) * exp(-ND[current_sources[,2],]^2/(2*sigma^2)) +
                          0.25 * 1/sqrt(2*pi*sigma^2) * exp(-ND[current_sources[,3],]^2/(2*sigma^2)) +
                          0.25 * 1/sqrt(2*pi*sigma^2) * exp(-ND[current_sources[,4],]^2/(2*sigma^2))) )
    
  } else {
    
    coeff <- rep(0, length(x))
    
    for(i in 1:length(x)){
      
      coeff[i] <- weight * ( 0.25 * 1/sqrt(2*pi*sigma^2) * exp(-ND[current_sources[i,1],i]^2/(2*sigma^2)) + 
                               0.25 * 1/sqrt(2*pi*sigma^2) * exp(-ND[current_sources[i,2],i]^2/(2*sigma^2)) +
                               0.25 * 1/sqrt(2*pi*sigma^2) * exp(-ND[current_sources[i,3],i]^2/(2*sigma^2)) +
                               0.25 * 1/sqrt(2*pi*sigma^2) * exp(-ND[current_sources[i,4],i]^2/(2*sigma^2)))
      
    }
    
    return ( coeff )
    
  }
  
}

## HOMOGENEOUS SPATIO-TEMPORAL INTENSITY FUNCTION [SIMULATION 3] -----
homogeneous_intensity3 <- function(x, y, t, seg, tp,
                                   n_obs = N, tlim = c(0,1),
                                   L = spatstat.linnet){
  
  domain_volume <- volume(L)*diff(tlim) # equivalent to: domain_area * (Tf-Ti)
  coeff <- rep(n_obs/domain_volume, length(x))
  
  return ( coeff )
  
}

## INHOMOGENEOUS SPATIO-TEMPORAL INTENSITY FUNCTION [SIMULATION 4] -----
inhomogeneous_intensity4 <- function(x, y, t, seg, tp,
                                     weight = w, tlim = c(0,1),
                                     L = spatstat.linnet){
  
  coeff <- weight * 1.5 * exp(sin(t * (x + y)))
  
  return ( coeff )
  
}

## INHOMOGENEOUS SPATIO-TEMPORAL INTENSITY FUNCTION [SIMULATION 5] -----
inhomogeneous_intensity5 <- function(x, y, t, seg, tp,
                                     weight = w, sigma = 0.3, tlim = c(0,1),
                                     L = spatstat.linnet){
  
  source <- as.matrix(rep(258, length(t)))
  
  nodes.lpp <- ppp(x = L$vertices$x, y = L$vertices$y, window = L$window)
  
  PP <- ppp(x = x, y = y, window = nodes.lpp$window)
  ND <- crossdist.lpp(lpp(nodes.lpp, L), lpp(PP, L))
  
  coeff <- rep(0, length(t))
  
  if(length(t) == 1){
    
    return ( weight *
               ( 1/sqrt(2*pi*sigma^2) * exp(-ND[source[,1],]^2/(2*sigma^2)) ) *
               ( exp(-t) ) )
    
  } else {
    
    coeff <- rep(0, length(x))
    
    for(i in 1:length(t)){
      
      coeff[i] <- ( weight *
                      ( 1/sqrt(2*pi*sigma^2) * exp(-ND[source[i,1],i]^2/(2*sigma^2)) ) *
                      ( exp(-(t[i])) ) )
        
    }
    
    return ( coeff )
    
  }
  
}

## COMPUTE THE KERNEL DENSITY ESTIMATE ON A LINEAR NETWORK USING THE HEAT EQUATION -----
density1Dkernel <- function(X,at=c("points","pixels"),leaveoneout=FALSE,
                            kernel=NULL,...){
  
  if(!inherits(X, "stlpp")) stop("X should an object of class stlpp")
  
  if(missing(at)) at <- "pixels"
  
  n <- npoints(X)
  
  Xs <- as.lpp.stlpp(X)
  
  Xt <- lpp(X=cbind(X$data$t,rep(0,n)), 
            L=linnet_interval(startp=0, endp=1))
  
  bw <- bw.lppl(X=Xs)
  if(is.null(kernel)){
    IntEstL <- densityHeat(x=Xs,sigma=as.numeric(bw),...)
  } else {
    IntEstL <- densityEqualSplit(x=Xs,sigma=as.numeric(bw),kernel=kernel,...)
  } 

  bw <- bw.lppl(X=Xt)
  if(is.null(kernel)){
    IntEstT <- densityHeat(x=Xt,sigma=as.numeric(bw),...)
  } else {
    IntEstT <- densityEqualSplit(x=Xt,sigma=as.numeric(bw),kernel=kernel,...)
  } 

  tgrid <- IntEstT$xcol
    
  IntEstTv <- IntEstT$v[!is.na(IntEstT$v)]/n
    
  out <- lapply(X=(1:length(IntEstTv)), FUN=function(j){IntEstL*IntEstTv[j]}) 
    
  if(at=="points"){
    t <- X$data$t
    id <- findInterval(t,tgrid)
    out1 <- c()
    for (i in 1:n){
      out1[i] <- out[[id[i]]][Xs[i]]
    }
    out <- out1
  }
  
  if(at=="points") class(out) <- c("numeric","stlppint")
  if(at=="pixels") class(out) <- c("list","stlppint")
  
  attr(out,"lint") <- IntEstL
  attr(out,"tint") <- IntEstT$v[!is.na(IntEstT$v)]

  attr(out,"tgrid") <- tgrid
  attr(out,"time") <- X$data$t
  attr(out,"stlpp") <- X
  
  return(out)
  
}
  
linnet_interval <- function(startp=0, endp=1,...){
  
  Wt <- c(startp, endp)
  vertices <- ppp(x=c(Wt[1],Wt[2]), y=c(0,0), 
                  window=owin(Wt,c(-0.5,0.5)))
  m <- matrix(data=c(FALSE,TRUE,TRUE,FALSE), nrow=2, ncol=2, 
              byrow=TRUE)
  out <- linnet(vertices=vertices, m=m,...)
  
  return(out)
  
}

## COMPUTE ERRORS -----
compute_err <- function(FEM, eval, true){
  
  num <- sum( fdaPDE:::CPP_get.FEM.Mass.Matrix(FEM) %*% ( ( eval - true )^2) )
  den <- sum( fdaPDE:::CPP_get.FEM.Mass.Matrix(FEM) %*% ( ( true )^2 ) )
  return ( num / den )
  
}
