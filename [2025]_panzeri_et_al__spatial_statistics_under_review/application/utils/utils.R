# CONVERT: FROM mesh.1.5D TO SpatialLinesDataFrame --------------------------
setGeneric("as.SpatialLinesDataFrame", function(x) standardGeneric("as.SpatialLinesDataFrame"))

# as.SpatialLinesDataFrame function overloaded for mesh.1.5D class.
# as.SpatialLinesDataFrame.mesh.1.5D <- function(x){
#   nodes = mesh$nodes
#   edges = mesh$edges
#   
#   x0 = vector(mode="numeric",length=nrow(edges))
#   y0 = x0
#   x1 = x0
#   y1 = x0
#   
#   for(e in 1:nrow(edges)){
#     x0[e] = nodes[edges[e,1],1]
#     y0[e] = nodes[edges[e,1],2]
#     x1[e] = nodes[edges[e,2],1]
#     y1[e] = nodes[edges[e,2],2]
#   }
#   
#   psp_ =psp(x0,y0,x1,y1, window=owin(xrange=c(min(x0,x1),max(x0,x1)), 
#                                      yrange=c(min(y0,y1),max(y0,y1))) ) 
#   # maptools 
#   SpatialLines_ = as.SpatialLines.psp(psp_)
#   
#   # rgeos
#   df <- data.frame(len = sapply(1:length(SpatialLines_), function(i) gLength(SpatialLines_[i, ])))
#   rownames(df) <- sapply(1:length(SpatialLines_), function(i) SpatialLines_@lines[[i]]@ID)
#   
#   # object to return
#   SpatialLinesDf <- SpatialLinesDataFrame(SpatialLines_, data = df)
#   
#   return(SpatialLinesDf)
# }

# CONVERT: TO mesh.1.5D --------------------------------------------------------
setGeneric("as.mesh.1.5D", function(x) standardGeneric("as.mesh.1.5D"))

# as.mesh.1.5D function overloading for SpatialLinesDataFrames
setMethod("as.mesh.1.5D", signature="SpatialLinesDataFrame", function(x){
  spat.stat.linnet = maptools::as.linnet.SpatialLines(x)
  nodes = cbind(spat.stat.linnet$vertices$x, spat.stat.linnet$vertices$y)
  edges = cbind(spat.stat.linnet$from, spat.stat.linnet$to)
  mesh = create.mesh.1.5D(nodes, edges)
})

# as.mesh.1.5D function overloading for linnet
setMethod("as.mesh.1.5D", signature="linnet", function(x){
  nodes = cbind(x$vertices$x, x$vertices$y)
  edges = cbind(x$from, x$to)
  mesh = create.mesh.1.5D(nodes, edges)
  return(mesh)
})

# as.mesh.1.5D function overloading for sfnetwork
setMethod("as.mesh.1.5D", signature="sfnetwork", function(x){
  nodes = st_coordinates(x, "nodes")
  edges = st_as_sf(x, "edges")
  edges = cbind(edges$from, edges$to) 
  mesh = create.mesh.1.5D(nodes, edges)
  return(mesh)
})

# CONVERT: FROM mesh.1.5D ------------------------------------------------------

# Converts fdaPDE Linear Network mesh into lattice objects. 
# returns a list containing a matrix of lattice nodes 
# sand a sparse adjacency matrix.
as.lattice <- function(mesh){
  
  # diffusionMaps nb. z coords are set to 0. 
  nodes.lattice = cbind(mesh$nodes, rep(0,times=nnodes))
  
  # Sparse matrix (spam package)
  adj_matrix = spam(0,nrow=nnodes, ncol=nnodes)
  nedges = nrow(mesh$edges)
  for( e in 1:nedges){
    sx = mesh$edges[e,1]
    dx = mesh$edges[e,2]
    
    adj_matrix[sx, dx] = 1
    adj_matrix[dx, sx] = 1
  }
  adj_matrix = spam::cleanup(adj_matrix)
  
  out = list(nodes.lattice = nodes.lattice,
             adj_matrix = adj_matrix)
  
  return(out)
}

# as.linnet function overloading for mesh.1.5D class 
as.linnet.mesh.1.5D <- function(X,...){
  vertices = ppp(x = X$nodes[,1], 
                 y = X$nodes[,2], 
                 window = owin(xrange = c(min(X$nodes[,1]),max(X$nodes[,1])),
                               yrange = c(min(X$nodes[,2]),max(X$nodes[,2]))))
  spat.stat.linnet = linnet(vertices = vertices, edges = X$edges, sparse = T)
  
  return(spat.stat.linnet)
}

# computing network distance matrix --------------------------------------------
as.LPP <- function(points_,L){
  return(as.lpp(x=points_[,1],y=points_[,2],L=L))
}

compute_dist_matrix <- function(points1, points2, L){
  lpp1 <- as.LPP(points1, L)
  lpp2 <- as.LPP(points2, L)
  ND <- crossdist.lpp(lpp1,lpp2)
  return( ND )
}

# MESH NORMALIZATION -----------------------------------------------------------
normalize.mesh <- function(mesh){
  
  x.m = mean(mesh$nodes[,1])
  y.m = mean(mesh$nodes[,2])
  
  x.sd = sd(mesh$nodes[,1])
  y.sd = sd(mesh$nodes[,1])
  
  x.norm = (mesh$nodes[,1] - x.m)/x.sd
  y.norm = (mesh$nodes[,2] - y.m)/y.sd
  
  mesh.norm = create.mesh.1.5D(nodes=cbind(x.norm,y.norm), edges = mesh$edges)  
  return(mesh.norm)
}

normalize.mesh.unit <- function(mesh){
  
  x.min = min(mesh$nodes[,1])
  y.min = min(mesh$nodes[,2])
  
  x.max = max(mesh$nodes[,1])
  y.max = max(mesh$nodes[,2])
  
  x.norm = (mesh$nodes[,1] - x.min)/(x.max-x.min)
  y.norm = (mesh$nodes[,2] - y.min)/(y.max-y.min)
  
  mesh.norm = create.mesh.1.5D(nodes=cbind(x.norm,y.norm), edges = mesh$edges)  
  
  ret = list(mesh=mesh.norm,
             x.min=x.min, y.min = y.min, x.max= x.max, y.max= y.max)
  
  return(ret)
}


# post processing --------------------------------------------------------------
# boxplot_RMSE <- function(RMSE,
#                          methods,
#                          methods.names,
#                          nsim,
#                          title.size=20,
#                          begin=0.95, #color
#                          end=0.25,   #color
#                          width =0.75,
#                          n = c(50, 100, 150, 250),
#                          title="RMSE")
# {
#   RMSE = RMSE[, methods]
#   
#   METHODS = rep(methods.names[methods], each = nrow(RMSE)) 
#   N = as.character(rep(rep(n, each = nsim), times = sum(methods)))
#   
#   RMSE =  as.vector(RMSE)
#   dataFrame = data.frame(RMSE=RMSE, METHODS = METHODS, N = N)
#   
#   MyTheme <- theme(
#     axis.text = element_text(size=title.size-5),
#     axis.title = element_text(size=title.size),
#     title = element_text(size=title.size),
#     plot.title = element_text(hjust = 0.5),
#     legend.text = element_text(size=title.size-5),
#     legend.key.size = unit(1,"cm"),
#     legend.key.height = unit(1,"cm"),
#     legend.title = element_blank(),
#     legend.background = element_rect(fill="white", color="black",
#                                      linewidth =c(1,0.5))
#   )
#   
#   border_col = darken(viridis(length(methods), begin=end,end=begin), amount=0.25)
#   fill_col = viridis(length(methods), begin=end, end=begin)
#   
#   BORDER = c()
#   FILL = c()
#   for(i in 1:length(methods)){
#     if(methods[i]){ 
#       FILL = append(FILL, fill_col[i])
#       BORDER = append(BORDER, border_col[i])
#     }
#   }
#   
#   dataFrame$METHODS = factor(dataFrame$METHODS, 
#                              levels=methods.names) 
#   
#   p<-ggplot(dataFrame)+
#     geom_boxplot(aes(x=N,
#                      y=RMSE, group=interaction(METHODS,N),
#                      fill=METHODS, color = METHODS))+
#     scale_x_discrete(limits=as.character(n))+
#     labs(x="", y="",
#          title=title)+
#     scale_fill_manual(values = FILL) +
#     scale_color_manual(values= BORDER) + 
#     MyTheme + 
#     theme(#plot.title=element_blank(),
#       axis.ticks.x = element_blank(),
#       legend.position = c(0.85,0.85) )
#   return(p)  
# }
# 
# # FEM, oggetto FEM da plottare (fem = FEM(coeff, FEMbasis))
# # title, titolo plot
# # legend.pos, posizione legenda
# # color.min, valore del colore più piccolo da inserire nella scala colore
# # color.max, valore del colore più grande da inserire nella scala colore
# # ratio, ratio fra asse x e asse y (consiglio lasciare valore di default)
# # palette, color palette (es jet.col, viridis, magma, inferno, plasma)
# # num.col, numero di colori della palette da generare
# # title.size, dimensione titolo plot
# R_plot_graph.ggplot2.2<-function(FEM,
#                                  title="", 
#                                  line.size=0.5,
#                                  legend.pos="right",
#                                  color.min=min(FEM$coeff), 
#                                  color.max=max(FEM$coeff),
#                                  ratio = diff(range(mesh$nodes[,1]))/diff(range(mesh$nodes[,2]) ),
#                                  palette = jet.col,
#                                  num.color = 28,
#                                  title.size = 26,
#                                  color.gradient = T){
#   
#   mesh=FEM$FEMbasis$mesh
#   x=vector(mode="double")
#   y=vector(mode="double")
#   coef=vector(mode="double")
#   grp.nodes=vector(mode="integer")
#   
#   num_edges= dim(mesh$edges)[1]
#   for(e in 1:num_edges){
#     x = append(x, c(mesh$nodes[ mesh$edges[e,1], 1], mesh$nodes[ mesh$edges[e,2], 1]))
#     y = append(y, c(mesh$nodes[ mesh$edges[e,1], 2], mesh$nodes[ mesh$edges[e,2], 2]))
#     coef=append(coef, rep((FEM$coeff[mesh$edges[e,1]] + FEM$coeff[mesh$edges[e,2]])/2,times=2) )  
#     grp.nodes = append(grp.nodes, rep(e,times=2))
#   }
#   
#   p <- palette(n=num.color,alpha=1)
#   
#   MyTheme <- theme(
#     axis.text = element_text(size=(title.size-2)),
#     axis.title = element_text(size=title.size),
#     title = element_text(size=title.size),
#     legend.text = element_text(size=(title.size-6)),
#     legend.key.size = unit(1,"cm"),
#     legend.position = legend.pos,
#   )
#   
#   data=data.frame(x,y,grp.nodes,coef)
#   MyTheme <- theme(
#     axis.text = element_text(size=(title.size-2)),
#     axis.title = element_text(size=title.size),
#     title = element_text(size=title.size),
#     legend.text = element_text(size=(title.size-6)),
#     legend.key.size = unit(1,"cm"),
#     legend.position = legend.pos
#   )
#   
#   gplot <- ggplot(data=data, aes(x=x,y=y,group=grp.nodes)) + 
#     geom_point(alpha=0.0) + 
#     geom_line(aes(color=coef), size=line.size)+
#     scale_color_gradientn(colours=p, limits = c(color.min, color.max))+
#     labs(x="",y="",color="", title=title) +  
#     coord_fixed(ratio=ratio) + 
#     theme_void() +
#     MyTheme +
#     theme(plot.title = element_text(hjust=0.5),
#           legend.title = element_blank(),
#           axis.title = element_blank(),
#           axis.text.x=element_blank(),
#           axis.text.y=element_blank(),
#           legend.key.width = unit(0.5,"cm"),
#           legend.position = legend.pos)
#   
#   if(!color.gradient){
#     gplot <- gplot + theme(legend.position = "none") 
#   }
#   
#   return (gplot)
# }
# 
# plot_density<-function(true.density,
#                        max.col,
#                        min.col,
#                        main =  "Density", 
#                        palette = NULL, #palette = NULL (ggplot default colors), "viridis", "magma"
#                        line.size=1)
# {
#   
#   if(!is.null(palette)){
#     if(palette=="viridis")
#       p=viridis
#     else if(palette=="magma")
#       p=magma
#     else if(palette=="inferno")
#       p=inferno
#   }else{
#     p=jet.col
#   }
#   
#   max.col = max(max.col, true.density$coeff)
#   min.col = min(min.col, true.density$coeff)
#   
#   density.plot <- R_plot_graph.ggplot2.2( true.density, # FEM object
#                                           line.size = line.size,
#                                           color.min = min.col,
#                                           color.max = max.col,
#                                           title = main,
#                                           palette=p,
#                                           legend.pos = "right")
#   
#   ret = list(max.col = max.col, min.col = min.col, density.plot = density.plot)  
#   return(ret)
# }
# 
# plot_estimates <-function(estimates, # list of estimates
#                           true.density,
#                           methods,
#                           methods.names =  c(rep("",times=length(estimates))), 
#                           true.density.main = "Density",
#                           palette = NULL, #palette = NULL (ggplot default colors), "viridis", "magma"
#                           line.size=1)
# {
#   
#   if(!is.null(palette)){
#     if(palette=="viridis")
#       p=viridis
#     else if(palette=="magma")
#       p=magma
#     else if(palette=="inferno")
#       p=inferno
#   }else{
#     p=jet.col
#   }
#   
#   num_edges= dim(mesh$edges)[1]
#   coef=matrix(0, nrow= num_edges, ncol=length(estimates) )
#   
#   for(i in 1:(length(estimates))){
#     for(e in 1:num_edges){
#       
#       coef[e,i]= (estimates[[i]]$coeff[mesh$edges[e,1]] + estimates[[i]]$coeff[mesh$edges[e,2]])/2  
#       
#     }
#   }
#   
#   max.col = max(coef)
#   min.col = min(coef)
#   
#   estimates.plot = list()
#   
#   for(i in 1:(length(estimates))){
#     
#     estimates.plot[[i]] = R_plot_graph.ggplot2.2( estimates[[i]], # FEM object
#                                                   line.size = line.size,
#                                                   color.min = min.col,
#                                                   color.max = max.col,
#                                                   title = methods.names[[i]],
#                                                   palette=p,
#                                                   legend.pos = "right")
#     
#   }
#   
#   ret.list = plot_density(true.density, max.col, min.col, main =  true.density.main, 
#                           palette = palette, 
#                           line.size=line.size)
#   
#   
#   ret = list(estimates.plot = estimates.plot, density.plot = ret.list$density.plot) 
#   return(ret)
#   
# }
