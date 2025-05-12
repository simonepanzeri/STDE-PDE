if(!require(pacman)) install.packages("pacman")
pacman::p_load("ggplot2", "ggforce", "viridis", "plotrix")

# Overload plot function for class mesh.1.5D
plot.mesh.1.5D <- function(x, ...){
  mesh = x
  num_edges = dim(mesh$edges)[1]
  
  x=vector(mode="double", length=2*num_edges)
  y=vector(mode="double", length=2*num_edges)
  grp.nodes=vector(mode="integer", length=2*num_edges)
  
  for(e in 1:num_edges){
    x[(2*(e-1)+1):(2*(e-1)+2)] = c(mesh$nodes[mesh$edges[e,1], 1], mesh$nodes[mesh$edges[e,2],1])
    y[(2*(e-1)+1):(2*(e-1)+2)] = c(mesh$nodes[mesh$edges[e,1], 2], mesh$nodes[mesh$edges[e,2],2])
    grp.nodes[(2*(e-1)+1):(2*(e-1)+2)] = rep(e,times=2)
  }
  
  data_plot <- data.frame(x=x, y=y, grp.nodes)
  
  margin_size <- 0.1
  plot_height <- 2
  
  ggplot(data_plot) +
    geom_link2(aes(x = x, y = y, group = grp.nodes),
               lineend = 'round', n = 1, ...) +
    labs(x="",y="",color="", title="") +  
    coord_fixed(ratio=1) + theme_void() +
    theme(plot.margin = unit(c(margin_size, margin_size, margin_size, margin_size), "cm"))
  
}

# Overloaded plot function for class FEM
plot.FEM <- function(x, m = NULL, M = NULL, colormap = "heat.colors",
                     filename = NULL, showLegend = FALSE, ...){
  
  FEM <- x
  mesh <- FEM$FEMbasis$mesh
  num_edges= dim(mesh$edges)[1]
  
  if(is.null(m)){
    m <- min(FEM$coeff, na.rm = TRUE)
  }
  
  if(is.null(M)){
    M <- max(FEM$coeff, na.rm = TRUE)
  }
  
  x=vector(mode="double", length=2*num_edges)
  y=vector(mode="double", length=2*num_edges)
  coeff=vector(mode="double", length=2*num_edges)
  grp.nodes=vector(mode="integer", length=2*num_edges)
  
  for(e in 1:num_edges){
    x[(2*(e-1)+1):(2*(e-1)+2)] = c(mesh$nodes[mesh$edges[e,1], 1], mesh$nodes[mesh$edges[e,2],1])
    y[(2*(e-1)+1):(2*(e-1)+2)]=  c(mesh$nodes[mesh$edges[e,1], 2], mesh$nodes[ mesh$edges[e,2],2])
    coeff[(2*(e-1)+1):(2*(e-1)+2)]= c(FEM$coeff[mesh$edges[e,1]],
                                      FEM$coeff[mesh$edges[e,2]])  
    grp.nodes[(2*(e-1)+1):(2*(e-1)+2)] = rep(e,times=2)
  }
  
  data_plot <- data.frame(x=x, y=y,
                          coeff=coeff, grp.nodes)
  
  margin_size <- 0.1
  plot_height <- 2
  
  color_palette <- match.fun(colormap)
  color_palette <- color_palette(100)
  
  plot <- ggplot(data_plot) +
          geom_link2(aes(x = x, y = y, colour = coeff, group = grp.nodes),
                     lineend = 'round', n = 10, ...) +
          labs(x="",y="",color="", title="") +  
          coord_fixed(ratio=1) + theme_void() +
          theme(plot.margin = unit(c(margin_size, margin_size, margin_size, margin_size), "cm")) +
          scale_color_gradientn(colors = color_palette, limits = c(m, M))
  
  if(showLegend){
    
    plot <- plot + theme(legend.position = "right",
                         legend.direction = "vertical",
                         legend.box = "vertical",
                         legend.key.width = unit(0.25, "in"),
                         legend.key.height = unit(1, "in"))
    
    w <- 7 # width
    h <- 6 # height
                   
  } else {
    
    plot <- plot + theme(legend.position = "none")
    
    w <- 7 # width
    h <- 7 # height
    
  }
  
  if (!is.null(filename)){
    
    #x11(width = w, height = h)
    pdf(paste0(filename, ".pdf"), family = "serif", width = w, height = h)
    print(plot)
    dev.off()
    
  } else {
    
    plot
    
  }
}

# Plot the trajectories followed by 4 sources of intensity
plot.trajectories <- function(mesh, sources, colors = NULL, linewidth = NULL,...){
  
  if(is.null(colors)) colors = brewer.pal(length(sources), "Spectral")[1:length(sources)]
  if(is.null(linewidth)) linewidth = 1.25
  linewidth <- c(linewidth, 1.5)
  
  num_edges = dim(mesh$edges)[1]
  
  x=vector(mode="double", length=2*num_edges)
  y=vector(mode="double", length=2*num_edges)
  grp.nodes=vector(mode="integer", length=2*num_edges)
  col.nodes=vector(mode="integer", length=2*num_edges)
  lw.nodes=vector(mode="integer", length=2*num_edges)
  
  for(e in 1:num_edges){
    x[(2*(e-1)+1):(2*(e-1)+2)] = c(mesh$nodes[mesh$edges[e,1], 1], mesh$nodes[mesh$edges[e,2],1])
    y[(2*(e-1)+1):(2*(e-1)+2)] = c(mesh$nodes[mesh$edges[e,1], 2], mesh$nodes[mesh$edges[e,2],2])
    grp.nodes[(2*(e-1)+1):(2*(e-1)+2)] = rep(e,times=2)
    col.nodes[(2*(e-1)+1):(2*(e-1)+2)] = rep(0,times=2)
    lw.nodes[(2*(e-1)+1):(2*(e-1)+2)] = rep(1,times=2)
    for(s in 1:length(sources)){
      if(mesh$edges[e,1] %in% sources[[s]] & mesh$edges[e,2] %in% sources[[s]]){
        col.nodes[(2*(e-1)+1):(2*(e-1)+2)] = rep(s,times=2)
        lw.nodes[(2*(e-1)+1):(2*(e-1)+2)] = rep(2,times=2)
      }
    }
  }
  
  data_plot <- data.frame(x=x, y=y, grp.nodes)
  
  margin_size <- 0.1
  plot_height <- 2
  
  palette <- c("black", colors)
  
  start <- vector(mode="integer", length=length(s))
  end <- vector(mode="integer", length=length(s))
  labels <- LETTERS[1:(2*length(sources))]
  for(s in 1:length(sources)){
    start[s] <- sources[[s]][1]
    end[s] <- sources[[s]][length(sources[[s]])]
  }
  
  ggplot(data_plot) +
    geom_link2(aes(x = x, y = y, group = grp.nodes), color = palette[col.nodes+1],
               lineend = 'square', n = 1, linewidth = linewidth[lw.nodes]) +
    labs(x="",y="",color="", title="") +  
    coord_fixed(ratio=1) + theme_void() +
    theme(plot.margin = unit(c(margin_size, margin_size, margin_size, margin_size), "cm")) +
    geom_label(data = data.frame(x = mesh$nodes[start,1], y = mesh$nodes[start,2]),
               aes(x = x, y = y, label = labels[seq(1,2*length(sources)-1,by=2)]),
               color = palette[2:length(palette)], family = "serif", fontface = "bold",
               alpha = 0.9, nudge_x = 0, nudge_y = 0.03, size = 6) +
    geom_label(data = data.frame(x = mesh$nodes[end,1], y = mesh$nodes[end,2]),
              aes(x = x, y = y, label = labels[seq(2,2*length(sources),by=2)]),
              color = palette[2:length(palette)], family = "serif", fontface = "bold",
              alpha = 0.9, nudge_x = 0.02, nudge_y = -0.03, size = 6)
  
}
