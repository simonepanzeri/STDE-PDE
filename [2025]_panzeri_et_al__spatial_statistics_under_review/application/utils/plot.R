if(!require(pacman)) install.packages("pacman")
pacman::p_load("ggplot2", "ggforce", "viridis", "plotrix", "chromote")

# chrome_path <- "C:/Users/simon/Downloads/chrome-win64/chrome-win64/chrome.exe"
# m <- Chromote$new(browser = Chrome$new(path = chrome_path))
# b <- m$new_session()

ses <- chromote::ChromoteSession$new()
ses$default_timeout <- 10*60

### Graphical parameters -------------------------------------------------------
#map.type = "Esri.WorldTopoMap"
map.type = "CartoDB.Positron"
col_mesh = rgb(60, 60, 60, max = 255)
col_bdy = rgb(128, 128, 128, max = 255)
col_nodes = rgb(0, 0, 0, max = 255)
col_locations = rgb(190, 0, 0, max = 255)

### Overload plot function for class mesh.1.5D ---------------------------------
plot.mesh.1.5D <- function(obj, map = TRUE, bdy = NULL, center = NULL,
                           foldername = NULL, filename = NULL, ...){
  
  if(!map){
    mesh = obj
    num_edges = dim(mesh$edges)[1]
    
    x = vector(mode = "double", length = 2*num_edges)
    y = vector(mode = "double", length = 2*num_edges)
    grp.nodes = vector(mode = "integer", length = 2*num_edges)
    
    for(e in 1:num_edges){
      x[(2*(e-1)+1):(2*(e-1)+2)] = c(mesh$nodes[mesh$edges[e,1], 1], mesh$nodes[mesh$edges[e,2],1])
      y[(2*(e-1)+1):(2*(e-1)+2)] = c(mesh$nodes[mesh$edges[e,1], 2], mesh$nodes[mesh$edges[e,2],2])
      grp.nodes[(2*(e-1)+1):(2*(e-1)+2)] = rep(e,times=2)
    }
    
    data_plot = data.frame(x = x, y = y, grp.nodes)
    
    margin_size = 0.1
    plot_height = 2
    
    plot = ggplot(data_plot) + geom_link2(aes(x = x, y = y, group = grp.nodes),
                                          lineend = 'round', n = 1, ...) +
      labs(x = "", y = "", color = "", title = "") +  
      coord_fixed(ratio = 1/0.6984) + theme_void() +
      theme(plot.margin = unit(c(margin_size, margin_size, margin_size, margin_size), "cm"))
    
    if(!is.null(bdy)){
      bdy_sp = as(bdy, "Spatial")
      bdy_df = fortify(bdy_sp)
      plot = plot + geom_polygon(data = bdy_df, aes(x = long, y = lat, group = group),
                                 fill = col_bdy, color = col_bdy, alpha = 0.3)
    }
    
    if (!is.null(filename) & !is.null(foldername)){
      
      w = 7 # width
      h = 7 # height
      
      #x11(width = w, height = h)
      pdf(paste0(foldername, "/pdf/", filename, ".pdf"), family = "serif", width = w, height = h)
      print(plot)
      dev.off()
      
    } else {
      
      plot
      
    }

  } else {
    
    mesh = obj
    mesh_df = data.frame(lon = mesh$nodes[,1], lat = mesh$nodes[,2])
    mesh_df = st_as_sf(mesh_df, coords = c("lon", "lat"), crs = 4326)
    mesh_linnet = as.linnet(mesh)
    mesh_sfnetwork = as_sfnetwork(mesh_linnet, directed = FALSE, edges_as_lines = TRUE)
    st_crs(mesh_sfnetwork) = 4326
    mesh_sf = st_as_sf(mesh_sfnetwork, "edges")
    
    map = mapview(mesh_sf, legend = FALSE, map.type = map.type, color = col_mesh,
                  layer.name = "road-network", lwd = 1.5, alpha.regions = 1) +
      mapview(mesh_df, legend = FALSE, col.region = col_nodes, layer.name = "mesh-nodes",
              alpha.regions = 1, cex = 1)
    
    if(!is.null(bdy)){
      map = map + mapview(st_sf(data.frame(x = 1), geoemtry = bdy), legend = FALSE,
                          col.region = col_bdy, alpha.regions = 0.3, layer.name = "boundary", lwd = 0.3)
    }
    
    if(!is.null(foldername) & !is.null(filename)){
      
      map@map = removeMapJunk(map, junk = c("zoomControl", "homeButton", "layersControl", "scaleBar",
                                            "drawToolbar", "easyButton", "control"))
      
      if(!is.null(center)){
        map@map = map@map %>% setView(center[1], center[2], zoom = 12)  
      }
      html_filename = paste0(foldername, "/html/mesh.html")
      pdf_filename = paste0(foldername, "/pdf/mesh.pdf")
      htmlwidgets::saveWidget(map@map, html_filename)
      webshot2::webshot(url = html_filename, file = pdf_filename, delay = 10)
      file.rename(paste0(foldername, "/html/mesh.html"),
                  paste0(foldername, "/html/", filename, ".html"))
      file.rename(paste0(foldername, "/pdf/mesh.pdf"),
                  paste0(foldername, "/pdf/", filename, ".pdf"))
      
      # Alternative:
      # html_filename = paste0(foldername, "/html/", filename, ".html")
      # pdf_filename = paste0(foldername, "/pdf/", filename, ".pdf")
      # mapshot2(map@map, url = html_filename, file = pdf_filename, delay = 10)
      
    } else {
      
      map
      
    }
    
  }
}

### Function to plot the sample data -------------------------------------------
plot_sample <- function(obj, mesh, map = TRUE, bdy = NULL, center = NULL,
                        foldername = NULL, filename = NULL,
                        zoom = NULL, cex = NULL, stroke = NULL,
                        box = FALSE, box1 = FALSE, box2 = FALSE, box3 = FALSE, box4 = FALSE, ...){
  
  if(is.null(zoom)) {zoom = 13}
  if(is.null(cex)) {cex = 2.5}
  if(is.null(stroke)) {stroke = 0.75}
  
  if(!map){
    
    locations_df = obj
    
    plot = plot(mesh, map = FALSE, bdy = bdy, linewidth = 0.5)
    plot = plot + geom_point(data = data.frame(x = locations_df[,1], y = locations_df[,2]),
                                               aes(x = x, y = y), color = "red3", size = 1.5)
    plot = plot + theme(legend.position = "none")
    
    if (!is.null(filename) & !is.null(foldername)){

      w = 7 # width
      h = 7 # height

      #x11(width = w, height = h)
      pdf(paste0(foldername, "/pdf/", filename, ".pdf"), family = "serif", width = w, height = h)
      print(plot)
      dev.off()
        
    } else {
        
      plot
        
    }
    
  } else {
    
    locations_df = st_as_sf(locations_df, coords = c("lon", "lat"), crs = 4326)
    
    mesh_df = data.frame(lon = mesh$nodes[,1], lat = mesh$nodes[,2])
    mesh_df = st_as_sf(mesh_df, coords = c("lon", "lat"), crs = 4326)
    mesh_linnet = as.linnet(mesh)
    mesh_sfnetwork = as_sfnetwork(mesh_linnet, directed = FALSE, edges_as_lines = TRUE)
    st_crs(mesh_sfnetwork) = 4326
    mesh_sf = st_as_sf(mesh_sfnetwork, "edges")
    
    map = mapview(mesh_sf, legend = FALSE, map.type = map.type, color = col_mesh,
                  layer.name = "road-network", lwd = 1, alpha.regions = 1) +
      mapview(locations_df, legend = FALSE, col.region = col_locations, alpha.regions = 1,
              layer.name = "road-accidents", cex = cex, stroke = stroke)
    
    if(!is.null(bdy)){
      map = map + mapview(st_sf(data.frame(x = 1), geoemtry = bdy), legend = FALSE,
                          col.region = col_bdy, alpha.regions = 0.3, layer.name = "boundary", lwd = 0.3)
    }
    
    if(box == TRUE){
      neighborhood = neighborhood_
      bbox = st_bbox(c(xmin = neighborhood$xmin, ymin = neighborhood$ymin, xmax = neighborhood$xmax, ymax = neighborhood$ymax), crs = 4326)
      rectangle = st_as_sfc(bbox)
      map = map + mapview(st_sf(data.frame(x = 1), geometry = rectangle),
                          legend = FALSE, color = "forestgreen", alpha.regions = 0, layer.name = "Subnetwork", lwd = 5)
    }
    
    if(box1 == TRUE){
      neigh_ = neighborhood_1
      center_ = st_sfc(st_point(c(neigh_$lon, neigh_$lat)), crs = 4326)
      circle_ = st_buffer(st_transform(center_, crs = 3857), dist = 900, nQuadSegs = 100)
      circle_ = st_transform(circle_, crs = 4326)
      map = map + mapview(st_sf(data.frame(x = 1), geometry = circle_),
                          legend = FALSE, color = "royalblue", alpha.regions = 0, layer.name = "Rondò delle Valli", lwd = 5)
    }
    
    if(box2 == TRUE){
      neigh_ = neighborhood_2
      center_ = st_sfc(st_point(c(neigh_$lon, neigh_$lat)), crs = 4326)
      circle_ = st_buffer(st_transform(center_, crs = 3857), dist = 900, nQuadSegs = 100)
      circle_ = st_transform(circle_, crs = 4326)
      map = map + mapview(st_sf(data.frame(x = 1), geometry = circle_),
                          legend = FALSE, color = "orangered", alpha.regions = 0, layer.name = "Bergamo Città Bassa", lwd = 5)
    }
    
    if(box3 == TRUE){
      neigh_ = neighborhood_3
      center_ = st_sfc(st_point(c(neigh_$lon, neigh_$lat)), crs = 4326)
      circle_ = st_buffer(st_transform(center_, crs = 3857), dist = 900, nQuadSegs = 100)
      circle_ = st_transform(circle_, crs = 4326)
      map = map + mapview(st_sf(data.frame(x = 1), geometry = circle_),
                          legend = FALSE, color = "orangered", alpha.regions = 0, layer.name = "Bergamo Città Alta", lwd = 5)
    }
    
    if(box4 == TRUE){
      neigh_ = neighborhood_4
      center_ = st_sfc(st_point(c(neigh_$lon, neigh_$lat)), crs = 4326)
      circle_ = st_buffer(st_transform(center_, crs = 3857), dist = 900, nQuadSegs = 100)
      circle_ = st_transform(circle_, crs = 4326)
      map = map + mapview(st_sf(data.frame(x = 1), geometry = circle_),
                          legend = FALSE, color = "royalblue", alpha.regions = 0, layer.name = "Raccordo Autostradale", lwd = 5)
    }
    
    if (!is.null(filename) & !is.null(foldername)){
      
      map@map = removeMapJunk(map, junk = c("zoomControl", "homeButton", "layersControl", "scaleBar",
                                        "drawToolbar", "easyButton", "control"))
      
      if(!is.null(center)){
        map@map = map@map %>% setView(center[1], center[2], zoom = zoom)  
      }
      html_filename = paste0(foldername, "/html/sample.html")
      pdf_filename = paste0(foldername, "/pdf/sample.pdf")
      htmlwidgets::saveWidget(map@map, html_filename)
      webshot2::webshot(url = html_filename, file = pdf_filename, delay = 10)
      file.rename(paste0(foldername, "/html/sample.html"),
                  paste0(foldername, "/html/", filename, ".html"))
      file.rename(paste0(foldername, "/pdf/sample.pdf"),
                  paste0(foldername, "/pdf/", filename, ".pdf"))
      
      # Alternative:
      # html_filename = paste0(foldername, "/html/", filename, ".html")
      # pdf_filename = paste0(foldername, "/pdf/", filename, ".pdf")
      # mapshot2(map@map, url = html_filename, file = pdf_filename, delay = 10)
      
    } else {
      
      map
      
    }
  }
  
}

### Overloaded plot function for class FEM -------------------------------------
plot.FEM <- function(obj, map = TRUE, bdy = NULL, center = NULL,
                     colormap = "heat.colors", m = NULL, M = NULL, showLegend = FALSE,
                     foldername = NULL, filename = NULL, zoom = NULL, ...){
  
  if(is.null(zoom)) {zoom = 13}
  
  color_palette = match.fun(colormap)
  color_palette = color_palette(100)
  
  if(is.null(m)){
    m = min(obj$coeff, na.rm = TRUE)
  }
  
  if(is.null(M)){
    M = max(obj$coeff, na.rm = TRUE)
  }
  
  if(!map){
    
    FEM = obj
    mesh = FEM$FEMbasis$mesh
    num_edges = dim(mesh$edges)[1]
    
    x = vector(mode = "double", length = 2*num_edges)
    y = vector(mode = "double", length = 2*num_edges)
    coeff = vector(mode = "double", length = 2*num_edges)
    grp.nodes = vector(mode = "integer", length = 2*num_edges)
    
    for(e in 1:num_edges){
      x[(2*(e-1)+1):(2*(e-1)+2)] = c(mesh$nodes[mesh$edges[e,1], 1], mesh$nodes[mesh$edges[e,2],1])
      y[(2*(e-1)+1):(2*(e-1)+2)] =  c(mesh$nodes[mesh$edges[e,1], 2], mesh$nodes[ mesh$edges[e,2],2])
      coeff[(2*(e-1)+1):(2*(e-1)+2)] = c(FEM$coeff[mesh$edges[e,1]], FEM$coeff[mesh$edges[e,2]])  
      grp.nodes[(2*(e-1)+1):(2*(e-1)+2)] = rep(e, times = 2)
    }
    
    data_plot = data.frame(x = x, y = y, coeff = coeff, grp.nodes)
    
    margin_size = 0.1
    plot_height = 2
    
    plot = ggplot(data_plot) +
      geom_link2(aes(x = x, y = y, colour = coeff, group = grp.nodes),
                 lineend = 'round', n = 10, linewidth = 1.5) +
      labs(x = "", y = "", color = "", title = "") +  
      coord_fixed(ratio = 1/0.6984) + theme_void() +
      coord_cartesian(xlim = c(neighborhood_4$lon - 0.012838, neighborhood_4$lon + 0.012838),
                      ylim = c(neighborhood_4$lat - 0.0025 - 0.008997, neighborhood_4$lat - 0.0025 + 0.008997)) +
      theme(plot.margin = unit(c(margin_size, margin_size, margin_size, margin_size), "cm")) +
      scale_color_gradientn(colors = color_palette, limits = c(m, M))
  
    if(!is.null(bdy)){
      bdy_sp = as(bdy, "Spatial")
      bdy_df = fortify(bdy_sp)
      plot = plot + geom_polygon(data = bdy_df, aes(x = long, y = lat, group = group),
                                 fill = col_bdy, color = col_bdy, alpha = 0.3)
    }
    
    if(showLegend){
      
      plot = plot + theme(legend.position = "right",
                           legend.direction = "vertical",
                           legend.box = "vertical",
                           legend.key.width = unit(0.25, "in"),
                           legend.key.height = unit(1, "in"))
      
      w = 7 # width
      h = 6 # height
      
    } else {
      
      plot = plot + theme(legend.position = "none")
      
      w = 7 # width
      h = 7 # height
      
    }
    
    if (!is.null(filename) & !is.null(foldername)){
      
      #x11(width = w, height = h)
      pdf(paste0(foldername, "/pdf/", filename, ".pdf"), family = "serif", width = w, height = h)
      print(plot)
      dev.off()
      
    } else {
      
      plot
      
    }
    
  } else {
    
    mesh.eval = obj$FEMbasis$mesh
    mesh.eval_linnet = as.linnet(mesh.eval)
    mesh.eval_sfnetwork = as_sfnetwork(mesh_linnet, directed = FALSE, edges_as_lines = TRUE)
    st_crs(mesh.eval_sfnetwork) = 4326
    mesh.eval_sf = st_sf(as.data.frame(st_as_sf(mesh.eval_sfnetwork, "edges")))
    
    pb = txtProgressBar(min = 0, max = nrow(mesh_sf) - 1, style = 3)
    for(e in 1:nrow(mesh_sf)){
    # pb = txtProgressBar(min = 0, max = 500 - 1, style = 3)
    # for(e in 1:500){  
      endpoints_coords = st_coordinates(mesh_sf$geom[e])[,1:2]
      endpoints_coords = endpoints_coords[c(1:nrow(endpoints_coords)),]
      
      dx = (endpoints_coords[2,1] - endpoints_coords[1,1])
      dy = (endpoints_coords[2,2] - endpoints_coords[1,2])
      
      x = endpoints_coords[1,1] + seq(0, 1, length.out = 11) * dx
      y = endpoints_coords[1,2] + seq(0, 1, length.out = 11) * dy
      
      # Split each segment in 10 subsegments (11 points, endpoints included)
      inner_coords = cbind(x, y)
      colors = interpolate.color(z_start = obj$coeff[mesh.eval$edges[e,1]],
                                 z_end = obj$coeff[mesh.eval$edges[e,2]],
                                 m = m, M = M, col = color_palette)
      
      # Create sf linestring with interpolated colors
      for (i in 1:(nrow(inner_coords) - 1)) {
        line = st_linestring(matrix(inner_coords[i:(i + 1), ], ncol = 2, byrow = FALSE))
        line = st_sfc(line)
        if(e == 1 & i == 1){
          all_lines = st_sf(geometry = st_sfc(line), color = colors[i])
        } else {
          all_lines = rbind(all_lines, st_sf(geometry = st_sfc(line), color = colors[i]))
        }
      }
      setTxtProgressBar(pb, e)
    }
    
    #mapview(all_lines, col.regions = unique(all_lines$color), legend = FALSE, lwd = 3)
    
    map = mapview(mesh.eval_sf, legend = FALSE, map.type = map.type, color = col_mesh,
                   layer.name = "road-network", lwd = 1.5, alpha.regions = 1) +
      mapview(all_lines, color = all_lines$color, alpha = 1, lwd = 3,
              layer.name = "road-estimates", map.type = map.type,
              legend = FALSE, alpha.regions = 1)
      # mapview(tmp, zcol = "coeff", alpha = 1, lwd = 3,
      #         layer.name = "road-estimates", map.type = map.type,
      #         legend = FALSE, color = jet.col(100), alpha.regions = 1)
      
    if(!is.null(bdy)){
      map = map + mapview(st_sf(data.frame(x = 1), geoemtry = bdy), legend = FALSE,
                          col.region = col_bdy, alpha.regions = 0.3, layer.name = "boundary", lwd = 0.3)
    }
    
    # tmp = as_Spatial(st_as_sf(mesh.eval_sfnetwork, "edges"))
    # 
    # num_edges = nrow(tmp)
    # 
    # coeff = matrix(nrow = num_edges, ncol = 1)
    # for(e in 1:num_edges){
    #   coeff[e]= (df[mesh.eval$edges[e,1],time_index] + df[mesh.eval$edges[e,2],time_index]) / 2
    # }
    # 
    # coeff[which(coeff > quantile(coeff, prob = 0.95))] = quantile(coeff, prob = 0.975)
    # 
    # tmp = st_as_sf(tmp) %>% mutate(coeff = as.vector(coeff))
    #
    # map = mapview(mesh.eval_sf, legend = FALSE, map.type = map.type, color = col_mesh,
    #                layer.name = "road-network", lwd = 1.5, alpha.regions = 1) +
    #   mapview(st_sf(data.frame(x = 1), geoemtry = bdy), legend = FALSE,
    #           col.region = col_bdy, alpha.regions = 0.3, layer.name = "boundary", lwd = 0.3) +
    #   mapview(tmp, zcol = "coeff", alpha = 1, lwd = 4,
    #           layer.name = "road-estimates", map.type = map.type,
    #           legend = FALSE, color = jet.col(100), alpha.regions = 1)
    
    if (!is.null(filename) & !is.null(foldername)){

      if(!is.null(center)){
        map@map = map@map %>% setView(center[1], center[2], zoom = zoom)  
      }
      html_filename = paste0(foldername, "/html/", filename, ".html")
      pdf_filename = paste0(foldername, "/pdf/", filename, ".pdf")
      mapshot2(map, url = html_filename, file = pdf_filename, delay = 10)
      
    } else {
      
      map
      
    }
    
  }
  
}

### Function to interpolate colors ---------------------------------------------
interpolate.color <- function(z_start, z_end, m, M, col, num_points = 11) {
  
  ncolor = length(col)
  diffrange = M - m
  
  col_start = (z_start - m)/diffrange*(ncolor-1)+1
  col_start = col[col_start]
  if(is.na(col_start)){
    col_start = col[1]
  }
  
  col_end = (z_end - m)/diffrange*(ncolor-1)+1
  col_end = col[col_end]
  if(is.na(col_end)){
    col_end = col[length(col)]
  }
  
  col_interp = colorRamp(c(col_start, col_end))
  z_interp = seq(z_start, z_end, length.out = num_points)

  if(col_start != col_end){
    colors = rgb(col_interp((z_interp - min(z_interp)) / diff(range(z_interp))), maxColorValue = 255)
  } else {
    colors = rep(col_start, num_points)
  }
  
  return(colors)
}

### Plot the trajectories followed by 4 sources of intensity -------------------
plot.trajectories <- function(mesh, sources, colors = NULL, linewidth = NULL,...){
  
  if(is.null(colors)) colors = brewer.pal(length(sources), "Spectral")[1:length(sources)]
  if(is.null(linewidth)) linewidth = 1
  linewidth = c(linewidth, 1.25)
  
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
  
  data_plot = data.frame(x=x, y=y, grp.nodes)
  
  margin_size = 0.1
  plot_height = 2
  
  palette = c("black", colors)
  
  start = vector(mode="integer", length=length(s))
  end = vector(mode="integer", length=length(s))
  labels = LETTERS[1:(2*length(sources))]
  for(s in 1:length(sources)){
    start[s] = sources[[s]][1]
    end[s] = sources[[s]][length(sources[[s]])]
  }
  
  ggplot(data_plot) +
    geom_link2(aes(x = x, y = y, group = grp.nodes), color = palette[col.nodes+1],
               lineend = 'square', n = 1, linewidth = linewidth[lw.nodes])+#, ...) +
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

### Plot the data over time at one specific neighborhood -----------------------
plot_neighborhood <- function(df, mean_sol_STDEPDE, mean_sol_STLNPP, FEMbasis, center, filename){
  
  STDEPDE_estimates = vector(mode = "numeric", length = ncol(mean_sol_STDEPDE))
  STLNPP_estimates = vector(mode = "numeric", length = ncol(mean_sol_STLNPP))
  
  # h = 1/8
  # t = seq(from = 0+h/2, to = 1-h/2, by = h)
  
  t = seq(from = 0, to = 24, by = 1)/24
  
  for(j in 1:ncol(mean_sol_STDEPDE)){
    STDEPDE_estimates[j] = eval.FEM(FEM(mean_sol_STDEPDE[,j], FEMbasis), locations = center)
    STLNPP_estimates[j] = eval.FEM(FEM(mean_sol_STLNPP[,j], FEMbasis), locations = center)
  }
  
  R = 6371000
  
  df_ = df
  df_$delta_lon <- (df_$lon - center$lon) * (pi / 180) * cos(center$lat * pi / 180) * R
  df_$delta_lat <- (df$lat - center$lat) * (pi / 180) * R
  df_$distance <- sqrt(df_$delta_lon^2 + df_$delta_lat^2)
  
  df_ <- df[df_$distance <= 500, ]
  
  blue_polimi <- rgb(114, 143, 165, max = 255, alpha = 65, names = "blue_polimi25") # #728FA5
  
  x11(width = 8, height = 8)
  #pdf(filename, family = "serif", width = 11.75, height = 5.25)
  par(mai = c(1,1,1,0.5))
  
  plot(x = df_$times, y = rep(0, nrow(df_)), xlim = c(0,1), ylim = c(0, max(STDEPDE_estimates, STLNPP_estimates, na.rm = TRUE)),
       pch = "|", lwd = 2, col = "red3", xaxt = 'n', xlab = "Hour", ylab = "Intensity estimate",
       cex = 2, cex.lab = 2, cex.axis = 2, cex.main = 2)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = blue_polimi)
  grid(nx = NULL, ny = NULL, col = "white", lty = 2, lwd = 1)
  par(new = TRUE)
  plot(x = df_$times, y = rep(0, nrow(df_)), xlim = c(0,1), ylim = c(0, max(STDEPDE_estimates, STLNPP_estimates, na.rm = TRUE)),
       pch = "|", lwd = 2, col = "red3", xaxt = 'n', xlab = "Hour", ylab = "Intensity estimate",
       cex = 2, cex.lab = 2, cex.axis = 2, cex.main = 2, main = "Workdays")
  
  points(x = t, y = STDEPDE_estimates, type='l', lwd=2, col=brewer.pal(6, "YlGnBu")[1])
  points(x = t, y = STLNPP_estimates, type='l', lwd=2, col=brewer.pal(6, "YlGnBu")[4])
  
  axis(side = 1, at = seq(from = 0, to = 1, by = 1/8),
       labels = paste0(as.character(24*seq(from = 0, to = 1, by = 1/8)), rep(":00", 9)),
       cex.axis = 2)
  
  #dev.off()
}
