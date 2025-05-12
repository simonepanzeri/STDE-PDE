## FUNCION TO GENERATE THE ROAD NETWORK OVER DOMAIN OF INTEREST ----------------
generateDomain <- function(district = "Rome"){
  
  foldername <- paste0(getwd(), "/utils/domain/", district, "/")
    
  if(!file.exists(paste0(foldername, district, "_osm.RData"))){
    
    if(!dir.exists(foldername)){
      dir.create(foldername, showWarnings = FALSE)
    }
    
    if(district == "Rome"){
      
      bounds <- c(12.432014, 41.849761, 12.548286, 41.943804)
      
    } else if(district == "Bergamo") {
      
      bounds <- c(9.626203, 45.672331, 9.709679, 45.722368)
      
    } else if(district == "Milan") {
      
      bounds <- c(9.041083, 9.277905, 45.387103, 45.535973)
      
    }
    
    xbounds <- c(bounds[1], bounds[3])
    ybounds <- c(bounds[2], bounds[4])
    
    its_bbox <- st_bbox(c(xmin = xbounds[1], ymin = ybounds[1],
                          xmax = xbounds[2], ymax = ybounds[2]), crs = 4326) %>% st_as_sfc()
    
    # Complete Network
    # osm_data <- oe_get(district, boundary = its_bbox, download_directory = foldername)
    
    # unique(osm_data[[3]])
    # [1] "motorway_link"  "trunk_link"     "trunk"          "primary"        "secondary"      "tertiary"       "unclassified"   "living_street" 
    # [9] "residential"    "pedestrian"     "tertiary_link"  "cycleway"       "service"        "secondary_link" "primary_link"   "steps"         
    # [17] "footway"        NA               "track"          "path"           "services"       "construction"   "rest_area"      "elevator"      
    # [25] "proposed"
    
    # mapview(osm_data[which(osm_data$highway == "unclassified"),]) +
    #   mapview(st_sf(data.frame(x = 1), geoemtry = adm_bdy_Bergamo), legend = FALSE,
    #           alpha.regions = 0.3, layer.name = "boundary", lwd = 0.1) 
    
    # Simplified Network
    my_vectortranslate <- c(
      
      "-select", "highway",
      
      "-where", "highway IN ('motorway_link', 'trunk_link', 'trunk', 'primary',
                             'secondary', 'tertiary', 'unclassified', 'residential',
                             'tertiary_link', 'secondary_link', 'primary_link',
                             'track')"
      
    )
    
    osm_data <- oe_get(district, vectortranslate_options = my_vectortranslate,
                       boundary = its_bbox, download_directory = foldername)

    mapview(osm_data)
    save(osm_data, bounds, file = paste0(foldername, district, "_osm.RData"))
  
  } else {
    
    load(paste0(foldername, district, "_osm.RData"))
    
  }
  
  return (list(osm_data, bounds))
  
}

## FUNCION TO GENERATE THE MESH OVER THE DOMAIN OF INTEREST --------------------
generateMesh <- function(district = "Rome"){
  
  foldername <- paste0(getwd(), "/utils/domain/", district, "/")
  
  if(!file.exists(paste0(foldername, district, "_mesh.RData"))){
    
    domain <- generateDomain(district = district)
    osm_data <- domain[[1]]
    xbounds <- domain[[2]][c(1,3)]
    ybounds <- domain[[2]][c(2,4)]
    
    adm_bdy_IT <- st_read("utils/map/adm_bdy_IT/Com01012023/Com01012023_WGS84.shx")
    adm_bdy_ <- adm_bdy_IT$geometry[which(adm_bdy_IT$COMUNE == district)]
    adm_bdy_ <- st_transform(adm_bdy_, crs = 4326)
    
    filter <- st_within(osm_data, adm_bdy_)
    osm_data <- osm_data[which(lengths(filter) != 0), ]
  
    sf_network <- as_sfnetwork(osm_data, directed = FALSE)
    
    # Simplify the sfnetwork
    sf_network <- sf_network %>%
      activate("edges") %>%
      filter(!edge_is_multiple()) %>%
      filter(!edge_is_loop())
    
    # Clean the sfnetwork
    sf_network <- sf_network %>% 
      convert(to_spatial_subdivision, .clean = TRUE)
    
    # Select a Full Connected Graph
    sf_network <- sf_network %>% 
      convert(to_components, .clean = TRUE, .select = 1L)

    filtered <- sf_network %>%
      activate("edges") %>%
      filter(edge_intersects(st_cast(adm_bdy_, to = "POLYGON"))) %>%
      activate("nodes") %>%
      filter(!node_is_isolated())
    
    mesh <- as.mesh.1.5D(filtered)
    
    save(mesh, file = paste0(foldername, district, "_mesh.RData"))
  
  } else {
    
    load(paste0(foldername, district, "_mesh.RData"))
    
  }
  
  return (mesh)
}

generateSubMesh <- function(mesh_sfnetwork,
                            xmin = NULL, ymin = NULL, xmax = NULL, ymax =  NULL){
  
  foldername <- paste0(getwd(), "/utils/domain/", district, "/")
  
  if(!file.exists(paste0(foldername, district, "_submesh.RData"))){
    
    if(is.null(xmin) | is.null(ymin) | is.null(xmax) | is.null(ymax)){
      xylim <- data.frame(xmin = 9.649472, ymin = 45.685616,
                          xmax = 9.680371, ymax = 45.705100)
    } else {
      xylim <- data.frame(xmin = xmin, ymin = ymin,
                          xmax = xmax, ymax = ymax)
    }
    
    p1 <- st_point(c(xylim$xmin, xylim$ymin))
    p2 <- st_point(c(xylim$xmax, xylim$ymin))
    p3 <- st_point(c(xylim$xmax, xylim$ymax))
    p4 <- st_point(c(xylim$xmin, xylim$ymax))

    poly <- st_multipoint(c(p1, p2, p3, p4)) %>%
      st_cast("POLYGON") %>%
      st_sfc(crs = 4326)

    filtered <- mesh_sfnetwork %>%
      activate("edges") %>%
      filter(edge_intersects(poly)) %>%
      activate("nodes") %>%
      filter(!node_is_isolated())
    
    mesh <- as.mesh.1.5D(filtered)
    
    save(mesh, file = paste0(foldername, district, "_submesh.RData"))
  
  } else {
    
    load(paste0(foldername, district, "_submesh.RData"))
    
  }
  
  return (mesh)
}

# Subnetwork Region
neighborhood_ <- data.frame(xmin = 9.6465, ymin = 45.6715, xmax = 9.7035, ymax = 45.72)

# Rondò delle Valli
neighborhood_1 <- data.frame(lon = 9.69469, lat = 45.70444)
neighborhood_1_box <- data.frame(xmin = 9.68702, ymin = 45.69725, xmax = 9.70143, ymax = 45.71305)

# Bergamo Centro (Città Bassa)
neighborhood_2 <- data.frame(lon = 9.67101, lat = 45.69301)
neighborhood_2_box <- data.frame(xmin = 9.65257, ymin = 45.68385, xmax = 9.67864, ymax = 45.69773)

# Bergamo Centro (Città Alta)
neighborhood_3 <- data.frame(lon = 9.66531, lat = 45.70454)
neighborhood_3_box <- data.frame(xmin = 9.65510, ymin = 45.69944, xmax = 9.67479, ymax = 45.71274)

# Raccordo Autostradale
neighborhood_4 <- data.frame(lon = 9.66949, lat = 45.67778)
neighborhood_4_box <- data.frame(xmin = 9.65660, ymin = 45.67213, xmax = 9.68704, ymax = 45.68228)
