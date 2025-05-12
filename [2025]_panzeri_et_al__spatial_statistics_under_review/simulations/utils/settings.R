setting <- function(network = "" ){
  
  if(network == "ontario"){
  
    data("ORN", package = "shp2graph")
  
    mesh = as.mesh.1.5D(ORN.nt)
    mesh = normalize_mesh_unit(mesh)$mesh
    mesh = refine.mesh.1.5D(mesh, delta=0.0125)
   
    FEMbasis = create.FEM.basis(mesh)
    
    nnodes = nrow(mesh$nodes)
    spatstat.linnet = as.linnet(mesh)
    
    res = list(mesh = mesh, FEMbasis=FEMbasis, nnodes = nnodes,
               spatstat.linnet = spatstat.linnet)
    return(res)
    
  }else if( network == "estevan"){
      
    data("ERN_OSM_correct", package = "shp2graph")
      
    mesh = as.mesh.1.5D(ERN_OSM_cor.nt)
    mesh = normalize_mesh_unit(mesh)$mesh
    mesh = refine.mesh.1.5D(mesh, delta=0.0125)
        
    FEMbasis = create.FEM.basis(mesh)
      
    nnodes = nrow(mesh$nodes)
    spatstat.linnet = as.linnet(mesh)
      
    res = list(mesh = mesh, FEMbasis=FEMbasis, nnodes = nnodes,
               spatstat.linnet = spatstat.linnet)
  
  }else if(network == "london"){
      
    data("LNNT", package = "shp2graph")
      
    network <- st_as_sf(LN.nt) 
      
    # setting the coordinate reference system
    st_crs(network) <- 27700
      
    # lat/lon reference system 
    network <- st_transform(network, 4326)
      
    # building sfnetwork
    sfnetwork <- as_sfnetwork(network, directed = FALSE, edges_as_lines = TRUE)
    
    # simplifying 
    sfnetwork <- sfnetwork %>%
      activate("edges") %>%
      filter(!edge_is_multiple()) %>%
      filter(!edge_is_loop())
      
    # cleaning
    sfnetwork <- sfnetwork %>% 
      convert(to_spatial_subdivision, .clean = TRUE)
      
    # selecting full connected graph
    sfnetwork <- sfnetwork %>% 
      convert(to_components, .clean = TRUE, .select = 1L)
      
    mesh <- as.mesh.1.5D(sfnetwork)
    FEMbasis <- create.FEM.basis(mesh)
    nnodes = nrow(mesh$nodes)
    spatstat.linnet = as.linnet(mesh)
      
    res = list(mesh = mesh, FEMbasis=FEMbasis, nnodes = nnodes,
               spatstat.linnet = spatstat.linnet)
    
  }else if(network == "munster"){
    
    data("roxel", package = "sfnetworks")
    
    mesh = as.mesh.1.5D(as_sfnetwork(roxel$geometry))
    
    FEMbasis = create.FEM.basis(mesh)
    
    nnodes = nrow(mesh$nodes)
    spatstat.linnet = as.linnet(mesh)
    
    res = list(mesh = mesh, FEMbasis=FEMbasis, nnodes = nnodes,
               spatstat.linnet = spatstat.linnet)
    
  }else if(network == "chicago"){
    
    data("chicagonet", package = "stopp")
    
    mesh = as.mesh.1.5D(chicagonet)
    mesh = normalize_mesh_unit(mesh)$mesh
    
    FEMbasis = create.FEM.basis(mesh)
    
    nnodes = nrow(mesh$nodes)
    spatstat.linnet = as.linnet(mesh)
    
    res = list(mesh = mesh, FEMbasis=FEMbasis, nnodes = nnodes,
               spatstat.linnet = spatstat.linnet)
  
  }else if(network == "valencia"){
    
    data("valencianet", package = "stopp")
    
    mesh = as.mesh.1.5D(valencianet)
    mesh = normalize_mesh_unit(mesh)$mesh
    
    FEMbasis = create.FEM.basis(mesh)
    
    nnodes = nrow(mesh$nodes)
    spatstat.linnet = as.linnet(mesh)
    
    res = list(mesh = mesh, FEMbasis=FEMbasis, nnodes = nnodes,
               spatstat.linnet = spatstat.linnet)
    
  }else if(network == "eastbourne"){
    
    data("Eastbourne", package = "stlnpp")
    
    # building sfnetwork
    sfnetwork <- as_sfnetwork(Eastbourne$domain, directed = FALSE, edges_as_lines = TRUE)
    
    # simplifying 
    sfnetwork <- sfnetwork %>%
      activate("edges") %>%
      filter(!edge_is_multiple()) %>%
      filter(!edge_is_loop())
    
    # cleaning
    sfnetwork <- sfnetwork %>% 
      convert(to_spatial_subdivision, .clean = TRUE)
    
    # selecting full connected graph
    sfnetwork <- sfnetwork %>% 
      convert(to_components, .clean = TRUE, .select = 1L)
    
    mesh = as.mesh.1.5D(sfnetwork)
    mesh = normalize_mesh_unit(mesh)$mesh
    
    FEMbasis = create.FEM.basis(mesh)
    
    nnodes = nrow(mesh$nodes)
    spatstat.linnet = as.linnet(mesh)
    
    res = list(mesh = mesh, FEMbasis=FEMbasis, nnodes = nnodes,
               spatstat.linnet = spatstat.linnet)
    
  }else if(network == "medellin"){
    
    data("Medellin", package = "stlnpp")
    
    mesh = as.mesh.1.5D(Medellin$domain)
    mesh = normalize_mesh_unit(mesh)$mesh
    
    FEMbasis = create.FEM.basis(mesh)
    
    nnodes = nrow(mesh$nodes)
    spatstat.linnet = as.linnet(mesh)
    
    res = list(mesh = mesh, FEMbasis=FEMbasis, nnodes = nnodes,
               spatstat.linnet = spatstat.linnet)
    
  }else if(network == "easynet"){
    
    data("easynet", package = "stlnpp")
    
    mesh = as.mesh.1.5D(easynet)
    mesh = normalize_mesh_unit(mesh)$mesh
    
    FEMbasis = create.FEM.basis(mesh)
    
    nnodes = nrow(mesh$nodes)
    spatstat.linnet = as.linnet(mesh)
    
    res = list(mesh = mesh, FEMbasis=FEMbasis, nnodes = nnodes,
               spatstat.linnet = spatstat.linnet)
    
  }else if(network == "simplenet"){
    
    data("simplenet", package = "spatstat.data")
    mesh = as.mesh.1.5D(simplenet)
    mesh = refine.mesh.1.5D(mesh,0.025)
      
    FEMbasis = create.FEM.basis(mesh)
    nnodes = nrow(mesh$nodes)
    spatstat.linnet = as.linnet(mesh)
      
    res = list(mesh=mesh, FEMbasis = FEMbasis, nnodes=nnodes,
               spatstat.linnet = spatstat.linnet)  
  }
}
