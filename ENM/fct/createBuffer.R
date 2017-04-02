createBuffer <- function(coord_ = coord, sp_ = sp, occs_ = occs, seed_ = seed, n.back_ = n.back, 
  buffer.type_ = "mean") {
  # Transformando em spatial points
  coordinates(coord_) <- ~lon + lat
  
  if (buffer.type_ == "mean") 
    dist.buf <- mean(spDists(x = coord_, longlat = FALSE, segments = TRUE))
  if (buffer.type_ == "max") 
    dist.buf <- max(spDists(x = coord_, longlat = FALSE, segments = TRUE))
  
  buffer <- raster::buffer(coord_, width = dist.buf, dissolve = TRUE)
  
  # Transformando coords de novo em matriz para rodar resto script
  coord_ <- occs_[occs_$sp == sp_, c("lon", "lat")]
  
  # Transformando em spatial polygon data frame
  buffer <- SpatialPolygonsDataFrame(buffer, data = as.data.frame(buffer@plotOrder), 
    match.ID = FALSE)
  crs(buffer) <- crs(predictors)
  
  ######### TENHO CERTEZA DE QUE ISTO PODE FICAR MENOS PESADO Reference raster com mesmo
  ######### extent e resolution que predictors
  r_buffer <- crop(predictors, buffer)
  r_buffer <- mask(r_buffer, buffer)
  
  # r_buffer <- raster(ext=extent(predictors), resolution=res(predictors))
  
  # Rasterizando o buffer p/ geração dos ptos aleatorios r_buffer <-
  # rasterize(buffer, r_buffer, field=buffer@plotOrder) Limitando a mascara
  # ambiental r_buffer <- r_buffer*(predictors[[1]]!=0)
  
  
  # Gerando pontos aleatorios no buffer
  set.seed(seed_ + 2)
  backgr <- randomPoints(r_buffer, n.back_)
  rm(buffer)
  gc()
  return(backgr)
}
