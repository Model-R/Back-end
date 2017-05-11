# coordinates: coordenadas para a espécie
# coordinates <- occs[occs$sp == sp, c("lon", "lat")]

# partitions: quantidade de partições
do_bioclim <- function(sp,
		       coordinates,
		       partitions,
		       buffer = FALSE,
		       seed = 512,
		       predictors,
		       models.dir,
		       project.model,
		       projections,
		       mask,
		       n.back = 500) {
  cat(paste("Bioclim", "\n"))

  if (file.exists(paste0(models.dir)) == FALSE)
    dir.create(paste0(models.dir))
  if (file.exists(paste0(models.dir, "/", sp)) == FALSE) 
    dir.create(paste0(models.dir, "/", sp))
  if (project.model == T) {
    for (proj in projections) {
      if (file.exists(paste0(models.dir, "/", sp, "/", proj)) == FALSE) 
        dir.create(paste0(models.dir, "/", sp, "/", proj))
    }
  }

  # tabela de valores
  presvals <- raster::extract(predictors, coordinates)
  
  if (buffer %in% c("mean", "max")) {
    source("../../fct/createBuffer.R")
    backgr <- createBuffer(coord = coordinates, n.back = n.back, buffer.type = buffer, 
      occs = coordinates, sp = sp, seed = seed, predictors = predictors)
  } else {
    set.seed(seed + 2)
    backgr <- randomPoints(predictors, n.back)
  }

  colnames(backgr) <- c("lon", "lat")
  
  # Extraindo dados ambientais dos bckgr
  backvals <- raster::extract(predictors, backgr)
  pa <- c(rep(1, nrow(presvals)), rep(0, nrow(backvals)))
  
  # Data partition
  if (nrow(coordinates) < 11) 
    partitions <- nrow(coordinates)
  set.seed(seed)  #reproducibility
  group <- kfold(coordinates, partitions)
  set.seed(seed + 1)
  bg.grp <- kfold(backgr, partitions)
  group.all <- c(group, bg.grp)
  
  pres <- cbind(coordinates, presvals)
  back <- cbind(backgr, backvals)
  rbind_1 <- rbind(pres, back)
  sdmdata <- data.frame(cbind(group.all, pa, rbind_1))
  rm(rbind_1)
  rm(pres)
  rm(back)
  gc()
  write.table(sdmdata, file = paste0(models.dir, "/", sp, "/sdmdata.txt"))
  
#  if (! file.exists(file = paste0(models.dir, "/", sp, "/evaluate", sp, "_", i, ".txt"))) {
#    write.table(data.frame(kappa = numeric(), spec_sens = numeric(), no_omission = numeric(), prevalence = numeric(), 
#			         equal_sens_spec = numeric(), sensitivity = numeric(), AUC = numeric(), TSS = numeric(), algoritmo = character(), 
#				 partition = numeric()), file = paste0(models.dir, "/", sp, "/evaluate", sp, "_", i, ".txt"))
#  }
  
  ##### Hace los modelos
  for (i in unique(group)) {
    cat(paste(sp, "partition number", i, "\n"))
    pres_train <- coordinates[group != i, ]
    if (nrow(coordinates) == 1) 
      pres_train <- coordinates[group == i, ]
    pres_test <- coordinates[group == i, ]
    
    backg_test <- backgr[bg.grp == i, ]  #new
    
    bc <- bioclim(predictors, pres_train)
    ebc <- evaluate(pres_test, backg_test, bc, predictors)
    thresholdbc <- ebc@t[which.max(ebc@TPR + ebc@TNR)]
    thbc <- threshold(ebc)
    bc_TSS <- max(ebc@TPR + ebc@TNR) - 1
    bc_cont <- predict(predictors, bc, progress = "text")
    bc_bin <- bc_cont > thresholdbc
    bc_cut <- bc_cont * bc_bin
    thbc$AUC <- ebc@auc
    thbc$TSS <- bc_TSS
    thbc$algoritmo <- "BioClim"
    thbc$partition <- i
    row.names(thbc) <- paste(sp, i, "BioClim")

    write.table(thbc, file = paste0(models.dir, "/", sp, "/evaluate", 
      sp, "_", i, "_bioclim.txt"))
    
    if (class(mask) == "SpatialPolygonsDataFrame") {
    source("../../fct/cropModel.R")
    bc_cont <- cropModel(bc_cont, mask)
    bc_bin <- cropModel(bc_bin, mask)
    bc_cut <- cropModel(bc_cut, mask)
    }
    writeRaster(x = bc_cont, filename = paste0(models.dir, "/", sp, "/BioClim_cont_", 
      sp, "_", i, ".tif"), overwrite = T)
    writeRaster(x = bc_bin, filename = paste0(models.dir, "/", sp, "/BioClim_bin_", 
      sp, "_", i, ".tif"), overwrite = T)
    writeRaster(x = bc_cut, filename = paste0(models.dir, "/", sp, "/BioClim_cut_", 
      sp, "_", i, ".tif"), overwrite = T)
  
    png(filename = paste0(models.dir, "/", sp, "/Bioclim", sp, "_", i, "%03d.png"))
    plot(bc_cont, main = paste("Bioclim raw", "\n", "AUC =", round(ebc@auc, 2), "-", 
      "TSS =", round(bc_TSS, 2)))
    plot(bc_bin, main = paste("Bioclim P/A", "\n", "AUC =", round(ebc@auc, 2), "-", 
      "TSS =", round(bc_TSS, 2)))
    plot(bc_cut, main = paste("Bioclim cut", "\n", "AUC =", round(ebc@auc, 2), "-", 
      "TSS =", round(bc_TSS, 2)))
    dev.off()
  
  
    if (project.model == T) {
      for (proj in projections) {
        data <- list.files(paste0("./env/", proj), pattern = proj)
        data2 <- stack(data)
        bc_proj <- predict(data2, do, progress = "text")
        bc_proj_bin <- bc_proj > thresholdbc
        bc_proj_cut <- bc_proj_bin * bc_proj
        # Normaliza o modelo cut do_proj_cut <- do_proj_cut/maxValue(do_proj_cut)
        if (class(mask) == "SpatialPolygonsDataFrame") {
          source("./fct/cropModel.R")
          bc_proj <- cropModel(bc_proj, mask)
          bc_proj_bin <- cropModel(bc_proj_bin, mask)
          bc_proj_cut <- cropModel(bc_proj_cut, mask)
        }
        writeRaster(x = bc_proj, filename = paste0(models.dir, "/", sp, "/", 
          proj, "/Domain_cont_", sp, "_", i, ".tif"), overwrite = T)
        writeRaster(x = bc_proj_bin, filename = paste0(models.dir, "/", sp, "/", 
          proj, "/Domain_bin_", sp, "_", i, ".tif"), overwrite = T)
        writeRaster(x = bc_proj_cut, filename = paste0(models.dir, "/", sp, "/", 
          proj, "/Domain_cut_", sp, "_", i, ".tif"), overwrite = T)
        rm(data2)
      }
    }
  }
  return(thbc)
}
