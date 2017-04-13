# coordinates: coordenadas para a espécie
# coordinates <- occs[occs$sp == sp, c("lon", "lat")]

# partitions: quantidade de partições
do_maxent <- function(sp,
		      coordinates,
		      partitions,
		      buffer = FALSE,
		      seed = 512,
		      predictors = predictors,
		      models.dir = model.dir,
		      project.model = project.model,
		      projections = projections,
		      mask = mask,
		      n.back = 500) {
  cat(paste("Maxent", "\n"))

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
  
  if (buffer == FALSE) {
    set.seed(seed + 2)
    backgr <- randomPoints(predictors, n.back)
  } else {
    source("../fct/createBuffer.R")
    backgr <- createBuffer(coord_ = coordinates, n.back_ = n.back, buffer.type_ = buffer.type, 
      occs_ = occs, sp_ = sp, seed_ = seed)
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
    
    mx <- maxent(predictors, pres_train)
    emx <- evaluate(pres_test, backg_test, mx, predictors)
    thresholdmx <- emx@t[which.max(emx@TPR + emx@TNR)]
    thmx <- threshold(emx)
    mx_TSS <- max(emx@TPR + emx@TNR) - 1
    mx_cont <- predict(predictors, mx, progress = "text")
    mx_bin <- mx_cont > thresholdmx
    mx_cut <- mx_cont * mx_bin
    thmx$AUC <- emx@auc
    thmx$TSS <- mx_TSS
    thmx$algoritmo <- "maxent"
    thmx$partition <- i
    row.names(thmx) <- paste(sp, i, "maxent")

    write.table(thmx, file = paste0(models.dir, "/", sp, "/evaluate", 
      sp, "_", i, "_maxent.txt"))
    
    if (class(mask) == "SpatialPolygonsDataFrame") {
      source("../fct/cropModel.R")
      mx_cont <- cropModel(mx_cont, mask)
      mx_bin <- cropModel(mx_bin, mask)
      mx_cut <- cropModel(mx_cut, mask)
    }
    writeRaster(x = mx_cont, filename = paste0(models.dir, "/", sp, "/maxent_cont_", 
      sp, "_", i, ".tif"), overwrite = T)
    writeRaster(x = mx_bin, filename = paste0(models.dir, "/", sp, "/maxent_bin_", 
      sp, "_", i, ".tif"), overwrite = T)
    writeRaster(x = mx_cut, filename = paste0(models.dir, "/", sp, "/maxent_cut_", 
      sp, "_", i, ".tif"), overwrite = T)
    
    png(filename = paste0(models.dir, "/", sp, "/maxent", sp, "_", i, "%03d.png"))
    plot(mx_cont, main = paste("maxent raw", "\n", "AUC =", round(emx@auc, 2), "-", 
      "TSS =", round(mx_TSS, 2)))
    plot(mx_bin, main = paste("maxent P/A", "\n", "AUC =", round(emx@auc, 2), "-", 
      "TSS =", round(mx_TSS, 2)))
    plot(mx_cut, main = paste("maxent cut", "\n", "AUC =", round(emx@auc, 2), "-", 
      "TSS =", round(mx_TSS, 2)))
    dev.off()
  
  
    if (project.model == T) {
      for (proj in projections) {
        data <- list.files(paste0("./env/", proj), pattern = proj)
        data2 <- stack(data)
        mx_proj <- predict(data2, mx, progress = "text")
        mx_proj_bin <- mx_proj > thresholdmx
        mx_proj_cut <- mx_proj_bin * mx_proj
        # Normaliza o modelo cut do_proj_cut <- do_proj_cut/maxValue(do_proj_cut)
        if (class(mask) == "SpatialPolygonsDataFrame") {
          source("../fct/cropModel.R")
          mx_proj <- cropModel(mx_proj, mask)
          mx_proj_bin <- cropModel(mx_proj_bin, mask)
          mx_proj_cut <- cropModel(mx_proj_cut, mask)
        }
        writeRaster(x = mx_proj, filename = paste0(models.dir, "/", sp, "/", 
          proj, "/maxent_cont_", sp, "_", i, ".tif"), overwrite = T)
        writeRaster(x = mx_proj_bin, filename = paste0(models.dir, "/", sp, "/", 
          proj, "/maxent_bin_", sp, "_", i, ".tif"), overwrite = T)
        writeRaster(x = mx_proj_cut, filename = paste0(models.dir, "/", sp, "/", 
          proj, "/maxent_cut_", sp, "_", i, ".tif"), overwrite = T)
        rm(data2)
      }
    }
  }
  return(thmx)
}
#    eval_df <- data.frame(kappa = 1, spec_sens = 1, no_omission = 1, prevalence = 1, 
#      equal_sens_spec = 1, sensitivity = 1, AUC = 1, TSS = 1, algoritmo = "foo", 
#      partition = 1)
#      eval_df <- rbind(eval_df, thdo)
