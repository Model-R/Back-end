write_Raster <- function(x,mask=mask){
  if (class(mask) == "RasterLayer"){
    x <- mask(x , mask)
    x <- crop(x , mask)
  writeRaster(x,filename=paste0("./",output.folder,"/",sp,"/",proj,"/",algo,"_cont_",sp,"_",i,".tif"),overwrite=T,datatype="INT1U")
  }
}