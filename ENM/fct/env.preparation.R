#Reading aspect
library(raster)
aspect1k <- raster(list.files("./1km", full.names=TRUE, pattern='aspect.asc$'))
recs<-which(values(aspect1k==-1))
values(aspect1k)[recs]<-trunc(runif(length(recs),0,36000))
aspect1k<-aspect1k*pi/18000
#plot(aspect1k)

writeRaster(aspect1k, filename = './1km/h_aspect.asc',format='ascii',overwrite=TRUE)
sin.aspect1k<-sin(aspect1k)
writeRaster(sin.aspect1k, filename = './1km/sin.aspect.asc',format='ascii',overwrite=TRUE)
cos.aspect1k<-cos(aspect1k)
writeRaster(cos.aspect1k, filename = './1km/cos.aspect.asc',format='ascii',overwrite=TRUE)
# Variables not used:
#Alt;
#h.aspect;
#h.flowacc;
#h.flowdir
