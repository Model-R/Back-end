ensemble <- function(sp,
                     input.folder1="models",##onde estao os modelos em geral ("models")
                     input.folder2="presfinal", ##onde estao os modelos finais ("presfinal")
                     output.folder="ensemble",
                     occs=spp.filt,
                     which.models = c("Final.bin.mean3","Final.mean.bin7"),
					 consensus = F,
                     consensus.level=0.5#cuanto de los modelos sea retenido: 0.5majority
                     ){

    ##pasta de output
    if (file.exists(paste0("./",input.folder1,"/",sp,"/",output.folder,"/"))==FALSE){
    dir.create(paste0("./",input.folder1,"/",sp,"/",output.folder,"/"))
     }#
    library(raster)
    library(scales)
    library(maps)

    ## para cada tipo de modelo
    for (whi in which.models){

    cat(paste(whi,"-",sp,"\n"))        #lê os arquivos
        tif.files <- list.files(paste0("./",input.folder1,"/",sp,"/",input.folder2),full.names=T,pattern=paste0(whi,'.*tif$'))

        if(length(tif.files)==0) {
            cat(paste("No models to ensemble from for",sp,"\n"))
        } else {
            cat(paste(length(tif.files),"models to ensemble from for",sp,"\n"))
            mod2 <- raster::stack(tif.files)
            if(length(tif.files)==1){
                ensemble.m <- mod2
            } else {
        #plot(mod2)
        ensemble.m <- mean(mod2)
        ensemble.sd <- raster::overlay(mod2,fun=function(x){return(sd(x,na.rm=T))})
            }
        coord <- occs[occs$sp==sp,c('lon','lat')]

        #par(mar=c(5,4,1,1))
        #plot(ensemble.m,main=paste(sp,whi,"ensemble"),font.main=3)
        #points(coord,pch=19,col=alpha("grey60",0.6))
        #map('world',c('',"South America"),add=T,col="grey")

        png(filename=paste0("./",input.folder1,"/",sp,"/",output.folder,"/",sp,"_",whi,"_ensemble.png"),res=300,width=410*300/72,height=480*300/72)
        par(mfrow=c(1,1),mar=c(3,4,4,0))
        plot(ensemble.m,main=paste(sp,whi),legend=F,
        cex.main=1,font.main=3)
        map('world',c('',"South America"),add=T,col="grey")
        points(coord,pch=21,cex=0.6,bg=scales::alpha("cyan",0.6))
        dev.off()

        # o ensemble cru
        writeRaster(ensemble.m,filename=paste0("./",input.folder1,"/",sp,"/",output.folder,"/",sp,"_",whi,"_ensemble.tif"),overwrite=T)

        ####Consensus models
       if (consensus == TRUE){
         ensemble.consensus <- ensemble.m >= consensus.level
        writeRaster(ensemble.consensus,filename=paste0("./",input.folder1,"/",sp,"/",output.folder,"/",sp,"_",whi,"_ensemble",consensus.level*100,".tif"),overwrite=T)


        png(filename=paste0("./",input.folder1,"/",sp,"/",output.folder,"/",sp,"_",whi,"_ensemble",consensus.level*100,".png"),res=300,width=410*300/72,height=480*300/72)
        par(mfrow=c(1,1),mar=c(3,4,4,0))
        plot(ensemble.consensus,main=paste(whi,consensus.level*100),legend=F,cex.main=1,font.main=3)
        map('world',c('',"South America"),add=T,col="grey")
        points(coord,pch=19,cex=0.3,col=scales::alpha("cyan",0.6))
        dev.off()
		}
    }
}
    #return(ensemble.m)
}
