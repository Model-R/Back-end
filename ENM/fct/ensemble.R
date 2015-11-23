

ensemble <- function(sp,
                     input.folder1="models",##onde estao os modelos em geral ("models")
                     input.folder2="presfinal", ##onde estao os modelos finais ("presfinal")
                     occs=spp.filt,
                     which.models=c("Final.bin.mean3","Final.mean.bin7"),
                     output.folder="ensemble"){
    
    ##pasta de output
    if (file.exists(paste0("./",input.folder1,"/",sp,"/",output.folder,"/"))==FALSE) dir.create(paste0("./",input.folder1,"/",sp,"/",output.folder,"/")) # 
    library(raster)
    library(scales)
    library(maps) 

    ## para cada tipo de modelo    
    for (whi in unique(which.models)){
        cat(paste(whi,"\n"))
        #lê os arquivos 
        files <- list.files(paste0("./",input.folder1,"/",sp,"/",input.folder2),full.names=T,pattern=paste0(whi,sp))
        
        tiffiles<-files[grep(".tif",files,fixed=T)]##aqui eu pego só os tif
        
        if(length(tiffiles)==0) cat(paste("No models to ensemble from for",sp,"\n"))
        if(length(tiffiles)!=0) {
            cat(paste("ensemble:",tiffiles,"/n"))
            mod2<-stack(tiffiles)
        }
        plot(mod2)
        
        ensemble.m <- mean(mod2)
        coord <- occs[occs$sp==sp,c('lon','lat')]
        
        par(mar=c(5,4,1,1))
        plot(ensemble.m,main=paste(sp,whi,"ensemble"),font.main=3)
        points(coord,pch=19,col=alpha("grey60",0.6))
        map('world',c('',"South America"),add=T,col="grey")
        
        
       
        
        
        png(filename=paste0("./",input.folder1,"/",sp,"/",output.folder,"/",sp,"_",whi,"_ensemble.png"),res=300,width=410*300/72,height=480*300/72)
        par(mfrow=c(1,1),mar=c(3,4,4,0))
        plot(ensemble.m,main=paste(sp,whi),legend=F,
        cex.main=1,font.main=3)
        map('world',c('',"South America"),add=T,col="grey")
        points(coord,pch=19,cex=0.8,col=alpha("grey60",0.6))
        dev.off()
        
        # o ensemble cru
        writeRaster(ensemble.m,filename=paste0("./",input.folder1,"/",sp,"/",output.folder,"/",sp,"_",whi,"_ensemble.tif"),overwrite=T)
        ####aqui filtra 30% dos modelos e cria um binário de novo
        writeRaster(ensemble.m>=0.3,filename=paste0("./",input.folder1,"/",sp,"/",output.folder,"/",sp,"_",whi,"_ensemble30.tif"),overwrite=T)
    }

    return(ensemble.m)    
}