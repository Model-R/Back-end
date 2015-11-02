

ensemble <- function(sp,
                     input.folder1,##onde estao os modelos em geral ("models")
                     input.folder2, ##onde estao os modelos finais ("presfinal")
                     occs=spp.filt,
                     which.models=c("Final.bin.mean3","Final.mean.bin7"),
                     output.folder="ensemble"){
    
    ##aqui muda sempre entre projetos e provavelmente neste também
    #if (file.exists(paste0("./",output.folder))==FALSE) dir.create(paste0("./",output.folder)) ###tirei este aqui, para os ensemble ficarem dentro das pastas de espécies
    dir.create(paste0("./",input.folder1,"/",sp,"/",output.folder,"/")) # adicionado por Guilherme -- ok!
    
    library(raster)
## para cada tipo de modelo    
    for (whi in unique(which.models)){
        cat(paste(whi,"\n"))
        #lê os arquivos -- dependendo da esturua de pastas vai mudar:
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
        points(coord,pch=21,bg="grey70")
        
        
        # o ensemble cru
        writeRaster(ensemble.m,filename=paste0("./",input.folder1,"/",sp,"/",output.folder,"/",sp,"_",whi,"_ensemble.tif"),overwrite=T)
        ####aqui filtra 30% dos modelos e cria um binário de novo
        writeRaster(ensemble.m>=0.3,filename=paste0("./",input.folder1,"/",sp,"/",output.folder,"/",sp,"_",whi,"_ensemble30.tif"),overwrite=T)
        
        
        library(maps) 
        png(filename=paste0("./",input.folder1,"/",sp,"/",output.folder,"/",sp,"_",whi,"_ensemble.png"),res=300,width=410*300/72,height=480*300/72)
        par(mfrow=c(1,1),mar=c(3,4,4,0))
        plot(ensemble.m,main=paste(sp,whi),legend=F,
             #xlim=c(-57,-35),ylim=c(-31,-3),
             cex.main=2,font.main=3)
        map('world',c('',"South America"),
            #xlim=c(-58,-35),ylim=c(-33,-5),
            add=T)
        points(coord,pch=19,cex=0.8)
        dev.off()
        
    }

    return(ensemble.m)    
}