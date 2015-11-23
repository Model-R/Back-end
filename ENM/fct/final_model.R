######### Final modeling: one model per algorithm ----

final.model <- function(sp,
                        select.partitions=T,
                        weight.partitions=F,
                        threshold=c("spec_sens"),
                        TSS.value=0.2,
                        #para weight partitions
                        weight.par=c("TSS","AUC"),
                        input.folder="models",
                        output.folder="presfinal"){
  
  if (file.exists(paste0("./",input.folder,"/",sp,"/",output.folder))==FALSE) dir.create(paste0("./",input.folder,"/",sp,"/",output.folder), recursive = TRUE)
  
  cat(sp)
  library("raster")
  library("data.table")
  cat(paste("Reading the evaluation files","\n"))
  evall<- list.files(path = paste0(input.folder,"/",sp), pattern=paste0("evaluate",sp),full.names = T)
  lista<-list()
  for (i in 1:length(evall)) {
    lista[[i]]<-read.table(file = evall[i],header=T,row.names=1)
  }
  stats<-rbindlist(lista)
  stats<-as.data.frame(stats)
  
  #Extracts only for the selected algorithm
  algoritmos <- unique(stats$algoritmo)
  for (algo in algoritmos){
    stats2 <- stats[stats$algoritmo==algo,]
    if(algo=="Mahal"){
       stats2$spec_sens[stats2$spec_sens<0]<-0
    }
    
    part <- nrow(stats2)#How many partitions were there
    
    cat(paste("Reading models from .tif files","\n"))
    modelos.cont <- list.files(path = paste0(input.folder,"/",sp),full.names=T,pattern=paste0(algo,"_cont_",sp))
    modelos.cut <- list.files(path = paste0(input.folder,"/",sp),full.names=T,pattern=paste0(algo,"_cut_",sp))
    modelos.bin <- list.files(path = paste0(input.folder,"/",sp),full.names=T,pattern=paste0(algo,"_bin_",sp))
    mod.cont<-stack(modelos.cont)#(0)
    mod.cut <- stack(modelos.cut)
    mod.bin <- stack(modelos.bin)
    
    names(mod.cont) <- paste0(sp,algo,"Partition",1:part)
    names(mod.cut) <- names(mod.cont)
    names(mod.bin) <- names(mod.cont)
    
  
    if (select.partitions==T){
      sel.index<- which(stats2[,"TSS"]>=TSS.value)
      cont.sel<- mod.cont[[sel.index]]
      bin.sel<- mod.bin[[sel.index]]#
      cut.sel<- mod.cut[[sel.index]]#
      
      th.mean<-mean(stats2[,names(stats2)==threshold][sel.index])
      
      if (length(sel.index)==0) cat(paste("No partition was selected for",sp,"\n"))
      
      #en caso de que sea solo uno varios modelos son el mismo
      if (length(sel.index)==1){
        cat(paste(length(sel.index), "partitions was selected for",sp))
        #cont.sel#(1)(2)
        #bin.sel#(5)(3)(7) (8)
        #cut.sel#(4)(6)(9)(10)
        
        final <- stack(cont.sel,#[2]
                       bin.sel,#[3],
                       cut.sel,#[4]
                       bin.sel,#[7]
                       bin.sel,#[8]
                       cut.sel,#[9]
                       cut.sel)#[10]
        
      }
      
      #en caso de que sean más aplica el mapa          
      
      if (length(sel.index)>1){
        cat(paste(length(sel.index), "partitions were selected for",sp))
        
        final.cont.mean <- mean(cont.sel)#(2)
        final.bin.mean <- (final.cont.mean>th.mean)#(3)
        final.cut.mean <- final.bin.mean*final.cont.mean #(4)
        
        final.sel.bin <- mean(bin.sel)#(7)
        final.inter<-prod(bin.sel)#(8)
        
        mean.cut.sel <- mean(cut.sel)#(9)
        inter.cut.sel <- prod(cut.sel)#(10)
        
        final <- stack(final.cont.mean,final.bin.mean,final.cut.mean,final.sel.bin,final.inter,mean.cut.sel,inter.cut.sel)
      }
        names(final) <- c("Final.cont.mean2","Final.bin.mean3","Final.cut.mean4","Final.mean.bin7","Final.inter.bin8","Mean.cut.sel9","inter.cut.sel10")
      
        if(exists("final")) {
        plot(final)
        
            #Escribe final
        writeRaster(x=final,filename=paste0("./",input.folder,"/",sp,"/",output.folder,"/",names(final),sp,algo),bylayer=T,overwrite=T,format="GTiff")
       
         for (i in 1:dim(final)[[3]]){
          png(filename=paste0(input.folder,"/",sp,"/",output.folder,"/",names(final)[i],sp,algo,".png"))
          plot(final[[i]],main=names(final)[i])
          dev.off()
        }
      
        }
      
    }
    
    if (weight.partitions==TRUE){
      for (par in unique(weight.par)){
        pond.stats<-stats2[,par]
        if (par=="TSS") pond.stats<-(pond.stats+1)/2
        pond <- mod.cont[[1:part]]*pond.stats
        final.pond <- mean(pond)    
        names(final.pond)<-paste0(par,"-Weighted")
        
        png(filename=paste0("./",input.folder,"/",sp,"/",output.folder,"/Final",par,"weighted",sp,algo,".png"))
        par(mar=c(4,4,3,3),mfrow=c(1,1))
        plot(final.pond,main=paste0(par,"-weighted ",sp," ",algo))
        dev.off()
        ##
        writeRaster(final.pond,filename=paste0("./",input.folder,"/",sp,"/",output.folder,"/Final",par,"weighted",sp,algo),overwrite=T,bylayer=T,format="GTiff")
        
        if(exists("final")){
            final <-  addLayer(final,final.pond)
            names(final)[length(names(final))]<-paste0(par,"-Weighted")
        }
        else {
            final <- final.pond
            names(final)
      }
    }
    plot(final)
    if(exists("final")) return(final)
    
  }
}
}
