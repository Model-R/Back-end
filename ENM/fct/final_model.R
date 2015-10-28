######### Final modeling: one model per algorithm ----

final.model <- function(sp,
                        select.partitions=T,
                        weight.partitions=F,                        
                        threshold=c("spec_sens"),
                        TSS.value=0.2,
                        #para select.partitions
                        mean=TRUE,
                        intersection=TRUE,
                        #para weigh partitions
                        weight.par=c("TSS","AUC"),
                        input.folder="2406models",
                        output.folder="presfinal"){
 
    if (file.exists(paste0("./",input.folder,"/",sp,"/",output.folder))==FALSE) dir.create(paste0("./",input.folder,"/",sp,"/",output.folder))
    cat(sp,"\n")
 library("raster")
 library("data.table")
cat(paste("Reading the evaluation files","\n"))
evall<- list.files(path = paste0(input.folder,"/",sp),pattern=paste0("evaluate",sp),full.names = T)
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
#if(algoritmo=="Mahal"){
 #   stats2$spec_sens[stats2$spec_sens<0]<-0
#}
#stats2 <- stats2[order(stats2$partition),]
part <- nrow(stats2)#How many partitions were there
    
cat(paste("Reading models from .tif files","\n"))
modelos <- list.files(path = paste0(input.folder,"/",sp),full.names=T,pattern=paste0(algo,"_cont_",sp))
mod<-stack(modelos)#(0)
names(mod)<-paste0(sp,"Partition",1:part)

#Binary by TSSth and Cut    
bin <- mod>stats2[,names(stats2)=="spec_sens"]#stack
cut <- bin * mod#stack


    if (select.partitions==T){
        sel.index<- which(stats2[,"TSS"]>=TSS.value)
        mod.sel<- mod[[sel.index]]
        if (length(sel.index)==0) cat(paste("No partition was selected for",sp,"\n"))
        if (length(sel.index)>0){
        
        mod.sel<- mod[[sel.index]] #(1)
        bin.sel<-mod.sel>stats2[,names(stats2)=="spec_sens"][sel.index] #(5)
        cut.sel<-bin.sel*mod.sel#(8)
        
        th.mean<-mean(stats2[,names(stats2)=="spec_sens"][sel.index])
        
        }
        #en caso de que sea solo uno varios modelos son el mismo
        if (length(sel.index)==1){
          cat(paste(length(sel.index), "partitions was selected for",sp))
          
        final.sel.cont<-mod.sel#(1)(2)
        final.sel.bin<-bin.sel#(5)(3)(7) (8)
        final.sel.cut<-cut.sel#(4)(6)(9)(10)
        
        final <- stack(mod.sel,bin.sel,cut.sel,bin.sel,bin.sel,cut.sel,cut.sel)
        names(final) <- c("Final.cont.mean2","Final.bin.mean3","Final.cut.mean4","Final.mean.bin7","Final.inter.bin8","Mean.cut.sel9","inter.cut.sel10")
        
        }
            
        #en caso de que sean más aplica el mapa          
        if (length(sel.index)>1){
          cat(paste(length(sel.index), "partitions were selected for",sp))
        final.cont.mean <- mean(mod.sel)#(2)
        final.bin.mean <- (final.cont.mean>th.mean)#(3)
        final.cut.mean <- final.bin.mean*final.cont.mean #(4)
        
        final.sel.bin <- mean(bin.sel)#(7)
        final.inter<-prod(bin.sel)#(8)
        
        mean.cut.sel <- mean(cut.sel)#(9)
        inter.cut.sel <- prod(cut.sel)#(10)
        final <- stack(final.cont.mean,final.bin.mean,final.cut.mean,final.sel.bin,final.inter,mean.cut.sel,inter.cut.sel)
        names(final) <- c("Final.cont.mean2","Final.bin.mean3","Final.cut.mean4","Final.mean.bin7","Final.inter.bin8","Mean.cut.sel9","inter.cut.sel10")
        
        }
        if(exists("final")) {
            plot(final)
        #Escribe mean binary de los seleccionados
        writeRaster(x=final,filename=paste0("./",input.folder,"/",sp,"/",output.folder,"/",names(final),sp,algo),bylayer=T,overwrite=T,format="GTiff")
        for (i in 1:dim(final)[[3]]){
        png(filename=paste0(input.folder,"/",sp,"/",output.folder,"/",names(final)[i],sp,algo,".png"))
            plot(final[[i]],main=names(final)[i])
        dev.off()
        }
        }
#         #.sel.bin,filename=paste0("./",output.folder,"/Mean_binary",threshold,"_sel_",sp,"_",algoritmo,".grd"),overwrite=T)
#               png(filename=paste0("./",output.folder,"/Mean_binary",threshold,"_sel_",sp,"_",algoritmo,".png"),width=500 ,height=500)
#               par(mar=c(4,4,3,3),mfrow=c(1,1))
#               plot(final.sel.bin,main=paste("Mean binary model",threshold,">",TSS.value,"\n",sp,algoritmo))
#               dev.off()
#         #Escribe mean continous cortado por el TSS medio, como maxent:
#         writeRaster(x=final.mean.threshold,filename=paste0("./",output.folder,"/Mean_continuous_meanTSS",threshold,"_sel_",sp,"_",algoritmo,".grd"),overwrite=T)
#         png(filename=paste0("./",output.folder,"/Mean_continuous_meanTSS",threshold,"_sel_",sp,"_",algoritmo,".png"),width=500 ,height=500)
#         par(mar=c(4,4,3,3),mfrow=c(1,1))
#         plot(final.mean.threshold,main=paste("Mean continuous model cut by mean TSS",threshold,">",TSS.value,"\n",sp,algoritmo))
#         dev.off()
#           }
          #if (length(sel.index)>1 & intersection==TRUE) { 
              #final.sel <- prod(bin.sel)
           #   writeRaster(x=final.sel,filename=paste0("./",output.folder,"/Intersection",threshold,"_sel_",sp,"_",algoritmo,".grd"),overwrite=T)
          #png(filename=paste0("./",output.folder,"/Intersection",threshold,"_sel_",sp,"_",algoritmo,".png"),width=500 ,height=500)
          #par(mar=c(4,4,3,3),mfrow=c(1,1))
          #plot(final.sel,main=paste("Intersection model",threshold,">",TSS.value,"\n",sp,algoritmo))
          #dev.off()
          #}
        }
    
    if (weight.partitions==TRUE){
        for (par in unique(weight.par)){
        pond.stats<-stats2[,par]
        #if (par=="TSS") pond.stats<-(pond.stats+1)/2 será?
        pond <- mod[[1:part]]*pond.stats
        final <- mean(pond)    
    
        png(filename=paste0("./",input.folder,"/",sp,"/",output.folder,"/Final",par,"Wmean",sp,algo,".png"))
        par(mar=c(4,4,3,3),mfrow=c(1,1))
        plot(final,main=paste(par,"Wmean",sp,algo))
        dev.off()
        ##
        writeRaster(final,filename=paste0("./",input.folder,"/",sp,"/",output.folder,"/Final",par,"Wmean",sp,algo),overwrite=T,bylayer=T,format="GTiff")
        
        #final<-addLayer(final,mod)
        #names(final)[length(names(final))]<-paste("Weighted",par)
            }
    }
 if(exists("final")) return(final)

    }
}