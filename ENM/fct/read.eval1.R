#####-----
#Function to organize the evaluate output from SDM to be used with xonsolidate.data().
read.eval<-function(sp,
                    input.folder1="FLORA_buffermax",
                    input.folder2="FLORA_sem_buffer",
                    output.folder="TSS_Evaluate"){

library("data.table")
  
  evall.buffer <- list.files(path = paste0(input.folder1,"/",sp), pattern=paste0("evaluate",sp),full.names = T)
  evall <- list.files(path = paste0(input.folder2,"/",sp), pattern=paste0("evaluate",sp),full.names = T)
  
  # Resultados com Buffer ----
  lista.buffer<-list()
  lista.buffer.final<-list()
  for (i in 1:length(evall.buffer)) {
    lista.buffer[[i]]<-read.table(file = evall.buffer[i],header=T,row.names=1)
    lista.buffer.final[[i]]<-lista.buffer[[i]]['TSS']
    lista.buffer.final[[i]][2]<-lista.buffer[[i]]['algoritmo']
  }
  
  mean.table.buffer<-matrix('NULL', nrow=length(lista.buffer.final)+2, ncol=dim(lista.buffer.final[[1]][1])[1]+2)
  for (a in 1:length(lista.buffer.final)){
    for (b in 1:dim(lista.buffer.final[[1]][1])[1]){
      mean.table.buffer[a,b] <- as.matrix(lista.buffer.final[[a]][1])[b]
    }
  }
  
  for (c in 1:length(mean.table.buffer[1,])){
    mean.table.buffer[nrow(mean.table.buffer)-1,c]<- mean(as.numeric((mean.table.buffer[,c])),na.rm=TRUE)
    mean.table.buffer[nrow(mean.table.buffer),c]<- sd(as.numeric((mean.table.buffer[,c])),na.rm=TRUE)
  }
  colnames(mean.table.buffer)<- c(as.list(paste(lista.buffer.final[[1]][,2])),'stat','tratamento')
  mean.table.buffer[,ncol(mean.table.buffer)-1]<-c(paste0('Com Buffer ',1:(nrow(mean.table.buffer)-2)),'mean','sd')
  mean.table.buffer[,ncol(mean.table.buffer)]<-rep('Com Buffer')
  
  #Resultados sem buffer----
  lista<-list()
  lista.final<-list()
  for (i in 1:length(evall)) {
    lista[[i]]<-read.table(file = evall[i],header=T,row.names=1)
    lista.final[[i]]<-lista[[i]]['TSS']
    lista.final[[i]][2]<-lista[[i]]['algoritmo']
  }
  mean.table<-matrix('NULL', nrow=length(lista.final)+2, ncol=dim(lista[[1]][1])[1]+2)
  for (a in 1:length(lista.final)){
    for (b in 1:dim(lista.final[[1]][1])[1]){
      mean.table[a,b] <- as.matrix(lista.final[[a]][1])[b]
    }
  }
  colnames(mean.table)<-c(as.list(paste(lista.final[[1]][,2])), 'stat','tratamento')
  for (c in 1:length(mean.table[1,])){
    mean.table[nrow(mean.table)-1,c]<- mean(as.numeric((mean.table[,c])),na.rm=TRUE)
    mean.table[nrow(mean.table),c]<- sd(as.numeric((mean.table[,c])),na.rm=TRUE)
  }
  
  
  mean.table[,ncol(mean.table)-1]<-c(paste0('sem Buffer ',1:(nrow(mean.table)-2)),'mean','sd')
  mean.table[,ncol(mean.table)]<-rep('sem Buffer')
  
  # Criando data.frame final ----
  #Finalizando mean
  intermed_mean <- rbind(mean.table[nrow(mean.table)-1,], mean.table.buffer[nrow(mean.table.buffer)-1,])
  final_mean <- data.frame()
  for (d in 1:(ncol(intermed_mean)-2)){
    #Sem bufer----
    final_mean[d,1]<-intermed_mean[1,d]
    final_mean[d,2]<- colnames(intermed_mean)[d]
    final_mean[d,3]<- intermed_mean[1,(ncol(intermed_mean))]
    #com buffer -----
    final_mean[(ncol(intermed_mean)-2)+d,1]<-intermed_mean[2,d]
    final_mean[(ncol(intermed_mean)-2)+d,2]<- colnames(intermed_mean)[d]
    final_mean[(ncol(intermed_mean)-2)+d,3]<- intermed_mean[2,(ncol(intermed_mean))]
  }
  
  #Finalizando SD
  intermed_sd <- rbind(mean.table[nrow(mean.table),], mean.table.buffer[nrow(mean.table.buffer),])
  final_sd <- data.frame()
  for (d in 1:(ncol(intermed_sd)-2)){
    #Sem bufer----
    final_sd[d,1]<-intermed_sd[1,d]
    final_sd[d,2]<- colnames(intermed_sd)[d]
    final_sd[d,3]<- intermed_sd[1,(ncol(intermed_sd))]
    #com buffer -----
    final_sd[(ncol(intermed_sd)-2)+d,1]<-intermed_sd[2,d]
    final_sd[(ncol(intermed_sd)-2)+d,2]<- colnames(intermed_sd)[d]
    final_sd[(ncol(intermed_sd)-2)+d,3]<- intermed_mean[2,(ncol(intermed_sd))]
  }
  colnames(final_mean)<-c('mean','Algorithm','Tratamento')
  colnames(final_sd)<-c('mean','Algorithm','Tratamento')
  
  #removendo dados já utilizados:
  #Buffer
  rm(evall.buffer)
  rm(lista.buffer)
  rm(lista.buffer.final)
  rm(mean.table.buffer)
  #s/buffer
  rm(evall)
  rm(lista)
  rm(lista.final)
  rm(mean.table)
  #demais
  rm(intermed_mean)
  rm(intermed_sd)
  
#Salvando em HD
  if (file.exists(paste0("./",output.folder))==FALSE) dir.create(paste0("./",output.folder))
  
  write.table(final_mean,file=paste0('./',output.folder,'/',sp,'_mean_results_tss.csv'), col.names=TRUE,row.names=FALSE)
  
  write.table(final_sd,file=paste0('./',output.folder,'/',sp,'_sd_results_tss.csv'), col.names=TRUE,row.names=FALSE)
  #Nao remover
  rm(final_mean)
  rm(final_sd)
  rm(input.folder1)
  rm(input.folder2)
  rm(output.folder)
}

#####-----
#Function to analyse(compare) the evaluate output from SDM, producing TSS Boxplot for each algorithm
consolidate.data <- function(
  input.folder1="TSS_Evaluate",
  stat='mean',
  name='TSS_analysis'
){
  
  library("data.table")
  evall.list <- list.files(path = input.folder1, pattern=stat,full.names = T) #listing all evaluate files previusly organized. see read.eval()
  lista<-list()
  for (i in 1:length(evall.list)) {
    lista[[i]]<-read.table(file = evall.list[i],header=T)
  }
  stats<-rbindlist(lista) #combining list
  stats<-as.data.frame(stats)
  #return(stats)
  #Saving the boxplot graph in PNG
  library(ggplot2)
  ggsave(paste0("./",input.folder1,"/",name,".png"))
  ggplot(stats, aes(x=Tratamento, y=mean, fill=Tratamento)) + geom_boxplot(notch=TRUE) + facet_grid(~Algorithm) + 
    labs(y='') +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_text(size=rel(1.5)),
          strip.text = element_text(face="bold", size=rel(1.25)),
          legend.title=element_text(face='bold',size=rel(1.5)),
          legend.text=element_text(size=rel(1.5)),
          legend.position="bottom")
  dev.off()
}

#####-----
#Small script to save the species name and final occurrence.
species.table <- function(
  output.folder="TSS_Evaluate",
  table.name="n.sp.modelado"){
n.sp<-table(occs$sp)
n.sp<-as.data.frame(n.sp)
#head(n.sp)
#dim(n.sp)
n.sp <- subset(n.sp, Freq > 10) 
#dim(n.sp)
n.sp.modelado <- subset(n.sp, Var1 %in% (names.sp[N])) 
#head(n.sp.modelado)
#dim(n.sp.modelado)
write.csv(n.sp.modelado, paste0(output.folder,'/',table.name,'.csv'))
}
