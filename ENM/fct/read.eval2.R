#Function to analyse(compare) the evaluate output from SDM.and analyse the TSS results
#####-----
consolidate.data <- function(
                    input.folder1="Evaluate_Flora",
                    stat='mean'
                    ){
library("data.table")
  evall.list <- list.files(path = input.folder1, pattern=stat,full.names = T)
  lista<-list()
  for (i in 1:length(evall.list)) {
    lista[[i]]<-read.table(file = evall.list[i],header=T)
  }
  stats<-rbindlist(lista)
  stats<-as.data.frame(stats)
  return(stats)}
