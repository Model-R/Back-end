
library(raster)
source("https://raw.githubusercontent.com/diogosbr/modelagem/master/clean.R")
abio=raster::stack(list.files("./env", pattern = '.gri',full.names = T))



pts=read.table("./data/FLORA_occs_final.csv", h=T, sep=",")
pts


#sp1
eug=pts[pts$sp==levels(pts$sp)[654],]
dim(eug)

raster::plot(abio[[1]])
points(eug[,3:4])

eug.c=clean(eug[,3:4], abio[[1]])
eug.c$sp=eug$sp[1]

#sp2
sp2=pts[pts$sp==levels(pts$sp)[5],]
dim(sp2)

raster::plot(abio[[1]])
points(sp2[,3:4])

sp2.c=clean(sp2[,3:4], abio[[1]])
sp2.c$sp=sp2$sp[1]

#sp3
#Não
sp3=pts[pts$sp==levels(pts$sp)[2046],]
dim(sp3)

raster::plot(abio[[1]])
points(sp3[,3:4])

sp3.c=clean(sp3[,3:4], abio[[1]])
sp3.c$sp=sp3$sp[1]

#sp4
sp4=pts[pts$sp==levels(pts$sp)[1431],]
dim(sp4)

raster::plot(abio[[1]])
points(sp4[,3:4])

sp4.c=clean(sp4[,3:4], abio[[1]])
sp4.c$sp=sp4$sp[1]

#sp5
sp5=pts[pts$sp==levels(pts$sp)[969],]
dim(sp5)

raster::plot(abio[[1]])
points(sp5[,3:4])

sp5.c=clean(sp5[,3:4], abio[[1]])
sp5.c$sp=sp5$sp[1]


#concatenando

dados=rbind(eug.c,sp2.c,sp4.c,sp5.c)
dados=dados[,c(3,1,2)]

dim(dados)

write.table(dados,"./data/dados_clean.csv", sep=",", row.names = F)








