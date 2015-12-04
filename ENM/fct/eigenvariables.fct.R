eigenvariables.fct <- function(vars, name){
library("raster")
if (file.exists("./pca")==FALSE) dir.create("./pca")
#Running PCA: 
# Counts not NA cels
non.na <- sum(!is.na(values(vars[[1]])))
# Sample the study area with n-non.na and creates an environmental table
sr <- sampleRandom(vars, non.na)
# faz o PCA dessa tabela padronizada
pca <- prcomp(scale(sr))

#Saving results:
capture.output(pca, file=sprintf('./pca/%s.pca.txt',name))

#saving summary
capture.output(summary(pca), file=sprintf('./pca/%s.summary.pca.txt',name))

#Plotting results
#GGPLOT
#####
#library(ggplot2)
# create data frame with scores
#scores = as.data.frame(pca$x)
# plot of observations
#ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores))) +
#  geom_hline(yintercept = 0, colour = "gray65") +
#  geom_vline(xintercept = 0, colour = "gray65") +
#  geom_text(colour = "tomato", alpha = 0.8, size = 4) +
#  ggtitle("PCA plot of USA States - Crime Rates")
#####
?png
png(filename = sprintf('./pca/%s.pca.biplot.png',name),
  bg = "white")
  biplot(pca)
  dev.off()

# Creating eigenvariable in space
eigenvariables <- predict(vars, pca, index=1:length(pca))
if (file.exists("./env")==FALSE) dir.create("./env")
writeRaster(eigenvariables,sprintf('./env/%s.eigenvariables',name),overwrite=T)
}