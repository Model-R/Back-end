cropModel <- function(modelo, mascara){
    modelo <- crop(modelo, mascara)
    modelo <- mask(modelo, mascara)
return(modelo)
    }