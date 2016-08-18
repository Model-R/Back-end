# Cada função de modelagem do_ALGORITMO retorna o threshold e escreve os arquivos
# de saída no diretório informado pelo parâmetro models.dir

do_bioclim <- function(sp,
                       predictors = predictors,
                       pres_train = pres_train, #NEW
                       pres_test = pres_test, #NEW
                       backg_test = backg_test, #NEW
                       i = i, #NEW
                       models.dir = model.dir,
                       project.model = project.model,
                       projections = projections,
                       projdata = projdata,#um vector con nombres
                       mask = mask# a SpatialPolygonsDataFrame layer to mask and crop the predicted model
                       ){
    cat(paste("Bioclim",'\n'))
    bc <- bioclim (predictors, pres_train)
    ebc <- evaluate(pres_test,backg_test,bc,predictors)
    thresholdbc <- ebc@t[which.max(ebc@TPR + ebc@TNR)]
    thbc <- threshold(ebc)
    bc_TSS <- max(ebc@TPR + ebc@TNR)-1

    bc_cont <- predict(predictors, bc, progress='text')
    bc_bin <- bc_cont > thresholdbc
    bc_cut <- bc_cont * bc_bin
    thbc$AUC <- ebc@auc
    thbc$TSS <- bc_TSS
    thbc$algoritmo <- "BioClim"
    thbc$partition <- i
    row.names(thbc) <- paste(sp,i,"BioClim")

    if (class(mask) == "SpatialPolygonsDataFrame"){
        source("./fct/cropModel.R")
        bc_cont <- cropModel(bc_cont , mask)
        bc_bin <- cropModel(bc_bin , mask)
        bc_cut <- cropModel(bc_cut , mask)

    }
    writeRaster(x=bc_cont,filename=paste0(models.dir,"/",sp,"/BioClim_cont_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
    writeRaster(x=bc_bin,filename=paste0(models.dir,"/",sp,"/BioClim_bin_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
    writeRaster(x=bc_cut,filename=paste0(models.dir,"/",sp,"/BioClim_cut_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")

    png(filename=paste0(models.dir,"/",sp,"/Bioclim",sp,"_",i,"%03d.png"))
    plot(bc_cont,main=paste("Bioclim raw","\n","AUC =", round(ebc@auc,2),'-',"TSS =",round(bc_TSS,2)))
    plot(bc_bin,main=paste("Bioclim P/A","\n","AUC =", round(ebc@auc,2),'-',"TSS =",round(bc_TSS,2)))
    plot(bc_cut,main=paste("Bioclim cut","\n","AUC =", round(ebc@auc,2),'-',"TSS =",round(bc_TSS,2)))
    dev.off()

    if (project.model == T){
        for (proj in projections){
            #data <- list.files(paste0("./env/",proj),pattern=proj)
            #data2 <- stack(data)
            stopifnot(names(projdata) == names(predictors))
            bc_proj <- predict(projdata,bc,progress='text')
            bc_proj_bin <- bc_proj > thresholdbc
            bc_proj_cut <- bc_proj_bin * bc_proj
            # Normaliza o modelo cut
            #bc_proj_cut <- bc_proj_cut/maxValue(bc_proj_cut)
            if (class(mask) == "SpatialPolygonsDataFrame"){
                source("./fct/cropModel.R")
                bc_proj <- cropModel(bc_proj , mask)
                bc_proj_bin <- cropModel(bc_proj_bin , mask)
                bc_proj_cut <- cropModel(bc_proj_cut , mask)
            }
            writeRaster(x=bc_proj,filename=paste0(models.dir,"/",sp,"/",proj,"/BioClim_cont_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
            writeRaster(x=bc_proj_bin,filename=paste0(models.dir,"/",sp,"/",proj,"/BioClim_bin_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
            writeRaster(x=bc_proj_cut,filename=paste0(models.dir,"/",sp,"/",proj,"/BioClim_cut_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
            rm(data2)
            png(filename=paste0(models.dir,"/",sp,"/",proj,"/Bioclim",sp,"_",i,"%03d.png"))
            plot(bc_proj,main=paste("Bioclim raw","\n","AUC =", round(ebc@auc,2),'-',"TSS =",round(bc_TSS,2)))
            plot(bc_proj_bin,main=paste("Bioclim P/A","\n","AUC =", round(ebc@auc,2),'-',"TSS =",round(bc_TSS,2)))
            plot(bc_proj_cut,main=paste("Bioclim cut","\n","AUC =", round(ebc@auc,2),'-',"TSS =",round(bc_TSS,2)))
            dev.off()
        }
    }
    return(thbc)
}

do_randomForest <- function(sp,
                            predictors = predictors,
                            sdmdata_train = sdmdata_train,
                            envtest_pre = envtest_pre,
                            envtest_back = envtest_back,
                            i = i,
                            envtrain = envtrain,
                            models.dir = model.dir,
                            project.model = project.model,
                            projections = projections,
                            projdata = projdata,#um vector con nombres
                            mask = mask# a SpatialPolygonsDataFrame layer to mask and crop the predicted model
                            ){
    library(randomForest)
    cat(paste("Random Forests",'\n'))
    #rf1 <- tuneRF(x=envtrain,y=sdmdata_train$pa,stepFactor = 0.5)
    rf <- randomForest (sdmdata_train$pa~.,data=envtrain)
    #rf <- randomForest (x =envtrain ,y=factor(sdmdata_train$pa),xtest=envtest,ytest = factor(sdmdata_teste$pa))#fazendo teste interno a funcao evaluate nao serve :(

    erf <- evaluate(envtest_pre,envtest_back,rf)
    rf_TSS <- max(erf@TPR + erf@TNR)-1

    thresholdrf <- erf@t[which.max(erf@TPR + erf@TNR)]
    thrf <- threshold(erf)
    thrf$AUC <- erf@auc
    thrf$TSS <- rf_TSS#raro
    thrf$algoritmo <- "rf"
    thrf$partition <- i
    row.names(thrf) <- paste(sp,i,"rf")

    rf_cont <- predict(predictors,rf,progress='text',type="response")
    rf_bin <- rf_cont>thresholdrf
    rf_cut <- rf_bin * rf_cont
    #rf1_cut <- rf1_cut/maxValue(rf1_cut)
    if (class(mask) == "SpatialPolygonsDataFrame"){
        source("./fct/cropModel.R")
        rf_cont <- cropModel(rf_cont , mask)
        rf_bin <- cropModel(rf_bin , mask)
        rf_cut <- cropModel(rf_cut , mask)
    }
    writeRaster(x=rf_cont,filename=paste0(models.dir,"/",sp,"/rf_cont_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
    writeRaster(x=rf_bin,filename=paste0(models.dir,"/",sp,"/rf_bin_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
    writeRaster(x=rf_cut,filename=paste0(models.dir,"/",sp,"/rf_cut_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")

    png(filename=paste0(models.dir,"/",sp,"/rf",sp,"_",i,"%03d.png"))
    plot(rf_cont,main=paste("RF raw","\n","AUC =", round(erf@auc,2),'-',"TSS =",round(rf_TSS,2)))
    plot(rf_bin,main=paste("RF P/A","\n","AUC =", round(erf@auc,2),'-',"TSS =",round(rf_TSS,2)))
    plot(rf_cut,main=paste("RF cut","\n","AUC =", round(erf@auc,2),'-',"TSS =",round(rf_TSS,2)))
    dev.off()


    if (project.model == T){
        for (proj in projections){
            data <- list.files(paste0("./env/",proj),pattern=proj)
            data2 <- stack(data)
            rf_proj <- predict(data2,rf,progress='text')
            rf_proj_bin <- rf_proj > thresholdrf
            rf_proj_cut <- rf_proj_bin * rf_proj
            # Normaliza o modelo cut
            #rf_proj_cut <- rf_proj_cut/maxValue(rf_proj_cut)
            if (class(mask) == "SpatialPolygonsDataFrame"){
                source("./fct/cropModel.R")
                rf_proj <- cropModel(rf_proj , mask)
                rf_proj_bin <- cropModel(rf_proj_bin , mask)
                rf_proj_cut <- cropModel(rf_proj_cut , mask)
            }
            writeRaster(x=rf_proj,filename=paste0(models.dir,"/",sp,"/",proj,"/rf_cont_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
            writeRaster(x=rf_proj_bin,filename=paste0(models.dir,"/",sp,"/",proj,"/rf_bin_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
            writeRaster(x=rf_proj_cut,filename=paste0(models.dir,"/",sp,"/",proj,"/rf_cut_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
            rm(data2)
        }
    }
    return(thrf)
}

do_SVM <- function(sp,
                   predictors = predictors,
                   sdmdata_train = sdmdata_train,
                   envtest_pre = envtest_pre,
                   envtest_back = envtest_back,
                   i = i,
                   envtrain = envtrain,
                   part = part,
                   models.dir = model.dir,
                   project.model = project.model,
                   projections = projections,
                   projdata = projdata,#um vector con nombres
                   mask = mask# a SpatialPolygonsDataFrame layer to mask and crop the predicted model
                   ){
    cat(paste("SVM",'\n'))
    library(kernlab)
    svm <- ksvm(sdmdata_train$pa~.,data=envtrain,cross=part)##svm deve ser com a variável resposta binária ou contínua, eu acho que binária
    esvm <- evaluate(envtest_pre,envtest_back,svm)
    #esvm <- evaluate(pres_test,backg_test,model = svm,x = predictors)
    svm_TSS <- max(esvm@TPR + esvm@TNR)-1
    thresholdsvm <- esvm@t[which.max(esvm@TPR + esvm@TNR)]
    thsvm <- threshold (esvm)
    thsvm$AUC <- esvm@auc
    thsvm$TSS <- svm_TSS
    thsvm$algoritmo <- "svm"
    thsvm$partition <- i
    row.names(thsvm) <- paste(sp,i,"svm")
    svm_cont <- predict(predictors,svm,progress='text')
    svm_bin <- svm_cont>thresholdsvm
    svm_cut <- svm_bin * svm_cont

    #TRANSFORMA 0 A 1
    svm_cont <- svm_cont/maxValue(svm_cont)
    svm_cut <- svm_cut/maxValue(svm_cut)
    if (class(mask) == "SpatialPolygonsDataFrame"){
        source("./fct/cropModel.R")
        svm_cont <- cropModel(svm_cont , mask)
        svm_bin <- cropModel(svm_bin , mask)
        svm_cut <- cropModel(svm_cut , mask)
    }

    writeRaster(x=svm_cont,filename=paste0(models.dir,"/",sp,"/svm_cont_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
    writeRaster(x=svm_bin,filename=paste0(models.dir,"/",sp,"/svm_bin_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
    writeRaster(x=svm_cut,filename=paste0(models.dir,"/",sp,"/svm_cut_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")

    png(filename=paste0(models.dir,"/",sp,"/svm",sp,"_",i,"%03d.png"))
    plot(svm_cont,main=paste("SVM raw","\n","AUC =", round(esvm@auc,2),'-',"TSS =",round(svm_TSS,2)))
    plot(svm_bin,main=paste("SVM P/A","\n","AUC =", round(esvm@auc,2),'-',"TSS =",round(svm_TSS,2)))
    plot(svm_cut,main=paste("SVM cut","\n","AUC =", round(esvm@auc,2),'-',"TSS =",round(svm_TSS,2)))
    dev.off()

    if (project.model == T){
        for (proj in projections){
            data <- list.files(paste0("./env/",proj),pattern=proj)
            data2 <- stack(data)
            svm_proj <- predict(data2,svm,progress='text')
            svm_proj_bin <- svm_proj > thresholdsvm
            svm_proj_cut <- svm_proj_bin * svm_proj

            # Normaliza o modelo cut
            #svm_proj_cut <- svm_proj_cut/maxValue(svm_proj_cut)
            if (class(mask) == "SpatialPolygonsDataFrame"){
                source("./fct/cropModel.R")
                svm_proj <- cropModel(svm_proj , mask)
                svm_proj_bin <- cropModel(svm_proj_bin , mask)
                svm_proj_cut <- cropModel(svm_proj_cut , mask)
            }
            writeRaster(x=svm_proj,filename=paste0(models.dir,"/",sp,"/",proj,"/svm_cont_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
            writeRaster(x=svm_proj_bin,filename=paste0(models.dir,"/",sp,"/",proj,"/svm_bin_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
            writeRaster(x=svm_proj_cut,filename=paste0(models.dir,"/",sp,"/",proj,"/svm_cut_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
            rm(data2)
        }
    }
    return(thsvm)
}

do_maxent <- function(sp,
                      predictors = predictors,
                      pres_train = pres_train, #NEW
                      pres_test = pres_test, #NEW
                      backg_test = backg_test, #NEW
                      i = i, #NEW
                      models.dir = model.dir,
                      project.model = project.model,
                      projections = projections,
                      projdata = projdata,#um vector con nombres
                      mask = mask# a SpatialPolygonsDataFrame layer to mask and crop the predicted model
                      ){
    cat(paste("maxent",'\n'))
    Sys.setenv(NOAWT=TRUE)#descomentei para ver
    library(rJava)
    mx <- maxent (predictors, pres_train)
    png(filename = paste0(models.dir,"/",sp,"/maxent_variable_contribution_",sp,"_",i,".png"))
    plot(mx)
    dev.off()

    png(filename = paste0(models.dir,"/",sp,"/maxent_variable_response_",sp,"_",i,".png"))
    response(mx)
    dev.off()

    emx <- evaluate(pres_test,backg_test,mx,predictors)
    thresholdmx <- emx@t[which.max(emx@TPR + emx@TNR)]
    thmx <- threshold(emx)
    mx_TSS <- max(emx@TPR + emx@TNR)-1
    mx_cont <- predict(mx, predictors,progress='text')
    mx_bin <- mx_cont > thresholdmx
    mx_cut <- mx_cont * mx_bin
    thmx$AUC <- emx@auc
    thmx$TSS <- mx_TSS
    thmx$algoritmo <- "maxent"
    thmx$partition <- i
    row.names(thmx) <- paste(sp,i,"maxent")
    if (class(mask) == "SpatialPolygonsDataFrame"){
        source("./fct/cropModel.R")
        mx_cont <- cropModel(mx_cont , mask)
        mx_bin <- cropModel(mx_bin , mask)
        mx_cut <- cropModel(mx_cut , mask)
    }
    writeRaster(x=mx_cont,filename=paste0(models.dir,"/",sp,"/maxent_cont_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
    writeRaster(x=mx_bin,filename=paste0(models.dir,"/",sp,"/maxent_bin_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
    writeRaster(x=mx_cut,filename=paste0(models.dir,"/",sp,"/maxent_cut_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")

    png(filename=paste0(models.dir,"/",sp,"/maxent",sp,"_",i,"%03d.png"))
    plot(mx_cont,main=paste("Maxent raw","\n","AUC =", round(emx@auc,2),'-',"TSS =",round(mx_TSS,2)))
    plot(mx_bin,main=paste("Maxent P/A","\n","AUC =", round(emx@auc,2),'-',"TSS =",round(mx_TSS,2)))
    plot(mx_cut,main=paste("Maxent cut","\n","AUC =", round(emx@auc,2),'-',"TSS =",round(mx_TSS,2)))
    dev.off()

    if (project.model == T){
        for (proj in projections){
            data <- list.files(paste0("./env/",proj),pattern=proj)
            data2 <- stack(data)
            mx_proj <- predict(data2,mx,progress='text')
            mx_proj_bin <- mx_proj > thresholdmx
            mx_proj_cut <- mx_proj_bin * mx_proj
            # Normaliza o modelo cut
            #mx_proj_cut <- mx_proj_cut/maxValue(mx_proj_cut)
            if (class(mask) == "SpatialPolygonsDataFrame"){
                source("./fct/cropModel.R")
                mx_proj <- cropModel(mx_proj , mask)
                mx_proj_bin <- cropModel(mx_proj_bin , mask)
                mx_proj_cut <- cropModel(mx_proj_cut , mask)
            }
            writeRaster(x=mx_proj,filename=paste0(models.dir,"/",sp,"/",proj,"/maxent_cont_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
            writeRaster(x=mx_proj_bin,filename=paste0(models.dir,"/",sp,"/",proj,"/maxent_bin_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
            writeRaster(x=mx_proj_cut,filename=paste0(models.dir,"/",sp,"/",proj,"/maxent_cut_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
            rm(data2)
        }
    }
    return(thmx)
}

do_GLM <- function(sp,
                   predictors = predictors,
                   sdmdata_train = sdmdata_train,
                   envtest_pre = envtest_pre,
                   envtest_back = envtest_back,
                   i = i,
                   envtrain = envtrain,
                   models.dir = model.dir,
                   project.model = project.model,
                   projections = projections,
                   projdata = projdata,#um vector con nombres
                   mask = mask# a SpatialPolygonsDataFrame layer to mask and crop the predicted model
                   ){##
    cat(paste("GLM",'\n'))
    null.model <- glm(sdmdata_train$pa~1,data=envtrain,family="binomial")
    full.model <- glm(sdmdata_train$pa~.,data=envtrain,family="binomial")
    glm <- step(object = null.model,scope = formula(full.model),direction = "both",trace=F)
    eglm <- evaluate(envtest_pre,envtest_back,model=glm,type="response")#####
    #eglm <- evaluate(pres_test,backg_test,glm,predictors,type="response")
    glm_TSS <- max(eglm@TPR + eglm@TNR)-1
    thresholdglm <- eglm@t[which.max(eglm@TPR + eglm@TNR)]
    thglm <- threshold (eglm)
    thglm$AUC <- eglm@auc
    thglm$TSS <- glm_TSS
    thglm$algoritmo <- "glm"
    thglm$partition <- i
    row.names(thglm) <- paste(sp,i,"glm")

    glm_cont <- predict(predictors,glm,progress='text',type="response")
    glm_bin <- glm_cont>thresholdglm
    glm_cut <- glm_bin * glm_cont
    # Normaliza o modelo cut
    #glm_cut <- glm_cut/maxValue(glm_cut)
    if (class(mask) == "SpatialPolygonsDataFrame"){
        source("./fct/cropModel.R")
        glm_cont <- cropModel(glm_cont , mask)
        glm_bin <- cropModel(glm_bin , mask)
        glm_cut <- cropModel(glm_cut , mask)
    }
    writeRaster(x=glm_cont,filename=paste0(models.dir,"/",sp,"/glm_cont_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
    writeRaster(x=glm_bin,filename=paste0(models.dir,"/",sp,"/glm_bin_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
    writeRaster(x=glm_cut,filename=paste0(models.dir,"/",sp,"/glm_cut_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")

    png(filename=paste0(models.dir,"/",sp,"/glm",sp,"_",i,"%03d.png"))
    plot(glm_cont,main=paste("GLM raw","\n","AUC =", round(eglm@auc,2),'-',"TSS =",round(glm_TSS,2)))
    plot(glm_bin,main=paste("GLM P/A","\n","AUC =", round(eglm@auc,2),'-',"TSS =",round(glm_TSS,2)))
    plot(glm_cut,main=paste("GLM cut","\n","AUC =", round(eglm@auc,2),'-',"TSS =",round(glm_TSS,2)))
    dev.off()

    if (project.model == T){
        for (proj in projections){
            data <- list.files(paste0("./env/",proj),pattern=proj)
            data2 <- stack(data)
            glm_proj <- predict(data2,glm,progress='text')
            glm_proj_bin <- glm_proj > thresholdglm
            glm_proj_cut <- glm_proj_bin * glm_proj
            # Normaliza o modelo cut
            #glm_proj_cut <- glm_proj_cut/maxValue(glm_proj_cut)

            if (class(mask) == "SpatialPolygonsDataFrame"){
                source("./fct/cropModel.R")
                glm_proj <- cropModel(glm_proj , mask)
                glm_proj_bin <- cropModel(glm_proj_bin , mask)
                glm_proj_cut <- cropModel(glm_proj_cut , mask)
            }
            writeRaster(x=glm_proj,filename=paste0(models.dir,"/",sp,"/",proj,"/glm_cont_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
            writeRaster(x=glm_proj_bin,filename=paste0(models.dir,"/",sp,"/",proj,"/glm_bin_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
            writeRaster(x=glm_proj_cut,filename=paste0(models.dir,"/",sp,"/",proj,"/glm_cut_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
            rm(data2)
        }
    }
    return(thglm)
}

do_domain <- function(sp,
                      predictors = predictors,
                      pres_train = pres_train, #NEW
                      pres_test = pres_test, #NEW
                      backg_test = backg_test, #NEW
                      i = i, #NEW
                      models.dir = model.dir,
                      project.model = project.model,
                      projections = projections,
                      projdata = projdata,#um vector con nombres
                      mask = mask# a SpatialPolygonsDataFrame layer to mask and crop the predicted model
                      ){
    cat(paste("Domain",'\n'))
    do <- domain (predictors, pres_train)
    edo <- evaluate(pres_test,backg_test,do,predictors)
    thresholddo <- edo@t[which.max(edo@TPR + edo@TNR)]
    thdo <- threshold(edo)
    do_TSS <- max(edo@TPR + edo@TNR)-1
    do_cont <- predict(predictors, do, progress='text')
    do_bin <- do_cont > thresholddo
    do_cut <- do_cont * do_bin
    thdo$AUC <- edo@auc
    thdo$TSS <- do_TSS
    thdo$algoritmo <- "Domain"
    thdo$partition <- i
    row.names(thdo) <- paste(sp,i,"Domain")

    if (class(mask) == "SpatialPolygonsDataFrame"){
        source("./fct/cropModel.R")
        do_cont <- cropModel(do_cont , mask)
        do_bin <- cropModel(do_bin , mask)
        do_cut <- cropModel(do_cut , mask)
    }
    writeRaster(x=do_cont,filename=paste0(models.dir,"/",sp,"/Domain_cont_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
    writeRaster(x=do_bin,filename=paste0(models.dir,"/",sp,"/Domain_bin_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
    writeRaster(x=do_cut,filename=paste0(models.dir,"/",sp,"/Domain_cut_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")

    png(filename=paste0(models.dir,"/",sp,"/Domain",sp,"_",i,"%03d.png"))
    plot(do_cont,main=paste("Domain raw","\n","AUC =", round(edo@auc,2),'-',"TSS =",round(do_TSS,2)))
    plot(do_bin,main=paste("Domain P/A","\n","AUC =", round(edo@auc,2),'-',"TSS =",round(do_TSS,2)))
    plot(do_cut,main=paste("Domain cut","\n","AUC =", round(edo@auc,2),'-',"TSS =",round(do_TSS,2)))
    dev.off()


    if (project.model == T){
        for (proj in projections){
            data <- list.files(paste0("./env/",proj),pattern=proj)
            data2 <- stack(data)
            do_proj <- predict(data2,do,progress='text')
            do_proj_bin <- do_proj > thresholddo
            do_proj_cut <- do_proj_bin * do_proj
            # Normaliza o modelo cut
            #do_proj_cut <- do_proj_cut/maxValue(do_proj_cut)
            if (class(mask) == "SpatialPolygonsDataFrame"){
                source("./fct/cropModel.R")
                do_proj <- cropModel(do_proj , mask)
                do_proj_bin <- cropModel(do_proj_bin , mask)
                do_proj_cut <- cropModel(do_proj_cut , mask)
            }
            writeRaster(x=do_proj,filename=paste0(models.dir,"/",sp,"/",proj,"/Domain_cont_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
            writeRaster(x=do_proj_bin,filename=paste0(models.dir,"/",sp,"/",proj,"/Domain_bin_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
            writeRaster(x=do_proj_cut,filename=paste0(models.dir,"/",sp,"/",proj,"/Domain_cut_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
            rm(data2)
        }
    }
    return(thdo)
}

do_mahal <- function(sp,
                     predictors = predictors,
                     pres_train = pres_train, #NEW
                     pres_test = pres_test, #NEW
                     backg_test = backg_test, #NEW
                     i = i, #NEW
                     models.dir = model.dir,
                     project.model = project.model,
                     projections = projections,
                     projdata = projdata,#um vector con nombres
                     mask = mask# a SpatialPolygonsDataFrame layer to mask and crop the predicted model
                     ){
    cat(paste("Mahalanobis distance",'\n'))
    ma <- mahal (predictors, pres_train)
    if (exists("ma")){
        ema <- evaluate(pres_test,backg_test,ma,predictors)
        thresholdma <- ema@t[which.max(ema@TPR + ema@TNR)]
        thma <- threshold(ema)
        ma_TSS <- max(ema@TPR + ema@TNR)-1
        ma_cont <- predict(ma, predictors,progress='text')
        ma_cont[ma_cont < threshold(ema,'no_omission')] <- threshold(ema,'no_omission')
        ma_bin <- ma_cont > thresholdma
        ma_cut <- ma_cont
        ma_cut[ma_cut < thresholdma] <- thresholdma
        if(minValue(ma_cut)<0) {
            ma_cut<-(ma_cut-minValue(ma_cut))/maxValue(ma_cut-minValue(ma_cut))}

        thma$AUC <- ema@auc
        thma$TSS <- ma_TSS
        thma$algoritmo <- "Mahal"
        thma$partition <- i
        row.names(thma) <- paste(sp,i,"Mahal")

        if (class(mask) == "SpatialPolygonsDataFrame"){
            source("./fct/cropModel.R")
            ma_cont <- cropModel(ma_cont , mask)
            ma_bin <- cropModel(ma_bin , mask)
            ma_cut <- cropModel(ma_cut , mask)
        }
        writeRaster(x=ma_cont,filename=paste0(models.dir,"/",sp,"/Mahal_cont_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
        writeRaster(x=ma_bin,filename=paste0(models.dir,"/",sp,"/Mahal_bin_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
        writeRaster(x=ma_cut,filename=paste0(models.dir,"/",sp,"/Mahal_cut_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")

        png(filename=paste0(models.dir,"/",sp,"/Mahal",sp,"_",i,"%03d.png"))
        plot(ma_cont,main=paste("Mahal raw","\n","AUC =", round(ema@auc,2),'-',"TSS =",round(ma_TSS,2)))
        plot(ma_bin,main=paste("Mahal P/A","\n","AUC =", round(ema@auc,2),'-',"TSS =",round(ma_TSS,2)))
        plot(ma_cut,main=paste("Mahal cut","\n","AUC =", round(ema@auc,2),'-',"TSS =",round(ma_TSS,2)))
        dev.off()

        if (project.model == T){
            for (proj in projections){
                data <- list.files(paste0("./env/",proj),pattern=proj)
                data2 <- stack(data)
                ma_proj <- predict(data2,ma,progress='text')
                ma_proj_bin <- ma_proj > thresholdma
                ma_proj_cut <- ma_proj_bin * ma_proj
                # Normaliza o modelo cut
                #ma_proj_cut <- ma_proj_cut/maxValue(ma_proj_cut)

                if (class(mask) == "SpatialPolygonsDataFrame"){
                    source("./fct/cropModel.R")
                    ma_proj <- cropModel(ma_proj , mask)
                    ma_proj_bin <- cropModel(ma_proj_bin , mask)
                    ma_proj_cut <- cropModel(ma_proj_cut , mask)
                }
                writeRaster(x=ma_proj,filename=paste0(models.dir,"/",sp,"/",proj,"/mahal_cont_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
                writeRaster(x=ma_proj_bin,filename=paste0(models.dir,"/",sp,"/",proj,"/mahal_bin_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
                writeRaster(x=ma_proj_cut,filename=paste0(models.dir,"/",sp,"/",proj,"/mahal_cut_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
                rm(data2)
            }
        }
    }
    else cat("Mahalanobis distance did not run")
    return(thma)
}

do_SVM2 <- function(sp,
                    predictors = predictors,
                    sdmdata_train = sdmdata_train,
                    envtest_pre = envtest_pre,
                    envtest_back = envtest_back,
                    i = i,
                    envtrain = envtrain,
                    models.dir = model.dir,
                    project.model = project.model,
                    projections = projections,
                    projdata = projdata,#um vector con nombres
                    mask = mask# a SpatialPolygonsDataFrame layer to mask and crop the predicted model
                    ){
    cat(paste("SVM2",'\n'))
    library(e1071)
    svm2 <- best.tune("svm",envtrain,sdmdata_train$pa,data=envtrain)##svm deve ser com a variável resposta binária ou contínua, eu acho que binária
    esvm2 <- evaluate(envtest_pre,envtest_back,svm2)
    #esvm <- evaluate(pres_test,backg_test,model = svm,x = predictors)
    svm2_TSS <- max(esvm2@TPR + esvm2@TNR)-1
    thresholdsvm2 <- esvm2@t[which.max(esvm2@TPR + esvm2@TNR)]
    thsvm2 <- threshold (esvm2)
    thsvm2$AUC <- esvm2@auc
    thsvm2$TSS <- svm2_TSS
    thsvm2$algoritmo <- "svm2"
    thsvm2$partition <- i
    row.names(thsvm2) <- paste(sp,i,"svm2")
    svm2_cont <- predict(predictors,svm2,progress='text')
    svm2_bin <- svm2_cont>thresholdsvm2
    svm2_cut <- svm2_bin * svm2_cont

    #TRANSFORMA 0 A 1
    svm2_cont <- svm2_cont/maxValue(svm2_cont)
    svm2_cut <- svm2_cut/maxValue(svm2_cut)
    if (class(mask) == "SpatialPolygonsDataFrame"){
        source("./fct/cropModel.R")
	    svm2_cont <- cropModel(svm2_cont , mask)
	    svm2_bin <- cropModel(svm2_bin , mask)
	    svm2_cut <- cropModel(svm2_cut , mask)
    }
    writeRaster(x=svm2_cont,filename=paste0(models.dir,"/",sp,"/svm2_cont_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
    writeRaster(x=svm2_bin,filename=paste0(models.dir,"/",sp,"/svm2_bin_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
    writeRaster(x=svm2_cut,filename=paste0(models.dir,"/",sp,"/svm2_cut_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")

    png(filename=paste0(models.dir,"/",sp,"/svm2",sp,"_",i,"%03d.png"))
    plot(svm2_cont,main=paste("SVM2 raw","\n","AUC =", round(esvm2@auc,2),'-',"TSS =",round(svm2_TSS,2)))
    plot(svm2_bin,main=paste("SVM2 P/A","\n","AUC =", round(esvm2@auc,2),'-',"TSS =",round(svm2_TSS,2)))
    plot(svm2_cut,main=paste("SVM2 cut","\n","AUC =", round(esvm2@auc,2),'-',"TSS =",round(svm2_TSS,2)))
    dev.off()

    if (project.model == T){
        for (proj in projections){
            data <- list.files(paste0("./env/",proj),pattern=proj)
            data2 <- stack(data)
            svm2_proj <- predict(data2,svm2,progress='text')
            svm2_proj_bin <- svm2_proj > thresholdsvm2
            svm2_proj_cut <- svm2_proj_bin * svm2_proj
            # Normaliza o modelo cut
            #svm2_proj_cut <- svm2_proj_cut/maxValue(svm2_proj_cut)
            if (class(mask) == "SpatialPolygonsDataFrame"){
                source("./fct/cropModel.R")
                svm2_proj <- cropModel(svm2_proj , mask)
		    svm2_proj_bin <- cropModel(svm2_proj_bin , mask)
		    svm2_proj_cut <- cropModel(svm2_proj_cut , mask)
            }
            writeRaster(x=svm2_proj,filename=paste0(models.dir,"/",sp,"/",proj,"/svm2_cont_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
            writeRaster(x=svm2_proj_bin,filename=paste0(models.dir,"/",sp,"/",proj,"/svm2_bin_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
            writeRaster(x=svm2_proj_cut,filename=paste0(models.dir,"/",sp,"/",proj,"/svm2_cut_",sp,"_",i,".tif"),overwrite=T, datatype="INT1U")
            rm(data2)
        }
    }
    return(thsvm2)
}

dismo.mod <- function(sp,
                      occs = spp.filt,#complete occurrence table
                      predictors = predictors,
                      buffer = TRUE,
                      buffer.type = "max",#"mean"
                      maxent = F,
                      Bioclim = T,
                      Domain = F,
                      Mahal = F,
                      GLM = F,
                      RF = F,
                      SVM = F,
                      SVM2 = F,
                      part = 3,
                      seed = NULL,#for reproducibility purposes
                      models.dir = "./models",
                      project.model = F,
                      projections = NULL,
                      projdata = NULL,#um vector con nombres
                      #stack_gcms = "future_vars", # Lista dos stacks de cada GCM. Ex: stack1 <- stack(variaveis_HADGEM); stack2<-stack(variaveis_CANESM); stack_gcms<-c(stack1,stack2)
                      mask = mask,# a SpatialPolygonsDataFrame layer to mask and crop the predicted model
                      n.back = 500){

  if (file.exists(paste0(models.dir)) == FALSE) dir.create(paste0(models.dir))
  if (file.exists(paste0(models.dir,"/",sp)) == FALSE) dir.create(paste0(models.dir,"/",sp))
  if (project.model == T) {
    for (proj in projections){
      if (file.exists(paste0(models.dir,"/",sp,"/",proj)) == FALSE) dir.create(paste0(models.dir,"/",sp,"/",proj))
    }
  }

  #Teste se pacotes estao instalados. Caso negativo, instala os pacotes e dependencias
  pacotes <- c('dismo','XML','raster','rgdal','maps','rgeos')
 # if (!pacotes %in% installed.packages()) install.packages(pacotes, dependencies = TRUE)
  #Loading packages
  lapply(pacotes, require, character.only = TRUE)


  print(date())

  cat(paste("Modeling",sp,"...",'\n'))
  #extrai as coordenadas de cada especie
  coord <- occs[occs$sp == sp,c('lon','lat')]
  n <- nrow(coord)
  n.back <- n * 2 #descomentado
  #tabela de valores
  presvals <- raster::extract(predictors, coord)

 if (buffer == FALSE) {
   set.seed(seed+2)
 backgr <- randomPoints(predictors, n.back)
 } else {
   source("./fct/createBuffer.R")
backgr <- createBuffer(coord, n.back, buffer.type)
  }

  colnames(backgr) <- c('lon', 'lat')

  #Extraindo dados ambientais dos bckgr
  backvals <- raster::extract(predictors, backgr)
  pa <- c(rep(1, nrow(presvals)), rep(0, nrow(backvals)))

  #Data partition
  if (n<11) part <- n else part  <- part
  set.seed(seed)#reproducibility
  group <- kfold(coord,part)
  set.seed(seed+1)
  bg.grp <- kfold(backgr,part)
  group.all <- c(group,bg.grp)

  pres <- cbind(coord,presvals)
  back <- cbind(backgr,backvals)
  rbind_1 <- rbind(pres,back)
  sdmdata <- data.frame(cbind(group.all,pa,rbind_1))
   rm(rbind_1)
   rm(pres)
   rm(back)
   gc()
  write.table(sdmdata,file = paste0(models.dir,"/",sp,"/sdmdata.txt"))


  #####Hace los modelos
  for (i in unique(group)){
    cat(paste(sp,"partition number",i,'\n'))
    pres_train <- coord[group != i, ]
    if(n == 1) pres_train <- coord[group == i,]
    pres_test  <-  coord[group  ==  i, ]

    backg_train <- backgr[bg.grp != i,]#not used?
    backg_test <- backgr[bg.grp == i,]#new

    sdmdata_train <- subset(sdmdata,group!=i)#new
    sdmdata_test <- subset(sdmdata,group == i)#new

    envtrain <- subset(sdmdata_train,select= c(-group,-lon,-lat))#new
    envtest <- subset(sdmdata_test,select=c(-group,-lon,-lat))
    envtest_pre <- subset(sdmdata_test,pa == 1,select= c(-group,-lon,-lat,-pa))#new
    envtest_back <- subset(sdmdata_test,pa == 0,select= c(-group,-lon,-lat,-pa))#new

    ##### Creates a .png plot of the initial dataset
    cat(paste("Plotting the dataset...",'\n'))
    extent <- extent(predictors)
    png(filename=paste0(models.dir,"/",sp,"/",i,sp,"dataset.png"))
    par(mfrow=c(1,1),mar=c(5,4,3,0))
    plot(predictors[[1]]!=0,col="grey95",main=paste(sp,"part.",i),legend=F)
    map('world',c('',"South America"),xlim=c(extent@xmin,extent@xmax),ylim=c(extent@ymin,extent@ymax),add=T)
    points(pres_train,pch=21,bg="red")
    points(pres_test,pch=21,bg="blue")
    points(backg_train,pch=3,col="red",cex=0.3)
    points(backg_test,pch=3,col="blue",cex=0.3)
    legend("topright",pch=c(21,21,3,3),pt.bg=c("red","blue"),legend=c("PresTrain","PresTest","BackTrain","BackTest"),col=c("black","black","red","blue"))
    dev.off()

    #SSB <- ssb(pres_test,backg_train,pres_train)
    #cat("Spatial sorting bias (approaches 0 with bias)")
    #print(SSB[,1]/SSB[,2])
    #faz os modelos

    cat(paste("Modeling",sp,"Partition",i,'\n'))
    eval_df <- data.frame("kappa"=1,"spec_sens"=1,"no_omission"=1,"prevalence"=1,"equal_sens_spec"=1,"sensitivity"=1,"AUC"=1,"TSS"=1,"algoritmo"="foo","partition"=1)

    if (Bioclim == T){
        thbc <- do_bioclim(sp = sp,
                           predictors = predictors,
                           pres_train = pres_train,
                           pres_test = pres_test,
                           backg_test = backg_test,
                           i = i,
                           models.dir = models.dir,
                           project.model = F,
                           projections = NULL,
                           projdata = NULL,
                           mask = NULL# a SpatialPolygonsDataFrame layer to mask and crop the predicted model
        )
        eval_df <- rbind(eval_df, thbc)
    }

    if (Domain == T){
        thdo <- do_domain(sp = sp,
                          predictors = predictors,
                          pres_train = pres_train,
                          pres_test = pres_test,
                          backg_test = backg_test,
                          i = i,
                          models.dir = models.dir,
                          project.model = F,
                          projections = NULL,
                          mask = NULL# a SpatialPolygonsDataFrame layer to mask and crop the predicted model
        )

        eval_df <- rbind(eval_df, thdo)
    }

    if (maxent == T){
        thmx <- do_maxent(sp = sp,
                          predictors = predictors,
                          pres_train = pres_train,
                          pres_test = pres_test,
                          backg_test = backg_test,
                          i = i,
                          models.dir = models.dir,
                          project.model = F,
                          projections = NULL,
                          mask = NULL# a SpatialPolygonsDataFrame layer to mask and crop the predicted model
        )
        eval_df <- rbind(eval_df, thmx)
    }

    if (Mahal == T){
        thma <- do_mahal(sp = sp,
                         predictors = predictors,
                         pres_train = pres_train,
                         pres_test = pres_test,
                         backg_test = backg_test,
                         i = i,
                         models.dir = models.dir,
                         project.model = F,
                         projections = NULL,
                         mask = NULL# a SpatialPolygonsDataFrame layer to mask and crop the predicted model
        )
        eval_df <- rbind(eval_df, thma)
    }

    if (GLM == T){
        thglm <- do_GLM (sp,
                         predictors = predictors,
                         sdmdata_train = sdmdata_train,
                         envtest_pre = envtest_pre,
                         envtest_back = envtest_back,
                         i = i,
                         envtrain = envtrain,
                         models.dir = models.dir,
                         project.model = F,
                         projections = NULL,
                         mask = NULL# a SpatialPolygonsDataFrame layer to mask and crop the predicted model
        )
        eval_df <- rbind(eval_df, thglm)
    }

    if (RF == T){
        thrf <- do_randomForest(sp,
                                predictors = predictors,
                                sdmdata_train = sdmdata_train,
                                envtest_pre = envtest_pre,
                                envtest_back = envtest_back,
                                i = i,
                                envtrain = envtrain,
                                models.dir = models.dir,
                                project.model = F,
                                projections = NULL,
                                mask = NULL# a SpatialPolygonsDataFrame layer to mask and crop the predicted model
        )
        eval_df <- rbind(eval_df, thrf)
    }

    if (SVM == T){
        thsvm <- do_SVM (sp,
                         predictors = predictors,
                         sdmdata_train = sdmdata_train,
                         envtest_pre = envtest_pre,
                         envtest_back = envtest_back,
                         i = i,
                         envtrain = envtrain,
                         part = part,
                         models.dir = models.dir,
                         project.model = F,
                         projections = NULL,
                         mask = NULL# a SpatialPolygonsDataFrame layer to mask and crop the predicted model
        )
        eval_df <- rbind(eval_df, thsvm)
    }

    if (SVM2 == T){
        thsvm2 <- do_SVM2(sp,
                          predictors = predictors,
                          sdmdata_train = sdmdata_train,
                          envtest_pre = envtest_pre,
                          envtest_back = envtest_back,
                          i = i,
                          envtrain = envtrain,
                          models.dir = models.dir,
                          project.model = F,
                          projections = NULL,
                          mask = NULL# a SpatialPolygonsDataFrame layer to mask and crop the predicted model
        )
        eval_df <- rbind(eval_df, thsvm2)
    }

    cat(paste("Saving the evaluation file...",sp,i,'\n'))
    write.table(eval_df[-1,],file = paste0(models.dir,"/",sp,"/evaluate",sp,"_",i,".txt"))
    }

   #cat(paste("Saving the evaluation file...",sp,i,'\n'))
   #write.table(eval,file = paste0(models.dir,"/",sp,"/evaluate",sp,"_",i,".txt"))
    cat(as.character(sp), "DONE on", date(), '\n')
    print(date())
  }
