Climate‐driven shifts in the distribution of koala browse species. doi: [10.1111/ecog.04530].

Input Layer Preparation:


setwd('C:/Users/....')

library('rgdal')
library('rgeos')
library('maptools')
library('raster')
library('sp')
library('dismo')

Aus <- readShapeSpatial('layers/Climatic/Input_variables/Shape/Australia_Dissolved.shp')
proj4string(Aus)<- CRS("++proj=longlat +datum=WGS84")

#------Background points-----------

bkgr <- randomPoints(current, 10000)
bkgr <- data.frame(bkgr)
bkgr <- data.frame(Longitude=bkgr$x, Latitude=bkgr$y)
write.csv(bkgr,'layers/Climatic/Input_variables/CSV/bkgr.csv')
bkgr <- SpatialPointsDataFrame(bkgr[,1:2], bkgr, proj4string = CRS("++proj=longlat +datum=WGS84"))
writeOGR(bkgr, dsn='layers/Climatic/Input_variables/Shape', layer="bkgr", driver="ESRI Shapefile")

#------Current-----------

current = list.files(pattern='asc$', path='layers/Climatic/Input_variables/Current', full.name=T)
current = stack(current)

crop <- crop(current, extent(Aus), snap="out")
plot(crop,1)
raster <- rasterize(Aus, crop)   
raster <- mask(x=crop, mask=raster)
plot(raster,1)
writeRaster(raster, filename=names(crop), bylayer=TRUE, format="GTiff")

#---------Past----------

cclgmbi = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Past/cclgmbi/', full.name=T)
cclgmbi = stack(cclgmbi)
crop1 <- crop(cclgmbi, extent(Aus), snap="out")
raster1 <- mask(x=crop1, mask=raster)
plot(raster1,1)
writeRaster(raster1, filename=names(crop1), bylayer=TRUE, format="GTiff")

melgmbi = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Past/melgmbi', full.name=T)
melgmbi = stack(melgmbi)
crop2 <- crop(melgmbi, extent(Aus), snap="out")
raster2 <- mask(x=crop2, mask=raster)
plot(raster2,1)
writeRaster(raster2, filename=names(crop2), bylayer=TRUE, format="GTiff")

mrlgmbi = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Past/mrlgmbi', full.name=T)
mrlgmbi = stack(mrlgmbi)
crop3 <- crop(mrlgmbi, extent(Aus), snap="out")
raster3 <- mask(x=crop3, mask=raster)
plot(raster3,1)
writeRaster(raster3, filename=names(crop3), bylayer=TRUE, format="GTiff")

ccmidbi = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Past/ccmidbi', full.name=T)
ccmidbi = stack(ccmidbi)
crop4 <- crop(ccmidbi, extent(Aus), snap="out")
raster4 <- mask(x=crop4, mask=raster)
plot(raster4,1)
writeRaster(raster4, filename=names(crop4), bylayer=TRUE, format="GTiff")

memidbi = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Past/memidbi', full.name=T)
memidbi = stack(memidbi)
crop5 <- crop(memidbi, extent(Aus), snap="out")
raster5 <- mask(x=crop5, mask=raster)
plot(raster5,1)
writeRaster(raster5, filename=names(crop5), bylayer=TRUE, format="GTiff")

mrmidbi = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Past/mrmidbi', full.name=T)
mrmidbi = stack(mrmidbi)
crop6 <- crop(mrmidbi, extent(Aus), snap="out")
raster6 <- mask(x=crop6, mask=raster)
plot(raster6,1)
writeRaster(raster6, filename=names(crop6), bylayer=TRUE, format="GTiff")

lig = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Past/lig_30s_bio', full.name=T)
lig = stack(lig)
crop7 <- crop(lig, extent(Aus), snap="out")
raster7 <- mask(x=crop7, mask=raster)
plot(crop7,1)
writeRaster(raster7, filename=names(crop7), bylayer=TRUE, format="GTiff")

bio5 <- raster('layers/Climatic/Input_variables/Past/lig_30s_bio/lig_30s_bio_5.bil')
writeRaster(bio5, 'layers/Climatic/Input_variables/Past/lig_30s_bio/bio5.tif')


#---------Furure 2050----------

cc2650 = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Future/2050/cc26bi50', full.name=T)
cc2650 = stack(cc2650)
crop8 <- crop(cc2650, extent(Aus), snap="out")
raster8 <- mask(x=crop8, mask=raster)
plot(raster8,1)
writeRaster(raster8, filename=names(crop8), bylayer=TRUE, format="GTiff")

cc4550 = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Future/2050/cc45bi50', full.name=T)
cc4550 = stack(cc4550)
crop9 <- crop(cc4550, extent(Aus), snap="out")
raster9 <- mask(x=crop9, mask=raster)
plot(raster9,1)
writeRaster(raster9, filename=names(crop9), bylayer=TRUE, format="GTiff")

cc6050 = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Future/2050/cc60bi50', full.name=T)
cc6050 = stack(cc6050)
crop10 <- crop(cc6050, extent(Aus), snap="out")
raster10 <- mask(x=crop10, mask=raster)
plot(raster10,1)
writeRaster(raster10, filename=names(crop10), bylayer=TRUE, format="GTiff")

cc8550 = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Future/2050/cc85bi50', full.name=T)
cc8550 = stack(cc8550)
crop11 <- crop(cc8550, extent(Aus), snap="out")
raster11 <- mask(x=crop11, mask=raster)
plot(raster11,1)
writeRaster(raster11, filename=names(crop11), bylayer=TRUE, format="GTiff")

mr2650 = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Future/2050/mr26bi50', full.name=T)
mr2650 = stack(mr2650)
crop12 <- crop(mr2650, extent(Aus), snap="out")
raster12 <- mask(x=crop12, mask=raster)
plot(raster12,1)
writeRaster(raster12, filename=names(crop12), bylayer=TRUE, format="GTiff")

mr4550 = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Future/2050/mr45bi50', full.name=T)
mr4550 = stack(mr4550)
crop13 <- crop(mr4550, extent(Aus), snap="out")
raster13 <- mask(x=crop13, mask=raster)
plot(raster13,1)
writeRaster(raster13, filename=names(crop13), bylayer=TRUE, format="GTiff")

mr6050 = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Future/2050/mr60bi50', full.name=T)
mr6050 = stack(mr6050)
crop14 <- crop(mr6050, extent(Aus), snap="out")
raster14 <- mask(x=crop14, mask=raster)
plot(raster14,1)
writeRaster(raster14, filename=names(crop14), bylayer=TRUE, format="GTiff")

mr8550 = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Future/2050/mr85bi50', full.name=T)
mr8550 = stack(mr8550)
crop15 <- crop(mr8550, extent(Aus), snap="out")
raster15 <- mask(x=crop15, mask=raster)
plot(raster15,1)
writeRaster(raster15, filename=names(crop15), bylayer=TRUE, format="GTiff")

#---------Furure 2070----------

cc2670 = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Future/2070/cc26bi70', full.name=T)
cc2670 = stack(cc2670)
crop16 <- crop(cc2670, extent(Aus), snap="out")
raster16 <- mask(x=crop16, mask=raster)
plot(raster16,1)
writeRaster(raster16, filename=names(crop16), bylayer=TRUE, format="GTiff")

cc4570 = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Future/2070/cc45bi70', full.name=T)
cc4570 = stack(cc4570)
crop17 <- crop(cc4570, extent(Aus), snap="out")
raster17 <- mask(x=crop17, mask=raster)
plot(raster17,1)
writeRaster(raster17, filename=names(crop17), bylayer=TRUE, format="GTiff")

cc6070 = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Future/2070/cc60bi70', full.name=T)
cc6070 = stack(cc6070)
crop18 <- crop(cc6070, extent(Aus), snap="out")
raster18 <- mask(x=crop18, mask=raster)
plot(raster18,1)
writeRaster(raster18, filename=names(crop18), bylayer=TRUE, format="GTiff")

cc8570 = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Future/2070/cc85bi70', full.name=T)
cc8570 = stack(cc8570)
crop19 <- crop(cc8570, extent(Aus), snap="out")
raster19 <- mask(x=crop19, mask=raster)
plot(raster19,1)
writeRaster(raster19, filename=names(crop19), bylayer=TRUE, format="GTiff")

mr2670 = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Future/2070/mr26bi70', full.name=T)
mr2670 = stack(mr2670)
crop20 <- crop(mr2670, extent(Aus), snap="out")
raster20 <- mask(x=crop20, mask=raster)
plot(raster20,1)
writeRaster(raster20, filename=names(crop20), bylayer=TRUE, format="GTiff")

mr4570 = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Future/2070/mr45bi70', full.name=T)
mr4570 = stack(mr4570)
crop21 <- crop(mr4570, extent(Aus), snap="out")
raster21 <- mask(x=crop21, mask=raster)
plot(raster21,1)
writeRaster(raster21, filename=names(crop21), bylayer=TRUE, format="GTiff")

mr6070 = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Future/2070/mr60bi70', full.name=T)
mr6070 = stack(mr6070)
crop22 <- crop(mr6070, extent(Aus), snap="out")
raster22 <- mask(x=crop22, mask=raster)
plot(raster22,1)
writeRaster(raster22, filename=names(crop22), bylayer=TRUE, format="GTiff")

mr8570 = list.files(pattern='tif$', path='layers/Climatic/Input_variables/Future/2070/mr85bi70', full.name=T)
mr8570 = stack(mr8570)
crop23<- crop(mr8570, extent(Aus), snap="out")
raster23 <- mask(x=crop23, mask=raster)
plot(raster23,1)
writeRaster(raster23, filename=names(crop23), bylayer=TRUE, format="GTiff")

#---------SDMs----------

setwd('C:/Users/....')

library('rgdal')
library('rgeos')
library('maptools')
library('raster')
library('sp')
library('dismo')
library('biomod2')


sp <- read.csv('SP1.csv')

sp <- data.frame(Longitude=sp$decimalLongitude, Latitude=sp$decimalLatitude, Species="SP1")
sp1 <- na.omit(sp)
sp2 <- SpatialPointsDataFrame(sp1[,1:2],sp1, proj4string = CRS("++proj=longlat +datum=WGS84"))
writeOGR(sp2, dsn='Shape', layer="SP1", driver="ESRI Shapefile")

Aus <- readShapeSpatial('Shape/Australia_Dissolved.shp')
proj4string(Aus)<- CRS("++proj=longlat +datum=WGS84")

Prs <- readShapeSpatial('Shape/SP1.shp')
proj4string(Prs)<- CRS("++proj=longlat +datum=WGS84")

Prs_Aus <- crop(Prs, extent(Aus))

SP_rmvdup <- remove.duplicates(Prs_Aus, zero = 8)
writeOGR(SP_rmvdup, dsn='Shape', layer="SP1_rmvdup", driver="ESRI Shapefile")
plot(SP_rmvdup)

bkgr <- readShapeSpatial('Shape/bkgr.shp')
proj4string(Prs)<- CRS("++proj=longlat +datum=WGS84")

df_Prs <- cbind(Longitude=Prs_Aus$Longitude, Latitude=Prs_Aus$Latitude, "SP1"=rep(1, nrow(Prs_Aus)))
df_bkgr <- cbind(Longitude=bkgr$Longitude, Latitude=bkgr$Latitude, "SP1"=rep(0, nrow(bkgr)))
sdm <- rbind(df_Prs, df_bkgr)
head(sdm)

#------Layers-----------

current = list.files(pattern='tif$', path='Current/AUS', full.name=T)
current = stack(current)
current <- dropLayer(current, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(current) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

cclgmbi = list.files(pattern='tif$', path='Past/cclgmbi/AUS', full.name=T)
cclgmbi = stack(cclgmbi)
cclgmbi <- dropLayer(cclgmbi, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(cclgmbi) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

melgmbi = list.files(pattern='tif$', path='Past/melgmbi/AUS', full.name=T)
melgmbi = stack(melgmbi)
melgmbi <- dropLayer(melgmbi, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(melgmbi) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

mrlgmbi = list.files(pattern='tif$', path='Past/mrlgmbi/AUS', full.name=T)
mrlgmbi = stack(mrlgmbi)
mrlgmbi <- dropLayer(mrlgmbi, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(mrlgmbi) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

ccmidbi = list.files(pattern='tif$', path='Past/ccmidbi/AUS/', full.name=T)
ccmidbi = stack(ccmidbi)
ccmidbi <- dropLayer(ccmidbi, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(ccmidbi) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

memidbi = list.files(pattern='tif$', path='Past/memidbi/AUS', full.name=T)
memidbi = stack(memidbi)
memidbi <- dropLayer(memidbi, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(memidbi) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

mrmidbi = list.files(pattern='tif$', path='Past/mrmidbi/AUS', full.name=T)
mrmidbi = stack(mrmidbi)
mrmidbi <- dropLayer(mrmidbi, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(mrmidbi) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

lig = list.files(pattern='tif$', path='Past/lig_30s_bio/AUS', full.name=T)
lig = stack(lig)
lig <- dropLayer(lig, c(2,3,5,7,8,9,10,12,13))
names(lig) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

#---------Furure 2050----------

cc2650 = list.files(pattern='tif$', path='Future/2050/cc26bi50/AUS', full.name=T)
cc2650 = stack(cc2650)
cc2650 <- dropLayer(cc2650, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(cc2650) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

cc4550 = list.files(pattern='tif$', path='Future/2050/cc45bi50/AUS', full.name=T)
cc4550 = stack(cc4550)
cc4550 <- dropLayer(cc4550, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(cc4550) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

cc6050 = list.files(pattern='tif$', path='Future/2050/cc60bi50/AUS', full.name=T)
cc6050 = stack(cc6050)
cc6050 <- dropLayer(cc6050, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(cc6050) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

cc8550 = list.files(pattern='tif$', path='Future/2050/cc85bi50/AUS', full.name=T)
cc8550 = stack(cc8550)
cc8550 <- dropLayer(cc8550, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(cc8550) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

mr2650 = list.files(pattern='tif$', path='Future/2050/mr26bi50/AUS', full.name=T)
mr2650 = stack(mr2650)
mr2650 <- dropLayer(mr2650, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(mr2650) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

mr4550 = list.files(pattern='tif$', path='Future/2050/mr45bi50/AUS', full.name=T)
mr4550 = stack(mr4550)
mr4550 <- dropLayer(mr4550, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(mr4550) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

mr6050 = list.files(pattern='tif$', path='Future/2050/mr60bi50/AUS', full.name=T)
mr6050 = stack(mr6050)
mr6050 <- dropLayer(mr6050, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(mr6050) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

mr8550 = list.files(pattern='tif$', path='Future/2050/mr85bi50/AUS', full.name=T)
mr8550 = stack(mr8550)
mr8550 <- dropLayer(mr8550, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(mr8550) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

#---------Furure 2070----------

cc2670 = list.files(pattern='tif$', path='Future/2070/cc26bi70/AUS', full.name=T)
cc2670 = stack(cc2670)
cc2670 <- dropLayer(cc2670, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(cc2670) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

cc4570 = list.files(pattern='tif$', path='Future/2070/cc45bi70/AUS', full.name=T)
cc4570 = stack(cc4570)
cc4570 <- dropLayer(cc4570, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(cc4570) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

cc6070 = list.files(pattern='tif$', path='Future/2070/cc60bi70/AUS', full.name=T)
cc6070 = stack(cc6070)
cc6070 <- dropLayer(cc6070, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(cc6070) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

cc8570 = list.files(pattern='tif$', path='Future/2070/cc85bi70/AUS', full.name=T)
cc8570 = stack(cc8570)
cc8570 <- dropLayer(cc8570, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(cc8570) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

mr2670 = list.files(pattern='tif$', path='Future/2070/mr26bi70/AUS', full.name=T)
mr2670 = stack(mr2670)
mr2670 <- dropLayer(mr2670, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(mr2670) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

mr4570 = list.files(pattern='tif$', path='Future/2070/mr45bi70/AUS', full.name=T)
mr4570 = stack(mr4570)
mr4570 <- dropLayer(mr4570, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(mr4570) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

mr6070 = list.files(pattern='tif$', path='Future/2070/mr60bi70/AUS', full.name=T)
mr6070 = stack(mr6070)
mr6070 <- dropLayer(mr6070, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(mr6070) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')

mr8570 = list.files(pattern='tif$', path='Future/2070/mr85bi70/AUS', full.name=T)
mr8570 = stack(mr8570)
mr8570 <- dropLayer(mr8570, c(2,3,5,7,8,9,10,12,13,16,17,18,19))
names(mr8570) <- c('bio1','bio12','bio14','bio19','bio4','bio5','tclay')


###---------Biomod------
myBiomodOption <- BIOMOD_ModelingOptions()
myRespName <- 'SP1'
DataSpecies <- sdm
head(sdm)
myResp <- as.numeric(DataSpecies[,myRespName])
myRespXY <- DataSpecies[,c("Longitude","Latitude")]
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp, expl.var = current, resp.xy = myRespXY, resp.name = myRespName)

myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,models = c('X','Y','Z'),  models.options = myBiomodOption,  NbRunEval=XX,  DataSplit=75,  Prevalence=0.5,  VarImport=3,  models.eval.meth = c('TSS','ROC'),  SaveObj = F,  rescal.all.models= T,  do.full.models = F,  modeling.id = paste(myRespName,"FirstModeling",sep=""))
myBiomodModelEval <- getModelsEvaluations(myBiomodModelOut)
dimnames(myBiomodModelEval)
get_evaluations(myBiomodModelOut)
getModelsVarImport(myBiomodModelOut)
write.csv(get_evaluations(myBiomodModelOut), 'SP1/AUC.TSS_SP1.csv')
write.csv(getModelsVarImport(myBiomodModelOut), 'SP1/Var.imp_SP1.csv')


Projection_current <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                        new.env = current,
                                        proj.name = 'current',
                                        selected.models = 'all',
                                        binary.meth = 'TSS',
                                        compress = FALSE,
                                        build.clamping.mask = FALSE,
                                        output.format = '.img')

Projection_cclgmbi <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                         new.env = cclgmbi,
                                         proj.name = 'cclgmbi',
                                         selected.models = 'all',
                                         binary.meth = 'TSS',
                                         compress = FALSE,
                                         build.clamping.mask = FALSE,
                                         output.format = '.img')

Projection_melgmbi <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                        new.env = melgmbi,
                                        proj.name = 'melgmbi',
                                        selected.models = 'all',
                                        binary.meth = 'TSS',
                                        compress = FALSE,
                                        build.clamping.mask = FALSE,
                                        output.format = '.img')

Projection_mrlgmbi <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                        new.env = mrlgmbi,
                                        proj.name = 'mrlgmbi',
                                        selected.models = 'all',
                                        binary.meth = 'TSS',
                                        compress = FALSE,
                                        build.clamping.mask = FALSE,
                                        output.format = '.img')

Projection_ccmidbi <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                        new.env = ccmidbi,
                                        proj.name = 'ccmidbi',
                                        selected.models = 'all',
                                        binary.meth = 'TSS',
                                        compress = FALSE,
                                        build.clamping.mask = FALSE,
                                        output.format = '.img')

Projection_memidbi <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                        new.env = memidbi,
                                        proj.name = 'memidbi',
                                        selected.models = 'all',
                                        binary.meth = 'TSS',
                                        compress = FALSE,
                                        build.clamping.mask = FALSE,
                                        output.format = '.img')

Projection_mrmidbi <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                        new.env = mrmidbi,
                                        proj.name = 'mrmidbi',
                                        selected.models = 'all',
                                        binary.meth = 'TSS',
                                        compress = FALSE,
                                        build.clamping.mask = FALSE,
                                        output.format = '.img')

Projection_lig <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                        new.env = lig,
                                        proj.name = 'lig',
                                        selected.models = 'all',
                                        binary.meth = 'TSS',
                                        compress = FALSE,
                                        build.clamping.mask = FALSE,
                                        output.format = '.img')

Projection_cc2650 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                    new.env = cc2650,
                                    proj.name = 'cc2650',
                                    selected.models = 'all',
                                    binary.meth = 'TSS',
                                    compress = FALSE,
                                    build.clamping.mask = FALSE,
                                    output.format = '.img')

Projection_cc4550 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                       new.env = cc4550,
                                       proj.name = 'cc4550',
                                       selected.models = 'all',
                                       binary.meth = 'TSS',
                                       compress = FALSE,
                                       build.clamping.mask = FALSE,
                                       output.format = '.img')

Projection_cc6050 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                       new.env = cc6050,
                                       proj.name = 'cc6050',
                                       selected.models = 'all',
                                       binary.meth = 'TSS',
                                       compress = FALSE,
                                       build.clamping.mask = FALSE,
                                       output.format = '.img')

Projection_cc8550 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                       new.env = cc8550,
                                       proj.name = 'cc8550',
                                       selected.models = 'all',
                                       binary.meth = 'TSS',
                                       compress = FALSE,
                                       build.clamping.mask = FALSE,
                                       output.format = '.img')

Projection_mr2650 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                       new.env = mr2650,
                                       proj.name = 'mr2650',
                                       selected.models = 'all',
                                       binary.meth = 'TSS',
                                       compress = FALSE,
                                       build.clamping.mask = FALSE,
                                       output.format = '.img')

Projection_mr4550 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                       new.env = mr4550,
                                       proj.name = 'mr4550',
                                       selected.models = 'all',
                                       binary.meth = 'TSS',
                                       compress = FALSE,
                                       build.clamping.mask = FALSE,
                                       output.format = '.img')

Projection_mr6050 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                       new.env = mr6050,
                                       proj.name = 'mr6050',
                                       selected.models = 'all',
                                       binary.meth = 'TSS',
                                       compress = FALSE,
                                       build.clamping.mask = FALSE,
                                       output.format = '.img')

Projection_mr8550 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                       new.env = mr8550,
                                       proj.name = 'mr8550',
                                       selected.models = 'all',
                                       binary.meth = 'TSS',
                                       compress = FALSE,
                                       build.clamping.mask = FALSE,
                                       output.format = '.img')

Projection_cc2670 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                       new.env = cc2670,
                                       proj.name = 'cc2670',
                                       selected.models = 'all',
                                       binary.meth = 'TSS',
                                       compress = FALSE,
                                       build.clamping.mask = FALSE,
                                       output.format = '.img')

Projection_cc4570 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                       new.env = cc4570,
                                       proj.name = 'cc4570',
                                       selected.models = 'all',
                                       binary.meth = 'TSS',
                                       compress = FALSE,
                                       build.clamping.mask = FALSE,
                                       output.format = '.img')

Projection_cc6070 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                       new.env = cc6070,
                                       proj.name = 'cc6070',
                                       selected.models = 'all',
                                       binary.meth = 'TSS',
                                       compress = FALSE,
                                       build.clamping.mask = FALSE,
                                       output.format = '.img')

Projection_cc8570 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                       new.env = cc8570,
                                       proj.name = 'cc8570',
                                       selected.models = 'all',
                                       binary.meth = 'TSS',
                                       compress = FALSE,
                                       build.clamping.mask = FALSE,
                                       output.format = '.img')

Projection_mr2670 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                       new.env = mr2670,
                                       proj.name = 'mr2670',
                                       selected.models = 'all',
                                       binary.meth = 'TSS',
                                       compress = FALSE,
                                       build.clamping.mask = FALSE,
                                       output.format = '.img')

Projection_mr4570 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                       new.env = mr4570,
                                       proj.name = 'mr4570',
                                       selected.models = 'all',
                                       binary.meth = 'TSS',
                                       compress = FALSE,
                                       build.clamping.mask = FALSE,
                                       output.format = '.img')

Projection_mr6070 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                       new.env = mr6070,
                                       proj.name = 'mr6070',
                                       selected.models = 'all',
                                       binary.meth = 'TSS',
                                       compress = FALSE,
                                       build.clamping.mask = FALSE,
                                       output.format = '.img')

Projection_mr8570 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                       new.env = mr8570,
                                       proj.name = 'mr8570',
                                       selected.models = 'all',
                                       binary.meth = 'TSS',
                                       compress = FALSE,
                                       build.clamping.mask = FALSE,
                                       output.format = '.img')


myBiomodEM <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut,
                                       chosen.models = 'all',
                                       em.by = 'all',
                                       eval.metric = c('ROC'),
                                       eval.metric.quality.threshold = NULL,
                                       models.eval.meth = c('ROC'),
                                       prob.mean = TRUE,
                                       prob.cv = FALSE,
                                       prob.ci = FALSE,
                                       prob.ci.alpha = 0.05,
                                       prob.median = FALSE,
                                       committee.averaging = FALSE,
                                       prob.mean.weight = FALSE,
                                       prob.mean.weight.decay = 'proportional')

ensm <- BIOMOD_EnsembleForecasting( projection.output = Projection_current,
                                    EM.output = myBiomodEM, output.format = '.img')

ensm_cclgmbi <- BIOMOD_EnsembleForecasting( projection.output = Projection_cclgmbi,
                                            EM.output = myBiomodEM, output.format = '.img')

ensm_melgmbi <- BIOMOD_EnsembleForecasting( projection.output = Projection_melgmbi,
                                            EM.output = myBiomodEM, output.format = '.img')

ensm_mrlgmbi <- BIOMOD_EnsembleForecasting( projection.output = Projection_mrlgmbi,
                                            EM.output = myBiomodEM, output.format = '.img')

ensm_ccmidbi <- BIOMOD_EnsembleForecasting( projection.output = Projection_ccmidbi,
                                            EM.output = myBiomodEM, output.format = '.img')

ensm_memidbi <- BIOMOD_EnsembleForecasting( projection.output = Projection_memidbi,
                                            EM.output = myBiomodEM, output.format = '.img')

ensm_mrmidbi <- BIOMOD_EnsembleForecasting( projection.output = Projection_mrmidbi,
                                            EM.output = myBiomodEM, output.format = '.img')

ensm_lig <- BIOMOD_EnsembleForecasting( projection.output = Projection_lig,
                                        EM.output = myBiomodEM, output.format = '.img')

ensm_cc2650 <- BIOMOD_EnsembleForecasting( projection.output = Projection_cc2650,
                                           EM.output = myBiomodEM, output.format = '.img')

ensm_cc4550 <- BIOMOD_EnsembleForecasting( projection.output = Projection_cc4550,
                                           EM.output = myBiomodEM, output.format = '.img')

ensm_cc6050 <- BIOMOD_EnsembleForecasting( projection.output = Projection_cc6050,
                                           EM.output = myBiomodEM, output.format = '.img')

ensm_cc8550 <- BIOMOD_EnsembleForecasting( projection.output = Projection_cc8550,
                                           EM.output = myBiomodEM, output.format = '.img')

ensm_cc2670 <- BIOMOD_EnsembleForecasting( projection.output = Projection_cc2670,
                                           EM.output = myBiomodEM, output.format = '.img')

ensm_cc4570 <- BIOMOD_EnsembleForecasting( projection.output = Projection_cc4570,
                                           EM.output = myBiomodEM, output.format = '.img')

ensm_cc6070 <- BIOMOD_EnsembleForecasting( projection.output = Projection_cc6070,
                                           EM.output = myBiomodEM, output.format = '.img')

ensm_cc8570 <- BIOMOD_EnsembleForecasting( projection.output = Projection_cc8570,
                                           EM.output = myBiomodEM, output.format = '.img')

ensm_mr2650 <- BIOMOD_EnsembleForecasting( projection.output = Projection_mr2650,
                                           EM.output = myBiomodEM, output.format = '.img')

ensm_mr4550 <- BIOMOD_EnsembleForecasting( projection.output = Projection_mr4550,
                                           EM.output = myBiomodEM, output.format = '.img')

ensm_mr6050 <- BIOMOD_EnsembleForecasting( projection.output = Projection_mr6050,
                                           EM.output = myBiomodEM, output.format = '.img')

ensm_mr8550 <- BIOMOD_EnsembleForecasting( projection.output = Projection_mr8550,
                                           EM.output = myBiomodEM, output.format = '.img')

ensm_mr2670 <- BIOMOD_EnsembleForecasting( projection.output = Projection_mr2670,
                                           EM.output = myBiomodEM, output.format = '.img')

ensm_mr4570 <- BIOMOD_EnsembleForecasting( projection.output = Projection_mr4570,
                                           EM.output = myBiomodEM, output.format = '.img')

ensm_mr6070 <- BIOMOD_EnsembleForecasting( projection.output = Projection_mr6070,
                                           EM.output = myBiomodEM, output.format = '.img')

ensm_mr8570 <- BIOMOD_EnsembleForecasting( projection.output = Projection_mr8570,
                                           EM.output = myBiomodEM, output.format = '.img')


current <- raster('SP1/proj_current/proj_current_SP1_ensemble.img')
cur_val <- extract(current, Prs)
threshold10 <- quantile(cur_val, probs=0.1, na.rm=T)
current_bin <- reclassify(current, c(0,threshold10,0,  threshold10,Inf, 1))
writeRaster(current_bin, 'SP1/Current.bin_SP1.tif')

cclgm <- raster('SP1/proj_cclgmbi/proj_cclgmbi_SP1_Ensemble.img')
melgm <- raster('SP1/proj_melgmbi/proj_melgmbi_SP1_Ensemble.img')
mrlgm <- raster('SP1/proj_mrlgmbi/proj_mrlgmbi_SP1_Ensemble.img')
Ens_lgm <- (cclgm + melgm + mrlgm)/3
writeRaster(Ens_lgm, 'SP1/Ens_lgm_SP1.tif')
Ens_lgm_bin <- reclassify(Ens_lgm, c(0,threshold10,0,  threshold10,Inf, 1))
writeRaster(Ens_lgm_bin, 'SP1/LGM_bin_SP1.tif')


ccmid <- raster('SP1/proj_ccmidbi/proj_ccmidbi_SP1_Ensemble.img')
memid <- raster('SP1/proj_memidbi/proj_memidbi_SP1_Ensemble.img')
mrmid <- raster('SP1/proj_mrmidbi/proj_mrmidbi_SP1_Ensemble.img')
Ens_mid <- (ccmid + memid + mrmid)/3
writeRaster(Ens_mid, 'SP1/Ens_mid_E_camaldulensis.tif')
Ens_mid_bin <- reclassify(Ens_mid, c(0,threshold10,0,  threshold10,Inf, 1))
writeRaster(Ens_mid_bin, 'SP1/MID_bin_SP1.tif')

lig <- raster('SP1/proj_lig/proj_lig_SP1_Ensemble.img')
lig_bin <- reclassify(lig, c(0,threshold10,0,  threshold10,Inf, 1))
writeRaster(lig_bin, 'SP1/LIG_bin_SP1.tif')

cc26.50 <- raster('SP1/proj_cc2650/proj_cc2650_SP1_Ensemble.img')
mr26.50 <- raster('SP1/proj_mr2650/proj_mr2650_SP1_Ensemble.img')
Ens26.50 <- (cc26.50 + mr26.50)/2
writeRaster(Ens26.50, 'SP1/Ens26.50_SP1.tif')
Ens26.50_bin <- reclassify(Ens26.50, c(0,threshold10,0,  threshold10,Inf, 1))
writeRaster(Ens26.50_bin, 'SP1/Ens26.50_bin_SP1.tif')

cc45.50 <- raster('SP1/proj_cc4550/proj_cc4550_SP1_Ensemble.img')
mr45.50 <- raster('SP1/proj_mr4550/proj_mr4550_SP1_Ensemble.img')
Ens45.50 <- (cc45.50 + mr45.50)/2
writeRaster(Ens45.50, 'SP1/Ens45.50_SP1.tif')
Ens45.50_bin <- reclassify(Ens45.50, c(0,threshold10,0,  threshold10,Inf, 1))
writeRaster(Ens45.50_bin, 'SP1/Ens45.50_bin_SP1.tif')

cc60.50 <- raster('SP1/proj_cc6050/proj_cc6050_SP1_Ensemble.img')
mr60.50 <- raster('SP1/proj_mr6050/proj_mr6050_SP1_Ensemble.img')
Ens60.50 <- (cc60.50 + mr60.50)/2
writeRaster(Ens60.50, 'SP1/Ens60.50_SP1.tif')
Ens60.50_bin <- reclassify(Ens60.50, c(0,threshold10,0,  threshold10,Inf, 1))
writeRaster(Ens60.50_bin, 'SP1/Ens60.50_bin_SP1.tif')

cc85.50 <- raster('SP1/proj_cc8550/proj_cc8550_SP1_Ensemble.img')
mr85.50 <- raster('SP1/proj_mr8550/proj_mr8550_SP1_Ensemble.img')
Ens85.50 <- (cc85.50 + mr85.50)/2
writeRaster(Ens85.50, 'SP1/Ens85.50_SP1.tif')
Ens85.50_bin <- reclassify(Ens85.50, c(0,threshold10,0,  threshold10,Inf, 1))
writeRaster(Ens85.50_bin, 'SP1/Ens85.50_bin_SP1.tif')

cc26.70 <- raster('SP1/proj_cc2670/proj_cc2670_SP1_Ensemble.img')
mr26.70 <- raster('SP1/proj_mr2670/proj_mr2670_SP1_Ensemble.img')
Ens26.70 <- (cc26.70 + mr26.70)/2
writeRaster(Ens26.70, 'SP1/Ens26.70_SP1.tif')
Ens26.70_bin <- reclassify(Ens26.70, c(0,threshold10,0,  threshold10,Inf, 1))
writeRaster(Ens26.70_bin, 'SP1/Ens26.70_bin_SP1.tif')

cc45.70 <- raster('SP1/proj_cc4570/proj_cc4570_SP1_Ensemble.img')
mr45.70 <- raster('SP1/proj_mr4570/proj_mr4570_SP1_Ensemble.img')
Ens45.70 <- (cc45.70 + mr45.70)/2
writeRaster(Ens45.70, 'SP1/Ens45.70_SP1.tif')
Ens45.70_bin <- reclassify(Ens45.70, c(0,threshold10,0,  threshold10,Inf, 1))
writeRaster(Ens45.70_bin, 'SP1/Ens45.70_bin_SP1.tif')

cc60.70 <- raster('SP1/proj_cc6070/proj_cc6070_SP1_Ensemble.img')
mr60.70 <- raster('SP1/proj_mr6070/proj_mr6070_SP1_Ensemble.img')
Ens60.70 <- (cc60.70 + mr60.70)/2
writeRaster(Ens60.70, 'SP1/Ens60.70_SP1.tif')
Ens60.70_bin <- reclassify(Ens60.70, c(0,threshold10,0,  threshold10,Inf, 1))
writeRaster(Ens60.70_bin, 'SP1/Ens60.70_bin_SP1.tif')

cc85.70 <- raster('SP1/proj_cc8570/proj_cc8570_SP1_Ensemble.img')
mr85.70 <- raster('SP1/proj_mr8570/proj_mr8570_SP1_Ensemble.img')
Ens85.70 <- (cc85.70 + mr85.70)/2
writeRaster(Ens85.70, 'SP1/Ens85.70_SP1.tif')
Ens85.70_bin <- reclassify(Ens85.70, c(0,threshold10,0,  threshold10,Inf, 1))
writeRaster(Ens85.70_bin, 'SP1/Ens85.70_bin_SP1.tif')


#------Range size----

RangeSize.lgm <- BIOMOD_RangeSize(current_bin, Ens_lgm_bin)
RangeSize.mid <- BIOMOD_RangeSize(current_bin, Ens_mid_bin)
RangeSize.lig <- BIOMOD_RangeSize(current_bin, lig_bin)

RangeSize26.50 <- BIOMOD_RangeSize(current_bin, Ens26.50_bin)
RangeSize45.50 <- BIOMOD_RangeSize(current_bin, Ens45.50_bin)
RangeSize60.50 <- BIOMOD_RangeSize(current_bin, Ens60.50_bin)
RangeSize85.50 <- BIOMOD_RangeSize(current_bin, Ens85.50_bin)

RangeSize26.70 <- BIOMOD_RangeSize(current_bin, Ens26.70_bin)
RangeSize45.70 <- BIOMOD_RangeSize(current_bin, Ens45.70_bin)
RangeSize60.70 <- BIOMOD_RangeSize(current_bin, Ens60.70_bin)
RangeSize85.70 <- BIOMOD_RangeSize(current_bin, Ens85.70_bin)

RangeSize <- rbind(RangeSize.lgm$Compt.By.Models, RangeSize.mid$Compt.By.Models,
                   RangeSize.lig$Compt.By.Models, RangeSize26.50$Compt.By.Models,
                   RangeSize45.50$Compt.By.Models, RangeSize60.50$Compt.By.Models,
                   RangeSize85.50$Compt.By.Models, RangeSize26.70$Compt.By.Models,
                   RangeSize45.70$Compt.By.Models, RangeSize60.70$Compt.By.Models, RangeSize85.70$Compt.By.Models)

row.names(RangeSize) <- c("Current_LGM", "Current_Mid", "Current_Lig", "Current_26.50",
                          "Current_45.50", "Current_60.50", "Current_85.50", "Current_26.70",
                          "Current_45.70", "Current_60.70", "Current_85.70")
write.csv(RangeSize, 'SP1/RangeSize_SP1.csv')

writeRaster(RangeSize.lgm$Diff.By.Pixel, 'SP1/RangeSize.lgm.tif')
writeRaster(RangeSize.mid$Diff.By.Pixel, 'SP1/RangeSize.mid.tif')
writeRaster(RangeSize.lig$Diff.By.Pixel, 'SP1/RangeSize.lig.tif')
writeRaster(RangeSize26.50$Diff.By.Pixel, 'SP1/RangeSize26.50.tif')
writeRaster(RangeSize45.50$Diff.By.Pixel, 'SP1/RangeSize45.50.tif')
writeRaster(RangeSize60.50$Diff.By.Pixel, 'SP1/RangeSize60.50.tif')
writeRaster(RangeSize85.50$Diff.By.Pixel, 'SP1/RangeSize85.50.tif')
writeRaster(RangeSize26.70$Diff.By.Pixel, 'SP1/RangeSize26.70.tif')
writeRaster(RangeSize45.70$Diff.By.Pixel, 'SP1/RangeSize45.70.tif')
writeRaster(RangeSize60.70$Diff.By.Pixel, 'SP1/RangeSize60.70.tif')
writeRaster(RangeSize85.70$Diff.By.Pixel, 'SP1/RangeSize85.70.tif')
