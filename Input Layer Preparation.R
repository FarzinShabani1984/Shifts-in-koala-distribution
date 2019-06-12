
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

