
library(raster)

setwd("C:/Users/mmann/Google Drive/India Night Lights/Data/2015")
output2=list.files(getwd(),pattern='output2')
load(paste(getwd(),'/',output2,sep=''))


lapply(1:length(output2),function(x) if(class(output2[[x]])=='RasterLayer'){ 
    print(x)
    windows()
    image(output2[[x]])
    writeRaster(output2[[x]], filename=paste(getwd(),'//',
    substr(d[x],36,47),'_dnb_v2.tif',sep=""),format='GTiff',overwrite=TRUE)})

