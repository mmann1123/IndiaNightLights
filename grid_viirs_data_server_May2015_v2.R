##This script reprojects VIIRS DNB and Cloud Mask data to common
##30 arc-second grid centered on Maharashtra, India

#source("http://bioconductor.org/biocLite.R")
#biocLite("rhdf5")


# Run the following in bash before starting R
#module load proj.4/4.8.0
#module load gdal
#module load R/3.0.2
#module load gcc/4.9.0
# R


library('rhdf5')
library('raster')
#library('R.utils')
library('sp')
library('rgdal')
library('foreach')
library('iterators')
library('doParallel')

#Register the parallel backend
registerDoParallel(32)



setwd("/groups/manngroup/India VIIRS/May2015")


d <- list.files(path=getwd(),pattern=glob2rx("*h5"),full.names=T,include.dirs=T)

#For each 5-min swath file (DNB/CMASK) in directory (HDF5s)
 
output <- foreach(i = 1:length(d),.packages =c('rhdf5','raster')) %dopar% {
    print(i)  
   
    fname <- d[i]
   # h5ls(fname)
    
    ### Grid DNB Radiance ###
    #http://neondataskills.org/HDF5/Create-Raster-Stack-Spectroscopy-HDF5-In-R/
        
    #Read in HDF5 files
    dnb <- h5read(fname,'/Radiance')
    lon_dnb <- h5read(fname,'/Longitude')
    lat_dnb <- h5read(fname,'/Latitude')
    #sza_dnb <- h5read(fname3,'/SolarZenithAngle')
    H5close()
    
    xmin = min(lon_dnb)
    ymin = min(lat_dnb)
    xmax = max(lon_dnb)
    ymax = max(lat_dnb)
    
    example =raster(matrix(NA,nrow=dim(lat_dnb)[1],ncol=dim(lat_dnb)[2]), xmn=xmin, xmx=xmax, 
	ymn=ymin, ymx=ymax, crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 ") )
    data = data.frame(lon=as.numeric(lon_dnb), lat = as.numeric(lat_dnb),dnb=as.numeric(dnb))
    data[data==-1.5e-09,]=NA
    coordinates(data) =~lon+lat
    a = rasterize(x=data,y=example,field='dnb',fun='last',background=NA)
    a
    #plot(a)
   
}
stopImplicitCluster()

#windows()
#plot(log(output[[1]]*1e9))

lapply(1:length(output),function(x) writeRaster(output[[x]], filename=paste(getwd(),'//',
                 substr(d[x],39,49),'_dnb_v2.tif',sep=""),format='GTiff',overwrite=TRUE) )
   
remove(output)
    
registerDoParallel(32)

output2 <- foreach(i = 1:length(d),.packages =c('rhdf5','raster')) %dopar% {
    
    fname <- d[i]
    
    cmask <- h5read(fname,'/CloudMask')
    lon_cmask <- h5read(fname,'/Longitude2')
    lat_cmask <- h5read(fname,'/Latitude2')
    
    H5close()
    
    xmin = min(lon_cmask)
    ymin = min(lat_cmask)
    xmax = max(lon_cmask)
    ymax = max(lat_cmask)
    
    example =raster(matrix(NA,nrow=dim(lat_cmask)[1],ncol=dim(lat_cmask)[2]), xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax, crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 ") )
    data = data.frame(lon=as.numeric(lon_cmask), lat = as.numeric(lat_cmask),cloud=as.numeric(cmask))
    coordinates(data) =~lon+lat
    cloud = rasterize(x=data,y=example,field='cloud',fun='last',background=NA)
    cloud
}

stopImplicitCluster()

#windows()
#plot(log(output2[[1]]*1e9))

lapply(1:length(output2),function(x) writeRaster(output2[[x]], filename=paste(getwd(),'//',
               substr(d[x],39,49),'_cld_v2.tif',sep=""),format='GTiff',overwrite=TRUE))
remove(output2)
