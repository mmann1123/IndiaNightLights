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
registerDoParallel(4)



setwd("C://Users/mmann/Desktop/NightTimeData/")

# read in list of files and set up iteration groups
d <- list.files(path=getwd(),pattern=glob2rx("*h5"),full.names=T,include.dirs=T)
iterator = split(1:length(d), cut(1:length(d),10))

#iterate through smaller groups
for(j in 1:length(iterator)){

#get fresh list 
d <- list.files(path=getwd(),pattern=glob2rx("*h5"),full.names=T,include.dirs=T)
print(d[iterator[[j]]])

d = d[iterator[[j]]]

#For each 5-min swath file (DNB/CMASK) in directory (HDF5s)
 
output <- foreach(i = 1:length(d), .inorder=FALSE,.packages =c('rhdf5','raster')) %dopar% {
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
}

windows()
plot(log(output[[1]]*1e9))

save(output, file=paste(getwd(),'/job_output1_', j,'_v2.RData',sep=""))


#lapply(1:length(output),function(x) if( class(output[[x]])=='RasterLayer'){ 
#	print(x)
#	writeRaster(output[[x]], filename=paste(getwd(),'//',
#        substr(d[x],36,47),'_dnb_v2.tif',sep=""),format='GTiff',overwrite=TRUE)})
   
#remove(output)
    

output2 <- foreach(i = 1:length(d), .inorder=FALSE,.packages =c('rhdf5','raster')) %dopar% {
    print(i)
    fname <- d[i]
    
    cmask <- h5read(fname,'/CloudMask')
    lon_cmask <- h5read(fname,'/Longitude2')
    lat_cmask <- h5read(fname,'/Latitude2')    
    H5close()
    
    xmin = min(lon_cmask)
    ymin = min(lat_cmask)
    xmax = max(lon_cmask)
    ymax = max(lat_cmask)
    example =raster(matrix(NA,nrow=dim(lat_cmask)[1],ncol=dim(lat_cmask)[2]), xmn=xmin, 
	xmx=xmax, ymn=ymin, ymx=ymax, crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 ") )
    
     if(length(lon_cmask)!=length(lat_cmask)| length(lon_cmask) !=length(cmask)){
	#special case for non matching data lengths
	a=c("this file had different lenghts for lat, lon and cmask")
	write.csv(a,paste(getwd(),'//',
               substr(d[i],40,51),'_cld_v2_FAILED_DIF_LENGTHS.csv',sep=""))
	return(NA)
	}else{

    data = data.frame(lon=as.numeric(lon_cmask), lat = as.numeric(lat_cmask),cloud=as.numeric(cmask))
    coordinates(data) =~lon+lat
    cloud = rasterize(x=data,y=example,field='cloud',fun='last',background=NA)
    cloud
	}
}


save(output2, file=paste(getwd(),'/job_output2_', j,'_v2.RData',sep=""))

#stopImplicitCluster()

#windows()
#plot(log(output2[[1]]*1e9))

#lapply(1:length(output2),function(x) if(  class(output2[[x]])=='RasterLayer'){ 
#	print(x)
#	writeRaster(output2[[x]], filename=paste(getwd(),'//',
#        substr(d[x],36,47),'_cld_v2.tif',sep=""),format='GTiff',overwrite=TRUE)})
#remove(output2)


}
