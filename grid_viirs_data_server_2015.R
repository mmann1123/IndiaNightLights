##This script reprojects VIIRS DNB and Cloud Mask data to grid centered on Maharashtra, India
## it saves output at lists stored in RData format (to avoid errors with writing files) 

# to install rhdf5
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
library('sp')
library('rgdal')
library('foreach')
library('iterators')
library('doParallel')

registerDoParallel(32)

setwd("/groups/manngroup/India VIIRS/2015")
 
# read in list of files and set up iteration groups
d = list.files(path=getwd(),pattern=glob2rx("*h5"),full.names=T,include.dirs=T)
iterator = split(1:length(d), cut(1:length(d),10))

#iterate through smaller groups
for(j in 1:length(iterator)){
    
    #get fresh list 
    d <- list.files(path=getwd(),pattern=glob2rx("*h5"),full.names=T,include.dirs=T)
    print(d[iterator[[j]]])
    
    # limit to the first group
    d = d[iterator[[j]]]
    
    # store names for later 
    name_date_time = substr(d[i],36,47)
    
    #For each 5-min swath file (DNB/CMASK) in directory (HDF5s)
    # Write out DNB band data -------------------------------------------
    
    output <- foreach(i = 1:length(d), .inorder=FALSE,.packages =c('rhdf5','raster')) %dopar% {
        print(i)  
        fname = d[i]
           
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
        a = rasterize(x=data,y=example,field='dnb',fun='mean',background=NA)
        a@data@names = name_date_time
        a
    }
    
    save(output, file=paste(getwd(),'/job_dnb_', j,'_v3.RData',sep=""))
    remove(output)
    
    
    # read in  cloud band data -------------------------------------------
    
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
        example =raster(matrix(NA,nrow=dim(lat_cmask)[1],ncol=dim(lat_cmask)[2]),xmn=xmin,xmx=xmax,
    	ymn=ymin,ymx=ymax,crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
        
         if(length(lon_cmask)!=length(lat_cmask)| length(lon_cmask) !=length(cmask)){
    	#special case for non matching data lengths
    	a=c("this file had different lenghts for lat, lon and cmask")
    	write.csv(a,paste(getwd(),'//',
    	            name_date_time,'_cld_v2_FAILED_DIF_LENGTHS.csv',sep=""))
    	return(NA)
    	}else{
    
        data = data.frame(lon=as.numeric(lon_cmask), lat = as.numeric(lat_cmask),cloud=as.numeric(cmask))
        coordinates(data) =~lon+lat
        cloud = rasterize(x=data,y=example,field='cloud',fun='last',background=NA)
        cloud@data@names = name_date_time
        cloud
    	}
       }
    
    save(output2, file=paste(getwd(),'/job_cld_', j,'_v3.RData',sep=""))
    remove(output2)
    
    
    output3 <- foreach(i = 1:length(d), .inorder=FALSE,.packages =c('rhdf5','raster')) %dopar% {
        print(i)
        fname = d[i]
        
        zenith = h5read(fname,'/LunarZenith')
        lon_cmask = h5read(fname,'/Longitude')
        lat_cmask = h5read(fname,'/Latitude')    
        H5close()
        
        xmin = min(lon_cmask)
        ymin = min(lat_cmask)
        xmax = max(lon_cmask)
        ymax = max(lat_cmask)
        example = raster(matrix(NA,nrow=dim(lat_cmask)[1],ncol=dim(lat_cmask)[2]),xmn=xmin,xmx=xmax,
                        ymn=ymin,ymx=ymax,crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
        
        if(length(lon_cmask)!=length(lat_cmask)| length(lon_cmask) !=length(zenith)){
            #special case for non matching data lengths
            a=c("this file had different lenghts for lat, lon and zenith")
            write.csv(a,paste(getwd(),'//',
                    name_date_time,'_zen_v3_FAILED_DIF_LENGTHS.csv',sep=""))
            return(NA)
        }else{
            data = data.frame(lon=as.numeric(lon_cmask), lat = as.numeric(lat_cmask),zenith=as.numeric(zenith))
            coordinates(data) =~lon+lat
            zenith = rasterize(x=data,y=example,field='zenith',fun='last',background=NA)
            zenith@data@names = name_date_time
            zenith
        }
    }
    
    save(output3, file=paste(getwd(),'/job_zen_', j,'_v3.RData',sep=""))
    remove(output3)
    
    
    
    output4 <- foreach(i = 1:length(d), .inorder=FALSE,.packages =c('rhdf5','raster')) %dopar% {
        print(i)
        fname = d[i]
        
        azimuth = h5read(fname,'/LunarAzimuth')
        lon_cmask = h5read(fname,'/Longitude')
        lat_cmask = h5read(fname,'/Latitude')    
        H5close()
        
        xmin = min(lon_cmask)
        ymin = min(lat_cmask)
        xmax = max(lon_cmask)
        ymax = max(lat_cmask)
        example = raster(matrix(NA,nrow=dim(lat_cmask)[1],ncol=dim(lat_cmask)[2]),xmn=xmin,xmx=xmax,
                         ymn=ymin,ymx=ymax,crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
        
        if(length(lon_cmask)!=length(lat_cmask)| length(lon_cmask) !=length(azimuth)){
            #special case for non matching data lengths
            a=c("this file had different lenghts for lat, lon and azimuth")
            write.csv(a,paste(getwd(),'//',
                        name_date_time,'_azt_v3_FAILED_DIF_LENGTHS.csv',sep=""))
            return(NA)
        }else{
            data = data.frame(lon=as.numeric(lon_cmask), lat = as.numeric(lat_cmask),azimuth=as.numeric(azimuth))
            coordinates(data) =~lon+lat
            azimuth = rasterize(x=data,y=example,field='azimuth',fun='last',background=NA)
            azimuth@data@names = name_date_time
            azimuth
        }
    }
    
    save(output4, file=paste(getwd(),'/job_zen_', j,'_v3.RData',sep=""))
    remove(output4)

}
