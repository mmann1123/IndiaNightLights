##This script reprojects VIIRS DNB and Cloud Mask data to grid centered on Maharashtra, India
## it saves output at lists stored in RData format (to avoid errors with writing files) 

# to install rhdf5
#source("http://bioconductor.org/biocLite.R")
#biocLite("rhdf5")


# Run the following in bash before starting R
#module load proj.4/4.8.0
#module load gdal/gcc/1.11
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

setwd("/groups/manngroup/India VIIRS/2015")   # all May2015 .h5 files moved to /2015/
 
# read in list of files and set up iteration groups
d = list.files(path=getwd(),pattern=glob2rx("*h5"),full.names=T,include.dirs=T)
iterator = split(1:length(d), cut(1:length(d),13))   #11


################################################
# parameter setup 

# check output names
start_c = 36
end_c = 47
substr(d[1],start_c,end_c)  # check here

# version number of output name
version = 'v5'


################################################
# New rasterize technique
e1 = extent(72, 81.50, 15, 22.5)

#figure out number of rows and columns to get close to 0.00675 resolution
res = 0.00784                     # 0.00675
cols = round((e1[2]-e1[1])/res,0)
rows = round((e1[4]-e1[3])/res,0)

# set a raster to resample to with Maharastra extent
#   originally was nrows=1000, ncols=1000
example = raster(e1, nrows=rows, ncols=cols, crs=CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'))

source('/groups/manngroup/India VIIRS/IndiaNightLights/test_intersection.R')


#iterate through smaller groups
for(j in 1:length(iterator)){
    
    #get fresh list 
    d <- list.files(path=getwd(),pattern=glob2rx("*h5"),full.names=T,include.dirs=T)
    print(d[iterator[[j]]])
    
    # limit to the first group
    d = d[iterator[[j]]]
    
    #For each 5-min swath file (DNB/CMASK) in directory (HDF5s)
    # Write out DNB band data -------------------------------------------
    
    output <- foreach(i = 1:length(d), .inorder=FALSE,.packages =c('rhdf5','raster')) %dopar% {
        print(i)  
        fname = d[i]

        # store names for later
        name_date_time = substr(d[i],start_c,end_c)
           
        ### Grid DNB Radiance ###
        #http://neondataskills.org/HDF5/Create-Raster-Stack-Spectroscopy-HDF5-In-R/
            
        #Read in HDF5 files
        dnb <- h5read(fname,'/Radiance')
        dnb[,1:224] = NA  # remove high noise cells file:///C:/Users/mmann/Downloads/8-133-1-PB%20(4).pdf
	columns = dim(dnb)[2]
	dnb[,(columns-224):columns] = NA
        dnb[dnb==-1.5e-09]=NA  # remove missing data verified with eli

	lon_dnb <- h5read(fname,'/Longitude')
        lat_dnb <- h5read(fname,'/Latitude')
        H5close()
        
        xmin = min(lon_dnb)
        ymin = min(lat_dnb)
        xmax = max(lon_dnb)
        ymax = max(lat_dnb)
        
        #example =raster(matrix(NA,nrow=dim(lat_dnb)[1],ncol=dim(lat_dnb)[2]), xmn=xmin, xmx=xmax, 
        #	ymn=ymin, ymx=ymax, crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 ") )
        data = data.frame(lon=as.numeric(lon_dnb), lat = as.numeric(lat_dnb),dnb=as.numeric(dnb))
       
        coordinates(data) =~lon+lat
        #a = rasterize(x=data,y=example,field='dnb',fun=mean,background=NA,na.rm=T)
        #a@data@names = name_date_time
        #a
	
	if( test_intersection(data,example) ){
		a = rasterize(x=data,y=example,field='dnb',fun=mean,background=NA,na.rm=T)
        	a@data@names = name_date_time
        	a
	}


    }
    
    save(output, file=paste(getwd(),'/job_dnb_', j,'_',version,'.RData',sep=""))
    remove(output)
    
    
    # read in  cloud band data -------------------------------------------
    
  output2 <- foreach(i = 1:length(d), .inorder=FALSE,.packages =c('rhdf5','raster')) %dopar% {
        print(i)
        fname <- d[i]
        
        # store names for later
        name_date_time = substr(d[i],start_c,end_c)

        cmask <- h5read(fname,'/CloudMask')
        lon_cmask <- h5read(fname,'/Longitude2')
        lat_cmask <- h5read(fname,'/Latitude2')    
        H5close()
        
        xmin = min(lon_cmask)
        ymin = min(lat_cmask)
        xmax = max(lon_cmask)
        ymax = max(lat_cmask)
        #example =raster(matrix(NA,nrow=dim(lat_cmask)[1],ncol=dim(lat_cmask)[2]),xmn=xmin,xmx=xmax,
    	#ymn=ymin,ymx=ymax,crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
        
        if(length(lon_cmask)!=length(lat_cmask)| length(lon_cmask) !=length(cmask)){
    	   #special case for non matching data lengths
    	   a=c("this file had different lenghts for lat, lon and cmask")
    	   write.csv(a,paste(getwd(),'//',
    	            name_date_time,'_cld_v2_FAILED_DIF_LENGTHS.csv',sep=""))
    	   return(NA)
    	}else{
           data = data.frame(lon=as.numeric(lon_cmask), lat = as.numeric(lat_cmask),cloud=as.numeric(cmask))
           coordinates(data) =~lon+lat
        
	   #cloud = rasterize(x=data,y=example,field='cloud',fun='max',background=NA,na.rm=T)
           #cloud@data@names = name_date_time
           #cloud

             if( test_intersection(data,example) ){
                cloud = rasterize(x=data,y=example,field='cloud',fun='max',background=NA,na.rm=T)
	        # deal with courser resolution
		cloud = focal(cloud, w=matrix(1,3,3),fun=function(x){max(x,na.rm=T)})
		cloud@data@names = name_date_time
		cloud
             }

    	  }
       }
    
    save(output2, file=paste(getwd(),'/job_cld_', j,'_',version,'.RData',sep=""))
    remove(output2)
    
    

    output3 <- foreach(i = 1:length(d), .inorder=FALSE,.packages =c('rhdf5','raster')) %dopar% {
        print(i)
        fname = d[i]

        # store names for later
        name_date_time = substr(d[i],start_c,end_c)
        
        zenith = h5read(fname,'/LunarZenith')
        lon_cmask = h5read(fname,'/Longitude')
        lat_cmask = h5read(fname,'/Latitude')    
        H5close()
        
        xmin = min(lon_cmask)
        ymin = min(lat_cmask)
        xmax = max(lon_cmask)
        ymax = max(lat_cmask)
        #example = raster(matrix(NA,nrow=dim(lat_cmask)[1],ncol=dim(lat_cmask)[2]),xmn=xmin,xmx=xmax,
        #                ymn=ymin,ymx=ymax,crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
        
        if(length(lon_cmask)!=length(lat_cmask)| length(lon_cmask) !=length(zenith)){
            #special case for non matching data lengths
            a=c("this file had different lenghts for lat, lon and zenith")
            write.csv(a,paste(getwd(),'//',
                    name_date_time,'_zen_v3_FAILED_DIF_LENGTHS.csv',sep=""))
            return(NA)
        }else{
            data = data.frame(lon=as.numeric(lon_cmask), lat = as.numeric(lat_cmask),zenith=as.numeric(zenith))
            coordinates(data) =~lon+lat
            #zenith = rasterize(x=data,y=example,field='zenith',fun='last',background=NA,na.rm=T)
            #zenith@data@names = name_date_time
            #zenith
	     if( test_intersection(data,example) ){
                zenith = rasterize(x=data,y=example,field='zenith',fun='last',background=NA,na.rm=T)
                zenith@data@names = name_date_time
                zenith
             }

        }
    }
    
    save(output3, file=paste(getwd(),'/job_zen_', j,'_',version,'.RData',sep=""))
    remove(output3)
    
    
    
    output4 <- foreach(i = 1:length(d), .inorder=FALSE,.packages =c('rhdf5','raster')) %dopar% {
        print(i)
        fname = d[i]
       
	# store names for later
        name_date_time = substr(d[i],start_c,end_c)

        azimuth = h5read(fname,'/LunarAzimuth')
        lon_cmask = h5read(fname,'/Longitude')
        lat_cmask = h5read(fname,'/Latitude')    
        H5close()
        
        xmin = min(lon_cmask)
        ymin = min(lat_cmask)
        xmax = max(lon_cmask)
        ymax = max(lat_cmask)
        #example = raster(matrix(NA,nrow=dim(lat_cmask)[1],ncol=dim(lat_cmask)[2]),xmn=xmin,xmx=xmax,
        #                 ymn=ymin,ymx=ymax,crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
        
        if(length(lon_cmask)!=length(lat_cmask)| length(lon_cmask) !=length(azimuth)){
            #special case for non matching data lengths
            a=c("this file had different lenghts for lat, lon and azimuth")
            write.csv(a,paste(getwd(),'//',
                        name_date_time,'_azt_v3_FAILED_DIF_LENGTHS.csv',sep=""))
            return(NA)
        }else{
            data = data.frame(lon=as.numeric(lon_cmask), lat = as.numeric(lat_cmask),azimuth=as.numeric(azimuth))
            coordinates(data) =~lon+lat
            #azimuth = rasterize(x=data,y=example,field='azimuth',fun='last',background=NA,na.rm=T)
            #azimuth@data@names = name_date_time
            #azimuth

	    if( test_intersection(data,example) ){
                azimuth = rasterize(x=data,y=example,field='azimuth',fun='last',background=NA,na.rm=T)
                azimuth@data@names = name_date_time
                azimuth
             }

        }
    }
    
    save(output4, file=paste(getwd(),'/job_azt_', j,'_',version,'.RData',sep=""))
    remove(output4)
}


# use Write_dnb_out.R to write these saved files to .tif
