##This script reprojects VIIRS DNB and Cloud Mask data to common
##30 arc-second grid centered on Maharashtra, India

library('rhdf5')
library('raster')
#library('R.utils')
library('sp')
library('rgdal')
library('foreach')
library('iterators')
library('doParallel')

#Register the parallel backend
#registerDoParallel(4)


# boundary = readOGR('C://Users//mmann//Google Drive//India Night Lights//Data//Administrative Borders','IND_adm1')
# plot(boundary[boundary$NAME_1=='Maharashtra',])
# boundary = boundary[boundary$NAME_1=='Maharashtra',]

#Set up 30 arc second grid centered on Maharashtra
grid <- raster(nrow=3072,ncol=4064,xmn=72,xmx=81,ymn=15,ymx=23)
res(grid) <- c(30/3600,30/3600)
temp <- seq(1,ncell(grid),1)
grid <- setValues(grid,temp)
remove(temp)


setwd('R://Mann Research//India Night Time Lights//Test Case//')

d <- list.files(path=getwd(),pattern=glob2rx("*h5"),full.names=T,include.dirs=T)

#For each 5-min swath file (DNB/CMASK) in directory (HDF5s)

for (i in 34:length(d)){
    print(i)
  
    #skip   A2012133.1935.h5  not working
    #if (i ==34 |i==35){next}
  
    fname <- d[i]
    h5ls(fname)
    
    ### Grid DNB Radiance ###
    
    #Read in HDF5 files
    dnb <- h5read(fname,'/Radiance')
    lon_dnb <- h5read(fname,'/Longitude')
    lat_dnb <- h5read(fname,'/Latitude')
    #sza_dnb <- h5read(fname3,'/SolarZenithAngle')
    H5close()
    
    #Resize matrices to single column vectors
    dim(lon_dnb) <- c(dim(lon_dnb)[1]*dim(lon_dnb)[2],1)
    dim(lat_dnb) <- c(dim(lat_dnb)[1]*dim(lat_dnb)[2],1)
    dim(dnb) <- c(dim(dnb)[1]*dim(dnb)[2],1)
    
    #This is in case we want to look at a smaller subset than the grid
    #defined above. If not, then use lat>=15, lat<=23, lon>=72, lon<=81.
    # Jumda BBox westlimit=76.9279; southlimit=19.8964; eastlimit=77.1818; northlimit=20.1037
    w <- which(lat_dnb>=19.5 & lat_dnb<=20.5 & lon_dnb>=76.5 & lon_dnb<=77.5)
    #w <- which(lat_dnb>=18.5 & lat_dnb<=20.5 & lon_dnb>=76 & lon_dnb<=78)
    
    lon_dnb <- lon_dnb[w]
    lat_dnb <- lat_dnb[w]
    dnb <- dnb[w]
    
    # if image doesn't cover area of interest move to next image
    if(length(dnb)==0){
      print('skipping')
      next}
    
    #Set up data frame with latitude, longitude, and DNB radiance for each pixel in 
    #the subset
    pts=data.frame(cbind(lat_dnb,lon_dnb,dnb))
     
    #Loop through each pixel in the DNB matrix and extract the location of the 
    #nearest grid cell
    pts_extr <- foreach(j = 1:length(w), .combine = rbind) %dopar% {
      #for (j in 1:length(w)){
         print(j)
         pts <- data.frame(lat_dnb[j],lon_dnb[j],dnb[j])
         colnames(pts) <- c('Lat','Lon','Radiance')
         coordinates(pts) = ~ Lon+Lat
         proj4string(pts) = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 ")
          
         pt_extr <- extract(grid,pts)
    }
  
    #Generate new matrix and assign dnb values to designated grid cells (pts_extr)
    pts_grid <- matrix(-9999,ncell(grid),1)
    pts_grid[pts_extr,1] <- dnb
    
    #Create new DNB raster
    grid_new <- setValues(grid,pts_grid)
      
    writeRaster(grid_new,filename=paste(getwd(),'//',
      substr(d[i],51,64),'.tif',sep=""),format='GTiff',overwrite=TRUE)
    
    # Grid cloud mask ---------------------------------------------------------
    
    cmask <- h5read(fname,'/CloudMask')
    lon_cmask <- h5read(fname,'/Longitude2')
    lat_cmask <- h5read(fname,'/Latitude2')
    
    dim(lon_cmask) <- c(dim(lon_cmask)[1]*dim(lon_cmask)[2],1)
    dim(lat_cmask) <- c(dim(lat_cmask)[1]*dim(lat_cmask)[2],1)
    dim(cmask) <- c(dim(cmask)[1]*dim(cmask)[2],1)
    
    w <- which(lat_cmask>=19.5 & lat_cmask<=20.5 & lon_cmask>=76.5 & lon_cmask<=77.5)
    lon_cmask <- lon_cmask[w]
    lat_cmask <- lat_cmask[w]
    cmask <- cmask[w]
    
    pts_extr <- foreach(j = 1:length(w), .combine = rbind) %dopar% {
        #for ( j in 1:length(w)){ 
        print(j)
        
        pts <- data.frame(lat_cmask[j],lon_cmask[j],cmask[j])
        colnames(pts) <- c('Lat','Lon','QF')
        coordinates(pts) <- ~Lon+Lat
        proj4string(pts) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
        
        pt_extr <- extract(grid,pts)
    }
    
    pts_grid <- matrix(-9999,ncell(grid),1)
    pts_grid[pts_extr,1] <- cmask
      
    grid_new <- setValues(grid,pts_grid)
    
    writeRaster(grid_new,filename=paste(getwd(),
        substr(d[i],51,64),'_cmask.tif',sep=""),format='GTiff',overwrite=TRUE)
    
    remove(grid_new,grid,pts,lat_cmask,lat_dnb,lon_cmask, long_dnb)
    #save(pts_extr,cmask,file='pts_extr')
  
}





