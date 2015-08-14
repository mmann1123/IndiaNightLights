# Michael Mann
# this script reads in data created by grid_viirs_data_server_2015.R
# it takes RData files created on the server and writes out tif files
# Limit extent to Maharashstra

library(raster)


setwd("C:\\Users\\mmann\\Desktop\\NightTimeData\\May2015/")

# mararashtra extent westlimit=71.99; southlimit=15.03; eastlimit=81.51; northlimit=22.5
# nice ideas on extent and resolution http://stackoverflow.com/questions/20733555/how-to-create-a-raster-brick-with-rasters-of-different-extents
e1 = extent(72, 81.50, 15, 22.5)
# set a raster to resample to 
s<-raster(e1, nrows=1000, ncols=1000, crs=CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'))

test_intersection <- function(a,b){
    #reads in two rasters and tests for overlap T or F
    # if returns TRUE then there is overlap
    !class(try(intersect(a,b),T ))=='try-error'
}


# Write out cloud band data -------------------------------------------

output_list=list.files(getwd(),pattern='output2')

for(i in 1:length(output_list)){
    load(paste(getwd(),'/',output_list[i],sep=''))
    
    if(length(output2)==0 ){ # check forempty lists
        warning(paste('error in ',output_list[i])); next}
   
     lapply(1:length(output2),function(x) if(class(output2[[x]])=='RasterLayer' & test_intersection(output2[[x]],s) ){  # avoid empty and non overlapping images
        print(x)
        output2[[x]]=resample(output2[[x]], s, method="ngb") # limit to new extent and res
        if(x==1){windows(); image(output2[[x]])}
        writeRaster(output2[[x]], filename=paste(getwd(),'/',
                output2[[x]]@data@names,'_cld_v3.tif',sep=""),format='GTiff',overwrite=TRUE)})
}    

# download job_output2_5_v3.RData  and correct output there was an error


# Write out day night band data -------------------------------------------

output_list=list.files(getwd(),pattern='output1')

for(i in 1:length(output_list)){
    load(paste(getwd(),'/',output_list[i],sep=''))
    
    if(length(output)==0){warning(paste('error in ',output_list[i]))
        next}
    lapply(1:length(output),function(x) if(class(output[[x]])=='RasterLayer'& test_intersection(output[[x]],s) ){ 
        print(x)
        if(x==1){ windows(); image(output[[x]])}
        
        output[[x]]=resample(output[[x]], s, method="bilinear")
        
        writeRaster(output[[x]], filename=paste(getwd(),'/',
            output[[x]]@data@names,'_dnb_v3.tif',sep=""),format='GTiff',overwrite=TRUE)})
}    


# # code to change names etc
# output_list=list.files(getwd(),pattern='X')
# 
# for(i in 1:length(output_list)){
#     file.rename(output_list[i],substr(output_list[i],2,30) )  #paste('2',substr(output_list[i],3,24),sep='')
# }


