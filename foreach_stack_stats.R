
require(raster)
require(doParallel)


# register the parallel backend
registerDoParallel(6)

setwd('C:\\Users\\mmann\\Desktop\\NightTimeData\\')
load('dnb_stack_wo_cld.RData')

stack_stat_mean = function(stack_rows){
    stack_rows = log(stack_rows*1e9+1)
    unlist(lapply( 1:dim(stack_rows)[1], function(i)  mean(stack_rows[i,],na.rm=T) ))
} 

stack_stat_min = function(stack_rows){
    stack_rows = log(stack_rows*1e9+1)
    unlist(lapply( 1:dim(stack_rows)[1], function(i)  min(stack_rows[i,],na.rm=T) ))
} 

stack_stat_pct_lessthan_value = function(stack_rows,value){
    # calculate the % of non-na values that are below some level
    stack_rows = log(stack_rows*1e9+1)
    unlist(lapply(1:dim(stack_rows)[1], function(i) sum(stack_rows[i,]<value,na.rm=T)/sum(!is.na(stack_rows[i,]),na.rm=T) ))
} 

#Determine optimal block size for loading in MODIS stack data
stack_in = dnb_stack
block_width = 10
nrows = dim(stack_in)[1]
nblocks <- nrows%/%block_width
bs_rows <- seq(1,nblocks*block_width,block_width)
bs_nrows <- rbind(matrix(block_width,length(bs_rows)-1,1),nrows-bs_rows[length(bs_rows)]+1)
print('Working on the following rows')
print(paste(bs_rows))


result <- foreach(i = 1:length(bs_rows),.packages = c('raster'),.inorder=T, .combine = c) %dopar% {
    stack_values = getValues(stack_in, bs_rows[i], bs_nrows[i])
    #dim(stack_values)
    stat_values =  stack_stat_mean(stack_rows=stack_values)
    #length(stat_values)
    stat_values
}
 
mean = dnb_stack[[1]]
mean[]=result
plot(mean,main='mean')
writeRaster(mean,'stat_mean.tif')
names(mean)='mean value of raster'
KML(mean, file='mean.kml',blur=30,overwrite=T)


plot(mean>1.5,main='greater than1.5')
mean_gt = (mean>1.5)
names(mean_gt)='mean value greater than 1.5'
KML(mean_gt, file='mean_gt_1.5.kml',blur=30,overwrite=T)

mean_gt = (mean>1.2)
names(mean_gt)='mean value greater than 1.2'
KML(mean_gt, file='mean_gt_1.2.kml',blur=30,overwrite=T)


result <- foreach(i = 1:length(bs_rows),.packages = c('raster'),.inorder=T, .combine = c) %dopar% {
    stack_values = getValues(stack_in, bs_rows[i], bs_nrows[i])
    #dim(stack_values)
    stat_values =  stack_stat_min(stack_rows=stack_values)
    #length(stat_values)
    stat_values
}

min = dnb_stack[[1]]
min[]=result
plot(min,main='min')
writeRaster(min,'stat_min.tif')
names(min)='minimum values'
KML(min, file='min.kml',blur=30,overwrite=T)

result <- foreach(i = 1:length(bs_rows),.packages = c('raster'),.inorder=T, .combine = c) %dopar% {
    stack_values = getValues(stack_in, bs_rows[i], bs_nrows[i])
    #dim(stack_values)
    stat_values =  stack_stat_pct_lessthan_value(stack_rows=stack_values,value=1.5)
    #length(sstack_stat_pct_lessthan_value
    stat_values
}
stopImplicitCluster()

lessone = dnb_stack[[1]]
lessone[]=result
plot((lessone),main='% less than 1')
writeRaster(lessone,'stat_prc_less_1pnt5.tif')
names(lessone)='pct val less 1.5'
KML(lessone, file='stat_prc_less_1pnt5.kml',blur=30,overwrite=T)




# Outlier Detection -------------------------------------------------------


library(tsoutliers)
stack_values = getValues(mean, 120, 1)
summary(as.numeric(stack_values))
plot(1:length(stack_values),stack_values)
 
tso()