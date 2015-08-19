library(raster)
library(ggplot2)
library(scales)
library(reshape2)
library(rgeos)


# this scripts reads in raster files exported from grid_viirs_data (3).R
setwd('C:\\Users\\mmann\\Desktop\\NightTimeData\\')



# Time Series Plots for locations of interest -----------------------------


# read in data ------------------------------------------------------------

# pull available files
files = dir(pattern = '.tif')
cld = files[grep('cld',files)]
dnb = files[grep('dnb',files)]

# create raster stacks  & extract data ----------------------------------------------------

cld_stack = stack(cld)
dnb_stack = stack(dnb)

#   windows()
#   plot(dnb_stack)
#   windows()
#   plot(cld_stack)

# Extract dates -----------------------------------------------------------
time_stamp_dnb = gsub(x=names(dnb_stack),pattern = "(.*X)(.*)(.*_dnb_v3)",replacement = "\\2")
time_stamp_dnb = strptime(time_stamp_dnb,"%Y%j.%H%M")
time_stamp_cld = gsub(x=names(cld_stack),pattern = "(.*X)(.*)(.*_cld_v3)",replacement = "\\2")
time_stamp_cld = strptime(time_stamp_cld,"%Y%j.%H%M")

# not all stacks have same dates
all.equal(time_stamp_dnb,time_stamp_cld)
# limit stacks to common elements
common_dnb = (time_stamp_dnb %in% intersect(time_stamp_dnb,time_stamp_cld))
dnb_stack = dnb_stack[[ (1:length(common_dnb))[common_dnb] ]]
common_cld = (time_stamp_cld %in% intersect(time_stamp_dnb,time_stamp_cld))
cld_stack = cld_stack[[ (1:length(common_cld))[common_cld] ]]

# extract dates from limited set and test equality
time_stamp_dnb = gsub(x=names(dnb_stack),pattern = "(.*X)(.*)(.*_dnb_v3)",replacement = "\\2")
time_stamp_dnb = strptime(time_stamp_dnb,"%Y%j.%H%M")
time_stamp_cld = gsub(x=names(cld_stack),pattern = "(.*X)(.*)(.*_cld_v3)",replacement = "\\2")
time_stamp_cld = strptime(time_stamp_cld,"%Y%j.%H%M")

# Extract data for sites of interest --------------------------------------

# define locations of interest
locations = read.csv('C://Users/mmann/Google Drive/India Night Lights/Data/MH-ESMI-Locations-Lat-Long-Overpass-Cuts-May-2015-ag.csv')
jumba.df = data.frame(STATE='MH', DISTRICT.CITY="NA", LOCATION='Jumda',LAT=20.010094,LON=77.044271, Ag.Rural=T)
locations=rbind(locations,jumba.df)

coordinates(locations)= ~LON+LAT
proj4string(locations) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# extract data (and surrounding area)
dnb_values = extract(dnb_stack,locations,buffer=1.2e3,fun= function(x) mean(x,na.rm=T), df=T)
cld_values = extract(cld_stack,locations,buffer=1.2e3,fun= function(x) mean(x,na.rm=T), df=T)

# NA out values with clouds  (zero = confident no clouds, 1 = probably no clouds)
dnb_values[cld_values>1]=NA

# Create plots ------------------------------

# organize data for plots

# get time stamps of extracted data , convert to date format
time_stamp_dnb = strptime(gsub(x=names(dnb_values),pattern = "(.*X)(.*)(.*_dnb_v3)",replacement = "\\2"),"%Y%j.%H%M")
#time_stamp_cld = strptime(gsub(x=names(cld_values),pattern = "(.*X)(.*)(.*_cld_v3)",replacement = "\\2"),"%Y%j.%H%M")
#all.equal(time_stamp_dnb,time_stamp_cld)

# put all locations names and city/rural into data frame  
dnb_values$LOCATION = as.character(locations$LOCATION)
dnb_values$Ag.Rural = as.character(locations$Ag.Rural)


dnb_values$data_type = 'dnb'
#cld_values$LOCATION = as.character(locations$LOCATION)
#cld_values$data_type = 'cld'

# create long format data sets
dnb_values2= cbind(melt(dnb_values, id.vars=c("LOCATION",'data_type','Ag.Rural')),date_time=time_stamp_dnb,date=as.Date(time_stamp_dnb))
#cld_values2= cbind(melt(cld_values, id.vars=c("LOCATION",'data_type')),date_time=time_stamp_cld,date=as.Date(time_stamp_cld))
head(dnb_values2)

# Combind cloud and dnb long format data 
data = rbind(dnb_values2)
data$Ag.Rural[data$Ag.Rural==T]='Agricultural'
data$Ag.Rural[data$Ag.Rural==F]='Urban'

# make plots 
windows()
ggplot(data,aes(x=date,y=log(value*1e9+1)))+geom_point()+scale_x_date(labels = date_format("%m-%d"))+ 
    facet_wrap( ~ LOCATION+Ag.Rural, ncol = 2,scales = 'free_y')  


write.csv(data,'C:/Users/mmann/Google Drive/India Night Lights/Data/MH-ESMI-Locations-DNB-output-1kmneighbors.csv')







# Stable Lights Map -------------------------------------------------------

# read in data 
setwd('C:\\Users\\mmann\\Desktop\\NightTimeData\\')

# pull available files
files = dir(pattern = '.tif')
cld = files[grep('cld',files)]
dnb = files[grep('dnb',files)]

# create raster stacks  & extract data  
cld_stack = stack(cld)
dnb_stack = stack(dnb)

# Extract dates  
time_stamp_dnb = gsub(x=names(dnb_stack),pattern = "(.*X)(.*)(.*_dnb_v3)",replacement = "\\2")
time_stamp_dnb = strptime(time_stamp_dnb,"%Y%j.%H%M")
time_stamp_cld = gsub(x=names(cld_stack),pattern = "(.*X)(.*)(.*_cld_v3)",replacement = "\\2")
time_stamp_cld = strptime(time_stamp_cld,"%Y%j.%H%M")

# not all stacks have same dates
all.equal(time_stamp_dnb,time_stamp_cld)

# limit stacks to common elements
common_dnb = (time_stamp_dnb %in% intersect(time_stamp_dnb,time_stamp_cld)) # find common times
dnb_stack = dnb_stack[[ (1:length(common_dnb))[common_dnb] ]]       # limit rasters to common times
common_cld = (time_stamp_cld %in% intersect(time_stamp_dnb,time_stamp_cld))
cld_stack = cld_stack[[ (1:length(common_cld))[common_cld] ]]
time_stamp = time_stamp_dnb[(time_stamp_dnb %in% intersect(time_stamp_dnb,time_stamp_cld))]

# remove cloud cells
#dnb_stack[cld_stack>1]=NA                      # already ran this just load
#save(dnb_stack,file = 'dnb_stack_wo_cld.RData')
load('dnb_stack_wo_cld.RData')                  # loads a cloud free version of dnb_stack
setZ(dnb_stack,as.Date(time_stamp))

# remove outliers 
source('G://Faculty/Mann/Scripts/SplineAndOutlierRemoval.R')

row = 900
plot(log(getValues(dnb_stack,500,1)[row,]*1e9+1))
points(log(SplineAndOutlierRemoval(getValues(dnb_stack,500,1)[400,], dates=time_stamp, out_sigma=3, spline_spar=0.55, out_iterations=30,pred_dates=time_stamp)*1e9+1),col='red')
    
# calculate stats   USER foreach_stack_stats.R for all examples

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



# Unused code -------------------------------------------------------------
#   proj_locations = spTransform( locations, CRS( "+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=7 +lat_2=-32 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs " ) )  #http://spatialreference.org/ref/esri/102029/
#   locations_buffer = gBuffer(proj_locations,capStyle = 'round',width = 3000)
#   locations_buffer = spTransform( locations_buffer, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") )  #http://spatialreference.org/ref/esri/102029/
#   
#   dnb_buffer_values =  extract(dnb_stack,locations_buffer)
#   
