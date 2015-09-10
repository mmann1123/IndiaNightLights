
# Run the following in bash before starting R
# module load proj.4/4.8.0
# module load gdal
# module load R/3.0.2
# module load gcc/4.9.0
# R




library(raster)
library(ggplot2)
library(scales)
library(reshape2)
library(rgeos)


# this scripts reads in raster files exported from grid_viirs_data (3).R



# Time Series Plots for locations of interest -----------------------------


# read in data ------------------------------------------------------------
setwd('/groups/manngroup/India\ VIIRS/2015')

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
locations = read.csv('/groups/manngroup/India\ VIIRS/2015/MH-ESMI-Locations-Lat-Long-Overpass-Cuts-May-2015-ag.csv')
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


write.csv(data,'./MH-ESMI-Locations-DNB-output-1kmneighbors.csv')



# Lunar adjustments ------------------------------------------------------
# try to remove cyclical lunar signal from cells

# read in data
setwd('R:/Mann Research/India Night Time Lights/2015 v2/')

# pull available files
files = dir(pattern = '.tif')
cld = files[grep('cld',files)]
dnb = files[grep('dnb',files)]
zen = files[grep('zen',files)]
azt = files[grep('azt',files)]


# create raster stacks  & extract data
cld_stack = stack(cld)
dnb_stack = stack(dnb)
zen_stack = stack(zen)
azt_stack = stack(azt)


# Extract dates
time_stamp_dnb = gsub(x=names(dnb_stack),pattern = "(.*X)(.*)(.*_dnb_v3)",replacement = "\\2")
time_stamp_dnb = strptime(time_stamp_dnb,"%Y%j.%H%M")
time_stamp_cld = gsub(x=names(cld_stack),pattern = "(.*X)(.*)(.*_cld_v3)",replacement = "\\2")
time_stamp_cld = strptime(time_stamp_cld,"%Y%j.%H%M")
time_stamp_zen = gsub(x=names(zen_stack),pattern = "(.*X)(.*)(.*_zen_v3)",replacement = "\\2")
time_stamp_zen = strptime(time_stamp_zen,"%Y%j.%H%M")
time_stamp_azt = gsub(x=names(azt_stack),pattern = "(.*X)(.*)(.*_azt_v3)",replacement = "\\2")
time_stamp_azt = strptime(time_stamp_azt,"%Y%j.%H%M")


# not all stacks have same dates
all.equal(time_stamp_dnb,time_stamp_cld)
all.equal(time_stamp_dnb,time_stamp_zen)
all.equal(time_stamp_dnb,time_stamp_azt)


# limit stacks to common elements
common_dnb = (time_stamp_dnb %in% intersect(time_stamp_dnb,time_stamp_cld))
dnb_stack = dnb_stack[[ (1:length(common_dnb))[common_dnb] ]]
common_cld = (time_stamp_cld %in% intersect(time_stamp_dnb,time_stamp_cld))
cld_stack = cld_stack[[ (1:length(common_cld))[common_cld] ]]


# remove cloud cells multicore  returns NA but runs fast!
#library(foreach)
#library(doParallel)
#registerDoParallel(32)
# remove cloud cells multicore  returns NA but runs fast!
#foreach(i=1:dim(dnb_stack)[3]) %do% { dnb_stack[[i]][cld_stack[[i]]>1]=NA}
#save(dnb_stack,file = 'dnb_stack_wo_cld.RData')
load('dnb_stack_wo_cld.RData')



# Extract data for sites of interest --------------------------------------

# define locations of interest
locations = read.csv('../MH-ESMI-Locations-Lat-Long-Overpass-Cuts-May-2015-ag.csv')
jumba.df = data.frame(STATE='MH', DISTRICT.CITY="NA", LOCATION='Jumda',LAT=20.010094,LON=77.044271, Ag.Rural=T)
locations=rbind(locations,jumba.df)

coordinates(locations)= ~LON+LAT
proj4string(locations) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# extract data (and surrounding area)
dnb_values = extract(dnb_stack,locations,buffer=1.2e3,fun= function(x) mean(x,na.rm=T), df=T)
zen_values = extract(zen_stack,locations,buffer=1.2e3,fun= function(x) mean(x,na.rm=T), df=T)
azt_values = extract(azt_stack,locations,buffer=1.2e3,fun= function(x) mean(x,na.rm=T), df=T)

time_stamp_extract = gsub(x=colnames(dnb_values),pattern = "(.*X)(.*)(.*_dnb_v3)",replacement = "\\2")

# put into long form
names(dnb_values) = time_stamp_extract
dnb_values$location = locations$LOCATION
dnb_values = subset(dnb_values,select=-c(ID))
#dnb_values =t(dnb_values)
head((dnb_values[]))
dnb_values <- melt(dnb_values )
names(dnb_values)=c('location','date.time','dnb')
head(dnb_values)

# read in moon phase (year doy time moon_illum_frac moon_phase_angle)
phase = read.csv('moon_info.csv')
names(phase)=c('year','doy','time', 'illum', 'phase')
phase$date.time = paste(phase$year,sprintf('%03d',(phase$doy)),'.',phase$time,sep='')

# add moon phase
library(plyr)
dnb_values = join(dnb_values, phase)

#dnb_values = na.omit(dnb_values)      # to help plot raw data and modeled

# run regressions 
library(splines)
library(plm)

fixed <- plm(dnb~ns(zen,df=3)+ns(azt,df=3)+ns(phase,df=2)+zen:azt+zen:phase+azt:phase,data=dnb_values, index=c("location", "date.time"), model="within")
summary(fixed)


#######################################
#######################################

# look at dnb and lunar time series
ts_dnb  = (extract(dnb_stack,100,100))
ts_zen =  (extract(zen_stack,100,100))
ts_azt =  (extract(azt_stack,100,100))
time_stamp_extract = gsub(x=colnames(ts_dnb),pattern = "(.*X)(.*)(.*_dnb_v3)",replacement = "\\2")


# compile data
data = data.frame(date.time = time_stamp_extract,dnb = as.numeric(ts_dnb), zen = as.numeric(ts_zen), 
                  azt=as.numeric(ts_azt))

# read in moon phase (year doy time moon_illum_frac moon_phase_angle)
phase = read.csv('/groups/manngroup/India\ VIIRS/2015/moon_info.csv')
names(phase)=c('year','doy','time', 'illum', 'phase')
phase$date.time = paste(phase$year,sprintf('%03d',(phase$doy)),'.',phase$time,sep='')

# add moon phase
library(plyr)
data = join(data, phase)
data = na.omit(data)      # to help plot raw data and modeled

# compare time series and modeled
library(splines)


plot(1:length(data$dnb),data$dnb)
lm1 = lm(dnb~ns(zen,df=3)+ns(azt,df=3)+ns(phase,df=2)+zen:azt+zen:phase+azt:phase,data=(data))
resid = lm1$fitted.values #residuals  #fitted.values
summary(lm1)
points( 1:length(resid),resid, col='red')


### pull urban data


# look at dnb and lunar time series
urb = data.frame(lat=17.35875,lon=78.45525)
coordinates(urb)=~lon+lat

ts_dnb  = extract(dnb_stack,urb )
ts_zen =  extract(zen_stack,urb  )
ts_azt =  extract(azt_stack,urb)
time_stamp_extract = gsub(x=colnames(ts_dnb),pattern = "(.*X)(.*)(.*_dnb_v3)",replacement = "\\2")


# compile data
data_urban = data.frame(date.time = time_stamp_extract,dnb = as.numeric(ts_dnb), zen = as.numeric(ts_zen),
                        azt=as.numeric(ts_azt))

# add moon phase
data_urban = join(data_urban, phase)
data_urban = na.omit(data_urban)      # to help plot raw data and modeled

# plot urban 
plot(1:length(data_urban$dnb),data_urban$dnb)
lm_u = lm(dnb~ns(zen,df=3)+ns(azt,df=3)+ns(phase,df=2)+zen:azt+zen:phase+azt:phase,data=(data_urban))
resid_u = lm_u$fitted.values #residuals  #fitted.values
summary(lm_u)
lunar = predict(lm1,data_urban)
#points( 1:length(resid_u),resid_u, col='red')
points(1:length(resid_u),lunar,col='green')  # predictions based on rural 
#points(1:length(resid_u),(lm_u$residuals+lm_u$coefficients[1]),col='orange')
points(1:length(resid_u),data_urban$dnb-lunar,col='blue')  # predictions based on rural

# maybe create panel dataset of uninhabited areas? subtract prediction from that




# Outlier detection  ----------------------------------------------------------
# https://blog.twitter.com/2015/introducing-practical-and-robust-anomaly-detection-in-a-time-series

# install.packages("devtools")
# devtools::install_github("twitter/AnomalyDetection")
library(AnomalyDetection)

help(AnomalyDetectionTs)
data_urban$date.time = strptime(data_urban$date.time,"%Y%j.%H%M" )
data3= data_urban[,c('date.time','dnb')]
data3$dnb = data3$dnb*1e9
row.names(data3) = 1:dim(data3)[1] 
res = AnomalyDetectionTs(data3, max_anoms=0.03, direction='pos', plot=TRUE, 
                         piecewise_median_period_weeks=8,longterm=T)
res$plot


# impulse sateration https://cran.r-project.org/web/packages/gets/gets.pdf from felix
library(zoo)
library(gets)
data4 = zoo(data_urban[,'dnb'],data_urban[,'date.time'])

# add a time trend
t = zoo(1:length(data4[,1]),data_urban[,'date.time'])

sat = isat(data4, t.pval=min(0.01,(1/length(data4[,1]))),mxreg=t)
plot(data4)


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
common_dnb = (time_stamp_dnb %in% intersect(time_stamp_dnb,time_stamp_cld))
dnb_stack = dnb_stack[[ (1:length(common_dnb))[common_dnb] ]]
common_cld = (time_stamp_cld %in% intersect(time_stamp_dnb,time_stamp_cld))
cld_stack = cld_stack[[ (1:length(common_cld))[common_cld] ]]

# remove cloud cells
#dnb_stack[cld_stack>1]=NA
#save(dnb_stack,file = 'dnb_stack_wo_cld.RData')
load('dnb_stack_wo_cld.RData')


# calculate stats
beginCluster(type='SOCK') #http://rpackages.ianhowson.com/rforge/raster/man/cluster.html

dnb_stack_mean <- clusterR(dnb_stack, mean, args=list(na.rm=T))
save(dnb_stack_mean,file = 'dnb_stack_mean.RData')

dnb_stack_quan <- clusterR(dnb_stack, quantile, args=list(na.rm=T))
save(dnb_stack_quan,file = 'dnb_stack_quan.RData')

endCluster()

# Unused code -------------------------------------------------------------
#   proj_locations = spTransform( locations, CRS( "+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=7 +lat_2=-32 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs " ) )  #http://spatialreference.org/ref/esri/102029/
#   locations_buffer = gBuffer(proj_locations,capStyle = 'round',width = 3000)
#   locations_buffer = spTransform( locations_buffer, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") )  #http://spatialreference.org/ref/esri/102029/
#   
#   dnb_buffer_values =  extract(dnb_stack,locations_buffer)
#   
