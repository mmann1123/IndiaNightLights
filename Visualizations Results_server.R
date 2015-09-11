
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
library(splines)


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
setwd('/groups/manngroup/India\ VIIRS/2015')

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
library(foreach)
library(doParallel)
registerDoParallel(32)


# remove cloud cells multicore  returns NA but runs fast!
#foreach(i=1:dim(dnb_stack)[3]) %do% { dnb_stack[[i]][cld_stack[[i]]>1]=NA}
#save(dnb_stack,file = 'dnb_stack_wo_cld.RData')
#foreach(i=1:dim(zen_stack)[3]) %do% { zen_stack[[i]][cld_stack[[i]]>1]=NA}
#save(zen_stack,file = 'zen_stack_wo_cld.RData')
#foreach(i=1:dim(azt_stack)[3]) %do% { azt_stack[[i]][cld_stack[[i]]>1]=NA}
#save(azt_stack,file = 'azt_stack_wo_cld.RData')


setwd('/groups/manngroup/India\ VIIRS/2015')

load('dnb_stack_wo_cld.RData')    # start here cloud free images stored here. 
load('zen_stack_wo_cld.RData')
load('azt_stack_wo_cld.RData')



# Extract data for sites of interest --------------------------------------

# define locations of interest
locations = read.csv('MH-ESMI-Locations-Lat-Long-Overpass-Cuts-May-2015-ag.csv')
jumba.df = data.frame(STATE='MH', DISTRICT.CITY="NA", LOCATION='Jumda',LAT=20.010094,LON=77.044271, Ag.Rural=T)
locations=rbind(locations,jumba.df)

coordinates(locations)= ~LON+LAT
proj4string(locations) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# extract data (and surrounding area)
dnb_values = extract(dnb_stack,locations,fun= function(x) mean(x,na.rm=T), df=T)#buffer=1.2e3,
zen_values = extract(zen_stack,locations,fun= function(x) mean(x,na.rm=T), df=T)
azt_values = extract(azt_stack,locations,fun= function(x) mean(x,na.rm=T), df=T)

time_stamp_extract = gsub(x=colnames(dnb_values),pattern = "(.*X)(.*)(.*_dnb_v3)",replacement = "\\2")


# put into long form
put_in_long <- function(wide_data,abreviation){
	names(wide_data) = time_stamp_extract
	wide_data$location = locations$LOCATION
	wide_data = subset(wide_data,select=-c(ID))
	wide_data <- melt(wide_data )
	names(wide_data)=c('location','date.time',paste(abreviation))
	head(wide_data)
	return(wide_data)
}

dnb_values=put_in_long(dnb_values,'dnb')
zen_values= put_in_long(zen_values,'zen')
azt_values=put_in_long(azt_values,'azt')



# read in moon phase (year doy time moon_illum_frac moon_phase_angle)
phase = read.csv('moon_info.csv')
names(phase)=c('year','doy','time', 'illum', 'phase')
phase$date.time = paste(phase$year,sprintf('%03d',(phase$doy)),'.',phase$time,sep='')


# add moon phase
library(plyr)
dnb_values = join(dnb_values,zen_values)
dnb_values = join(dnb_values,azt_values)
dnb_values = join(dnb_values, phase)
head(dnb_values,20)


# run regressions 
library(splines)

fixed = lm(dnb~factor(location)+ns(zen,df=3)+ns(azt,df=3)+ns(phase,df=2)+zen:azt+zen:phase+azt:phase,data=dnb_values)
summary(fixed)


# compare actual and predicted

actual = na.omit(dnb_values)
actual$count = 1:dim(actual)[1]
actual$pred_actual = 'actual'
pred = na.omit(dnb_values)
pred$count = 1:dim(pred)[1]
pred$dnb = fixed$fitted.values
pred$pred_actual = 'predicted'

combined = rbind(actual,pred)

# plot actual vs predicted
ggplot(combined, aes(count, dnb,colour=pred_actual))+geom_point()+ 
  facet_wrap(~ location, scales="free_y")


# compare actual and resid plus constant for FE regression  # DON'T USE FE
resid = na.omit(dnb_values)
resid$count = 1:dim(resid)[1]
resid$pred_actual = 'resid+const'
# get intercepts for each village join to data
nam = names(fixed$coefficient)[grep('factor',names(fixed$coefficients))]
nam = substr(nam,17,31)
coeff = data.frame(location = nam, FE = fixed$coefficient[grep('factor',names(fixed$coefficients))]) 
resid = join(resid,coeff)
resid$dnb = fixed$residuals #+ fixed$coefficients[1] +resid$FE
head(resid)

combined2 = rbind(actual,subset(resid,select=-c(FE))   )

ggplot(combined2, aes(count, dnb,colour=pred_actual))+geom_point()+
  facet_wrap(~ location, scales="free_y")





# compare actual and resid plus constant for demeaned regression  

mean_dnb = aggregate(dnb~location,data=dnb_values,function(x) mean(x,na.rm=T))
names(mean_dnb)=c('location','mean_dnb')
dnb_values = join(dnb_values, mean_dnb,by='location')
dnb_values$demean_dnb = dnb_values$dnb - dnb_values$mean_dnb

mean_lm = lm(demean_dnb~0+ns(zen,df=3)+ns(azt,df=3)+ns(phase,df=2)+zen:azt+zen:phase+azt:phase,data=dnb_values) # omit intercept
summary(mean_lm)

resid3 = actual
resid3$pred_actual = 'resid+const'
resid3$count = 1:dim(resid3)[1]
resid3 = join(resid3,mean_dnb,by='location')
resid3$dnb = mean_lm$residuals  +resid3$mean_dnb  
head(resid3)

pred3 = actual
pred3$pred_actual = 'predicted'
pred3$count = 1:dim(pred3)[1]
pred3 = join(pred3,mean_dnb,by='location')
pred3$dnb = mean_lm$fitted.values+pred3$mean_dnb
head(pred3)


combined = rbind.fill(actual,resid3,pred3)

# plot actual vs predicted
ggplot(combined, aes(count, dnb,colour=pred_actual))+geom_point()+
  facet_wrap(~ location, scales="free_y")

ggplot(combined[combined$pred_actual != 'resid+const',], aes(count, dnb,colour=pred_actual))+geom_point()+
  facet_wrap(~ location, scales="free_y")

ggplot(combined[combined$pred_actual != 'resid+const',], aes(count, dnb,colour=pred_actual))+geom_point()+
  facet_wrap(~ location)

ggplot(combined[combined$pred_actual == 'resid+const',], aes(count, dnb,colour=pred_actual))+geom_point()+
  facet_wrap(~ location, scales="free_y")







# comare acutal and resid plus constant  pooled  # DONT USE POOLED
pool = lm(dnb~ns(zen,df=3)+ns(azt,df=3)+ns(phase,df=2)+zen:azt+zen:phase+azt:phase,data=dnb_values)
summary(pool)


# comare acutal and resid plus constant for linear regression
village = 'Tilak Nagar'
ind = lm(dnb~ns(zen,df=3)+ns(azt,df=3)+ns(phase,df=2)+zen:azt+zen:phase+azt:phase,
	data=dnb_values[dnb_values$location==village,])
summary(ind)
resid2 = actual[actual$location==village,]
resid2$pred_actual = 'resid+const'
resid2$dnb = ind$residuals + ind$coefficients[1]
head(resid2)


combined3 = rbind(actual[actual$location==village,],resid2   )

ggplot(combined3, aes(count, dnb,colour=pred_actual))+geom_point()







# Remove Lunar Cycle for raster stacks ------------------------------------------------



# plot urban
plot(1:length(data_urban$dnb),data_urban$dnb)
lm_u = lm(dnb~ns(zen,df=3)+ns(azt,df=3)+ns(phase,df=2)+zen:azt+zen:phase+azt:phase,data=(data_urban))
resid_u = lm_u$fitted.values #residuals  #fitted.values
summary(lm_u)
lunar = predict(lm1,data_urban)
#points( 1:length(resid_u),resid_u, col='red')
points(1:length(resid_u),lunar,col='green')  # predictions based on rural
#points(1:length(resid_u),(lm_u$residuals+lm_u$coefficients),col='orange')
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
