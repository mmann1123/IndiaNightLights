
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
library(flexclust)   # allow for multidimentional clustering with predict function
library(sampling)
library(plyr)

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

fixed = lm(dnb~factor(location)+ns(zen,df=3)+ns(azt,df=3)+ns(phase,df=2)+zen:azt+
    zen:phase+azt:phase,data=dnb_values)

fixed = lm(dnb~factor(location)+zen+azt+phase+zen:azt+
   zen:phase^2+azt:phase^2,data=dnb_values)

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


## compare actual and resid plus constant for FE regression  # DON'T USE FE
#resid = na.omit(dnb_values)
#resid$count = 1:dim(resid)[1]
#resid$pred_actual = 'resid+const'
## get intercepts for each village join to data
#nam = names(fixed$coefficient)[grep('factor',names(fixed$coefficients))]
#nam = substr(nam,17,31)
#coeff = data.frame(location = nam, FE = fixed$coefficient[grep('factor',names(fixed$coefficients))]) 
#resid = join(resid,coeff)
#resid$dnb = fixed$residuals #+ fixed$coefficients[1] +resid$FE
#head(resid)
#
#combined2 = rbind(actual,subset(resid,select=-c(FE))   )
#
#ggplot(combined2, aes(count, dnb,colour=pred_actual))+geom_point()+
#  facet_wrap(~ location, scales="free_y")





# compare actual and resid plus constant for demeaned regression  

mean_dnb = aggregate(dnb~location,data=dnb_values,function(x) mean(x,na.rm=T))
names(mean_dnb)=c('location','mean_dnb')
dnb_values = join(dnb_values, mean_dnb,by='location')
dnb_values$demean_dnb = dnb_values$dnb - dnb_values$mean_dnb

mean_lm = lm(demean_dnb~0+I(mean_dnb>1.5e-8)*(ns(zen,df=3)+ns(azt,df=3)+ns(phase,df=3)+
      zen:azt+poly(zen:phase,2)+azt:phase),data=dnb_values) # omit intercept

mean_lm = lm(demean_dnb~0+I(mean_dnb>1.5e-8)*(ns(zen,df=3)+ns(azt,df=3)+ns(phase,df=3)+
      zen*azt+ns(zen*phase,3)+ns(azt*phase,3)),data=dnb_values) # omit intercept

summary(mean_lm)
#I(mean_dnb>1e-8)


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


# test if residual is white noise  NOT NEEDED TOO STRICT
#H0: The data are independently distributed (i.e. the correlations in the population from which the sample is taken are 0, so that any observed correl$
#Ha: The data are not independently distributed; they exhibit serial correlation.
#
#library(zoo)
#for (location in as.character(unique(combined$location))){
#   test_data = combined[combined$location==location,]
#   test_data_ag = aggregate(dnb~date.time,data=test_data, FUN=mean)
#   test_data_ts = zoo(test_data_ag$dnb,test_data_ag$date.time)
#   print( Box.test(test_data_ts, lag = 1, type = c("Ljung-Box"), fitdf = 0))
#  }



## comare acutal and resid plus constant  pooled  # DONT USE POOLED
#pool = lm(dnb~ns(zen,df=3)+ns(azt,df=3)+ns(phase,df=2)+zen:azt+zen:phase+azt:phase,data=dnb_values)
#summary(pool)
#
#
## comare acutal and resid plus constant for linear regression
#village = 'Tilak Nagar'
#ind = lm(dnb~ns(zen,df=3)+ns(azt,df=3)+ns(phase,df=2)+zen:azt+zen:phase+azt:phase,
#	data=dnb_values[dnb_values$location==village,])
#summary(ind)
#resid2 = actual[actual$location==village,]
#resid2$pred_actual = 'resid+const'
#resid2$dnb = ind$residuals + ind$coefficients[1]
#head(resid2)
#
#
#combined3 = rbind(actual[actual$location==village,],resid2   )
#
#ggplot(combined3, aes(count, dnb,colour=pred_actual))+geom_point()





############################################################################################
############################################################################################




# Demeaned regression on 1000 point sample -------------------------------------------


# load data
setwd('/groups/manngroup/India\ VIIRS/2015')

load('dnb_stack_wo_cld.RData')    # start here cloud free images stored here.
load('zen_stack_wo_cld.RData')
load('azt_stack_wo_cld.RData')


# create a plygon to create samples within
source('/groups/manngroup/India VIIRS/IndiaNightLights/PolygonFromExtent.R')


# create sample of 10000 set seed for reproducibility
set.seed(2)
sampled = spsample(PolygonFromExtent(extent(72, 81.50, 15, 22.5),crs=CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ')),
     n=10000,type='random')


# extract data (and surrounding area)
dnb_values = extract(dnb_stack,sampled,fun= function(x) mean(x,na.rm=T), df=T)#buffer=1.2e3,
zen_values = extract(zen_stack,sampled,fun= function(x) mean(x,na.rm=T), df=T)
azt_values = extract(azt_stack,sampled,fun= function(x) mean(x,na.rm=T), df=T)
time_stamp_extract = gsub(x=colnames(dnb_values),pattern = "(.*X)(.*)(.*_dnb_v3)",replacement = "\\2")


# extract local testing data 
# define locations of interest
locations = read.csv('MH-ESMI-Locations-Lat-Long-Overpass-Cuts-May-2015-ag.csv')
jumba.df = data.frame(STATE='MH', DISTRICT.CITY="NA", LOCATION='Jumda',LAT=20.010094,LON=77.044271, Ag.Rural=T)
locations=rbind(locations,jumba.df)

coordinates(locations)= ~LON+LAT
proj4string(locations) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
dnb_values_loc = extract(dnb_stack,locations,fun= function(x) mean(x,na.rm=T), df=T)#buffer=1.2e3,
zen_values_loc = extract(zen_stack,locations,fun= function(x) mean(x,na.rm=T), df=T)
azt_values_loc  = extract(azt_stack,locations,fun= function(x) mean(x,na.rm=T), df=T)
time_stamp_extract_loc = gsub(x=colnames(dnb_values_loc),pattern = "(.*X)(.*)(.*_dnb_v3)",replacement = "\\2")



# put into long form
put_in_long2 <- function(wide_data,abreviation){
        names(wide_data) = time_stamp_extract
        wide_data$location = row.names(wide_data)
        wide_data = subset(wide_data,select=-c(ID))
        wide_data <- melt(wide_data )
        names(wide_data)=c('location','date.time',paste(abreviation))
        head(wide_data)
        return(wide_data)
}
put_in_long_loc <- function(wide_data,abreviation){
        names(wide_data) = time_stamp_extract_loc
        wide_data$location = locations$LOCATION
        wide_data = subset(wide_data,select=-c(ID)) # extra variable was added, remove
        wide_data <- melt(wide_data )
        names(wide_data)=c('location','date.time',paste(abreviation))
        head(wide_data)
        return(wide_data)
}

dnb_values=put_in_long2(dnb_values,'dnb')
zen_values= put_in_long2(zen_values,'zen')
azt_values=put_in_long2(azt_values,'azt')
dnb_values_loc=put_in_long_loc(dnb_values_loc,'dnb')
zen_values_loc= put_in_long_loc(zen_values_loc,'zen')
azt_values_loc=put_in_long_loc(azt_values_loc,'azt')


# add moon phase

# read in moon phase (year doy time moon_illum_frac moon_phase_angle)
phase = read.csv('moon_info.csv')
names(phase)=c('year','doy','time', 'illum', 'phase')
phase$date.time = paste(phase$year,sprintf('%03d',(phase$doy)),'.',phase$time,sep='')


# join in moon and dnb
dnb_values = join(dnb_values,zen_values) # add moon characteristics to dnb values
dnb_values = join(dnb_values,azt_values)
dnb_values = join(dnb_values, phase)

dnb_values_loc = join(dnb_values_loc,zen_values_loc) # add moon characteristics to dnb values
dnb_values_loc = join(dnb_values_loc,azt_values_loc)
dnb_values_loc = join(dnb_values_loc,phase)

head(dnb_values)
head(dnb_values_loc)


# save output to load quickly
setwd('/groups/manngroup/India\ VIIRS/2015')

#save(dnb_values,file='dnb_values_w_moon.RData')
#save(dnb_values_loc,file='dnb_values_w_moon_loc.RData')
load('dnb_values_w_moon.RData')
load('dnb_values_w_moon_loc.RData')


# merge training and testing data
dnb_values$type = 'training'
dnb_values_loc$type = 'testing'
dnb_values = rbind.fill(dnb_values,dnb_values_loc)


# create descriptive statistics
sd_dnb = aggregate(dnb~location,data=dnb_values,function(x) sd(x,na.rm=T))
mn_dnb = aggregate(dnb~location,data=dnb_values,function(x) mean(x,na.rm=T))
names(sd_dnb)=c('location','sd_dnb')
names(mn_dnb)=c('location','mn_dnb')
dnb_values = join(dnb_values, mn_dnb)
dnb_values = join(dnb_values, sd_dnb)


##########
## find k-means = 3 groups based on dnb values (based on mean  of values)
## NOT WORKING COMPARE WITH OLD KMEANS RESULTS.... 
#na_dnb = na.omit(dnb_values)
## cluster based on mean and sd of data 
#na_dnb_mn_sd = join(mn_dnb, sd_dnb,by='location')  # store mean and sd of dnb values for each location
#na_dnb = join(na_dnb, na_dnb_mn_sd,by='location')  # store mean and sd of dnb values for each location
## stratified sample to make sure we have adequate representation for each land use class
#na_dnb$quantile = cut(na_dnb$mn_dnb,quantile(na_dnb[,'mn_dnb']))
#set.seed(2)
##s=strata(na_dnb,'quantile',size=c(250,250,250,250,1), method="srswor") # create strata based on quantiles
##na_dnb_sample = getdata( na_dnb,s)
#na_dnb_sample = na_dnb[sample(nrow(na_dnb),10000),]
#cl1 = kcca(na_dnb_sample[,c('dnb')], k=3, kccaFamily("kmeans"))   # THIS MIGHT BE LETTING INDV OBS HAVE DIFFERENT MEMBERSHIP WITHIN 1 LOCATION
#image(cl1)
#na_dnb$kmn = predict(cl1, newdata=na_dnb[,c('dnb'),drop=F])
## add mean sd and cluster back into full data
#dnb_values = join(dnb_values, na_dnb, by=c('location','date.time'))
#table(na_dnb$kmn)
#table(dnb_values$kmn)

## OLD KMEANS METHOD
#dnb_num = data.frame(dnb=dnb_values$dnb)   # set up a dataframe to store row numbers so can be joined later
#kmn2 = kmeans(na.omit(dnb_num),3)
#kmn2 = data.frame(row = names(kmn2$cluster),kmn2 = kmn2$cluster)
#dnb_values$row = row.names(dnb_values)
#dnb_values = join(dnb_values, kmn2,by='row')
#head(dnb_values)
#table(dnb_values$kmn2)


## get summary of old kmeans breaks 
#tapply(dnb_values$dnb, dnb_values$kmn2, summary)


#########
# compare actual and resid plus constant for demeaned regression using kmean cluster for slope & intercept dummies

# demean the y variable to simulate FE
dnb_values$demean_dnb = dnb_values$dnb - dnb_values$mn_dnb

#mean_lm = lm(demean_dnb~0+factor(kmn)*(ns(zen,df=3)+ns(azt,df=3)+ns(phase,df=3)+zen*azt+
#      ns(zen*phase,3)+ns(azt*phase,3)),data=dnb_values) # omit intercept
#mean_lm = lm(demean_dnb~0+factor(kmn2)*(ns(zen*azt,2)+
#      ns(zen*phase,2)+ns(azt*phase,2)),data=na.omit(dnb_values)) # omit intercept

#mean_lm = lm(demean_dnb~0+I(dnb<1.099e-8)*(ns(zen*azt,2)+ns(zen*phase,2)+ns(azt*phase,2)),data=na.omit(dnb_values)
#	+I(dnb<4e-6)*(ns(zen*azt,2)+ns(zen*phase,2)+ns(azt*phase,2)),data=na.omit(dnb_values)) # omit intercept


mean_lm = lm(demean_dnb~0+I(mn_dnb<1.099e-8)*(ns(zen*azt,3)+
      ns(zen*phase,3)+ns(azt*phase,3)+ns(phase,4)),data=na.omit(dnb_values[dnb_values$type=='training',])) # omit intercept

summary(mean_lm)


# compare actual and predicted
remove(actual, pred, resid)
actual = na.omit(dnb_values)
actual$count = 1:dim(actual)[1]
actual$pred_actual = 'actual'

pred = actual
pred$pred_actual = 'predicted'
pred$count = 1:dim(pred)[1]
pred = join(pred,mn_dnb,by='location')
pred$dnb = mean_lm$fitted.values+pred$mn_dnb   # add back the mean from demeaning process
head(pred)

resid = actual
resid$pred_actual = 'resid+const'
resid$count = 1:dim(resid)[1]
resid = join(resid,mn_dnb,by='location')
resid$dnb = mean_lm$residuals  +resid$mn_dnb
head(resid)

combined = rbind.fill(actual,pred, resid)


# sample a portion for graphing

set.seed(11)  # 2 has declining values
actual$quantile = cut(actual$mn_dnb,quantile(actual[,'mn_dnb']))
s=strata(actual,'quantile',size=c(3,3,3,3,1), method="srswor") # create strata based on quantiles
combined_sample = getdata(combined,s)  # only get one observation for each location 
combined_sample = combined[combined$location %in% combined_sample$location,]   # limit to desired locations keeping all observations


# plot actual vs predicted

ggplot(combined_sample, aes(count, dnb,colour=pred_actual))+geom_point()+
  facet_wrap(~ location, scales="free_y")

ggplot(combined_sample[combined_sample$pred_actual != 'resid+const',], aes(count, dnb,colour=pred_actual,alpha=0.2))+geom_point()+
  facet_wrap(~ location, scales="free_y")

ggplot(combined_sample[combined_sample$pred_actual != 'predicted',], aes(count, dnb,colour=pred_actual))+geom_point()+
  facet_wrap(~ location, scales='free_y')

ggplot(combined_sample[combined_sample$pred_actual == 'resid+const',], aes(count, dnb,colour=pred_actual))+geom_point()+
  facet_wrap(~ location, scales="free_y")



# predict to testing locations
dnb_values_loc = dnb_values[dnb_values$type=='testing',]
mean_lm_loc = predict(mean_lm,dnb_values_loc)
summary(mean_lm_loc)
mean_lm_loc=data.frame(fitted.values=mean_lm_loc)
mean_lm_loc$residuals = dnb_values_loc$dnb - mean_lm_loc$fitted.values

# compare actual and predicted
actual_loc = dnb_values_loc
actual_loc$count = 1:dim(actual_loc)[1]
actual_loc$pred_actual = 'actual'

pred_loc = actual_loc
pred_loc$pred_actual = 'predicted'
pred_loc$count = 1:dim(pred_loc)[1]
pred_loc = join(pred_loc,mn_dnb,by='location')
pred_loc$dnb = mean_lm_loc$fitted.values+pred_loc$mn_dnb   # add back the mean from demeaning process
head(pred_loc)

resid_loc = actual_loc
resid_loc$pred_actual = 'resid+const'
resid_loc$count = 1:dim(resid_loc)[1]
resid_loc = join(resid_loc,mn_dnb,by='location')
resid_loc$dnb = mean_lm_loc$residuals # +resid_loc$mn_dnb
head(resid_loc)

combined_loc = rbind.fill(actual_loc,pred_loc, resid_loc)


ggplot(combined_loc[combined_loc$pred_actual != 'resid+const',], aes(count, dnb,colour=pred_actual,alpha=0.2))+geom_point()+
  facet_wrap(~ location, scales="free_y")

ggplot(combined_loc[combined_loc$pred_actual == 'resid+const',], aes(count, dnb,colour=pred_actual))+geom_point()+
  facet_wrap(~ location, scales="free_y")


ggplot(combined_loc[combined_loc$pred_actual != 'predicted',], aes(count, dnb,colour=pred_actual))+geom_point()+
  facet_wrap(~ location, scales="free_y")






# Compare to actual outage data  ----------------------------------------------------------

# read in each xls file and convert to long format and write to csv
# install.packages('readxl')
library(readxl)
file_list = list.files('..//TestingData//',pattern='.xlsx')
for(i in 1:length(file_list)){
	a_location = read_excel(paste("..//TestingData//",file_list[i],sep=''), na = "99")
	names(a_location) = c('location','date','hour',paste(1:60))	

	# put data into long form and sort 
	longs =  melt(a_location,id=c('location','date','hour') )
	names(longs) = c('location','date','hour','minute','voltage')
	longs = longs[with(longs,order(location,date,hour)),]
	#http://www.regexr.com/
	write.csv(longs,paste('..//TestingData//',gsub('(.[a-z]+$)','\\2',file_list[i]),'.csv',sep=''))
}

# stack all files together
file_list = list.files('..//TestingData//',pattern='.csv')
voltage = read.csv(paste('..//TestingData//',file_list[1],sep=''))
for(i in 2:length(file_list)){
	inner =  read.csv(paste('..//TestingData//',file_list[i],sep=''))
	voltage = rbind.fill(voltage,inner)
}

#save the output
#save(voltage,file='..//TestingData//Voltage.RData')
load('..//TestingData//Voltage.RData')


# deal with date time
OlsonNames() # full list of time zones (only works on unix?)
voltage$date.hour.minute = paste(voltage$date,sprintf('%02d',voltage$hour), sprintf('%02d',voltage$minute),sep='.')
voltage$date.hour.minute =  as.character(strptime(voltage$date.hour.minute, '%Y-%m-%d.%H.%M'))
voltage$date.hour.minute <- as.POSIXct(voltage$date.hour.minute, tz="Asia/Calcutta")
head(voltage$date.hour.minute)
attributes(voltage$date.hour.minute)$tzone = 'UTC'  # convert Calcutta time to UTC
voltage$date.time2 = voltage$date.hour.minute
head(voltage$date.hour.minute)

#same for all locations 
actual_loc$date.time2 = strptime(actual_loc$date.time,'%Y%j.%H%M')
actual_loc$date.time2 <- as.POSIXct(actual_loc$date.time2, tz="UTC")


# match to closest time in local data 
locales = unique(actual_loc$location)
locale = locales[11]
test_site = actual_loc[actual_loc$location==locale,]
test_voltage = voltage[voltage$location ==locale,]
test_join = cbind(test_site, test_voltage[ sapply(test_site$date.time2, 
                      function(x) which.min(abs(difftime(x, test_voltage$date.time2)))), ])
head(test_join,15)

test_join$count = 1:dim(test_join)[1]
test_join$lights_out = 'on'
test_join$lights_out[test_join$voltage==0]='off'
ggplot(test_join, aes(count, dnb,colour=lights_out))+geom_point() + 
     geom_vline(xintercept = test_join$count[test_join$voltage==0],colour='red',alpha=0.25)+
     geom_vline(xintercept = test_join$count[is.na(test_join$voltage)],colour='blue',alpha=0.25)










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



# look at distribution of values
#qplot(mean_dnb, data=mean_dnb, geom="histogram",binwidth=diff(range(mean_dnb$mean_dnb))/100)
#qplot(mean_dnb, data=mean_dnb, geom="histogram",binwidth=diff(range(mean_dnb$mean_dnb))/100)+coord_cartesian(xlim = c(0,5e-8))
#
#table(cut(mean_dnb$mean_dnb, breaks=c(0, 1e-9,,1e-8,1.5e-8,2e-8,10), include.lowest=TRUE))
#quantile(mean_dnb$mean_dnb)


# OLD KMEANS METHOD
#dnb_num = data.frame(dnb=dnb_values$dnb)   # set up a dataframe to store row numbers so can be joined later
#row.names(dnb_num)
#kmn= kmeans(na.omit(dnb_num),3)
#kmn = data.frame(row = names(kmn$cluster),kmn = kmn$cluster)
#dnb_values$row = row.names(dnb_values)
#dnb_values = join(dnb_values, kmn)



## create kmeans clusters on wide data
#na_dnb_wide_sample = dnb_values[sample(nrow(dnb_values),5000),]
#wide_kmn_class = kcca(na_dnb_wide_sample, k=3, kccaFamily("kmeans"))
#na_dnb$kmn = predict(cl1, newdata=na_dnb[,c('dnb'),drop=F])

