
# Run the following in bash before starting R
# module load proj.4/4.8.0
# module load gdal/gcc/1.11
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
library(e1071)
library(randomForest)
library(party)

# this scripts reads in raster files exported from grid_viirs_data (3).R



# Time Series Plots for locations of interest -----------------------------



# Lunar adjustments ------------------------------------------------------
# try to remove cyclical lunar signal from cells

# read in data
setwd('/groups/manngroup/India\ VIIRS/2015')

# pull available files
files = dir(pattern = '.tif')
cld = files[grep('cld_v5',files)]
dnb = files[grep('dnb_v5',files)]
zen = files[grep('zen_v5',files)]
azt = files[grep('azt_v5',files)]


# create raster stacks  & extract data
cld_stack = stack(cld)
dnb_stack = stack(dnb)
zen_stack = stack(zen)
azt_stack = stack(azt)


# Extract dates
time_stamp_dnb = gsub(x=names(dnb_stack),pattern = "(.*X)(.*)(.*_dnb_v5)",replacement = "\\2")
time_stamp_dnb = strptime(time_stamp_dnb,"%Y%j.%H%M", tz = 'UTC')
names(dnb_stack) = format(time_stamp_dnb,"%Y%j.%H%M",usetz=F)
time_stamp_cld = gsub(x=names(cld_stack),pattern = "(.*X)(.*)(.*_cld_v5)",replacement = "\\2")
time_stamp_cld = strptime(time_stamp_cld,"%Y%j.%H%M", tz = 'UTC')
names(cld_stack) = format(time_stamp_cld,"%Y%j.%H%M",usetz=F)
time_stamp_zen = gsub(x=names(zen_stack),pattern = "(.*X)(.*)(.*_zen_v5)",replacement = "\\2")
time_stamp_zen = strptime(time_stamp_zen,"%Y%j.%H%M", tz = 'UTC')
names(zen_stack) = format(time_stamp_zen,"%Y%j.%H%M",usetz=F)
time_stamp_azt = gsub(x=names(azt_stack),pattern = "(.*X)(.*)(.*_azt_v5)",replacement = "\\2")
time_stamp_azt = strptime(time_stamp_azt,"%Y%j.%H%M", tz = 'UTC')
names(azt_stack) = format(time_stamp_azt,"%Y%j.%H%M",usetz=F)

# TEST: not all stacks have same dates
all.equal(time_stamp_dnb,time_stamp_cld)
all.equal(time_stamp_dnb,time_stamp_zen)
all.equal(time_stamp_dnb,time_stamp_azt)


# limit stacks to common elements  NOT NEEDED IF ALL SAME TIMES
#common_dnb = (time_stamp_dnb %in% intersect(time_stamp_dnb,time_stamp_cld))
#dnb_stack = dnb_stack[[ (1:length(common_dnb))[common_dnb] ]]
#common_cld = (time_stamp_cld %in% intersect(time_stamp_dnb,time_stamp_cld))
#cld_stack = cld_stack[[ (1:length(common_cld))[common_cld] ]]


# remove cloud cells multicore  returns NA but runs fast!
library(foreach)
library(doParallel)
registerDoParallel(32)


# remove cloud cells multicore  returns NA but runs fast!
foreach(i=1:dim(dnb_stack)[3]) %do% { dnb_stack[[i]][cld_stack[[i]]>0]=NA}
save(dnb_stack,file = 'dnb_stack_wo_cld.RData')
foreach(i=1:dim(zen_stack)[3]) %do% { zen_stack[[i]][cld_stack[[i]]>0]=NA}
save(zen_stack,file = 'zen_stack_wo_cld.RData')
foreach(i=1:dim(azt_stack)[3]) %do% { azt_stack[[i]][cld_stack[[i]]>0]=NA}
save(azt_stack,file = 'azt_stack_wo_cld.RData')



# Demeaned regression on 1000 point sample -------------------------------------------


# load data
rm(list=ls())
setwd('/groups/manngroup/India\ VIIRS/2015')

load('dnb_stack_wo_cld.RData')    # start here cloud free images stored here.
load('zen_stack_wo_cld.RData')
load('azt_stack_wo_cld.RData')


# create a plygon to create samples within
source('/groups/manngroup/India VIIRS/IndiaNightLights/PolygonFromExtent.R')


# create sample of 10000 set seed for reproducibility
set.seed(2)
sampled = spsample(PolygonFromExtent(extent(72, 81.50, 15, 22.5),crs=CRS('+proj=longlat +ellps=WGS84 
	+datum=WGS84 +no_defs ')), n=10000,type='random')


# extract data (and surrounding area)
dnb_values = extract(dnb_stack,sampled,fun= function(x) mean(x,na.rm=T), df=T)#buffer=1.2e3,
zen_values = extract(zen_stack,sampled,fun= function(x) mean(x,na.rm=T), df=T)
azt_values = extract(azt_stack,sampled,fun= function(x) mean(x,na.rm=T), df=T)
time_stamp_extract = gsub(x=colnames(dnb_values),pattern = "(.*X)(.*)(.*)",replacement = "\\2")

# change time stamps 
names(dnb_values) = time_stamp_extract
names(zen_values) = time_stamp_extract
names(azt_values) = time_stamp_extract


# extract local testing data 
# define locations of interest
locations = read.csv('MH-ESMI-Locations-Lat-Long-Overpass-Cuts-May-2015-ag.csv')
locations2 = read.csv('/groups/manngroup/India VIIRS/TestingData/Other 28 ESMI MH Locations Coordinates.csv')
names(locations2) = c('STATE','DISTRICT.CITY','LOCATION','TYPE','LAT','LON','Ag.Rural')
locations2 = locations2[,!names(locations2) %in% 'TYPE']
jumba.df = data.frame(STATE='MH', DISTRICT.CITY="NA", LOCATION='Jumda',LAT=20.010094,LON=77.044271, Ag.Rural=T)
locations=rbind.fill(locations,jumba.df)
locations=rbind.fill(locations,locations2)
head(locations)


# convert to spatial objects and extract values
coordinates(locations)= ~LON+LAT
proj4string(locations) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
dnb_values_loc = extract(dnb_stack,locations,fun= function(x) mean(x,na.rm=T), df=T)#buffer=1.2e3,
zen_values_loc = extract(zen_stack,locations,fun= function(x) mean(x,na.rm=T), df=T)
azt_values_loc  = extract(azt_stack,locations,fun= function(x) mean(x,na.rm=T), df=T)
time_stamp_extract_loc = gsub(x=colnames(dnb_values_loc),pattern = "(.*X)(.*)(.*)",replacement = "\\2")


# change time stamps
names(dnb_values_loc) = time_stamp_extract_loc
names(zen_values_loc) = time_stamp_extract_loc
names(azt_values_loc) = time_stamp_extract_loc




# put into long form
put_in_long2 <- function(wide_data,abreviation){
        names(wide_data) = time_stamp_extract
        wide_data$location = row.names(wide_data)
        wide_data = subset(wide_data,select=-c(ID))
        wide_data <- melt(wide_data )
        names(wide_data)=c('location','date.time',paste(abreviation))
        head(wide_data)
	if(class(wide_data$date.time)=='factor'){wide_data$date.time=as.character(wide_data$date.time)}
        return(wide_data)
}
put_in_long_loc <- function(wide_data,abreviation){
        names(wide_data) = time_stamp_extract_loc
        wide_data$location = locations$LOCATION
        wide_data = subset(wide_data,select=-c(ID)) # extra variable was added, remove
        wide_data <- melt(wide_data )
        names(wide_data)=c('location','date.time',paste(abreviation))
        if(class(wide_data$date.time)=='factor'){wide_data$date.time=as.character(wide_data$date.time)}
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
head(dnb_values[20000:25000,],50)
head(dnb_values_loc[50:60,])



# LOAD DNB_VALUES --------------------------------------------------------------


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
head(dnb_values)


# create descriptive statistics
sd_dnb = aggregate(dnb~location,data=dnb_values,function(x) sd(x,na.rm=T))
mn_dnb = aggregate(dnb~location,data=dnb_values,function(x) mean(x,na.rm=T))
names(sd_dnb)=c('location','sd_dnb')
names(mn_dnb)=c('location','mn_dnb')
dnb_values = join(dnb_values, mn_dnb)
dnb_values = join(dnb_values, sd_dnb)


# New KMEANS out of sample   # NOTE K MEAN CAN BE ON MEAN OR ACTUAL VALUES
source('/groups/manngroup/India VIIRS/IndiaNightLights/clusters.R')
variable_to_cluster =c('sd_dnb')
dnb_num =  dnb_values[variable_to_cluster] # set up a dataframe to store row numbers so can be joined later
kmn2 = kmeans(na.omit(dnb_num),5)
# calculate cluster membership out of sample for large dataset
library(foreach)
library(doParallel)
registerDoParallel(32)
iterator_groups = split(1:length(dnb_values[,variable_to_cluster]), cut(1:length(dnb_values[,variable_to_cluster]),32))   
dnb_values$kmn2 = foreach(group=1:length(iterator_groups), .inorder=T, .combine = c) %dopar% { 
  x1 = data.frame(dnb_values[,variable_to_cluster][iterator_groups[[group]]])
  clusters(x1,centers = kmn2[["centers"]])
}
stopImplicitCluster()


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

#mean_lm = lm(demean_dnb~0+I(mn_dnb<1.099e-8)*(ns(zen*azt,3)+
#      ns(zen*phase,3)+ns(azt*phase,3)+ns(phase,4)),data=na.omit(dnb_values[dnb_values$type=='training',])) # omit intercept
#
#mean_lm = lm(demean_dnb~0+I(mn_dnb<1.099e-8)*(ns(zen*azt,3)+   
#      ns(zen*phase,3)+ns(azt*phase,3)+ns(phase,3)+ns(illum,3))+I(mn_dnb<1e-7)*(ns(zen*azt,3)+
#      ns(zen*phase,3)+ns(azt*phase,3)+ns(phase,3)+ns(illum,3)),data=na.omit(dnb_values[dnb_values$type=='training',])) # omit intercept

########### REGRESSIONS FOR V4 DATA ##################
mean_lm = lm(demean_dnb~0+I(factor(kmn2))*ns(zen*azt,3)+ns(zen*phase,3)+ns(azt*phase,3)+ns(phase,4),
           data=na.omit(dnb_values[dnb_values$type=='training',])) # omit intercept

summary(mean_lm)


# compare actual and predicted
remove(actual, pred, resid)
actual = na.omit(dnb_values[dnb_values$type=='training',])
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

set.seed(11)  # 2 has declining values   11 has a zero light pixel
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
resid_loc$dnb = mean_lm_loc$residuals  +resid_loc$mn_dnb
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
# only need to run this the first time... processing unformated files
# install.packages('readxl')
# library(readxl)
#file_list = list.files('..//TestingData//',pattern='2015.xlsx')
#for(i in 1:length(file_list)){
#	print(i)
#	a_location = read_excel(paste("..//TestingData//",file_list[i],sep=''), na = "NA")
#	names(a_location) = c('location','date','hour',paste(1:60))	
#	head(a_location)
#
#	if(file_list[i]=='Other 28 ESMI MH Voltage Data 2015.xlsx'){
#		a_location$date = 	# convert file to proper date time 
#		format(strptime(a_location$date,'%m/%d/%Y',tz='Asia/Calcutta'),'%Y-%m-%d',usetz=F)
#		}
#	# put data into long form and sort 
#	longs =  melt(a_location,id=c('location','date','hour') )
#	names(longs) = c('location','date','hour','minute','voltage')
#	longs = longs[with(longs,order(location,date,hour)),]
#	head(longs)
#	#http://www.regexr.com/
#	write.csv(longs,paste('..//TestingData//',gsub('(.[a-z]+$)','\\2',file_list[i]),'.csv',sep=''))
#}
#
#
#file_list = list.files('..//TestingData//',pattern='2015.csv')
#voltage = read.csv(paste('..//TestingData//',file_list[1],sep=''))
#voltage$hour = as.numeric(voltage$hour) + 1
#print(unique(voltage$hour))
#head(voltage)
#for(i in 2:length(file_list)){
#	inner =  read.csv(paste('..//TestingData//',file_list[i],sep=''))
#	print(class(inner$hour)) 
#	inner$hour = as.numeric(inner$hour) + 1
#	print(unique(inner$hour))
#	voltage = rbind.fill(voltage,inner)
#}


# save the output
#save(voltage,file='..//TestingData//Voltage.RData')
setwd('/groups/manngroup/India\ VIIRS/2015')
load('..//TestingData//Voltage.RData')


# deal with date time
#OlsonNames() # full list of time zones (only works on unix?)
voltage$date.hour.minute = paste(voltage$date,sprintf('%02d',voltage$hour), sprintf('%02d',voltage$minute),sep='.')
voltage$date.hour.minute = as.character(strptime(voltage$date.hour.minute, '%Y-%m-%d.%H.%M'))
voltage$date.hour.minute = as.POSIXct(voltage$date.hour.minute, tz="Asia/Calcutta")
head(voltage$date.hour.minute)
attributes(voltage$date.hour.minute)$tzone = 'UTC'  # convert Calcutta time to UTC
voltage$date.time2 = voltage$date.hour.minute
voltage$location = as.character(voltage$location)
head(voltage$date.hour.minute)


#same for all locations 
resid_loc$date.time2 = strptime(resid_loc$date.time,'%Y%j.%H%M')
resid_loc$date.time2 = as.POSIXct(resid_loc$date.time2, tz="UTC")
head(resid_loc$date.time2)


# match to closest time in local data 
locales = unique(resid_loc$location)
locales_v = unique(voltage$location)


# change names to match 
# find partial matches
for(i in 1:length(unique(voltage$location))){
	print(paste(locales[i],' -MATCH- ',locales_v[pmatch(locales[i], locales_v)]))}
# fill in holes create a dataframe as a lookup table
look_up = data.frame(locales =locales,
  	locales_v = apply(data.frame(locales),1,function(x) locales_v[pmatch(x, locales_v)]),
	stringsAsFactors=F)
look_up[4,2] = locales_v[37]
look_up[5,2] = locales_v[36]
look_up[27,2] = locales_v[23]
look_up


# switch names out use non-voltage data names 
for(i in 1:length(look_up$locales)){
	print(i)
	voltage$location[voltage$location == look_up$locales_v[i] ] = look_up$locales[i]
}
# double check that it worked
 sort(unique(resid_loc$location))
 sort(unique(voltage$location))


# create locale dnb and 
for(i in 1:length(locales)){
	locale = locales[i]
	# select a location
	test_site = resid_loc[resid_loc$location==locale,]
	test_voltage = voltage[voltage$location ==locale,]
	if(dim(test_voltage)[1]==0){next}   # avoid missing data
	# join dnb and voltage data based on name and lowest time difference
	test_join = join(test_site, test_voltage[ sapply(test_site$date.time2,
                      function(x) which.min(abs(difftime(x, test_voltage$date.time2)))), ])
	head(test_join,5)
	if(i == 1){test_join_holder = test_join    # store data for later 
		}else{test_join_holder=rbind.fill(test_join_holder,test_join)}
	test_join$count = 1:dim(test_join)[1]
	test_join$lights_out = 'on'
	test_join$lights_out[test_join$voltage<=100]='off'  #test_join$voltage==0,131 used by ngo
	test_join$lights_out[is.na(test_join$voltage)]='No Data'
	if(length(unique(test_join$lights_out))==3){color_codes=c("#CC00CC","#FF0000","#00CC00")}
        if(length(unique(test_join$lights_out))==2){color_codes=c("#CC00CC","#00CC00")}
        if(length(unique(test_join$lights_out))==1){color_codes=c("#CC00CC","#00CC00")}

	a = ggplot(test_join, aes(count, dnb,colour=lights_out)) + 
	     #geom_vline(xintercept = test_join$count[test_join$voltage==0],colour='red',alpha=0.25,, show_guide=T)+
	     geom_vline(xintercept = test_join$count[is.na(test_join$voltage)],colour='blue',alpha=0.25, show_guide=T)+
             geom_vline(xintercept = test_join$count[test_join$illum>50],colour='yellow',alpha=0.25, show_guide=T)+
	     geom_point()+ scale_fill_manual(values = color_codes)+
	ggtitle(locale)
	ggsave(a, file=paste('..//TestingData//plot_',locale,'.png',sep=''), width=6, height=4)
	remove(test_join) # USE: to avoid mismatch between locations
	print(locale)
}



# Train a classifier to find outages ----------------------------------------------
#save(.,file='..//TestingData//Voltage.RData')

test_join_holder=test_join_holder[!is.na(test_join_holder$dnb)&!is.na(test_join_holder$voltage),]
test_join_holder$lightsout = test_join_holder$voltage<100

model = randomForest(factor(lightsout)~dnb+zen+azt+illum+phase+I(dnb^2)+I(zen^2)+I(azt^2)+I(illum^2)+
	illum:zen+illum:azt+sd_dnb+I(sd_dnb*2)+I(sd_dnb*3)+I(sd_dnb*4)+I(sd_dnb*5)+I(sd_dnb*6), 
	data=test_join_holder, ntree=3000,importance=T,do.trace = 100) 
model$confusion  #91
importance(model)
varImpPlot(model)


# Tune the classifier ------------------------------------------------------------- 

set.seed(101)
alpha     = 0.8 # percentage of training set
inTrain   = sample(1:nrow(test_join_holder), alpha * nrow(test_join_holder))
train.set = test_join_holder[inTrain,]
test.set  = test_join_holder[-inTrain,]

form = factor(lightsout)~dnb+zen+azt+illum+phase+I(dnb^2)+I(zen^2)+I(azt^2)+I(illum^2)+
        illum:zen+illum:azt+sd_dnb+I(sd_dnb*2)+I(sd_dnb*3)+I(sd_dnb*4)+I(sd_dnb*5)+I(sd_dnb*6)

model = randomForest(form,data=train.set, ntree=3000,importance=T,do.trace = 100)
model$confusion  

# next function gives a graphical depiction of the marginal effect of a variable on the class probability (classification) or response (regression).
# partialPlot(model, train.set, dnb, "TRUE")
# partialPlot(model, train.set, illum, "TRUE")
# partialPlot(model, train.set, phase, "TRUE")

# histogram of treesize
hist(treesize(model))

# Choose tune parameters #ntree=number of trees, #mtry=# of features sampled for use at each node for splitting
rf_ranges = list(ntree=seq(1,1000,200),mtry=7:25)

# Tune the tree, multicore 
tuned.r = tune(randomForest, train.x = form, data = train.set, validation.x = test.set,
         tunecontrol = tune.control(sampling = "fix",fix = 9/10), ranges=rf_ranges)
tuned.r

# Get best parameters
tuned.r$best.parameters

# Store the best model 
best.model = tuned.r$best.model
predictions = predict(best.model, test.set)
table.random.forest = table(test.set$lightsout, predictions)
table.random.forest


# Calculate error rate
error.rate <- 1 - sum(diag(as.matrix(table.random.forest))) / sum(table.random.forest)
error.rate




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

