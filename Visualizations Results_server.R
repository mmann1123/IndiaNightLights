
# this scripts reads in raster files exported from grid_viirs_data (3).R
# it 1) reads in a dnb, zenith angle, moon phase etc data, 2) extracts this for all training sites (with voltage meters)
# 3) uses this info to train a ML model to predict outages 

# Run the following in bash before starting R
 module load proj.4/4.8.0
 module load gdal/gcc/1.11
 module load R/3.0.2
 module load gcc/4.9.0
 R




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
library(foreign)



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
#avoid problem image
azt = azt[azt!='2015226.2125_azt_v5.tif']


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
intersect_dates = intersect(intersect(time_stamp_dnb,time_stamp_cld),time_stamp_azt)
common_dnb = (time_stamp_dnb %in% intersect_dates)
dnb_stack = dnb_stack[[ (1:length(common_dnb))[common_dnb] ]]
common_cld = (time_stamp_cld %in% intersect_dates)
cld_stack = cld_stack[[ (1:length(common_cld))[common_cld] ]]
common_zen = (time_stamp_zen %in% intersect_dates)
zen_stack = zen_stack[[ (1:length(common_zen))[common_zen] ]]
common_azt = (time_stamp_azt %in% intersect_dates)
azt_stack = azt_stack[[ (1:length(common_azt))[common_azt] ]]



# remove cloud cells multicore  returns NA but runs fast!
library(foreach)
library(doParallel)
registerDoParallel(32)


# remove cloud cells multicore  returns NA but runs fast!
foreach(i=1:dim(dnb_stack)[3]) %dopar% { dnb_stack[[i]][cld_stack[[i]]>0]=NA}
save(dnb_stack,file = 'dnb_stack_wo_cld2.RData')
foreach(i=1:dim(zen_stack)[3]) %dopar% { zen_stack[[i]][cld_stack[[i]]>0]=NA}
save(zen_stack,file = 'zen_stack_wo_cld2.RData')
foreach(i=1:dim(azt_stack)[3]) %dopar% { azt_stack[[i]][cld_stack[[i]]>0]=NA}
save(azt_stack,file = 'azt_stack_wo_cld2.RData')



# Extract dnb values for all training data locations  -------------------------------------------

# load data
rm(list=ls())
setwd('/groups/manngroup/India\ VIIRS/2015')

load('dnb_stack_wo_cld2.RData')    # start here cloud free images stored here.
load('zen_stack_wo_cld2.RData')
load('azt_stack_wo_cld2.RData')


# extract local testing data 
# define locations of interest
locations = read.csv('/groups/manngroup/India VIIRS/Location Data/MH-ESMI-Locations-Lat-Long-Overpass-Cuts-May-2015-ag.csv')
locations2 = read.csv('/groups/manngroup/India VIIRS/Location Data/Other 28 ESMI MH Locations Coordinates.csv')
names(locations2) = c('STATE','DISTRICT.CITY','LOCATION','TYPE','LAT','LON','Ag.Rural')
locations2 = locations2[,!names(locations2) %in% 'TYPE']
jumba.df = data.frame(STATE='MH', DISTRICT.CITY="NA", LOCATION='Jumda',LAT=20.010094,LON=77.044271, Ag.Rural=T)
locations3 = read.dbf('/groups/manngroup/India VIIRS/Location Data/Emptyspaces.dbf')
names(locations3)=c('ID','LON','LAT')
locations3$STATE='MH'
locations3$LOCATION=paste('uninhabited',1:length(locations3$STATE),sep='')
locations=rbind.fill(locations,jumba.df)
locations=rbind.fill(locations,locations2)
locations=rbind.fill(locations,locations3)
head(locations,50)


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


dnb_values_loc=put_in_long_loc(dnb_values_loc,'dnb')
zen_values_loc= put_in_long_loc(zen_values_loc,'zen')
azt_values_loc=put_in_long_loc(azt_values_loc,'azt')


# add moon phase

# read in moon phase (year doy time moon_illum_frac moon_phase_angle)
phase = read.csv('moon_info.csv')
names(phase)=c('year','doy','time', 'illum', 'phase')
phase$date.time = paste(phase$year,sprintf('%03d',(phase$doy)),'.',phase$time,sep='')


# join in moon and dnb
dnb_values_loc = join(dnb_values_loc,zen_values_loc) # add moon characteristics to dnb values
dnb_values_loc = join(dnb_values_loc,azt_values_loc)
dnb_values_loc = join(dnb_values_loc,phase)
#head(dnb_values[20000:25000,],50)
dnb_values_loc[550:650,]





# LOAD DNB_VALUES --------------------------------------------------------------


# save output to load quickly
setwd('/groups/manngroup/India\ VIIRS/2015')

#save(dnb_values_loc,file='dnb_values_w_moon_loc2.RData')
load('dnb_values_w_moon_loc2.RData')


# merge training and testing data
dnb_values_loc$type = 'testing'
dnb_values  = dnb_values_loc
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
# demean the y variable to simulate FE
dnb_values$demean_dnb = dnb_values$dnb - dnb_values$mn_dnb


    

# Compare to actual outage data  ----------------------------------------------------------

# read in each xls file and convert to long format and write to csv
# only need to run this the first time... processing unformated files
# install.packages('readxl')
library(readxl)
file_list = list.files('..//Hourly Voltage Data//',pattern='Unified')
i=1

# read in hourly voltage data
a_location = read.csv(paste("..//Hourly Voltage Data//",file_list[i],sep=''), na = "NA",stringsAsFactors =F)
names(a_location) = c('anID','location','date','hour',paste(1:60))	
a_location = a_location[,c('location','date','hour',paste(1:60))]
head(a_location)

# two different time stamp types, put all into string version of '%Y-%m-%d'
new_date = format(strptime(a_location$date,'%m/%d/%Y',tz='Asia/Calcutta'),'%Y-%m-%d',usetz=F)
a_location$date[!is.na(new_date)] = new_date[!is.na(new_date)] 
	
# put data into long form and sort 
voltage =  melt(a_location,id=c('location','date','hour') )
names(voltage) = c('location','date','hour','minute','voltage')
voltage = voltage[with(voltage,order(location,date,hour)),]
head(voltage)

# go from 0-23 to 1-24 clock
voltage$hour =  as.numeric(voltage$hour+1)


# save the output
#save(voltage,file='..//Hourly Voltage Data//Voltage2.RData')
setwd('/groups/manngroup/India\ VIIRS/2015')
load('..//TestingData//Voltage2.RData')


# deal with date time 
# assign time stamps for voltage to Calcutta
#OlsonNames() # full list of time zones (only works on unix?)
voltage$date.hour.minute = paste(voltage$date,sprintf('%02d',voltage$hour), sprintf('%02d',voltage$minute),sep='.')
voltage$date.hour.minute = as.character(strptime(voltage$date.hour.minute, '%Y-%m-%d.%H.%M'))
voltage$date.hour.minute = as.POSIXct(voltage$date.hour.minute, tz="Asia/Calcutta")
head(voltage$date.hour.minute)
attributes(voltage$date.hour.minute)$tzone = 'UTC'  # convert Calcutta time to UTC
voltage$date.time2 = voltage$date.hour.minute
voltage$location = as.character(voltage$location)
head(voltage$date.hour.minute)


## assign UTC zone for DNB images
dnb_values$date.time2 = strptime(dnb_values$date.time,'%Y%j.%H%M')
dnb_values$date.time2 = as.POSIXct(dnb_values$date.time2, tz="UTC")
head(dnb_values$date.time2)


###################################################################
# match to closest time in local data 
# store location names for dnb and voltage data

locales = unique(dnb_values$location)
locales_v = unique(voltage$location)


# change names to match 
# find partial matches
for(i in 1:length(unique(voltage$location))){
	print(paste(locales[i],' -MATCH- ',locales_v[pmatch(locales[i], locales_v)]))}
# fill in holes create a dataframe as a lookup table
look_up = data.frame(locales =locales,
  	locales_v = apply(data.frame(locales),1,function(x) locales_v[pmatch(x, locales_v)]),
	stringsAsFactors=F)

head(look_up,40)

look_up[4,2] = locales_v[32]
look_up[5,2] = locales_v[31]
look_up[27,2] = locales_v[21]
head(look_up,40)


# switch names out use non-voltage data names 
for(i in 1:length(look_up$locales)){
	print(i)
	voltage$location[voltage$location == look_up$locales_v[i] ] = look_up$locales[i]
}
# double check that it worked
 sort(unique(dnb_values$location))
 sort(unique(voltage$location))


# create locale dnb and 
for(i in 1:length(locales)){
	locale = locales[i]
	# select a location and grab dnb and voltage data
	#test_site = resid_loc[resid_loc$location==locale,]
	test_site = dnb_values[dnb_values$location==locale,]
	test_voltage = voltage[voltage$location ==locale,]
	test_voltage = test_voltage[,c('date.time2','voltage')]
	# if locations are missing
	if( (dim(test_site)[1]==0 | dim(test_voltage)[1]==0) & !grepl('uninhabit',locale)   ){next}   # avoid missing data but not uninhabited  not finding voltage data for jumba
        # if locations are missing & the area is uninhabited
	# create an empty dataframe with dates that match dnb and fill voltage with 0s 
	if(dim(test_voltage)[1]==0 & grepl("uninhabit",locale) ){
		test_voltage = data.frame(date.time2=as.POSIXct(test_site$date.time2), voltage=0)
		test_voltage$date.time2 = as.POSIXct(test_site$date.time2, tz="UTC") # tz doesn't get assigned without this  	
		}   
	# join dnb and voltage data based on name and lowest time difference
	test_join = join(test_site, test_voltage[ sapply(test_site$date.time2,
                      function(x) which.min(abs(difftime(x, test_voltage$date.time2)))), ])
	head(test_join,5)
	if(i == 1){test_join_holder = test_join    # store data for later 
		}else{test_join_holder=rbind.fill(test_join_holder,test_join)}  # store data for use later
	test_join$count = 1:dim(test_join)[1]   # add an index for plotting
	test_join$lights_out = 'on'  # set what defines an outage 
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
	ggsave(a, file=paste('..//Voltage Plots//plot_',locale,'.png',sep=''), width=6, height=4)
	remove(test_join) # USE: to avoid mismatch between locations
	print(locale)
}



# save an image of all the data
#save.image(file='..//Hourly Voltage Data//AllData.RData')



# Train a classifier to find outages ----------------------------------------------
setwd('/groups/manngroup/India\ VIIRS/2015')
load('..//Hourly Voltage Data//AllData.RData')

# set up data for classifier (define what an outage is)
# read in joined dnb and voltage data where not missing
test_join_holder=test_join_holder[!is.na(test_join_holder$dnb)&!is.na(test_join_holder$voltage),]

# define classification classes,   voltage below 100v is considered an outage
hist(test_join_holder$voltage)
test_join_holder$lightsout = NA
test_join_holder$lightsout[ test_join_holder$voltage<100] = 2
test_join_holder$lightsout[ test_join_holder$voltage>=100] = 1
test_join_holder$lightsout[ grepl('uninhabit',test_join_holder$location) ] = 0

# 5 missing illum and zen
test_join_holder = na.omit(test_join_holder)

# store training and testing data
prop = 0.75 # proportion of subset data
set.seed(1234)
# training data set 
training.s = sample (1:nrow(test_join_holder), round(prop*nrow(test_join_holder),0))
test_join_holder.train = test_join_holder[training.s,]

# testing data set 
testing.s <- setdiff(1:nrow(test_join_holder), training.s)
test_join_holder.test <- test_join_holder[testing.s,]


# Tune the classifier ------------------------------------------------------------- 

# LOAD MULTICORE TUNE ----------------------------------------------------


# Choose tune parameters #ntree=number of trees, #mtry=# of features sampled for use at each node for splitting
rf_ranges = list(ntree=seq(1,200,10),mtry=10:30)

# features for random forest -  drop azth angle (angle from north, creates discontinuous predictions
form = factor(lightsout)~dnb+zen+illum+phase

source('..//IndiaNightLights//mctune.R')


# remove cloud cells multicore  returns NA but runs fast!
library(foreach)
library(doParallel)
registerDoParallel(32)


tuned.r3 = mctune(randomForest, train.x = form, data = test_join_holder.train,
         tunecontrol = tune.control(sampling = "cross",cross = 5), ranges=rf_ranges,
         mc.control=list(mc.cores=16, mc.preschedule=T),confusionmatrizes=T )

save(tuned.r3,file='..//Hourly Voltage Data//tunedTrees3.RData')


tuned.r5 = mctune(randomForest, train.x = form, data = test_join_holder.train,
         tunecontrol = tune.control(sampling = "cross",cross = 10), ranges=rf_ranges,
   	 mc.control=list(mc.cores=16, mc.preschedule=T),confusionmatrizes=T )

save(tuned.r5,file='..//Hourly Voltage Data//tunedTrees5.RData')





# LOAD REGRESSION TREE RESULTS ----------------------------------------------------


# read in stored regression trees look at performance 
tune.number = tuned.r5
# Get best parameters
tune.number$best.parameters
tune.number$best.confusionmatrizes



# Store the best model 
best.model = tune.number$best.model

predictions = predict(best.model, test_join_holder)
table.random.forest = table(test_join_holder$lightsout, predictions)
table.random.forest

# next function gives a graphical depiction of the marginal effect of a variable on the class probability (classification) or response (reg$
partialPlot(best.model, test_join_holder, dnb, "2")
partialPlot(best.model, test_join_holder, illum, "2")
partialPlot(best.model, test_join_holder, phase, "2")
partialPlot(best.model, test_join_holder, mn_dnb, "2")


# Calculate error rate
error.rate = 1 - sum(diag(as.matrix(table.random.forest))) / sum(table.random.forest)
error.rate



# PREDICT OUT OF SAMPLE ---------------------------------------------------------------

#remove(list=ls())
setwd('/groups/manngroup/India\ VIIRS/2015')

# load raster stacks
load('dnb_stack_wo_cld.RData')   
load('zen_stack_wo_cld.RData')
load('azt_stack_wo_cld.RData')
dnb_date_time = substr(names(dnb_stack),2,20)

# create sd_dnb and mn_dnb for each dnb time series
sd_dnb_layer = calc(dnb_stack,function(x){sd(x,na.rm=T)} )
mn_dnb_layer = calc(dnb_stack,function(x){mean(x,na.rm=T)} )

# read in moon phase (year doy time moon_illum_frac moon_phase_angle)
phase = read.csv('moon_info.csv')
names(phase)=c('year','doy','time', 'illum', 'phase')
phase$date_time = paste(phase$year,sprintf('%03d',(phase$doy)),'.',phase$time,sep='')

# load the tuned trees
load(file='..//TestingData//tunedTrees5.RData')
best.model = tuned.r5$best.model

# iterate through each date_time collect dnb, add sd_dnb illum phase
prediction_stack = stack()

# remove cloud cells multicore  returns NA but runs fast!
library(foreach)
library(doParallel)
registerDoParallel(32)


prediction_list = foreach(i = 1:length(dnb_date_time), .inorder=T,.packages=c('raster','e1071'),
	.errorhandling="remove") %dopar% {
	# packages must be loaded to each core, error handling removes output when 'next' is invoked, doesn't work in dopar
	#for (i in 1:1:length(dnb_date_time)){
	date_time = dnb_date_time[i]
	print(date_time)
	# find the dnb that matches date_time
	hold_stack_dnb = dnb_stack[[which(dnb_date_time %in% date_time)]]  
	hold_stack_zen = zen_stack[[which(dnb_date_time %in% date_time)]]  
        hold_stack_azt = azt_stack[[which(dnb_date_time %in% date_time)]]  
	# store illumination
	# avoid missing data
	if(length(phase[phase$date_time == date_time,'illum'])==0){next}
	hold_stack_illum = hold_stack_azt
	hold_stack_illum[] = phase[phase$date_time == date_time,'illum']
	# store phase
        hold_stack_phase = hold_stack_azt
        hold_stack_phase[] = phase[phase$date_time == date_time,'phase']
        # stack all data with appropriate names
	data_stack = stack(hold_stack_dnb,hold_stack_zen,hold_stack_azt,
		sd_dnb_layer,hold_stack_illum,hold_stack_phase,mn_dnb_layer)
	names(data_stack) = c('dnb','zen','azt','sd_dnb','illum','phase','mn_dnb')
	# predict to surface
	prediction = raster::predict(data_stack,best.model)
        names(prediction) = paste(date_time)
	writeRaster(prediction,paste('predictions//Prediction_',date_time,'.tif',sep=''),overwrite=T)
	return(prediction)
	}
prediction_list

# convert dopar list to raster stack
prediction_stack = stack(prediction_list)
save(prediction_stack, file='prediction_stack_w_cld.RData')


# remove clouds from prediction stack
# create raster stacks  & extract data
files = dir(pattern = '.tif')
cld = files[grep('cld_v5',files)]
cld_stack = stack(cld)
# find cloud masks that match prediction dates
time_stamp_cld = gsub(x=names(cld_stack),pattern = "(.*X)(.*)(.*_cld_v5)",replacement = "\\2")
time_stamp_prd = gsub(x=names(prediction_stack),pattern = "(.*X)(.*)",replacement = "\\2")
common_cld = (time_stamp_cld %in% intersect(time_stamp_cld,time_stamp_prd))
cld_stack = cld_stack[[ (1:length(common_cld))[common_cld] ]]


# Mask prediction stack
foreach(i=1:dim(prediction_stack)[3],.inorder=T,.packages=c('raster')) %dopar% { prediction_stack[[i]][cld_stack[[i]]>0]=NA}
save(prediction_stack, file='prediction_stack_wo_cld.RData')


#write out prediction mask wo clouds
#for (i in 1:dim(prediction_stack)[3]){
dumbout = foreach(i=1:dim(prediction_stack)[3],.inorder=F,.packages=c('raster'),.errorhandling='remove') %dopar% { 
   print(i)
   date_time = substr(names(prediction_stack)[i],2,15)
   writeRaster(prediction_stack[[i]],paste('predictions//Prediction_wo_cld',date_time,'.tif',sep=''),overwrite=T)
}
rm(dumbout)

# plot outages (2 = outage, 1 = normal, 0 = no electricity 
plot(prediction_stack[[11]])   # good example of outage in image #11


# IDEAS: GIVE CLASSIFIER TIME OF NIGHT TO SCALE ACORDING TO HOUR. 


# Estimate reliability ----------------------------------------------------

load(file='prediction_stack_wo_cld.RData')

outage_count = calc((prediction_stack==2),function(x){sum(x,na.rm=T)})
cloud_count = calc(is.na(prediction_stack),function(x){sum(x,na.rm=T)})
leng_days_raster = prediction_stack[[1]]
leng_days_raster[]=dim(prediction_stack)[3]
percent_outage = outage_count / (leng_days_raster-cloud_count)
writeRaster(percent_outage,'predictions/percent_outage.tif')




# Stable Lights Map -------------------------------------------------------

# read in data 
#rm(list=ls())
setwd('/groups/manngroup/India\ VIIRS/2015')
load('dnb_stack_wo_cld.RData')    # start here cloud free images stored here.

# calculate stats
# mean layer
mn_dnb_layer = calc(dnb_stack,function(x){mean(x,na.rm=T)} )
plot(log(mn_dnb_layer))
writeRaster(mn_dnb_layer,'predictions/mn_dnb_layer.tif')
# median layer
md_dnb_layer = calc(dnb_stack,function(x){median(x,na.rm=T)} )
plot(log(md_dnb_layer))
writeRaster(md_dnb_layer,'predictions/md_dnb_layer.tif', overwrite=T)
# quantile layer
qt_dnb_layer = calc(dnb_stack,function(x){quantile(x,na.rm=T)} )
names(qt_dnb_layer)=paste('Precentile',c(0,25,50,75,100))
writeRaster(qt_dnb_layer,'predictions/qt_dnb_layer.tif', overwrite=T)



# Write out files to KML ------------------------------------------------
# this works but is pretty low resolution
KML(md_dnb_layer,'predictions/md_dnb_layer',blur=50,overwrite=T)

KML(percent_outage,'predictions/percent_outage',blur=50,overwrite=T, col=rev(rainbow(25)))




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



#set.seed(101)
#alpha     = 0.8 # percentage of training set
#inTrain   = sample(1:nrow(test_join_holder), alpha * nrow(test_join_holder))
#train.set = test_join_holder[inTrain,]
#test.set  = test_join_holder[-inTrain,]



