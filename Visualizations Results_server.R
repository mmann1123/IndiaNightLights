
# this scripts reads in raster files exported from grid_viirs_data_server_2015.R
# it 1) reads in a dnb, zenith angle, moon phase etc data, 2) extracts this for all training sites (with voltage meters)
# 3) uses this info to train a ML model to predict outages 

# Run the following in bash before starting R
 module load proj.4/4.8.0
 module load gdal/gcc/1.11
 #module use /home/mmann1123/local/modulefiles
 #module load R/3.2.2
 #module load R/3.0.2
 module load R
 module load gcc/4.9.0
 R


library(raster)
library(ggplot2)
library(ggthemes)
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
library(caret)

# A: Stack data remove clouds, correct times ------------------------------------------------------
rm(list=ls())

# read in data
 #setwd('/groups/manngroup/India\ VIIRS/2015') #server
 setwd('R:\\Mann_Research\\India VIIRS SERVER/2015')

# pull available files
files = dir(pattern = '.tif')
cld = files[grep('cld_v5',files)]
dnb = files[grep('dnb_v5',files)]
#zen = files[grep('zen_v5',files)]
#azt = files[grep('azt_v5',files)]
#avoid problem image
#azt = azt[azt!='2015226.2125_azt_v5.tif']


# create raster stacks  & extract data
cld_stack = stack(cld)
dnb_stack = stack(dnb)
#zen_stack = stack(zen)
#azt_stack = stack(azt)


# Extract dates
time_stamp_dnb = gsub(x=names(dnb_stack),pattern = "(.*X)(.*)(.*_dnb_v5)",replacement = "\\2")
time_stamp_dnb = strptime(time_stamp_dnb,"%Y%j.%H%M", tz = 'UTC')
names(dnb_stack) = format(time_stamp_dnb,"%Y%j.%H%M",usetz=F)
time_stamp_cld = gsub(x=names(cld_stack),pattern = "(.*X)(.*)(.*_cld_v5)",replacement = "\\2")
time_stamp_cld = strptime(time_stamp_cld,"%Y%j.%H%M", tz = 'UTC')
names(cld_stack) = format(time_stamp_cld,"%Y%j.%H%M",usetz=F)
#time_stamp_zen = gsub(x=names(zen_stack),pattern = "(.*X)(.*)(.*_zen_v5)",replacement = "\\2")
#time_stamp_zen = strptime(time_stamp_zen,"%Y%j.%H%M", tz = 'UTC')
#names(zen_stack) = format(time_stamp_zen,"%Y%j.%H%M",usetz=F)
#time_stamp_azt = gsub(x=names(azt_stack),pattern = "(.*X)(.*)(.*_azt_v5)",replacement = "\\2")
#time_stamp_azt = strptime(time_stamp_azt,"%Y%j.%H%M", tz = 'UTC')
#names(azt_stack) = format(time_stamp_azt,"%Y%j.%H%M",usetz=F)

# TEST: not all stacks have same dates
all.equal(time_stamp_dnb,time_stamp_cld)
#all.equal(time_stamp_dnb,time_stamp_zen)
#all.equal(time_stamp_dnb,time_stamp_azt)


# limit stacks to common elements  NOT NEEDED IF ALL SAME TIMES
intersect_dates = intersect(time_stamp_dnb,time_stamp_cld)
#intersect_dates = intersect(intersect(time_stamp_dnb,time_stamp_cld),time_stamp_azt)
common_dnb = (time_stamp_dnb %in% intersect_dates)
dnb_stack = dnb_stack[[ (1:length(common_dnb))[common_dnb] ]]
common_cld = (time_stamp_cld %in% intersect_dates)
cld_stack = cld_stack[[ (1:length(common_cld))[common_cld] ]]
#common_zen = (time_stamp_zen %in% intersect_dates)
#zen_stack = zen_stack[[ (1:length(common_zen))[common_zen] ]]
#common_azt = (time_stamp_azt %in% intersect_dates)
#azt_stack = azt_stack[[ (1:length(common_azt))[common_azt] ]]



# remove cloud cells multicore  returns NA but runs fast!
library(foreach)
library(doParallel)
registerDoParallel(5)


# remove cloud cells multicore  returns NA but runs fast!
foreach(i=1:dim(dnb_stack)[3]) %dopar% { 
  
  # ASSIGNMENT IN A DOPAR LOOP DOESN'T WORK!!!!!!!!!!!!!!!!!
  
  dnb_stack[[i]][cld_stack[[i]]>0]=NA}
save(dnb_stack,file = 'dnb_stack_wo_cld2.RData')
#foreach(i=1:dim(zen_stack)[3]) %dopar% { zen_stack[[i]][cld_stack[[i]]>0]=NA}
#save(zen_stack,file = 'zen_stack_wo_cld2.RData')
#foreach(i=1:dim(azt_stack)[3]) %dopar% { azt_stack[[i]][cld_stack[[i]]>0]=NA}
#save(azt_stack,file = 'azt_stack_wo_cld2.RData')
save(cld_stack,file = 'cld_stack.RData')

# add global mean
Gmn_dnb_stack=dnb_stack
Gmn = foreach(i=1:dim(dnb_stack)[3], .combine='c') %dopar% {cellStats(dnb_stack[[i]],mean,rm.na=T)}
for(i in 1:dim(dnb_stack)[3]){Gmn_dnb_stack[[i]][]=Gmn[i]}
save(Gmn_dnb_stack,file = 'Gmn_dnb_stack_wo_cld2.RData')

# add global median
Gmd_dnb_stack=dnb_stack
Gmd = foreach(i=1:dim(dnb_stack)[3],.combine='c') %dopar% {cellStats(dnb_stack[[i]],median,rm.na=T) }
for(i in 1:dim(dnb_stack)[3]){Gmd_dnb_stack[[i]][]=Gmd[i]}
save(Gmd_dnb_stack,file = 'Gmd_dnb_stack_wo_cld2.RData')

save.image(file='..//AllData_A.RData')


# B: Extract dnb values for all training data locations  -------------------------------------------

# load data
rm(list=ls())
#setwd('/groups/manngroup/India\ VIIRS/2015') #server
setwd('R:\\Mann_Research\\India VIIRS SERVER/2015')

load('dnb_stack_wo_cld2.RData')    # start here cloud free images stored here.
#load('zen_stack_wo_cld2.RData')
#load('azt_stack_wo_cld2.RData')
load('cld_stack.RData')
load('Gmd_dnb_stack_wo_cld2.RData')
 

# extract local testing data 
# define locations of interest
locations = read.csv('../Location Data/MH-ESMI-Locations-Lat-Long-Overpass-Cuts-May-2015-ag.csv')
locations2 = read.csv('../Location Data/Other 28 ESMI MH Locations Coordinates.csv')
names(locations2) = c('STATE','DISTRICT.CITY','LOCATION','TYPE','LAT','LON','Ag.Rural')
locations2 = locations2[,!names(locations2) %in% 'TYPE']
jumba.df = data.frame(STATE='MH', DISTRICT.CITY="NA", LOCATION='Jumda',LAT=20.010094,LON=77.044271, Ag.Rural=T)
locations3 = read.dbf('../Location Data/Emptyspaces.dbf')
names(locations3)=c('ID','LON','LAT')
locations3$STATE='MH'
locations3$LOCATION=paste('uninhabited',1:length(locations3$STATE),sep='')
locations4 = read.csv('../Location Data/4 TS Hyderabad ESMI Locations with Coordinates.csv')
names(locations4) = c('STATE','DISTRICT.CITY','LOCATION','LAT','LON','Ag.Rural')

locations=rbind.fill(locations,jumba.df)
locations=rbind.fill(locations,locations2)
locations=rbind.fill(locations,locations3)
locations=rbind.fill(locations,locations4)
locations=subset(locations,select=-c(ID))
head(locations,50)


# convert to spatial objects and extract values
coordinates(locations)= ~LON+LAT
proj4string(locations) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
save(locations,file='locationsXY.RData')
dnb_values_loc = extract(dnb_stack,locations, df=T)
#zen_values_loc = extract(zen_stack,locations,  df=T)
#azt_values_loc  = extract(azt_stack,locations, df=T)
Gmd_dnb_values_loc = extract(Gmd_dnb_stack,locations, df=T)
Gmn_dnb_values_loc = extract(Gmn_dnb_stack,locations, df=T)     
cld_values_loc = extract(cld_stack,locations, df=T)
time_stamp_extract_loc = gsub(x=colnames(dnb_values_loc),pattern = "(.*X)(.*)(.*)",replacement = "\\2")


# change time stamps
names(dnb_values_loc) = time_stamp_extract_loc
#names(zen_values_loc) = time_stamp_extract_loc
#names(azt_values_loc) = time_stamp_extract_loc
names(Gmn_dnb_values_loc) = time_stamp_extract_loc
names(Gmd_dnb_values_loc) = time_stamp_extract_loc
names(cld_values_loc) = time_stamp_extract_loc


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
#zen_values_loc= put_in_long_loc(zen_values_loc,'zen')
#azt_values_loc=put_in_long_loc(azt_values_loc,'azt')
cld_values_loc=put_in_long_loc(cld_values_loc,'cld')
Gmd_dnb_values_loc=put_in_long_loc(Gmd_dnb_values_loc,'Gmd')
Gmn_dnb_values_loc=put_in_long_loc(Gmn_dnb_values_loc,'Gmn')


# add moon phase

# read in moon phase (year doy time moon_illum_frac moon_phase_angle)
phase = read.csv('moon_info.csv')
names(phase)=c('year','doy','time', 'illum', 'phase')
phase$date.time = paste(phase$year,sprintf('%03d',(phase$doy)),'.',phase$time,sep='')


# join in moon and dnb
#dnb_values_loc = join(dnb_values_loc,zen_values_loc) # add moon characteristics to dnb values
#dnb_values_loc = join(dnb_values_loc,azt_values_loc)
dnb_values_loc = join(dnb_values_loc,Gmd_dnb_values_loc)
dnb_values_loc = join(dnb_values_loc,Gmn_dnb_values_loc)
dnb_values_loc = join(dnb_values_loc,phase)
dnb_values_loc = join(dnb_values_loc,cld_values_loc)

dnb_values_loc[550:650,]

save.image(file='..//AllData_B.RData')


#C:  LOAD DNB_VALUES --------------------------------------------------------------


# save output to load quickly
#setwd('/groups/manngroup/India\ VIIRS/2015')
setwd('R:\\Mann_Research\\India VIIRS SERVER/2015')

#save(dnb_values_loc,file='dnb_values_w_moon_loc2.RData')
load('dnb_values_w_moon_loc2.RData')


# merge training and testing data
dnb_values_loc$type = 'testing'
dnb_values  = dnb_values_loc
head(dnb_values)


# create descriptive statistics
sd_dnb = aggregate(dnb~location,data=dnb_values,function(x) sd(x,na.rm=T))
mn_dnb = aggregate(dnb~location,data=dnb_values,function(x) mean(x,na.rm=T))
md_dnb = aggregate(dnb~location,data=dnb_values,function(x) median(x,na.rm=T))
names(sd_dnb)=c('location','sd_dnb')
names(mn_dnb)=c('location','mn_dnb')
names(md_dnb)=c('location','md_dnb')
dnb_values = join(dnb_values, mn_dnb)
dnb_values = join(dnb_values, sd_dnb)
dnb_values = join(dnb_values, md_dnb)


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


save.image(file='..//AllData_C.RData')
    

# D: Compare to actual outage data  ----------------------------------------------------------
#load('..//AllData_C.RData')

# read in each xls file and convert to long format and write to csv
# only need to run this the first time... processing unformated files
library(readxl)
file_list = list.files('..//Hourly Voltage Data//',pattern='Unified')
file_list

# read in hourly voltage data
a_location = read.csv(paste("..//Hourly Voltage Data//",file_list[1],sep=''), na = "NA",stringsAsFactors =F)
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

save.image(file='..//AllData_D.RData')



# E: CORRECT NAMES FOR NEW LOCATION ------------------------------------------
load('..//AllData_D.RData')

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
look_up

look_up[4,2] = locales_v[34]
look_up[5,2] = locales_v[33]
look_up[27,2] = locales_v[23]
look_up



# switch names out use non-voltage data names 
for(i in 1:length(look_up$locales)){
	print(i)
	voltage$location[voltage$location == look_up$locales_v[i] ] = look_up$locales[i]
}
# double check that spellings are same on both sheets
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
	# create an empty dataframe with dates that match dnb and fill voltage with os 
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
save.image(file='..//Hourly Voltage Data//AllData.RData')



# Train a classifier to find outages ----------------------------------------------
#setwd('/groups/manngroup/India\ VIIRS/2015')
setwd('R:\\Mann_Research\\India VIIRS SERVER/2015')

load('..//Hourly Voltage Data//AllData.RData')

# set up data for classifier (define what an outage is)
# read in joined dnb and voltage data where not missing
test_join_holder=test_join_holder[!is.na(test_join_holder$dnb)&!is.na(test_join_holder$voltage),]

# define classification classes,   voltage below 100v is considered an outage
hist(test_join_holder$voltage)
test_join_holder$lightsout = NA
test_join_holder$lightsout[ test_join_holder$voltage<100  ] = 2   # outage 
test_join_holder$lightsout[ test_join_holder$voltage>=100 & test_join_holder$voltage<1000  ] = 1 # normal
test_join_holder$lightsout[ grepl('uninhabit',test_join_holder$location) ] = 0   # uninhabited 
# remove two problem sites with extremely high outage rates (over 50%)
test_join_holder = test_join_holder[!(test_join_holder$location=='Chandikapur'|test_join_holder$location=='Kanheri Sarap'),]


# 5 missing illum and zen
test_join_holder = na.omit(test_join_holder)

# store training and testing data
prop = 0.75 # proportion of subset data
set.seed(1234)

# training data set 
training.s = sample (1:nrow(test_join_holder), round(prop*nrow(test_join_holder),0))
test_join_holder.train = test_join_holder[training.s,]

# testing data set 
testing.s = setdiff(1:nrow(test_join_holder), training.s)
test_join_holder.test = test_join_holder[testing.s,]

table(test_join_holder.train$lightsout)
table(test_join_holder.test$lightsout)



# Tune the classifier ------------------------------------------------------------- 

# Load multicore tuner
source('..//IndiaNightLights//mctune.R')

# Choose tune parameters #ntree=number of trees, #mtry=# of features sampled for use at each node for splitting
rf_ranges = list(ntree=seq(1,40,1),mtry=5:35)

# features for random forest -  drop azth angle (angle from north, creates discontinuous predictions

form = factor(lightsout)~dnb+illum+phase+md_dnb+I(dnb-md_dnb) #dont use: sd_dnb mn_dnb+I(dnb-mn_dnb)

# remove cloud cells multicore  returns NA but runs fast!
library(foreach)
library(doParallel)
registerDoParallel(4)

# custom error function
custom_error = function(y,pred){
  # this function is minimized if set in mc.control. compare actual rates to predicted outage rates
  # takes two parameters, actualvalues and predicted values
  # 2 = outage, 1 = normal, 0 = uninhabited

  # DOESN'T TAKE INTO ACCOUNT CLOUDS!
  y.outrate = sum(y==2,na.rm=T)/length(y)
  pred.outrate = sum(pred==2,na.rm=T)/length(pred)
  return((y.outrate-pred.outrate)^2)	
}

tuned.r10 = tune(randomForest, train.x = form, data = test_join_holder.train,
         tunecontrol = tune.control(sampling = "cross",cross = 5,error.fun=custom_error), ranges=rf_ranges,
         mc.control=list(mc.cores=16, mc.preschedule=T),confusionmatrizes=T )

save(tuned.r10,file='..//Hourly Voltage Data//tunedTrees10.RData')



# TEST ACCURACY OF TREES  ----------------------------------------------------
load('..//Hourly Voltage Data//tunedTrees10.RData')
#best parameters: ntree2, mtry 20

plot(tuned.r10)

# estimate variable importance
modeledRF = randomForest(form, data=test_join_holder.train, ntree=2,importance=T) 
varImpPlot(modeledRF)

#library(caret)
# read in stored regression trees look at performance 
tune.number = tuned.r10
# Get best parameters
tune.number$best.parameters
tune.number$best.confusionmatrizes
plot(tuned.r10)


# Store the best model 
best.model = tune.number$best.model
# Predict out of sample and create confusion matrix
predictions = predict(best.model, test_join_holder.test)
table.random.forest = table(test_join_holder.test$lightsout, predictions)
table.random.forest
#confusionMatrix(predictions, test_join_holder.test$lightsout)
    0     1      2
0 1741    0     0 
1    0   1877   52
2    0   48      4


# next function gives a graphical depiction of the marginal effect of a variable on the class probability (classification) or response (reg$
partialPlot(best.model, test_join_holder, dnb, "2")
partialPlot(best.model, test_join_holder, illum, "2")
partialPlot(best.model, test_join_holder, phase, "2")
partialPlot(best.model, test_join_holder, md_dnb, "2")


# Calculate error rate
error.rate = 1 - sum(diag(as.matrix(table.random.forest))) / sum(table.random.forest)
error.rate



# Sensativity Analysis ----------------------------------------------------
# test impact of reduced number of temporal observations on classifying outages

for(i in 1:15){
  set.seed(i)  
  rf_ranges = list(ntree=seq(1,5,1),mtry=10:30)

  tuned.model_temp_sense = tune(randomForest, train.x = form, data = test_join_holder.train,
                 tunecontrol = tune.control(sampling = "cross",cross = 5,error.fun=custom_error), ranges=rf_ranges,
                 mc.control=list(mc.cores=5, mc.preschedule=T),confusionmatrizes=T )

  # Store the best model 
  best.model = tuned.model_temp_sense$best.model
  # Predict out of sample and create confusion matrix
  predictions = predict(best.model, test_join_holder)
  table.random.forest = table(test_join_holder$lightsout, predictions)
  print(i)
  print(table.random.forest)
}

[1] 1
predictions
0    1    2
0 1741    0    0
1    0 1929    0
2    0    3   49
[1] 2
predictions
0    1    2
0 1741    0    0
1    0 1929    0
2    0   10   42
[1] 3
predictions
0    1    2
0 1741    0    0
1    0 1929    0
2    0    6   46
[1] 4
predictions
0    1    2
0 1741    0    0
1    0 1929    0
2    0    5   47
[1] 5
predictions
0    1    2
0 1741    0    0
1    0 1929    0
2    0    5   47
[1] 6
predictions
0    1    2
0 1741    0    0
1    0 1929    0
2    0    4   48
[1] 7
predictions
0    1    2
0 1741    0    0
1    0 1929    0
2    0    4   48
[1] 8
predictions
0    1    2
0 1741    0    0
1    0 1929    0
2    0    4   48
[1] 9
predictions
0    1    2
0 1741    0    0
1    0 1929    0
2    0    4   48
[1] 10
predictions
0    1    2
0 1741    0    0
1    0 1929    0
2    0    5   47
[1] 11
predictions
0    1    2
0 1741    0    0
1    0 1929    0
2    0    5   47
[1] 12
predictions
0    1    2
0 1741    0    0
1    0 1929    0
2    0    4   48


# PREDICT OUT OF SAMPLE ---------------------------------------------------------------

#remove(list=ls())
setwd('/groups/manngroup/India\ VIIRS/2015')

# load raster stacks
load('dnb_stack_wo_cld.RData')   
#load('zen_stack_wo_cld.RData')
#load('azt_stack_wo_cld.RData')
dnb_date_time = substr(names(dnb_stack),2,20)

# create sd_dnb and mn_dnb for each dnb time series
#sd_dnb_layer = calc(dnb_stack,function(x){sd(x,na.rm=T)} )
mn_dnb_layer = calc(dnb_stack,function(x){mean(x,na.rm=T)} )
md_dnb_layer = calc(dnb_stack,function(x){median(x,na.rm=T)} )


# read in moon phase (year doy time moon_illum_frac moon_phase_angle)
phase = read.csv('moon_info.csv')
names(phase)=c('year','doy','time', 'illum', 'phase')
phase$date_time = paste(phase$year,sprintf('%03d',(phase$doy)),'.',phase$time,sep='')

# load the tuned trees
load(file='..//Hourly Voltage Data//tunedTrees10.RData')
best.model = tuned.r10$best.model


# remove cloud cells multicore  returns NA but runs fast!
library(foreach)
library(doParallel)
registerDoParallel(32)

# iterate through each date_time collect dnb, add sd_dnb illum phase

prediction_list = foreach(i = 1:length(dnb_date_time), .inorder=T,.packages=c('raster','e1071'),
	.errorhandling="remove") %dopar% {
	# packages must be loaded to each core, error handling removes output when 'next' is invoked, doesn't work in dopar
	#for (i in 1:1:length(dnb_date_time)){
	date_time = dnb_date_time[i]
	print(date_time)
	# find the dnb that matches date_time
	hold_stack_dnb = dnb_stack[[which(dnb_date_time %in% date_time)]]  
#	hold_stack_zen = zen_stack[[which(dnb_date_time %in% date_time)]]  
#        hold_stack_azt = azt_stack[[which(dnb_date_time %in% date_time)]]  
	# store illumination
	# avoid missing data
	if(length(phase[phase$date_time == date_time,'illum'])==0){next}
	hold_stack_illum = hold_stack_dnb
	hold_stack_illum[] = phase[phase$date_time == date_time,'illum']
	# store phase
        hold_stack_phase = hold_stack_dnb
        hold_stack_phase[] = phase[phase$date_time == date_time,'phase']
        # stack all data with appropriate names
	data_stack = stack(hold_stack_dnb,
		sd_dnb_layer,hold_stack_illum,hold_stack_phase,mn_dnb_layer,md_dnb_layer)
	names(data_stack) = c('dnb','sd_dnb','illum','phase','mn_dnb','md_dnb')
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
foreach(i=1:dim(prediction_stack)[3],.inorder=T,.packages=c('raster')) %dopar% { 
	print(i)
  
  # ASSIGNMENT IN A DOPAR LOOP DOESN'T WORK!!!! THIS IS WORKING
  
	prediction_stack[[i]][cld_stack[[i]]>0]=NA}
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
setwd('/groups/manngroup/India\ VIIRS/2015')
load(file='prediction_stack_wo_cld.RData')

outage_count = calc((prediction_stack==2),function(x){sum(x,na.rm=T)})
cloud_count = calc(is.na(prediction_stack),function(x){sum(x,na.rm=T)})
leng_days_raster = prediction_stack[[1]]  # create a raster that holds the total # of images
leng_days_raster[]=dim(prediction_stack)[3]
percent_outage = outage_count / (leng_days_raster-cloud_count)
writeRaster(percent_outage,'predictions/percent_outage.tif', overwrite=T)
plot(percent_outage)

# Find uninhabited areas (sum of prediction stack should be ~0)
pred_sum_uninhab = sum(prediction_stack,na.rm=T)
writeRaster(pred_sum_uninhab,'predictions/pred_sum_uninhab.tif', overwrite=T)
plot(pred_sum_uninhab)



# compare to actual reliability 
load('..//Hourly Voltage Data//AllData.RData')
percent_outage = raster('predictions/percent_outage.tif')

# set up data for classifier (define what an outage is)
# read in joined dnb and voltage data where not missing
test_join_holder=test_join_holder[!is.na(test_join_holder$dnb)&!is.na(test_join_holder$voltage),]

# find # of days of outages  voltage below 100v is considered an outage
hist(test_join_holder$voltage)
test_join_holder$lightsout = NA
test_join_holder$lightsout[ test_join_holder$voltage<100  ] = 1   # outage
test_join_holder$lightsout[ grepl('uninhabit',test_join_holder$location)  ] = 0   # uninhabited areas assigned value of zero b/c no lights
head(test_join_holder,20)
count_out = aggregate(lightsout~location, data = test_join_holder, function(x){sum(x,na.rm=T)}) # count outages

#count number of days with obs
test_join_holder$obs = 1 
out_of = aggregate(obs~location, data = test_join_holder, function(x){sum(x,na.rm=T)})
reliability = join(count_out,out_of)
reliability$percent_outage = reliability$lightsout / reliability$obs

#extract dnb estimates of percent_outage to points
load('locationsXY.RData')
locations_reliability = locations
names(locations_reliability) = c('STATE','DISTRICT.CITY','location','Ag.Rural')
locations_reliability = locations_reliability[!(locations_reliability@data$location=='Chandikapur'|
	locations_reliability@data$location=='Kanheri Sarap'),]
locations_reliability@data = join(locations_reliability@data,reliability)
locations_reliability@data = cbind(locations_reliability@data ,extract(percent_outage,locations_reliability, df=T))
names(locations_reliability)=c("STATE","DISTRICT.CITY","location","Ag.Rural","lightsout",
	"obs","percent_outage","ID","percent_outage_estimate")


#visualize comparison
plot_data = locations_reliability@data
lm1 = lm(percent_outage~percent_outage_estimate,data=plot_data[-c(13),])

# test if intercept = 0 and slope = 1
(summary(lm1)$coefficients[1]-0)/summary(lm1)$coefficients[3]
(summary(lm1)$coefficients[2]-1)/summary(lm1)$coefficients[4]
#0.15 without sd or mn
ggplot(plot_data[-c(13),],aes(x=percent_outage_estimate, y=percent_outage))+geom_point(colour='grey30',size=3)+
        xlab('Predicted')+ylab('Actual')+ geom_abline(intercept = 0, slope = 1,size=1.1)+
         stat_smooth(method="lm", se=T,size=1.5,linetype = 2)+
	coord_cartesian(xlim = c(-.008,0.17), ylim = c(-.008,0.17))+ 
	annotate("text", x = 0.13, y = 0.02, label = 
	paste(round(summary(lm1)$coefficients[1],3),'+',round(summary(lm1)$coefficients[2],3)))+
 	annotate("text", x = 0.125, y = 0.013, label =paste('R2 =',round(summary(lm1)$adj.r.square,3)))+
        annotate("text", x = 0.118, y = 0.006, label =
		paste('(',round((summary(lm1)$coefficients[2]-1)/summary(lm1)$coefficients[4],2),')',sep=''))


# Plot time series for a location
unique(dnb_values$location)
loc2plot = 'Saharkar Nagar'
tilak = dnb_values[dnb_values$location==loc2plot & !is.na(dnb_values$dnb),]
ggplot(tilak,aes(x=as.Date(date.time2), y=dnb))+geom_point(colour='grey30',size=2)+
	xlab('2015')+ylab('Day Night Band Radiance')






# Stable Lights Map -------------------------------------------------------

# read in data 
#rm(list=ls())
setwd('/groups/manngroup/India\ VIIRS/2015')
load('dnb_stack_wo_cld.RData')    # start here cloud free images stored here.

# calculate stats
# mean layer
mn_dnb_layer = calc(dnb_stack,function(x){mean(x,na.rm=T)} )
plot(log(mn_dnb_layer))
writeRaster(mn_dnb_layer,'predictions/mn_dnb_layer.tif',overwrite=T)
# median layer
md_dnb_layer = calc(dnb_stack,function(x){median(x,na.rm=T)} )
plot(log(md_dnb_layer))
writeRaster(md_dnb_layer,'predictions/md_dnb_layer.tif', overwrite=T)
# quantile layer
qt_dnb_layer = calc(dnb_stack,function(x){quantile(x,na.rm=T)} )
names(qt_dnb_layer)=paste('Precentile',c(0,25,50,75,100))
writeRaster(qt_dnb_layer,'predictions/qt_dnb_layer.tif', overwrite=T)



# Extract reliability to shapefile ----------------------------------------

#mask out non-electrified
qplot(md, data=data.frame(md=getValues(md_dnb_layer)), geom="histogram",binwidth=1e-10)+coord_cartesian(xlim=c(0,5e-9))
# use 1e-9 as cutoff for non-electrified for now 

library(rgdal)
library(maptools)
proj = proj4string(readOGR('..//Boundary Data//Census_2011//Country',"2011_Dist"))
districts = readShapePoly("..//Boundary Data//Census_2011//Country//2011_Dist.shp", proj4string=CRS(proj))
plot(percent_outage)
plot(districts,add=T)

percent_outage[percent_outage<1e-9]=NA  # mask out non-electrified
mn_percent_outage = extract(percent_outage, districts, fun = mean, na.rm = T, df = T) # don't remove NA b/c of districts w few obs
mn_percent_outage[mn_percent_outage[,2]==0,2]=NA   #if mn_percent_outage is 0 then it doesn't have any data, na out
districts@data$mn_per_outage = mn_percent_outage[,2]
write.csv(districts@data,'.//predictions//Districts_Outage.csv')
writePolyShape(districts,'.//predictions//Districts Outage.shp')










# Write out files to KML ------------------------------------------------
# this works but is pretty low resolution
KML(md_dnb_layer,'predictions/md_dnb_layer',blur=50,overwrite=T)

KML(percent_outage,'predictions/percent_outage',blur=50,overwrite=T, col=rev(rainbow(25)))


# 
# 
# # Outlier detection  ----------------------------------------------------------
# # https://blog.twitter.com/2015/introducing-practical-and-robust-anomaly-detection-in-a-time-series
# 
# 
# 
# # install.packages("devtools")
# # devtools::install_github("twitter/AnomalyDetection")
# library(AnomalyDetection)
# 
# help(AnomalyDetectionTs)
# data_urban$date.time = strptime(data_urban$date.time,"%Y%j.%H%M" )
# data3= data_urban[,c('date.time','dnb')]
# data3$dnb = data3$dnb*1e9
# row.names(data3) = 1:dim(data3)[1]
# res = AnomalyDetectionTs(data3, max_anoms=0.03, direction='pos', plot=TRUE,
#         piecewise_median_period_weeks=8,longterm=T)
# res$plot
# 
# # impulse sateration https://cran.r-project.org/web/packages/gets/gets.pdf from felix
# library(zoo)
# library(gets)
# data4 = zoo(data_urban[,'dnb'],data_urban[,'date.time'])
# 
# # add a time trend
# t = zoo(1:length(data4[,1]),data_urban[,'date.time'])
# 
# sat = isat(data4, t.pval=min(0.01,(1/length(data4[,1]))),mxreg=t)
# plot(data4)
# 





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





