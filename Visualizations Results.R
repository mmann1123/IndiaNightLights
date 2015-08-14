library(raster)
library(ggplot2)
library(scales)
library(reshape2)
library(rgeos)


# this scripts reads in raster files exported from grid_viirs_data (3).R

# read in data ------------------------------------------------------------
  setwd('C:\\Users\\mmann\\Desktop\\NightTimeData\\')
  
  # pull available files
  files = dir(pattern = '.tif')
  cld = files[grep('cld',files)]
  dnb = files[grep('dnb',files)]
  
# create raster stacks  & extract data ----------------------------------------------------
  
  cld_stack = stack(cld)
  dnb_stack = stack(dnb)

  #e <- extent(76.5,77.5,19.5,20.5)
  #dnb_stack <- crop(dnb_stack,e)
  #cld_stack <- crop(cld_stack,e)
  windows()
  plot(dnb_stack)
  windows()
  plot(cld_stack)
  
  
# Extract dates -----------------------------------------------------------
  time_stamp_dnb = gsub(x=names(dnb_stack),pattern = "(.*X)(.*)(.*_dnb_v3)",replacement = "\\2")
  time_stamp_dnb = strptime(time_stamp_dnb,"%Y%j.%H%M")
  time_stamp_cld = gsub(x=names(cld_stack),pattern = "(.*X)(.*)(.*_cld_v3)",replacement = "\\2")
  time_stamp_cld = strptime(time_stamp_cld,"%Y%j.%H%M")

  
  # not all stacks have same dates
  all.equal(time_stamp_dnb,time_stamp_cld)
  # limit stacks to common elements
  stack2 = dnb_stack
  stack2[[c(T,rep(F,dim(stack2)[3]-1))]]
  stack2[[1:2]]
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
  
  locations = read.csv('C://Users/mmann/Google Drive/India Night Lights/Data/MH-ESMI-Locations-Lat-Long-Overpass-Cuts-May-2015-ag.csv')
  
   #jumda = data.frame(lat=20.010094,lon=77.044271)
  coordinates(locations)= ~LON+LAT
  proj4string(locations) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  
  dnb_values = extract(dnb_stack,locations,buffer=1.2e3,fun= function(x) mean(x,na.rm=T), df=T)
  cld_values = extract(cld_stack,locations,buffer=1.2e3,fun= function(x) mean(x,na.rm=T), df=T)
  
  # NA out values with clouds  (zero = confident no clouds, 1 = probably no clouds)
  dnb_values[cld_values>1]=NA
  
  # Create plots ------------------------------------------------------------

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
  ggplot(data,aes(x=date,y=log(value*1e9)))+geom_point()+scale_x_date(labels = date_format("%m-%d"))+ 
      facet_wrap( ~ LOCATION+Ag.Rural, ncol = 2,scales = 'free_y')  
  
  
  write.csv(data,'C:/Users/mmann/Google Drive/India Night Lights/Data/MH-ESMI-Locations-DNB-output.csv')
  
  
  
  #   proj_locations = spTransform( locations, CRS( "+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=7 +lat_2=-32 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs " ) )  #http://spatialreference.org/ref/esri/102029/
  #   locations_buffer = gBuffer(proj_locations,capStyle = 'round',width = 3000)
  #   locations_buffer = spTransform( locations_buffer, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") )  #http://spatialreference.org/ref/esri/102029/
  #   
  #   dnb_buffer_values =  extract(dnb_stack,locations_buffer)
  #   