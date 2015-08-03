# find list of swaths using:

# download viirs data       
# https://ladsweb.nascom.nasa.gov/data/search.html
# NPP_VDNE_L1 3144 (DNB data and associated Lat/Lon) 
# NPP_VMAE_L1 (Lat/Lon for cloud mask)
# NPP_CMIP_L2 (cloud mask)

# Bounding Box
# http://www.mapdevelopers.com/geocode_bounding_box.php
# North Latitude: 22.028131 South Latitude: 15.602412 East Longitude: 80.890924 West Longitude: 72.659363


# copy returned search to txt file

library(RCurl)
library(stringr)
library(lubridate)


# DNB Read in swath names -----------------------------------------------------
  setwd('C://Users//mmann//Google Drive//India Night Lights//')
  data = read.table('Data//VIIRS Swath//swath_june_2012_NPP_VDNE.txt',sep="\t", header=T,stringsAsFactors=F)
  
  # List files
  data$File.Name

  # Store all days of the year
  DOYs = unique(as.numeric(substr(lapply(strsplit(data$Time,'-'), function(x) x[2]),1,3)))
 
  # placeholder for download tracking
  data$downloaded = 'N'

# FTP  ---------------------------------------------------------------------
  for(DOY in DOYs){
    # iterate through days of the year and download data
    ftp1 = paste('ftp://ladsweb.nascom.nasa.gov//allData//3144//NPP_VDNE_L1//2012//',DOY,'//',sep='')
  
    # check files on FTP
    ftpfilelist = getURL(ftp1,dirlistonly=T)
    ftpfilelist = unlist(str_split(ftpfilelist,"\r\n"))
    ftpfilelist = ftpfilelist[!ftpfilelist == ""]
    head(ftpfilelist)
    
    # find matching date and time from data and store full name of matches
    data_DT   = unlist(lapply(str_split(data$File.Name,'.P1'),function(x) x[1]))
    ftp_DT    = unlist(lapply(str_split(ftpfilelist,'.P1'),function(x) x[1]))
    ftptopull = ftpfilelist[ftp_DT %in%  data_DT]
    print(paste("Downloading ",ftptopull))
    
    # create URLS
    dataURL = paste(ftp1,ftptopull,sep='')
    
    # loop over URLS and download files 
    for(i in 1:length(dataURL)){
      print(paste('working on file ',i))
      
      out_name = paste('Data//VIIRS Swath//',ftptopull[i],sep='')
      if(!file.exists(out_name) | url.exists(dataURL[i])){
        
        # download.file(dataURL[i],out_name)
        bin = getBinaryURL(url = dataURL[i])
        writeBin(bin, out_name)  
        # Mark off as downloaded
        namefind = str_sub(ftptopull[i],1, as.numeric(str_locate(ftptopull[i],'.P1')[,'start']-1))
        data$downloaded[data_DT %in%  namefind] = 'Yes'
        
      }else{if(file.exists(out_name)){
            print('File Exists - Skipping')
            # Mark off as skipped
            data$downloaded[data_DT %in%  namefind] = 'Skipped'}else{
            print(' URL doesnt Exist') 
            data$downloaded[data_DT %in%  namefind] = 'No URL'
            }}
    }  
  }
print(data)




# Day night band method 2  ------------------------------------------------
  # List of available files downloaded from http://peate.ssec.wisc.edu/flo/search
  
  # copy returned search to txt file
  
  library(RCurl)
  library(stringr)
  library(lubridate)
  
  # DNB Read in swath names -----------------------------------------------------
  setwd('C://Users//mmann//Google Drive//India Night Lights//')
  data = read.table('Data//VIIRS Swath//swath_june_2012_NPP_VDNE_PEATE.txt',sep="\t", header=F,stringsAsFactors=F)
  # remove header
  data=data[-1,]
  data=data[!data=="EOF"]
  data = data.frame(File.Name = data)

  # pull date and time
  data$date=  gsub(pattern = "(.*_d)(.*)(_t.*)", replacement = "\\2", data$File.Name)
  data$time=  gsub(pattern = "(.*_t)(.*)(_e.*)", replacement = "\\2", data$File.Name)
  data$datetime = paste(data$date,data$time,sep='-')
  data$datetime = strptime(data$datetime,"%Y%m%d-%H%M%S")

  # Store all days of the year
  DOYs=unique(as.numeric(strftime(data$datetime, format = "%j")))
  
 
  # Round up to nearest five minutes (to match 5 minute product on FTP)
  data$datetime_RND = as.POSIXlt(floor(as.double(data$datetime)/(60*5) )*60*5,origin=(as.POSIXlt('1970-01-01')))
 
  # Create time stamp to match files on ftp://ladsweb.nascom.nasa.gov
  data$YJ= as.numeric(strftime(data$datetime_RND, format = "%Y%j"))
  data$HM= paste(hour(data$datetime_RND),format(data$datetime_RND,"%M"),sep='')
  data$YJHM = paste('NPP_VDNE_L1.A',data$YJ,'.',data$HM,sep='') 
  
  # placeholder for download tracking
  data$downloaded = 'N'
    

  # FTP  ---------------------------------------------------------------------
  for(DOY in DOYs){
    # iterate through days of the year and download data
    ftp1 = paste('ftp://ladsweb.nascom.nasa.gov//allData//3110//NPP_VDNE_L1//2012//',DOY,'//',sep='')
    
    # check files on FTP
    ftpfilelist = getURL(ftp1,dirlistonly=T)
    ftpfilelist = unlist(str_split(ftpfilelist,"\r\n"))
    ftpfilelist = ftpfilelist[!ftpfilelist == ""]
    head(ftpfilelist)
    
    # find matching date and time from data and store full name of matches
    data_DT   = data$YJHM
    ftp_DT    = unlist(lapply(str_split(ftpfilelist,'.P1'),function(x) x[1]))
    ftptopull = ftpfilelist[ftp_DT %in%  data_DT]
    print(paste("Downloading ",ftptopull))
    
    # create URLS
    dataURL = paste(ftp1,ftptopull,sep='')
    
    # loop over URLS and download files 
    for(i in 1:length(dataURL)){
      print(paste('working on file ',i))
      
      out_name = paste('Data//VIIRS Swath//',ftptopull[i],sep='')
      if(!file.exists(out_name) | url.exists(dataURL[i])){
        
        # download.file(dataURL[i],out_name)
        bin = getBinaryURL(url = dataURL[i])
        writeBin(bin, out_name)  
        # Mark off as downloaded
        namefind = str_sub(ftptopull[i],1, as.numeric(str_locate(ftptopull[i],'.P1')[,'start']-1))
        data$downloaded[data_DT %in%  namefind] = 'Yes'
        
      }else{if(file.exists(out_name)){
        print('File Exists - Skipping')
        # Mark off as skipped
        data$downloaded[data_DT %in%  namefind] = 'Skipped'}else{
          print(' URL doesnt Exist') 
          data$downloaded[data_DT %in%  namefind] = 'No URL'
        }}
    }  
  }
  print(data)



#     strftime(strftime("2012-06-01"), format = "%j")
#     strftime(strftime("2012-06-30"), format = "%j")

ftp2 = 'ftp:/ladsweb.nascom.nasa.gov/allData/3110/NPP_VMAE_L1/2012/153/'
ftp3 = 'ftp:/ladsweb.nascom.nasa.gov/allData/3110/NPP_CMIP_L2/'