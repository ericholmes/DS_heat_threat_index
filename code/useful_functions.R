## Functions

downloadNWIS <- function(site_no, parameterCd, startDT, endDT){
  temp <- read.table(paste0("https://waterservices.usgs.gov/nwis/iv/?sites=", site_no,
                            "&parameterCd=", parameterCd,
                            "&startDT=", startDT,
                            "T00:00:00-07:00&endDT=", endDT,
                            "T00:00:00-07:00&siteStatus=all&format=rdb"),
                     skip = 30, stringsAsFactors = F, col.names = c("Agency", "Site_no", "Date", "Time", "tz", "Param_val", "qc_code"))
  temp$parameterCd <- parameterCd
  temp$Datetime <- as.POSIXct(paste(temp$Date, temp$Time), format = "%Y-%m-%d %H:%M")
  return(temp)
}

downloadCDEC <- function(site_no, parameterCd, startDT, endDT){
  temp <- read.table(paste("http://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=", site_no, "&SensorNums=",
                           parameterCd, "&dur_code=", "H", "&Start=", startDT, "&End=", endDT, sep=""),
                     header=FALSE, sep=",", skip=1, stringsAsFactors = F)
  
  temp <- temp[,c(5,7)]
  colnames(temp) <- c("Datetime", "Param_val")
  temp$site_no <- site_no
  temp$parameterCd <- parameterCd
  temp$Datetime <- as.POSIXct(temp$Datetime, format = "%Y%m%d %H%M")
  
  return(temp[,c("site_no", "Datetime", "parameterCd", "Param_val")])
}

