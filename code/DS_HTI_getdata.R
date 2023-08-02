library(WDLimportr)
library(tidyverse)

source("code/useful_functions.R")

# Download data from DWR water data library
continuous_station <- "B9156000" ##Lisbon station ID

test_parameters <- c("Water Temperature")

LIS <- WDLimportr::ImportContinuousRaw(stationNumber = "B9156000", parameters = test_parameters)
colnames(LIS) <- c("Datetime", "Param_val")
LIS$Site_no <- "LIS"
LIS$Param_val <- as.numeric(LIS$Param_val)

# Download data from CDEC
RVB <- downloadCDEC("RVB", "25", "1989-1-1","2023-9-30")
RVB$Param_val <- as.numeric(RVB$Param_val)
RVB <- RVB[is.na(RVB$Datetime) == F,]

RVB$Param_val <- (RVB$Param_val - 32)*5/9
RVB <- RVB[RVB$Param_val > 5 & RVB$Param_val < 32, ]
# Append data sources

wtempdata <- rbind(LIS[, c("Site_no", "Datetime", "Param_val")],
                   RVB[, c("Site_no", "Datetime", "Param_val")])
  
# Save as an Rdata file
save(wtempdata, file = "data/DS_HTI_Wtemp_data.Rdata")
