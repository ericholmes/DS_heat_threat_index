library(WDLimportr)
library(tidyverse)

source("code/useful_functions.R")

##Download data from DWR water data library
continuous_station <- "B9156000" ##Lisbon station ID

test_parameters <- c("Water Temperature")

LIS <- WDLimportr::ImportContinuousRaw(stationNumber = "B9156000", parameters = test_parameters)

##Download data from CDEC
RVB <- downloadCDEC("RVB", "25", "1989-1-1","2023-9-30")

##Save as an Rdata file
save(LIS, RVB, file = "data/DS_HTI_Wtemp_data.Rdata")
